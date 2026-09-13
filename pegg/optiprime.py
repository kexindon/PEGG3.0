"""
OptiPrime scoring for pegRNAs designed by pegg.

OptiPrime (Hsu et al., Nat Biotechnol 2026, doi:10.1038/s41587-026-03261-7) is a
mechanistic model of prime editing efficiency. Its score is the predicted
fraction of alleles carrying the intended edit, so a score of 0.30 means ~30%
editing -- unlike PEGG2_Score/RF_Score, which are arbitrary scales.

Two things about how it is used here:

* The model's own design script (DESIGN_PE.py) runs a combinatorial search over
  RTT/PBS lengths and silent edits, which would replace pegg's designer. That is
  not what we want: pegg has already chosen the pegRNAs. The model's actual
  input contract is only five columns -- spacer, rtt, pbs, full_unedited,
  full_edited (see _REQUIRED_COLUMNS in optiprime's scripts/pe/pe_utils.py) --
  so arbitrary pegRNAs can be scored directly, with nothing pruned.

* Scoring is done in HeLa (MMR-proficient) parameters, which is OptiPrime's own
  default. The paper's recommendation is explicit: "high-performing pegRNAs
  optimized using its HeLa cell parameters tend to also maintain their
  performance in HEK293T cells, while the converse is not true. Therefore, we
  recommend that researchers interested in designing pegRNAs with OptiPrime use
  its default settings for typical applications." Every prospective validation
  in the paper -- primary human T cells, patient fibroblasts, MEFs, in vivo
  mouse brain -- used these defaults. HEK293T is partially MMR deficient, and
  designing against it yields pegRNAs that need not transfer.

OptiPrime requires Python 3.11 and its own jax/flax stack, which cannot coexist
with pegg's pinned numpy<2 / scikit-learn 1.1.1 environment. It is therefore run
in a subprocess against a separate interpreter; see score() for configuration.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import Bio.Seq


#OptiPrime's unedited sequence starts 4 nt upstream of a 20 nt protospacer;
#scripts/pe/pe_constants.py calls this PS20_OFFSET.
PS20_OFFSET = 4

#Group factors in the released weights are named <lab>_<cell line>. Only these
#twelve exist, and a name that is not among them is NOT an error in OptiPrime:
#rx_model.load_rates() leaves the group factor at zero and logs at INFO level,
#so a typo silently yields a plausible-looking score from an undefined MMR
#context. Since cell type is how MMR status enters the model, that is a silent
#correctness failure, and it is validated against this list instead.
VALID_GROUPS = (
    'Kim_A549', 'Kim_DLD1', 'Kim_HCT116', 'Kim_HEK293T', 'Kim_HeLa',
    'Kim_MDA-MB-231', 'Kim_NIH3T3', 'Liu_HEK293T', 'Liu_HeLa',
    'Schwank_HEK293T', 'Schwank_K562', 'Schwank_U2OS',
)

#HeLa, MMR-proficient -- OptiPrime's default, see module docstring.
DEFAULT_GROUP = 'Liu_HeLa'

#Each worker process pays a fixed startup cost (import jax, build the model, JIT)
#of roughly 20 s, so there is no point splitting the input below this.
_MIN_SHARD = 500


class OptiPrimeError(RuntimeError):
    """Raised when OptiPrime is unavailable or fails to produce scores."""


def _env_or(value, var):
    if value is not None:
        return value
    return os.environ.get(var)


def optiprime_available(optiprime_src=None, optiprime_python=None):
    """
    Whether OptiPrime can be run, i.e. both its source tree and a suitable
    interpreter are configured and present. Lets callers switch the feature on
    only when the environment supports it, rather than failing mid-run.

    Parameters
    -----------
    optiprime_src
        *type = str or None*

        Path to a clone of github.com/alvin-hsu/optiprime-src. Falls back to the
        OPTIPRIME_SRC environment variable.

    optiprime_python
        *type = str or None*

        Python interpreter with OptiPrime's dependencies installed (jax, flax,
        rs3, ViennaRNA). Falls back to OPTIPRIME_PYTHON.
    """
    src = _env_or(optiprime_src, 'OPTIPRIME_SRC')
    py = _env_or(optiprime_python, 'OPTIPRIME_PYTHON')
    if not src or not py:
        return False
    return os.path.isdir(src) and (os.path.isfile(py) or shutil.which(py) is not None)


def _rc(seq):
    return str(Bio.Seq.Seq(seq).reverse_complement())


def _rows_from_pegdf(df):
    """
    Builds OptiPrime's five input columns from pegg's pegRNA table.

    Returns (rows, errors), each a list the same length as df. rows[i] is None
    where the pegRNA could not be expressed in OptiPrime's coordinates, with the
    reason in errors[i].

    Geometry. pegg stores, per pegRNA, the PAM-strand sequence implicitly via
    Protospacer_30 = seq_F[PAM_start-24 : PAM_start+6], so its first 24 nt are
    exactly OptiPrime's 4 nt pad + 20 nt protospacer. The nick sits at
    PAM_start-3, i.e. 21 nt into that window (PS20_OFFSET+17). RTT and PBS are
    stored reverse-complemented relative to the PAM strand, which is the same
    orientation OptiPrime expects, so they are passed through as-is.

    full_unedited/full_edited are the PAM-strand target either side of the edit.
    They are truncated to 25+RTT_length as DESIGN_PE.py does, because that sets
    hom_len, which feeds the MMR model; passing the full context changes the
    score.
    """
    rows = []
    errors = []

    has_bystanders = 'has_silent_bystander' in df.keys()

    #import here rather than at module scope to avoid a circular import:
    #prime imports this module.
    from .prime import _alt_context_from_RTT

    for _, val in df.iterrows():
        proto30 = val['Protospacer_30']
        rtt = val['RTT']
        pbs = val['PBS']

        if not isinstance(proto30, str) or len(proto30) < PS20_OFFSET + 20:
            rows.append(None)
            errors.append('protospacer context too short for OptiPrime')
            continue

        #--- the target sequence, in PAM-strand orientation -------------------
        #wt_w_context/alt_w_context are in the forward orientation of the input;
        #seq_F inside pegRNA_generator is the reverse complement of that when the
        #PAM is on the - strand, so flip to match.
        wt = val['wt_w_context']
        alt = val['alt_w_context']

        #A pegRNA carrying silent bystanders installs the intended edit PLUS the
        #bystanders, and the bystanders live in the RTT rather than in
        #alt_w_context. Rebuild the edited sequence from the RTT so the model
        #scores what the pegRNA actually makes -- this is what keeps the scoring
        #consistent with silent bystander designs.
        if has_bystanders and val['has_silent_bystander']:
            alt = _alt_context_from_RTT(val, wt, alt)
            if alt is None:
                rows.append(None)
                errors.append('could not place bystander RTT in context sequence')
                continue

        strand = val['PAM_strand']
        if strand == '-':
            wt_pam = _rc(wt)
            alt_pam = _rc(alt)
        else:
            wt_pam = wt
            alt_pam = alt

        #--- locate the protospacer window on the PAM-strand sequence ---------
        #Anchor on the 24 nt pad+protospacer, which is unchanged by the edit on
        #the unedited sequence. Require a unique hit so a repeated context does
        #not silently place the window in the wrong spot.
        anchor = proto30[:PS20_OFFSET + 20]
        hit = wt_pam.find(anchor)
        if hit < 0:
            rows.append(None)
            errors.append('could not locate protospacer in context sequence')
            continue
        if wt_pam.find(anchor, hit + 1) != -1:
            rows.append(None)
            errors.append('protospacer context is not unique')
            continue

        rtt_len = int(val['RTT_length'])
        #DESIGN_PE.py uses unedited[:25+rtt-len_delta] and edited[:25+rtt]
        len_delta = len(str(val['ALT'])) - len(str(val['REF']))
        end_u = 25 + rtt_len - len_delta
        end_e = 25 + rtt_len

        unedited = wt_pam[hit:hit + end_u]
        if len(unedited) < end_u:
            rows.append(None)
            errors.append('insufficient downstream context; increase context_size')
            continue

        #The edited sequence is anchored the same way. The pad+protospacer is 5'
        #of the nick, so on a PAM-disrupting edit it can itself differ; anchor on
        #the unedited coordinates instead by aligning the shared 5' flank.
        hit_e = alt_pam.find(anchor)
        if hit_e < 0:
            #the edit changed the protospacer window itself -- fall back to the
            #offset of the sequence 5' of it, which the edit cannot reach
            flank = wt_pam[max(0, hit - 20):hit]
            if not flank:
                rows.append(None)
                errors.append('could not anchor edited sequence')
                continue
            f = alt_pam.find(flank)
            if f < 0 or alt_pam.find(flank, f + 1) != -1:
                rows.append(None)
                errors.append('could not anchor edited sequence')
                continue
            hit_e = f + len(flank)

        edited = alt_pam[hit_e:hit_e + end_e]
        if len(edited) < end_e:
            rows.append(None)
            errors.append('insufficient downstream context; increase context_size')
            continue

        if unedited == edited:
            rows.append(None)
            errors.append('edited and unedited sequences identical')
            continue

        #OptiPrime works in RNA; format_pe_df does T->U itself, so DNA is fine.
        rows.append((proto30[PS20_OFFSET:PS20_OFFSET + 20], rtt, pbs,
                     unedited, edited))
        errors.append(None)

    return rows, errors


#Runs inside the OptiPrime interpreter. Kept as a string so pegg does not need
#to ship a file that only ever runs under a different Python.
_WORKER = r'''
import json, sys, tempfile
from functools import partial
from pathlib import Path

#Featurization is ~95% ViennaRNA partition-function folding: ScaffoldDefect
#folds the whole ~140 nt pegRNA five times per row. That is GIL-releasing C
#code, and RxDataset can pool it (multiprocess=True, which needs the 'fork'
#start method -- its worker looks up RxInput._instances, a registry populated by
#import side effects, so 'spawn' children start empty and raise KeyError).
#
#It is NOT used, deliberately. ViennaRNA's energy parameters are process-global
#and the inputs mutate them: HetBPP.pre_process() switches to DNA_Mathews2004
#and restores RNA_Turner2004 afterwards. Run serially, ScaffoldDefect (which is
#ordered before HetBPP) always folds under the RNA parameters; forked workers
#inherit whichever set was live when they were spawned, so the same pegRNA gets
#a different score depending on scheduling. Measured on 400 pegRNAs: max
#|difference| 0.057, and no row agreed to 1e-6. A 2.2x speedup is not worth
#scores that change between runs.

import jax
from jax import jit, vmap
import jax.numpy as jnp
from jax.random import PRNGKey
import numpy as np
import pandas as pd

from reaction.models import FlaxRateModel, LinearRateModel, SharedFlaxRateModel
from reaction.rx_dataset import RxDataset
from reaction.rx_graph import read_rxfile
from reaction.rx_module import RxModule
from reaction.rx_model import RxModel
from reaction.utils import get_param_sizes, get_param_idxs

from scripts.utils import deterministic_hash
from scripts.pe.models import EditModel, PegRNAMLP, SynModel
from scripts.pe.pe_utils import format_pe_df
from scripts.pe.pe_inputs import (PE_ON_INPUTS, SYN_INPUTS, EDIT_INPUTS,
                                  MMR_INPUTS, UneditedEncoding, EditedEncoding,
                                  HetBPP)

RATE_INPUTS = {'pe_on': PE_ON_INPUTS, 'syn': SYN_INPUTS, 'muts_off': EDIT_INPUTS,
               'rep_u': EDIT_INPUTS, 'rep_e': EDIT_INPUTS, 'mmr': MMR_INPUTS}
ALL_INPUTS = PE_ON_INPUTS + SYN_INPUTS + EDIT_INPUTS + MMR_INPUTS


def main():
    cfg = json.load(open(sys.argv[1]))
    df = pd.read_csv(cfg['in_csv'])
    weight_dirs = cfg['weight_dirs']
    group = cfg['group']
    time = cfg['time']

    rx_graph = read_rxfile(Path(cfg['graph_rx']))
    edited_idx = rx_graph.obs_names.index('edited')
    max_u = int(df['full_unedited'].str.len().max())
    max_e = int(df['full_edited'].str.len().max())

    def prep(_, d):
        d['scaffold_name'] = 'OG_F+E'
        d['cas9_type'] = 'PEmax-Cas9'
        d['cas9_pam'] = 'SpNGG'
        d['rt_name'] = 'PE2-RT'
        d['motif'] = 'tevoPreQ1'
        d['pe_type'] = 'PE2'
        d['pol3'] = True
        d['group'] = group
        d = format_pe_df(_, d)
        d['spacer_hash'] = d['spacer'].apply(deterministic_hash)
        d['pegrna_hash'] = d['pegrna'].apply(deterministic_hash)
        d['edit_hash'] = d['min_edit'].apply(deterministic_hash)
        return d

    with tempfile.TemporaryDirectory() as td:
        #multiprocess=False: see the note above -- parallel featurization races
        #on ViennaRNA's global energy parameters and changes the scores.
        ds = RxDataset(df=df.copy(), rate_plates=rx_graph.plate_names,
                       observables=rx_graph.obs_names,
                       json_path=Path(td) / 'ds.json',
                       preprocess_fn=prep, rx_inputs=ALL_INPUTS,
                       multiprocess=False)
        UneditedEncoding.set_max_len(ds._cache, max_u)
        EditedEncoding.set_max_len(ds._cache, max_e)
        HetBPP.set_max_lens(ds._cache, max_u, max_e)

        with (Path(weight_dirs[0]) / 'metadata.json').open() as f:
            meta = json.load(f)
        het = EditModel(**meta['edit_hparams'], max_u_len=max_u, max_e_len=max_e)
        rate_models = {
            'pe_on': FlaxRateModel(PegRNAMLP(4, 32), name='pe_on'),
            'syn': FlaxRateModel(SynModel(), name='syn'),
            'mmr': LinearRateModel(name='mmr'),
            'rep_u': SharedFlaxRateModel(het, name='rep_u', out_index=0),
            'rep_e': SharedFlaxRateModel(het, name='rep_u', out_index=1),
            'muts_off': SharedFlaxRateModel(het, name='rep_u', out_index=2)}
        for k, m in rate_models.items():
            m.init(inputs=ds.get_inputs(RATE_INPUTS[k], 0), key=PRNGKey(0))

        rx_module = RxModule(rx_graph=rx_graph,
                             param_sizes=get_param_sizes(rx_graph, ds),
                             init_name='unedited',
                             num_groups=len(ds.groups))
        model = RxModel(rx_graph=rx_graph, rx_module=rx_module, models=rate_models)
        model.init_rates()
        #jit(vmap(...)) rather than pmap: pmap maps over DEVICES, so on a CPU
        #box it scores one pegRNA per forward pass. vmap batches within the one
        #device, which is what PREDICT_PE.py does and is orders of magnitude
        #faster for library-scale input.
        apply_fn = jit(vmap(partial(model.full_apply, rngs=None, deterministic=True),
                            in_axes=(None, None, 0, 0, 0, 0, 0, 0)))

        bs = int(cfg.get('batch_size') or 256)
        full_rate_idxs = get_param_idxs(rx_graph, ds)
        rtt_lens_all = ds.df['rtt'].str.len().values
        gidx_all = ds.df['group_idx'].to_numpy(dtype=np.uint32)
        times = jnp.full(bs, time)

        n = len(ds)
        acc = np.zeros((len(weight_dirs), n))
        for w_i, wd in enumerate(weight_dirs):
            model.load_full(Path(wd), ds)
            for i in range(0, n, bs):
                sl = slice(i, i + bs)
                inputs = {k: ds.get_inputs(RATE_INPUTS[k], sl) for k in model.models}
                ridx = {k: v[sl] for k, v in full_rate_idxs.items()}
                rl = rtt_lens_all[sl]
                gi = gidx_all[sl]
                if i + bs > n:
                    pad = i + bs - n
                    inputs = {k: [jnp.pad(x, [(0, pad)] + (len(x.shape) - 1) * [(0, 0)])
                                  for x in v] for k, v in inputs.items()}
                    ridx = {k: jnp.pad(v, (0, pad)) for k, v in ridx.items()}
                    rl = jnp.pad(rl, (0, pad))
                    gi = jnp.pad(gi, (0, pad))
                preds, _ = apply_fn(model._rate_vars,
                                    {k: m.model_vars for k, m in model.models.items()},
                                    inputs, ridx, {}, {'syn': rl}, gi, times)
                p = np.asarray(preds[:, edited_idx])
                acc[w_i, i:min(i + bs, n)] = p[:min(bs, n - i)]

        #RxDataset drops rows (weight<=0 / NaN) and reindexes, so map scores back
        #by the row_id carried through rather than by position.
        out = pd.DataFrame({'row_id': ds.df['row_id'].values,
                            'OptiPrime_Score': acc.mean(axis=0)})
        out.to_csv(cfg['out_csv'], index=False)


main()
'''


def score(df, optiprime_src=None, optiprime_python=None, weight_dirs=None,
          graph_rx=None, group=DEFAULT_GROUP, time=4.0, batch_size=256,
          n_jobs=None, quiet=True):
    """
    Scores pegRNAs with OptiPrime, returning df with an OptiPrime_Score column.

    The score is the predicted fraction of alleles carrying the intended edit,
    in [0, 1]. Rows that could not be expressed in OptiPrime's coordinates get
    NaN, with the reason in OptiPrime_error.

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        pegRNAs generated by run().

    optiprime_src
        *type = str or None*

        Path to a clone of github.com/alvin-hsu/optiprime-src. Falls back to the
        OPTIPRIME_SRC environment variable.

    optiprime_python
        *type = str or None*

        Python interpreter with OptiPrime's dependencies. Falls back to
        OPTIPRIME_PYTHON.

    weight_dirs
        *type = list or None*

        Model weight directories. Default = all five of the released ensemble,
        found under <optiprime_src>/weights.

    graph_rx
        *type = str or None*

        Reaction graph. Default = <optiprime_src>/graphs/pe_model.rx.

    group
        *type = str*

        Training group whose parameters to score against, as <lab>_<cell line>.
        Default = 'Liu_HeLa', OptiPrime's own default and the setting its
        authors recommend for general use; see the module docstring. Must be one
        of VALID_GROUPS -- an unrecognised name would otherwise be accepted
        silently and scored with an undefined MMR context.

    time
        *type = float*

        Experiment duration passed to the model. Default = 4.0, as DESIGN_PE.py.

    batch_size
        *type = int*

        pegRNAs per forward pass. Default = 256. Larger is faster but uses more
        memory; the model pads every sequence to the longest in the whole input,
        so memory also grows with the longest RTT.

    n_jobs
        *type = int or None*

        How many worker processes to split the input across. Default = None,
        meaning one per CPU less one. Each worker pays ~20 s of startup, so the
        input is not split into shards smaller than _MIN_SHARD rows regardless.

        Separate processes rather than threads or a pool: featurization is
        single-threaded, and ViennaRNA's energy parameters are process-global
        state that the inputs mutate, so sharing an interpreter makes scores
        depend on scheduling. Sharding across processes is deterministic --
        agreement with single-process scoring is ~5e-7, i.e. float32 rounding.

    quiet
        *type = bool*

        Suppress OptiPrime's own stdout/stderr. Default = True.
    """
    src = _env_or(optiprime_src, 'OPTIPRIME_SRC')
    py = _env_or(optiprime_python, 'OPTIPRIME_PYTHON')

    if not src or not py:
        raise OptiPrimeError(
            'OptiPrime is not configured. Set optiprime_src/optiprime_python '
            '(or the OPTIPRIME_SRC/OPTIPRIME_PYTHON environment variables) to a '
            'clone of github.com/alvin-hsu/optiprime-src and a Python 3.11 '
            'interpreter with its requirements installed.')
    if not os.path.isdir(src):
        raise OptiPrimeError('OptiPrime source not found at %s' % src)

    if group not in VALID_GROUPS:
        raise OptiPrimeError(
            "unknown OptiPrime group %r. The released weights only carry group "
            "factors for: %s. An unrecognised name is not rejected by OptiPrime "
            "itself -- it silently scores with a zeroed group factor -- so it is "
            "checked here instead." % (group, ', '.join(VALID_GROUPS)))

    if weight_dirs is None:
        wroot = os.path.join(src, 'weights')
        weight_dirs = sorted(
            os.path.join(wroot, d) for d in os.listdir(wroot)
            if os.path.isdir(os.path.join(wroot, d)))
        if not weight_dirs:
            raise OptiPrimeError('no weight directories found under %s' % wroot)
    if graph_rx is None:
        graph_rx = os.path.join(src, 'graphs', 'pe_model.rx')

    rows, errors = _rows_from_pegdf(df)

    out = df.copy()
    out['OptiPrime_Score'] = np.nan
    out['OptiPrime_error'] = errors

    keep = [i for i, r in enumerate(rows) if r is not None]
    if not keep:
        return out

    in_df = pd.DataFrame(
        [rows[i] for i in keep],
        columns=['spacer', 'rtt', 'pbs', 'full_unedited', 'full_edited'])
    #carried through the model so scores can be mapped back to pegg's rows even
    #though RxDataset reindexes
    in_df['row_id'] = keep

    #Featurization is the bulk of the cost and is single-threaded, so split the
    #input across several worker processes. Separate processes rather than a
    #pool inside one: ViennaRNA's energy parameters are process-global and the
    #inputs mutate them (see the note above _WORKER), so sharing an interpreter
    #makes scores depend on scheduling. A fresh process per shard gets its own
    #copy of that state and stays deterministic -- measured against scoring the
    #whole set in one process, sharding agrees to 5e-7 (float32 rounding from
    #batch padding), versus 0.06 for a forked pool.
    if n_jobs is None or n_jobs <= 0:
        n_jobs = max(1, (os.cpu_count() or 2) - 1)
    #each worker loads the model and JIT-compiles, which is a fixed ~20 s, so
    #do not shard so finely that the overhead dominates
    n_jobs = max(1, min(int(n_jobs), (len(in_df) + _MIN_SHARD - 1) // _MIN_SHARD))

    with tempfile.TemporaryDirectory(prefix='pegg_optiprime_') as td:
        worker = os.path.join(td, 'worker.py')
        with open(worker, 'w') as f:
            f.write(_WORKER)

        env = dict(os.environ)
        env['PYTHONPATH'] = src + os.pathsep + env.get('PYTHONPATH', '')
        #JAX would otherwise grab every core in each worker, which oversubscribes
        #the machine once several are running; keep each to one thread and get
        #the parallelism from the shards instead.
        env.setdefault('JAX_PLATFORMS', 'cpu')
        if n_jobs > 1:
            env.setdefault('XLA_FLAGS', '--xla_force_host_platform_device_count=1')
            for var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                env.setdefault(var, '1')

        procs = []
        out_paths = []
        for k in range(n_jobs):
            shard = in_df.iloc[k::n_jobs]
            if not len(shard):
                continue
            in_csv = os.path.join(td, 'in_%d.csv' % k)
            out_csv = os.path.join(td, 'out_%d.csv' % k)
            cfg_path = os.path.join(td, 'cfg_%d.json' % k)
            shard.to_csv(in_csv, index=False)
            with open(cfg_path, 'w') as f:
                json.dump({'in_csv': in_csv, 'out_csv': out_csv,
                           'weight_dirs': list(weight_dirs), 'graph_rx': graph_rx,
                           'group': group, 'time': float(time),
                           'batch_size': int(batch_size)}, f)
            procs.append(subprocess.Popen(
                [py, worker, cfg_path], cwd=src, env=env,
                stdout=subprocess.DEVNULL if quiet else None,
                stderr=subprocess.PIPE if quiet else None))
            out_paths.append(out_csv)

        failures = []
        for p in procs:
            _, err = p.communicate()
            if p.returncode != 0:
                failures.append((p.returncode, (err or b'').decode('utf-8', 'replace')))

        if failures:
            code, msg = failures[0]
            raise OptiPrimeError(
                'OptiPrime scoring failed (%d of %d workers; first exit %d).\n%s'
                % (len(failures), len(procs), code, msg[-4000:]))

        missing = [p for p in out_paths if not os.path.isfile(p)]
        if missing:
            raise OptiPrimeError('OptiPrime produced no output for %d of %d shards'
                                 % (len(missing), len(out_paths)))

        scored = pd.concat([pd.read_csv(p) for p in out_paths], ignore_index=True)

    idx = out.index[scored['row_id'].to_numpy()]
    out.loc[idx, 'OptiPrime_Score'] = scored['OptiPrime_Score'].to_numpy()

    return out
