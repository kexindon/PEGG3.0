🎯 OptiPrime Scoring
=====================

New in version 3.1. PEGG can rank pegRNAs with `OptiPrime <https://github.com/alvin-hsu/optiprime-src>`_
(Hsu et al., Nature Biotechnology 2026) instead of its own heuristic score.

**We recommend using OptiPrime for library design.** PEGG's built-in ``PEGG2_Score`` is a hand-weighted
combination of design features -- GC content, PAM disruption, homology arm length, and so on. It is fast and needs
no extra installation, but the weights were chosen by hand rather than fitted to data. OptiPrime is a machine
learning model trained on large-scale measured prime editing outcomes, and it explicitly models mismatch repair
(MMR), which is often the dominant determinant of whether an edit is installed. For designing a library you intend
to order and screen, that difference matters.

**The only cost is time**, and it is bounded -- see :ref:`optiprime-runtime` below. There is no accuracy or
compatibility penalty: OptiPrime scoring is additive, and every column PEGG produced before is still produced.

.. note::
   OptiPrime is **optional**. PEGG installs and runs without it, and ``prime.run()`` behaves exactly as before
   unless you ask for OptiPrime explicitly.


Why a separate installation
****************************

OptiPrime cannot be a dependency of PEGG. PEGG pins ``numpy<2`` and runs on Python 3.9/3.10, because the
on-target scoring models are pickled with scikit-learn 1.1.1. OptiPrime is built on JAX and Flax, which require
Python 3.11+ and ``numpy>=2``. The two dependency sets cannot coexist in one environment.

PEGG therefore runs OptiPrime **in a subprocess against a separate interpreter**. You install OptiPrime once into
its own conda environment, tell PEGG where it lives, and PEGG shells out to it when scoring. JAX never enters
PEGG's own environment.


Installing OptiPrime
*********************

The upstream ``requirements.txt`` is unpinned and does not resolve. Use these versions, which are known to work:

.. code-block:: bash

   # 1. get the source
   git clone https://github.com/alvin-hsu/optiprime-src ~/optiprime/src

   # 2. a dedicated Python 3.11 environment
   conda create -p ~/optiprime/env python=3.11 -y
   conda activate ~/optiprime/env

   # 3. lightgbm and numpy from conda-forge -- pip builds fail without CMake
   conda install -c conda-forge lightgbm=3.3.5 numpy=1.26.4 -y

   # 4. the JAX stack
   pip install jax==0.4.35 jaxlib==0.4.35 flax==0.10.2 optax==0.2.4 pandas==2.2.3

   # 5. rs3 and sglearn pull in conflicting pins, so install without dependencies
   pip install --no-deps rs3 sglearn
   pip install seqfold

.. warning::
   Do not install ``jax>=0.10``. It requires ``numpy>=2``, which conflicts with the LightGBM path that ``rs3``
   depends on.

Two things to know about the upstream repository:

* It is missing files its own Dockerfile expects (``rs3.patch``, ``sglearn.patch``, ``graphs/pe_model_2.rx``).
  Only ``graphs/pe_model.rx`` ships, and that is the file PEGG uses.
* ``DESIGN_PE.py`` is the inference entry point. ``PREDICT_PE.py`` is a benchmarking script that needs labelled
  data and cannot be used for design. PEGG uses neither -- it calls the model directly.


Pointing PEGG at it
********************

Set two environment variables, or pass the equivalent parameters to ``prime.run()``:

.. code-block:: python

   import os
   os.environ['OPTIPRIME_SRC']    = os.path.expanduser('~/optiprime/src')
   os.environ['OPTIPRIME_PYTHON'] = os.path.expanduser('~/optiprime/env/bin/python')

Check that PEGG can find it:

.. code-block:: python

   from pegg import optiprime
   optiprime.optiprime_available()   # True when both paths resolve


Basic usage
************

Pass ``optiprime=True`` to ``prime.run()``:

.. code-block:: python

   import pandas as pd
   from pegg import prime

   mutations = pd.DataFrame({'SEQ': ['...ATGC(C/T)GCTA...']})

   pegRNAs = prime.run(
       mutations, 'PrimeDesign', None,
       pegRNAs_per_mut=10,
       optiprime=True,                 # score with OptiPrime
       rankby='OptiPrime_Score',       # and rank by it
   )

This adds two columns:

``OptiPrime_Score``
    Predicted editing efficiency, roughly 0--1. Higher is better. ``NaN`` where the pegRNA could not be scored.

``OptiPrime_error``
    Why a pegRNA was not scored, or ``NaN`` when it was. Scoring never fails silently: an unscored pegRNA always
    carries a reason.

To keep only pegRNAs above a threshold, use ``optiprime_cutoff``:

.. code-block:: python

   pegRNAs = prime.run(
       mutations, 'PrimeDesign', None,
       pegRNAs_per_mut=10,
       optiprime=True,
       optiprime_cutoff=0.1,           # drop anything below 0.1
   )

.. warning::
   ``run()`` scores **before** the downstream polyT and sensor filters that notebooks usually apply afterwards.
   A mutation that has 10 pegRNAs when ``run()`` returns can end up with none once those filters run. If coverage
   per mutation matters, check it after filtering, not before.


.. _optiprime-cell-line:

Choosing a cell line: this is an MMR decision
**********************************************

``optiprime_group`` selects which cell line's parameters the model scores with. **This is not a tuning knob --
it encodes mismatch repair status, and it changes scores substantially.** The same pegRNA can score 0.174 under
``Liu_HeLa`` and 0.407 under ``Liu_HEK293T``.

The reason is biological. HEK293T is partially MMR-deficient; HeLa is MMR-proficient. Prime edits that would be
reverted by mismatch repair survive in HEK293T and do not in HeLa.

Per Hsu et al., pegRNAs optimised under MMR-proficient (HeLa) parameters generally transfer to HEK293T, **but the
reverse is not true**. So:

.. tip::
   Design on ``Liu_HeLa`` (the default) unless you have a specific reason not to. It is the conservative choice
   and produces libraries that transfer across cell lines.

The released weights carry group factors for twelve groups:

.. code-block:: text

   Kim_A549       Kim_DLD1      Kim_HCT116   Kim_HEK293T
   Kim_HeLa       Kim_MDA-MB-231              Kim_NIH3T3
   Liu_HEK293T    Liu_HeLa
   Schwank_HEK293T               Schwank_K562  Schwank_U2OS

In practice only ``Liu_HeLa`` and ``Liu_HEK293T`` resolve, because OptiPrime's own ``DESIGN_PE.py`` hardcodes the
``Liu_`` prefix.

.. warning::
   **An unrecognised group name is not an error in OptiPrime.** Its rate loader leaves the group factor at zero and
   logs only at ``INFO`` level, so a typo returns a plausible-looking score that means nothing. PEGG validates the
   name against the twelve above and raises ``OptiPrimeError`` instead. Do not bypass this check.


.. _optiprime-runtime:

Runtime, and how to control it
*******************************

Scoring is the slow part of a PEGG run. The cost has two parts:

**Start-up (fixed).** Each worker process imports JAX, builds the model and JIT-compiles it. This costs roughly
**14 seconds** and does not depend on how many pegRNAs you score. It is paid once per worker.

**Per pegRNA (marginal).** After warm-up, roughly **0.1 seconds per pegRNA** on one CPU core.

Measured on an 8-core CPU laptop (no GPU):

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

   * - pegRNAs
     - Total
     - Per pegRNA (apparent)
     - Notes
   * - 10
     - 15 s
     - 1.50 s
     - almost entirely warm-up
   * - 50
     - 18 s
     - 0.37 s
     - warm-up still dominates
   * - 200
     - 33 s
     - 0.17 s
     - warm-up amortising
   * - 720
     - 82 s
     - 0.11 s
     - approaching the marginal rate
   * - 720
     - 53 s
     - 0.07 s
     - ``optiprime_jobs=3`` (2 workers used)

The per-pegRNA figure falls as the set grows because the fixed cost is spread over more rows. Fitting the larger
points gives **~14 s warm-up plus ~0.1 s per pegRNA**, so a rough estimate is:

.. code-block:: text

   total seconds  ~=  14  +  0.1 x (number of pegRNAs scored)

A 10,000-pegRNA library is therefore on the order of 15--20 minutes on one core.

.. note::
   These are order-of-magnitude figures from one machine, not a guarantee. On a GPU, scoring is substantially
   faster.

Two parameters control the total:

``optiprime_prefilter`` (default ``4``)
    Score only the top ``pegRNAs_per_mut x optiprime_prefilter`` pegRNAs per mutation, ranked by ``PEGG2_Score``,
    rather than every candidate. With the defaults this scores 4x what you keep -- generous enough that the
    shortlist rarely excludes a pegRNA OptiPrime would have ranked top, while cutting the work by a large factor.
    Set to ``None`` to score everything.

    **This is the first thing to lower if a run is taking too long**, and it usually has the largest effect.

``optiprime_jobs`` (default: CPU count minus one)
    How many worker processes to shard across.

.. important::
   **Sharding only engages above 500 pegRNAs per worker.** Each worker pays the full ~14 s warm-up, so PEGG refuses
   to create shards smaller than 500 rows: below 1,000 pegRNAs you get one worker no matter what ``optiprime_jobs``
   says, and raising it changes nothing.

   Above that threshold it does help -- 720 pegRNAs took 82 s on one worker and 53 s on two (1.5x). The speedup is
   sublinear because each extra worker repeats the warm-up, so it is worth setting only for large libraries.

.. note::
   Sharding is deterministic. Scores from a sharded run agree with single-process scoring to ~5e-7, which is
   float32 rounding. Workers are separate processes rather than threads because ViennaRNA's energy parameters are
   process-global and mutated during featurization -- sharing an interpreter would make scores depend on
   scheduling.


Silent bystanders and the per-mutation cap
*******************************************

When ``silent_bystander=True`` and ``optiprime=True`` are used together, PEGG scores the sequence each pegRNA
actually installs -- the intended edit plus its bystanders -- not just the intended edit.

``cap_total_per_mut`` (default ``True``) splits the per-mutation budget between plain and bystander designs rather
than taking the overall top N.

.. note::
   This matters more than it looks. Plain and bystander designs do not share a score scale: a bystander can disrupt
   the PAM and shift RTT GC content, both of which the model weights positively. Capping by overall rank lets
   bystander designs sweep the top, so a variant whose only viable designs are plain can lose every row. Splitting
   the budget prevents that.


Complete example
*****************

A runnable version of this, with outputs, is in
`optiprime_example.ipynb <https://github.com/kexindon/PEGG3.0/blob/main/optiprime_example.ipynb>`_.

.. code-block:: python

   import os
   import pandas as pd

   # 1. point PEGG at your OptiPrime installation
   os.environ['OPTIPRIME_SRC']    = os.path.expanduser('~/optiprime/src')
   os.environ['OPTIPRIME_PYTHON'] = os.path.expanduser('~/optiprime/env/bin/python')

   from pegg import prime, optiprime

   # 2. check it resolves before designing anything
   assert optiprime.optiprime_available(), 'OptiPrime not found -- check the two paths above'

   # 3. design and score
   mutations = pd.DataFrame({'SEQ': [
       'ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAG'
       'CAGCAGCAGCAGCAGCAGCAACAGCCACCTCCACCGCCGCCGCCGCCGCCGCCTCCTCCTCAGCTTCCTCAG'
       'CCGCCGCCG(C/T)AGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACC'
       'CGGCCCGGCTGTGGCTGAGGAGCCGCTGCACCGACCAAAAAAAGAGCTTTCTGCTACCAAGAAAGACCGTGT'
       'GAATCATTGTCTGACAATATGTGAAAACATAGTGGCACAGTCTGTCAGAAATTCTCCA'
   ]})

   pegRNAs = prime.run(
       mutations, 'PrimeDesign', None,
       pegRNAs_per_mut=10,
       optiprime=True,
       optiprime_group='Liu_HeLa',     # MMR-proficient; see above
       rankby='OptiPrime_Score',
   )

   # 4. every pegRNA is either scored or carries a reason why not
   print(pegRNAs['OptiPrime_Score'].describe())
   print(pegRNAs['OptiPrime_error'].value_counts(dropna=False))

   print(pegRNAs[['Protospacer_30', 'RTT', 'PBS',
                  'PEGG2_Score', 'OptiPrime_Score']].head(10))


Parameter reference
********************

All ``prime.run()`` parameters relating to OptiPrime:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Meaning
   * - ``optiprime``
     - ``False``
     - Score pegRNAs with OptiPrime. Implied by ``rankby='OptiPrime_Score'`` or ``optiprime_cutoff``.
   * - ``optiprime_cutoff``
     - ``None``
     - Drop pegRNAs scoring below this value.
   * - ``optiprime_group``
     - ``'Liu_HeLa'``
     - Cell line parameters. See :ref:`optiprime-cell-line`.
   * - ``optiprime_src``
     - ``None``
     - Path to the OptiPrime clone. Falls back to ``OPTIPRIME_SRC``.
   * - ``optiprime_python``
     - ``None``
     - Python 3.11 interpreter with OptiPrime's requirements. Falls back to ``OPTIPRIME_PYTHON``.
   * - ``optiprime_prefilter``
     - ``4``
     - Score the top ``pegRNAs_per_mut x`` this per mutation. ``None`` scores everything.
   * - ``optiprime_jobs``
     - ``None``
     - Worker processes. Defaults to CPU count minus one.
   * - ``cap_total_per_mut``
     - ``True``
     - Split the per-mutation budget between plain and bystander designs.

Full function documentation is on the :doc:`PEGG` page.


Citation
*********

If you use OptiPrime scoring, please cite both PEGG and OptiPrime:

Gould, S.I., Wuest, A.N., Dong, K. et al. High-throughput evaluation of genetic variants with prime editing sensor
libraries. *Nat Biotechnol* (2024). https://doi.org/10.1038/s41587-024-02172-9

Hsu, A. et al. Massively parallel measurement of prime editing outcomes enables predictive modelling of pegRNA
efficiency. *Nat Biotechnol* (2026). https://doi.org/10.1038/s41587-026-03261-7
