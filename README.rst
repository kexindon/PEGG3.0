|PEGG| PEGG: Prime Editing Guide Generator
======================================================
Version 3.0 (Original version released Sept. 2023; Version 2.0 updated Mar. 2024; Version 3.0 updated Aug. 2026)
**********************************************************
.. |PEGG| image:: docs/PEGG_3.png
   :width: 200px
   :height: 200px

`Full Documentation is available here (pegg.readthedocs.io) <https://pegg.readthedocs.io/en/latest/>`_

`Click here to read the Nature Biotechnology article <https://www.nature.com/articles/s41587-024-02172-9>`_ 

Installation
****************
PEGG is available through the python package index. To install, use pip: 

.. code-block:: python

   pip install pegg

Note
*****
PEGG has been tested with python versions 3.9 and 3.10. Python versions higher than 3.10 are not compatible with the scikit-learn package version needed to compute protospacer on-target scores.
To get it to install, you may need to use a `virtual environment <https://saturncloud.io/blog/how-to-install-python-39-with-conda-a-guide-for-data-scientists/>`_ :

.. code-block:: python

   conda create -n myenv python=3.9

Additionally, some users (particularly on Windows) have reported installation issues stemming from the cyvcf2 package used for translating ClinVar IDs to a format compatiable with PEGG.
A version without this package and its functionality is available for local pip installation in the following `dropbox link (pegg-2.0.92-py3-none-any.whl) <https://www.dropbox.com/sh/5xsdzyiyrjiu9pf/AADiFFA3BQ3vX7swja-i2NBqa?dl=0>`_ .


Usage
*******

PEGG is a python package that designs prime editing guide RNAs (pegRNAs) and base editing guide RNAs (gRNAs) for use in precision genome editing.
Unlike the existing, web-based programs for pegRNA design, PEGG is suitable for designing thousands of pegRNAs at once, giving users the ability to design entire libraries of pegRNAs
and gRNAs. Uniquely, PEGG can design paired pegRNA or gRNA-sensor cassettes that include a synthetic version of the target locus, allowing for 
the calibration of guide editing activity in pooled screens (see above bioRxiv preprint for more information).

PEGG's main functions are:

(1) Generating pegRNAs or gRNAs based on a list of input mutations. The input format is extremely flexible, allowing for users to input a list of genome coordinates or sequences with desired edits.

(2) Ranking and filtering these pegRNAs based on their properties, including On-Target (Azimuth) Scores.

(3) Automated oligo generation (with the option to include a synthetic "sensor" region).

(4) Automated pegRNA/gRNA library design with included safe-targeting mutations, non-targeting guides, and silent substitution controls.

(5) Visualization Tools for pegRNA and gRNA design.

Updates from Version 1.0 to Version 2.0:
PEGG has been updated to version 2.0, with new features including (1) increased input mutation format flexibility,
(2) Dynamically computed Azimuth on-target scores, (3) a new base editing module, (4) improved library design functionality, as well as some minor bug fixes with INS/DEL design.

Updates from Version 2.0 to Version 3.0:
PEGG has been updated to version 3.0, with new feature of adding synonymous bystander mutations nearby the editing site to increase prime editing efficiency.

The synonymous mutation logic is informed by the
silent bystander design approach in PRIDICT2.0 (https://github.com/uzh-dqbm-cmi/PRIDICT2).

Mathis, N., Marquart, K.F., Allam, A., Krauthammer, M. & Schwank, G. Systematic pegRNA design with PRIDICT2.0 and ePRIDICT for efficient prime editing. Nature Protocols (2025). https://doi.org/10.1038/s41596-025-01244-7

Usage
------

Switch the feature on with ``silent_bystander=True``. The output contains **both** the ordinary pegRNAs and the
bystander-carrying ones, distinguished by the ``has_silent_bystander`` column, so either set can be filtered out afterwards:

.. code-block:: python

   peg_df = prime.run(mutations, 'cBioPortal', chrom_dict=chrom_dict,
                      silent_bystander=True,
                      silent_per_mut=2,              # bystander designs per pegRNA
                      transcript_strand='-',         # strand the transcript is on
                      start_end_cds=start_end_cds,   # CDS blocks of the transcript
                      seed=0)                        # for a reproducible library

Sensors, oligos, and polyT/restriction-site filtration all reflect the bystander mutations. Additional output columns:
``n_bystander_muts``, ``bystander_positions``, ``bystander_dist_to_edit``, ``PAM_disrupted_edit``,
``PAM_disrupted_by_bystander``, and ``pegRNA_rank_within_group``.

A silent bystander can knock out the PAM as well as the intended edit can, so PAM disruption is reported three ways:

* ``PAM_disrupted`` -- disrupted by **either** cause. This is what the scoring functions read, since it is what matters
  biologically. Note this is a change in meaning from PEGG 2.0, where the column referred to the intended edit alone.
* ``PAM_disrupted_edit`` -- disrupted by the intended edit alone (the PEGG 2.0 meaning of ``PAM_disrupted``).
* ``PAM_disrupted_by_bystander`` -- disrupted by a silent bystander alone.

Other parameters: ``bystander_window_nt`` (how far from the edit silent mutations may be placed, default 5),
``max_bystander_muts`` (silent changes per pegRNA, default 2), and ``splice_buffer`` (minimum distance from an exon
boundary, default 3).

Building the CDS annotation
----------------------------

``start_end_cds`` can be extracted from a GTF/GFF3 annotation rather than typed by hand, and checked before use:

.. code-block:: python

   from pegg import bystander

   # one entry per gene: {'strand': ..., 'cds': [[start, end], ...]}
   gene_cds = bystander.cds_from_gtf('gencode.v19.annotation.gtf.gz',
                                     genes=list(mutations['Hugo_Symbol'].unique()))

   # check it before designing: frame validity, overlapping blocks, genes with no
   # annotation, and variants that fall outside their gene's CDS
   print(bystander.cds_annotation_report(gene_cds, mutations=mutations))

For a cBioPortal-format table the reading frame can be resolved straight from the standard ``Hugo_Symbol`` column,
with no extra annotation columns required:

.. code-block:: python

   import gffutils
   db = gffutils.FeatureDB('gencode.v19.db')      # or h2m.anno_loader(path)

   df, cds_lookup = bystander.cds_for_variants(mutations, db)
   # 4/6 variants have a usable reading frame (2 of 3 genes resolved)
   #   no canonical transcript known for: UNKNOWN
   #   CDS not a multiple of 3: GENE2
   #   those variants still get ordinary pegRNAs, just no bystanders

   parts = []
   for gene, sub in df.groupby('Hugo_Symbol'):
       ann = cds_lookup.get(gene)
       if ann is not None and ann['valid']:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal',
                           silent_bystander=True,
                           transcript_strand=ann['strand'],
                           start_end_cds=ann['cds'], seed=0, **run_params)
       else:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal', **run_params)
       parts.append(out)

   peg_df = pd.concat(parts, ignore_index=True)

``cds_for_variants()`` adds ``transcript_id``, ``transcript_strand``, ``cds_valid`` and ``cds_n_codons`` columns.

Each gene is designed against its **canonical transcript**, looked up from the curated table H2M uses (bundled with
pegg; see ``bystander.canonical_transcripts()``). A gene has many transcripts and the reading frame belongs to exactly
one of them, so this is not left to a heuristic: variant annotation such as ``HGVSp_Short`` is generally computed
against the canonical transcript, and designing against a different one can shift the reading frame and make those
protein-level labels disagree with the design.

Pass ``species='mouse'`` or ``genome_version=38`` to match your data, ``transcript_ids={'TP53': 'ENST00000269305'}``
to override particular genes, or ``tx_column='tx_id_h'`` to use ids already in the table (e.g. from H2M's
``get_tx_batch()``).

If transcript ids are already attached to the variant table -- for instance by H2M's ``get_tx_batch()`` -- the reading
frame can be resolved per row instead, which keeps each variant tied to the transcript it was annotated against:

.. code-block:: python

   import gffutils
   db = gffutils.FeatureDB('gencode.v19.db')      # or h2m.anno_loader(path)

   # df already has a tx_id_h column, e.g. from h2m.get_tx_batch(df, 'h', ver=37)
   df, cds_lookup = bystander.add_cds_to_variants(df, db, tx_column='tx_id_h')

   parts = []
   for tx, sub in df.groupby('tx_id_h'):
       ann = cds_lookup.get(tx)
       if ann is not None and ann['valid']:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal',
                           silent_bystander=True,
                           transcript_strand=ann['strand'],
                           start_end_cds=ann['cds'], seed=0, **run_params)
       else:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal', **run_params)
       parts.append(out)

   peg_df = pd.concat(parts, ignore_index=True)

``add_cds_to_variants()`` adds ``transcript_strand``, ``cds_valid`` and ``cds_n_codons`` columns and reports how many
variants have a usable reading frame, naming the transcripts that do not.

``cds_from_gtf()`` returns every transcript, keyed as ``"gene|transcript_id"``, unless ``transcript_ids={'TP53':
'ENST00000269305'}`` names the ones to keep. Ensembl, GENCODE and NCBI RefSeq files are all accepted, gzipped or not.

Designing a library across several genes then means one ``run()`` call per gene, since the reading frame belongs to a
single transcript:

.. code-block:: python

   parts = []
   for gene, sub in mutations.groupby('Hugo_Symbol'):
       ann = gene_cds.get(gene)
       if ann is not None and ann['valid']:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal',
                           silent_bystander=True, silent_per_mut=2,
                           transcript_strand=ann['strand'],
                           start_end_cds=ann['cds'], seed=0, **run_params)
       else:
           out = prime.run(sub.reset_index(drop=True), 'cBioPortal', **run_params)
       parts.append(out)

   peg_df = pd.concat(parts, ignore_index=True)

Genes with no annotation simply get no bystanders, so a mixed library is fine.

Reading frame requirements
---------------------------

Silent bystanders need reading frame information, which differs by input format:

* ``cBioPortal`` -- needs ``start_end_cds`` and ``transcript_strand``. The reading frame is looked up per genomic position,
  and bystanders are kept inside the CDS and at least ``splice_buffer`` nt away from splice sites.
* ``WT_ALT`` / ``PrimeDesign`` -- need ``ORF_start`` (0, 1 or 2), and **the input sequence must be in frame**, i.e. entirely
  coding sequence. CDS membership and splice sites cannot be checked for these formats.

All variant types are supported (SNP, ONP, INS, DEL, INDEL). Where a frameshifting indel would destroy the downstream
reading frame, bystanders are placed only upstream of the edit. If the reading frame cannot be established for a variant --
for example because it falls outside the supplied CDS -- that variant simply gets no bystanders, rather than a guess.

Ordinary pegRNA design is unaffected by any of this: with ``silent_bystander=False`` (the default) no reading frame
information is needed and the output is unchanged.

Version 3.0 change summary
****************************

Everything below is additive. With ``silent_bystander=False`` (the default) PEGG 3.0 produces byte-identical output to
version 2.0 for all three input formats and for the base editing module; this was verified column by column against the
previous release on real GRCh37 data.

**New capability**

pegRNAs can carry synonymous ("silent") bystander mutations alongside the intended edit. Both designs are returned in
the same table -- the ordinary pegRNA and up to ``silent_per_mut`` bystander-carrying variants of it -- distinguished by
the ``has_silent_bystander`` column, so either set can be filtered out afterwards. Sensors, oligos, and polyT /
restriction-site filtration all reflect the bystander mutations. All variant types are supported (SNP, ONP, INS, DEL,
INDEL).

**New parameters on** ``prime.run()`` (all keyword, all appended, all inert when the feature is off)

.. code-block:: text

   silent_bystander=False     switch the feature on
   silent_per_mut=2           bystander designs kept per pegRNA
   transcript_strand=None     '+' or '-'; strand the transcript is on
   start_end_cds=None         CDS blocks of the transcript, 1-based inclusive
   ORF_start=None             frame offset, for WT_ALT / PrimeDesign input
   bystander_window_nt=5      how far from the edit silent mutations may be placed
   max_bystander_muts=2       silent changes per pegRNA
   splice_buffer=3            minimum distance from an exon boundary
   seed=None                  makes the random choice of designs reproducible

**New output columns** (7 added, none removed, no dtype changes)

.. code-block:: text

   has_silent_bystander        bool      the column to filter on
   n_bystander_muts            int       silent changes carried by this pegRNA
   bystander_positions         str       ';'-separated offsets within the RTT
   bystander_dist_to_edit      float     nt from the nearest silent change to the edit
   PAM_disrupted_edit          bool      PAM disrupted by the intended edit alone
   PAM_disrupted_by_bystander  bool      PAM disrupted by a silent bystander alone
   pegRNA_rank_within_group    int       rank within each design type

**One changed meaning**

``PAM_disrupted`` now means "disrupted by the intended edit **or** by a silent bystander", since that is what matters
biologically and what the scoring functions should see. Its previous meaning is preserved exactly in
``PAM_disrupted_edit``. With the feature off the two columns are identical, so nothing changes for existing pipelines.

Similarly, ``pegRNA_rank`` still ranks all pegRNAs of a mutation together, which now includes bystander designs;
``pegRNA_rank_within_group`` ranks within each design type and reproduces the version 2.0 ``pegRNA_rank`` exactly. When
``pegRNAs_per_mut`` is set, the limit applies to each design type separately, so a mutation always keeps both.

**Bug fix:** ``sensor_viz()`` with an out-of-window 3' extension

``prime.sensor_viz()`` could fail with ``ValueError: setting an array element with a sequence. The requested array has
an inhomogeneous shape`` for some pegRNAs. This affected version 2.0 as well and is fixed here.

The cause: the plot draws the sensor, the protospacer and the 3' extension as rows of one array, padding each with
``'-'`` so they line up. Where the extension reached past the edge of the sensor window the padding offset went
negative, and ``['-'] * <negative>`` is an empty list rather than an error -- so that row silently came out longer than
the others and matplotlib refused to render the ragged array. It showed up on short sensors and in repetitive sequence,
where the protospacer match can land at an unexpected offset.

The row is now clipped to the part that is actually visible in the sensor window, and a protospacer that cannot be
located in the sensor at all raises a clear message naming the row instead of a matplotlib error. Plots that rendered
before are unchanged -- verified cell-for-cell against version 2.0 output.

**New module:** ``pegg.bystander``

Silent bystander generation plus reading frame utilities: ``cds_for_variants()`` (the batch entry point for
cBioPortal input), ``canonical_transcripts()``, ``cds_from_annotation_db()``, ``cds_from_gtf()``,
``cds_annotation_report()``, ``cds_frame_map()``.

**Backwards compatibility**

No function signature lost or reordered a parameter; no output column was removed or renamed; ``base.run_base()`` is
untouched. Existing code that reads sequence or score columns needs no change. Code that compares the full column list,
or that filters on ``PAM_disrupted`` / ``pegRNA_rank`` **with the feature switched on**, should move to
``PAM_disrupted_edit`` / ``pegRNA_rank_within_group``. Filtering ``df[~df['has_silent_bystander']]`` recovers the
version 2.0 design set at any point.


Updating the online documentation
***********************************

The documentation at https://pegg.readthedocs.io is built with Sphinx from the ``docs/`` directory of this repository.
The ``.rst`` sources are in the repo, so the edits below can be made here and pushed; whoever holds the Read the Docs
account only needs to trigger (or allow) a rebuild.

Four places need updating for version 3.0:

**1.** ``docs/index.rst`` -- *main function description*

The numbered list under "PEGG's main functions are:" (around line 21) ends at item (5). Add a sixth item:

.. code-block:: rst

   (6) Optional silent bystander mutations: synonymous edits placed alongside the intended edit to evade mismatch
       repair and, where they fall in the PAM, reduce re-nicking. Both the ordinary and bystander-carrying designs are
       returned, so they can be compared directly.

The paragraph immediately after that list still says "PEGG has recently been updated to version 2.0". Replace it with a
version 3.0 paragraph -- the "Version 3.0 change summary" section above can be reused verbatim.

**2.** ``docs/quickstart.rst`` -- *input formatting* (section starts line 24)

This section documents the three input formats. Add a note that silent bystanders need reading frame information, and
that the requirement differs by format:

.. code-block:: rst

   For silent bystander design, the reading frame must be known. With ``cBioPortal`` input it is looked up per genomic
   position from a CDS annotation, and CDS membership and splice sites are checked. With ``WT_ALT`` and ``PrimeDesign``
   input there are no genomic coordinates, so the frame is declared with ``ORF_start`` and **the input sequence must be
   in frame** -- entirely coding sequence. Ordinary pegRNA design is unaffected and needs no frame information.

**3.** ``docs/quickstart.rst`` -- *generating pegRNAs* (line 128) *and design options* (lines 152, 211)

Add ``silent_bystander`` and its companion parameters to the design options tables, and add a worked example after the
existing ``prime.run()`` example. The "Usage" and "Building the CDS annotation" sections of this README can be copied
across as-is.

**4.** ``docs/PEGG.rst`` -- *complete documentation*

This page autogenerates API docs with ``automodule``. Add a fourth block alongside ``pegg.prime``, ``pegg.base`` and
``pegg.library``:

.. code-block:: rst

   pegg.bystander
   ---------------

   .. automodule:: pegg.bystander
      :members:

Every function in ``pegg/bystander.py`` already carries a numpydoc-style docstring in the same format as the rest of
the package, so this block is all that is needed.

**Rebuilding locally to check**

.. code-block:: bash

   pip install -r docs/doc-requirements.txt
   sphinx-build -b html docs docs/_build/html

Read the Docs rebuilds automatically on push if the webhook is active; otherwise trigger a build from the project
dashboard.


PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S.I., Wuest, A.N., Dong, K. et al. High-throughput evaluation of genetic variants with prime editing sensor libraries. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02172-9