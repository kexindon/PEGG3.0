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

PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S.I., Wuest, A.N., Dong, K. et al. High-throughput evaluation of genetic variants with prime editing sensor libraries. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02172-9