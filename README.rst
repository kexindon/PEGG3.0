|PEGG| PEGG: Prime Editing Guide Generator
======================================================
`Link to PEGG2.0 <https://github.com/samgould2/PEGG2.0>`_
**************************************************************************************************************

Version 3.0 (Original version released Sept. 2023; Version 2.0 updated Mar. 2024; Version 3.0 updated Aug. 2026)
****************************************************************************************************************
.. |PEGG| image:: https://raw.githubusercontent.com/kexindon/PEGG3.0/main/docs/PEGG_3.png
   :width: 200px
   :height: 200px

`Full Documentation is available here (pegg30.readthedocs.io) <https://pegg30.readthedocs.io/en/latest/>`_

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

(6) Optional silent bystander mutations: synonymous edits placed alongside the intended edit to evade mismatch repair and, where they fall in the PAM, reduce re-nicking. Both the ordinary and bystander-carrying designs are returned, so they can be compared directly.

Updates from Version 1.0 to Version 2.0:
PEGG has been updated to version 2.0, with new features including (1) increased input mutation format flexibility,
(2) Dynamically computed Azimuth on-target scores, (3) a new base editing module, (4) improved library design functionality, as well as some minor bug fixes with INS/DEL design.

Updates from Version 2.0 to Version 3.0:
PEGG has been updated to version 3.0, with new feature of adding synonymous bystander mutations nearby the editing site to increase prime editing efficiency.

The synonymous mutation logic is informed by the
silent bystander design approach in PRIDICT2.0 (https://github.com/uzh-dqbm-cmi/PRIDICT2).

Mathis, N., Marquart, K.F., Allam, A., Krauthammer, M. & Schwank, G. Systematic pegRNA design with PRIDICT2.0 and ePRIDICT for efficient prime editing. Nature Protocols (2025). https://doi.org/10.1038/s41596-025-01244-7

Basic usage
------------

Unchanged from version 2.0. Input a table of mutations in one of three formats -- ``cBioPortal``, ``WT_ALT`` or
``PrimeDesign`` -- and ``prime.run()`` returns a dataframe of pegRNA-sensor designs:

.. code-block:: python

   from pegg import prime

   pegRNAs = prime.run(input_df, 'PrimeDesign')

The design parameters are all optional and all have defaults:

* **pegRNA:** ``PAM`` (default ``"NGG"``), ``rankby`` (``'PEGG2_Score'`` or ``'RF_Score'``), ``pegRNAs_per_mut``,
  ``RTT_lengths``, ``PBS_lengths``, ``min_RHA_size``, ``RE_sites``, ``polyT_threshold``, ``proto_size``,
  ``context_size``.
* **sensor:** ``sensor`` (True/False), ``sensor_length``, ``sensor_orientation``, ``before_proto_context``.
* ``chrom_dict`` -- required for ``cBioPortal`` input, from ``prime.genome_loader()``.

See the `full documentation <https://pegg30.readthedocs.io/en/latest/>`_ for what each one does.

Silent bystander mutations
---------------------------

New in version 3.0, and off by default. Switch it on with ``silent_bystander=True``. The output then contains **both**
the ordinary pegRNAs and the bystander-carrying ones, distinguished by the ``has_silent_bystander`` column, so either
set can be filtered out afterwards.

**For cBioPortal input, the mutation table must be prepared with** ``bystander.cds_for_variants()`` **first.** A
bystander is only silent with respect to a reading frame, and that frame comes from the gene's transcript; this step
resolves it and attaches it to each row. Passing an unprepared table raises an error naming the missing column, rather
than quietly designing without bystanders:

.. code-block:: python

   from pegg import bystander

   # attach each gene's canonical transcript, strand and CDS blocks to the table
   mutations, cds = bystander.cds_for_variants(mutations, db)

   peg_df = prime.run(mutations, 'cBioPortal', chrom_dict=chrom_dict,
                      silent_bystander=True,
                      silent_per_mut=2,      # bystander designs per pegRNA
                      seed=0)                # for a reproducible library

The reading frame is read from the table per variant, so one call covers a library of any number of genes -- each
variant is designed against its own transcript. Variants with no usable annotation simply get ordinary pegRNAs.

Further parameters: ``bystander_window_nt`` (how far from the edit silent mutations may be placed, default 5),
``max_bystander_muts`` (silent changes per pegRNA, default 2), and ``splice_buffer`` (minimum distance from an exon
boundary, default 3).

Sensors, oligos, and polyT/restriction-site filtration all reflect the bystander mutations. Additional output columns:
``n_bystander_muts``, ``bystander_positions``, ``bystander_dist_to_edit``, ``PAM_disrupted_edit``,
``PAM_disrupted_by_bystander``, and ``pegRNA_rank_within_group``.

A silent bystander can knock out the PAM as well as the intended edit can, so PAM disruption is reported three ways:
``PAM_disrupted`` (either cause -- this is what the scoring functions read, and a change in meaning from version 2.0),
``PAM_disrupted_edit`` (the intended edit alone, i.e. the version 2.0 meaning), and ``PAM_disrupted_by_bystander``.

**Reading frame.** Because bystanders must be synonymous, this feature -- and only this feature -- needs to know the
reading frame. For ``cBioPortal`` input it is read per row from the columns ``cds_for_variants()`` attaches, so
nothing has to be typed by hand. For ``WT_ALT`` and ``PrimeDesign`` input there are no genomic coordinates, so pass
``ORF_start`` (0, 1 or 2) and make sure the input sequence is entirely coding. Variants whose frame cannot be
established simply get no bystanders.

All variant types are supported (SNP, DNP, TNP, ONP, INS, DEL, INDEL).


Version 3.0 change summary
****************************

Everything is additive: with ``silent_bystander=False`` (the default), PEGG 3.0 produces byte-identical output to
version 2.0.

**New capability.** pegRNAs can carry synonymous ("silent") bystander mutations alongside the intended edit. Both
designs are returned in the same table, distinguished by the ``has_silent_bystander`` column. Sensors, oligos and
filtration all reflect the bystander mutations, and all variant types are supported (SNP, DNP, TNP, ONP, INS,
DEL, INDEL).

**New module** ``pegg.bystander`` -- silent bystander generation plus reading frame utilities, most usefully
``cds_for_variants()``.

**New parameters** on ``prime.run()``, all keyword and all inert when the feature is off: ``silent_bystander``,
``silent_per_mut``, ``ORF_start``, ``bystander_window_nt``, ``max_bystander_muts``, ``splice_buffer``, ``seed``.

**New output columns**, 7 added, none removed or renamed: ``has_silent_bystander``, ``n_bystander_muts``,
``bystander_positions``, ``bystander_dist_to_edit``, ``PAM_disrupted_edit``, ``PAM_disrupted_by_bystander``,
``pegRNA_rank_within_group``.

**One changed meaning.** ``PAM_disrupted`` now means "disrupted by the intended edit **or** by a silent bystander";
its previous meaning is preserved exactly in ``PAM_disrupted_edit``. Likewise ``pegRNA_rank`` now ranks bystander
designs together with ordinary ones, while ``pegRNA_rank_within_group`` reproduces the version 2.0 ranking. With the
feature off both pairs are identical.

**Bug fixes.** Two, both of which affected version 2.0 as well:

``prime.sensor_viz()`` could fail with a matplotlib ``ValueError`` when a pegRNA's 3' extension reached past the edge
of the sensor window. Plots that rendered before are unchanged.

A ``TNP`` variant was not recognised by ``df_formatter()``, which handled ``SNP``, ``DNP``, ``ONP``, ``INDEL``,
``INS`` and ``DEL`` but not ``TNP``. Because the flanking context is assigned inside a loop, an unrecognised type
silently kept the **previous row's** sequence and so described a different locus; the variant was then dropped by the
consistency check with a confusing "Error in mutant #N" message, or -- worse, had the lengths happened to agree --
designed against the wrong site. ``TNP`` is now handled, and any unrecognised ``Variant_Type`` raises a named error
rather than falling through.


PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S.I., Wuest, A.N., Dong, K. et al. High-throughput evaluation of genetic variants with prime editing sensor libraries. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02172-9
