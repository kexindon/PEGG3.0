.. PEGG documentation master file, created by
   sphinx-quickstart on Wed Sep 28 12:28:06 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|PEGG| PEGG: Prime Editing Guide Generator 
============================================

.. |PEGG| image:: PEGG_3.png
   :width: 200px
   :height: 200px

`This is PEGG 3.0. For the original PEGG 2.0 documentation, click here <https://pegg.readthedocs.io/en/latest/>`_
******************************************************************************************************************

`Click here to read the Nature Biotechnology article <https://www.nature.com/articles/s41587-024-02172-9>`_ 
******************************************************************************************************************

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

PEGG was updated to version 2.0 with new features including (1) increased input mutation format flexibility,
(2) Dynamically computed Azimuth on-target scores, (3) a new base editing module, (4) improved library design functionality, as well as some minor bug fixes with INS/DEL design.

Version 3.0 adds silent bystander mutations. pegRNAs can carry synonymous edits alongside the intended one, and both
designs are returned in the same table -- the ordinary pegRNA and up to ``silent_per_mut`` bystander-carrying variants
of it -- distinguished by the ``has_silent_bystander`` column, so either set can be filtered out afterwards. Sensors,
oligos, and polyT / restriction-site filtration all reflect the bystander mutations, and all variant types are
supported (SNP, ONP, INS, DEL, INDEL). The feature is off by default: with ``silent_bystander=False`` the output is
identical to version 2.0.

The synonymous mutation logic is informed by the silent bystander design approach in
`PRIDICT2.0 <https://github.com/uzh-dqbm-cmi/PRIDICT2>`_ (Mathis et al., Nature Protocols 2025).

Installation
**************
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


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   jupyter_tutorial
   PEGG
   about
   other_design_tools

PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S.I., Wuest, A.N., Dong, K. et al. High-throughput evaluation of genetic variants with prime editing sensor libraries. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02172-9

Function Index
***************
* :ref:`genindex`
