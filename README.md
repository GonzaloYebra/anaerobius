# anaerobius
Scripts used in the genomic analysis of <i>S. aureus</i> subsp. <i>anaerobius</i>.

<i><b>pseudogenes.py</i></b>

This script takes Roary's output from a pairwise comparison between the query genome and one reference genome. It requires as input the IDs of reference and query and the path to
a folder containing:
    - Roary's 'gene_presence_absence.csv' file.
    - annotated genomes (genbank format) of reference and query with names consistent with those in the Roary file. This is to calculate gene lengths from both genomes.

Usage:
    python pseudogenes.py queryID refID path-to-files

Dependencies: python 3, biopython, pandas.
