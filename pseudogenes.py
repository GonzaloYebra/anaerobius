#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes Roary's output from a pairwise comparison between reference
and query and identifies putative pseudogenes in the query.
It needs as input the Roary's 'gene_presence_absence.csv' file and the IDs
of reference and query as stated in the Roary file.

Usage:
    python pseudogenes.py queryID refID roaryFile

@author: gyebra
"""

import pandas as pd
import sys

queryID = sys.argv[1]
refID = sys.argv[2]
roaryFile = sys.argv[3]

# Count number of locus tags in each COG, for query and reference
df = pd.read_csv(roaryFile)
df["no._query"] = df[queryID].str.count("\t") + 1
df["no._ref"] = df[refID].str.count("\t") + 1

# Calculate proportional size difference between query and reference
df["Prop group size nuc"] = df["Min group size nuc"] / df["Max group size nuc"]


# Identify split pseudogenes: Filter those COGs with one locus tag in reference 
# but 1+ in query
df_pseudo_split = df[(df["no._query"] > 1) & (df["no._ref"] == 1) \
                     & (~df["Annotation"].str.contains("transposase")) \
                     & (df["Min group size nuc"] != df["Max group size nuc"])]

df_pseudo_split["Pseudogene type"] = "Split"    

# Identify shortened pseudogenes: Filter those COGs where gene length in reference
# and query differ in 80%
df_pseudo_short = df[(df["Prop group size nuc"] < 0.8) \
                     & (~df["Annotation"].str.contains("transposase")) \
                     & (df["no._query"] == 1) & (df["no._ref"] == 1)]

df_pseudo_short["Pseudogene type"] = "Truncated"

# Merge and export file
df_concat = pd.concat([df_pseudo_split, df_pseudo_short], axis = 0)
output_file = queryID + "_pseudogenes.csv"
df_concat.to_csv(output_file, sep=',')