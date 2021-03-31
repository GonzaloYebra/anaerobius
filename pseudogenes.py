#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes Roary's output from a pairwise comparison between reference
and query and identifies putative pseudogenes in the query.
It requires as input the IDs of reference and query and the path to
a folder containing:
    - the Roary's 'gene_presence_absence.csv' file
    - the annotated genomes (genbank format) of reference and query with 
      names consistent with those in the Roary file.

Usage:
    python pseudogenes.py queryID refID path-to-files

Dependencies: python3, pandas, biopython

@author: gyebra
"""

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from Bio import SeqIO
import sys

queryID = sys.argv[1]
refID = sys.argv[2]
path = sys.argv[3]

# Count number of locus tags in each COG, for query and reference
df = pd.read_csv(path + "/gene_presence_absence.csv")
df["no._query"] = df[queryID].str.count("\t") + 1
df["no._ref"] = df[refID].str.count("\t") + 1

# Calculate proportional size difference between query and reference
df["Prop group size nuc"] = df["Min group size nuc"] / df["Max group size nuc"]
df["Prop group size nuc"] = df["Prop group size nuc"].round(decimals=2)

# Disregard genes where max length is <100bp (<300aa)
df = df[df["Max group size nuc"] > 299]

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
df_pseudo_short["comment"] = ""

# Extract info from genbank files to see which is longer, the ref or the query
for i in range(0, len(df_pseudo_short[queryID])):
    tag1 = df_pseudo_short[refID].iloc[i]
    for record in SeqIO.parse(path + "/" + refID + ".gbk", "genbank"):
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                locus_tag = f.qualifiers["locus_tag"][0]
                if locus_tag  == tag1:
                    sizeRef = (f.location.end - f.location.start) / 3
    tag2 = df_pseudo_short[queryID].iloc[i]
    for record in SeqIO.parse(path + "/" + queryID + ".gbk", "genbank"):
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                locus_tag = f.qualifiers["locus_tag"][0]
                if locus_tag  == tag2:
                    sizeQuery = (f.location.end - f.location.start) / 3
    
    if int(sizeRef) > int(sizeQuery):
        df_pseudo_short["comment"].iloc[i] = "Ref is longer"

df_pseudo_short = df_pseudo_short[df_pseudo_short["comment"] == "Ref is longer"]

# Merge and export file
df_concat = pd.concat([df_pseudo_split, df_pseudo_short], axis = 0)

df_concat = df_concat[[
    "Gene",
    "Annotation",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    queryID,
    refID,
    "no._query",
    "no._ref",
    "Prop group size nuc",
    "Pseudogene type"
    ]]

output_file = queryID + "_pseudogenes.csv"
df_concat.to_csv(output_file, sep=',')
