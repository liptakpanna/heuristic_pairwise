# Introduction

Heuristic pairwise alignment algorithm for DNA sequences in MonetDB.

# Requirements

- MonetDB version 11.44.0. built from source
- Python at least version 3.5
- Sufficient RAM according to the sequence lengths

# Start database server with embedded Python

mserver5 --dbpath=\<path to dbfarm database\> --set embedded_py=true

# Example usage

## Create tables and functions

mclient -u \<user\> -d \<database\> create_tables.sql create_functions.sql create_functions_affine.sql

## Insert data from file
COPY INTO seq_data(seq) FROM '\<absolute path to data file\>'(seq);

## Insert data into lookup table

score_match \> 0, score_mismatch, score_gap, score_gap_open, score_gap_extend \< 0

k: length of k-mers \> 0

samplesize: how many sequence should be in the sample set \> 1

upperlimit: the largest id to choose the samples from \> 1

- For fix gap penalty:

  CALL update_lookup(score_match, score_mismatch, score_gap, k, samplesize, upperlimit);

- For affine gap penalty:

  CALL update_lookup_affine(score_match, score_mismatch, score_gap_open, score_gap_extend, k, samplesize, upperlimit);

## Use lookup table to get the alignment

If there is no record in the lookup table that could be used it will fall back to the Needleman-Wunsch algorithm.

If you want to log the alignments set the logfile in the create_functions.sql use_subalign and create_functions_affine.sql use_subalign_affine functions.

id1, id2: valid sequence ids from the seq_data table

- For fix gap penalty:

  SELECT use_lookup(id1, id2, score_match, score_mismatch, score_gap);

- For affine gap penalty:

  SELECT use_lookup_affine(id1, id2, score_match, score_mismatch, score_gap_open, score_gap_extend);

The above functions return the two aligned string and the score of the alignment.

# Measurements

For running the measurements the first steps are inserting the sequences and creating some additional tables.
In these tables the score and time of the Needleman-Wunsch algorithm for given pairs of sequences will be stored, because this data can be reused, it will not change during the executions.

mclient -u \<user\> -d \<database\> help_tables.sql

Before running, specify a logfile in the 'score_comp.R' and 'score_comp_affine.R' scripts.
