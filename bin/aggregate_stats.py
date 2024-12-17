#!/usr/bin/env python

# Load and analyze results

import pandas as pd
import os
from fitness_utils import *


result_dir = '.'
use_correct = True


# Get sample names
samps = list_samples(result_dir)

if len(samps) == 0:
	raise ValueError("No samples to analyze")


print("Loading stats...")
read_stats = load_read_stats(samps, 'read_stats.tsv')
barcode_stats = load_read_stats(samps, 'barcode_stats.tsv')
barcode_filt_stats = load_read_stats(samps, 'barcodes_filtered_stats.tsv')



# rename stats columns
bc_stats = barcode_stats[['sample', 'num_seqs']].rename(columns={'num_seqs': 'barcodes_extracted'}) 
bc_filt = barcode_filt_stats[['sample', 'num_seqs']].rename(columns={'num_seqs': 'barcodes_filtered'})


# format stats dataframe
stats = read_stats[['sample', 'num_seqs']] \
    .drop_duplicates() \
    .merge(bc_stats) \
    .merge(bc_filt) \
    .assign(pct_barcodes = lambda x: 100 * (x.barcodes_filtered / x.num_seqs)) \
    .rename(columns = 
			{'sample': 'Sample', 
             'num_seqs': 'Total reads', 
             'pct_barcodes': 'Perecent with barcodes',
             'barcodes_extracted': 'Barcodes extracted',
             'barcodes_filtered': 'Barcodes passed size filters'
            })

if use_correct:
    print("Loading corrected stats...")
    correct_stats = load_correct_stats(samps, 'correct_stats.csv')
    correct_stats = correct_stats[['sample', 'nUnique', 'nTrueBCs', 'nCorrected']] \
        .rename(columns = {'sample': 'Sample', 'nUnique': 'Unique Barcodes', 'nTrueBCs': 'Num True Barcodes', 'nCorrected': 'Barcodes Corrected'})
    stats = stats.merge(correct_stats)

# Write it all out
stats.to_csv('stats.csv', index = False)