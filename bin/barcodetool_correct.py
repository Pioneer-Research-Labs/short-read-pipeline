#!/usr/bin/env python

from DeletionCorrect import correct_bc_errors
import pandas as pd
import numpy as np
import sys


path = sys.argv[1]
min_centroid = 2
error_rate = 0.1
max_edits = 3
codes = pd.read_table(path, names = ['barcode', 'count']) \
    .sort_values('count', ascending = False)

corrector = correct_bc_errors(np.array(codes), min_counts_for_centroid = min_centroid, 
                              max_edits = max_edits, poisson_error_rate = error_rate)

codes['barcode_corrected'] = codes['barcode'].map(corrector)
removed_error = codes.query('barcode_corrected != "excluded error"')
low_count = removed_error.query('barcode_corrected == "excluded low count"')[['barcode', 'count']]
corrected = removed_error.query('barcode_corrected != "excluded low count"')[['barcode_corrected', 'count']] \
    .rename(columns={'barcode_corrected': 'barcode'})
all = pd.concat([corrected, low_count])

final_barcodes = (all
 .groupby('barcode', as_index=False)
 .sum()
 .rename(columns={'sum':'count'})
 .sort_values('count', ascending=False)
 )


final_barcodes.to_csv("barcodes_corrected.tsv", sep = "\t", index = False, header = False)


