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

not_excluded = codes[codes['barcode_corrected'].apply(lambda bc: 'xcluded' not in str(bc))]
excluded = codes[codes['barcode_corrected'].apply(lambda bc: 'xcluded' in str(bc))]

stats = {}
stats['nUnique'] = [len(set(codes['barcode']))]
stats['nReads'] = [np.sum(codes['count'])]
stats['nTrueBCs'] = [len(set(not_excluded['barcode_corrected']))]
stats['nCorrected'] = [len(not_excluded[not_excluded['barcode']!=not_excluded['barcode_corrected']])]
stats['readsUsed'] = [np.sum(not_excluded['count'])]
stats['readsTrueBCs'] = [np.sum(not_excluded[not_excluded['barcode']==not_excluded['barcode_corrected']]['count'])]
stats['readsCorrected'] = [np.sum(not_excluded[not_excluded['barcode']!=not_excluded['barcode_corrected']]['count'])]
stats['nExcludedLowCount'] = [len(excluded[excluded['barcode_corrected']=='excluded low count'])]
stats['readsExcludedLowCount'] = [np.sum(excluded[excluded['barcode_corrected']=='excluded low count']['count'])]
stats['nExcludedError'] = [len(excluded[excluded['barcode_corrected']=='excluded error'])]
stats['readsExcludedError'] = [np.sum(excluded[excluded['barcode_corrected']=='excluded error']['count'])]


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
pd.DataFrame.from_dict(stats).to_csv("correct_stats.csv", index = False)


