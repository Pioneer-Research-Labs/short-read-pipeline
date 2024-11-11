#!/usr/bin/env python

from DeletionCorrect import correct_bc_errors
import pandas as pd
import numpy as np
import click



@click.command()
@click.option('--path', type=click.Path(exists=True), help='Path to barcode count file.')
@click.option('--min-centroid', default=2, type=int, show_default=True, help='Minimum number of centroids.')
@click.option('--error-rate', default=0.1, type=float, show_default=True, help='Allowed error rate.')
@click.option('--max-edits', default=3, type=int, show_default=True, help='Maximum number of edits.')

def main(path, min_centroid, error_rate, max_edits):
    """
    Run barcode correction
    """
    
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


    final_barcodes.to_csv("barcodes_corrected_low_count.tsv", sep = "\t", index = False, header = False)
    corrected.to_csv("barcodes_corrected.tsv", sep = "\t", index = False, header = False)
    pd.DataFrame.from_dict(stats).to_csv("correct_stats.csv", index = False)


if __name__ == '__main__':
    main()
