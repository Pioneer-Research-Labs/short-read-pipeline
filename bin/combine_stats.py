#!/usr/bin/env python

import click
import pandas as pd

def load_read_stats(path, sample):
    x = pd.read_table(path)
    x['sample'] = sample
    return(x)

def load_correct_stats(path, sample):
    x = pd.read_csv(path)
    x.insert(0, 'sample', sample)
    return(x)

@click.command()
@click.argument('bc_stats_file')
@click.argument('bc_filt_stats_file')
@click.argument('read_stats_file')
@click.argument('correct_stats_file', required = False, default = None)
@click.option('--sample_name', required = True)
def main(bc_stats_file, bc_filt_stats_file, read_stats_file, correct_stats_file, sample_name):

    read_stats = pd.read_table(read_stats_file).assign()

    read_stats = load_read_stats(read_stats_file, sample_name)
    barcode_stats = load_read_stats(bc_stats_file, sample_name)
    barcode_filt_stats = load_read_stats(bc_filt_stats_file, sample_name)


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


    if correct_stats_file is not None:
        print("Loading corrected stats...")
        correct_stats = load_correct_stats(correct_stats_file, sample_name)
        correct_stats = correct_stats[['sample', 'nUnique', 'nTrueBCs', 'nCorrected']] \
            .rename(columns = {'sample': 'Sample', 'nUnique': 'Unique Barcodes', 'nTrueBCs': 'Num True Barcodes', 'nCorrected': 'Barcodes Corrected'})
        stats = stats.merge(correct_stats)

    # Write it all out
    stats.to_csv(f'{sample_name}_stats.csv', index = False)

if __name__ == '__main__':
    main()