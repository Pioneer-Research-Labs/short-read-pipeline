#!/usr/bin/env python

# --- process barcode results

import pandas as pd
import os
import click
from fitness_utils import *



@click.command()
@click.argument('path')
@click.option('--outfile', default='barcodes_counts.csv', help='Output file name.')
@click.option('--result_dir', default='.', help='Directory containing results.')
@click.option('--add_frequency', is_flag=True, default=True, help='Add frequency to the output.')
@click.option('--psi', default=0, help='Psi value for frequency calculation.')
@click.option('--cutoff', default=5, help='Cutoff value for filtering.')
@click.option('--overwrite', is_flag=True, default=False, help='Overwrite the output file if it exists.')
@click.help_option('--help', '-h', help='Show this message and exit.')
def main(result_dir, path, outfile, add_frequency, psi, cutoff, overwrite):
    """
    This script processes barcode files and generates a summary of unique barcode counts.

    PATH: The name of the per-sample barcode counts file genenerated by the pipeline.  
    This is 'barcode_counts.tsv' for uncorrected barcodes and 'barcodes_corrected.tsv' for the corrected ones.
    """
    samps = list_samples(result_dir)
    uniq_counts = process_barcodes_files(samps, path, outfile, add_frequency, psi, cutoff, overwrite)
    uniq_path = os.path.splitext(outfile)[0] + '_uniq_counts.csv'
    uniq_counts.to_csv(uniq_path, index=False)

if __name__ == '__main__':
    main()


