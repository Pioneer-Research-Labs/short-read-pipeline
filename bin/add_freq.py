#!/usr/bin/env python

import pandas as pd
import click
from pioneertools.fitness_utils import add_freq

@click.command()
@click.argument('path')
@click.argument('sample')
@click.option('--cutoff', type = int, required = True)
def main(path, sample, cutoff):

    x = pd.read_table(path, names=['barcode', 'count'])
    x.insert(0, 'sample', sample)
        
    x = add_freq(x, psi = 0)

    x.query('n > @cutoff', inplace = True)

    uniq = x.value_counts('sample').reset_index().rename(columns={'count': 'num_uniq_bc'})
    
    x.to_csv(f'{sample}_{cutoff}_barcodes_freq.csv', index = False)
    uniq.to_csv(f'{sample}_{cutoff}_uniq_barcodes.csv', index = False)
    

if __name__ == '__main__':
    main()