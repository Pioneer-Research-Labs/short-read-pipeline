#!/usr/bin/env python

import BarcodeErrorCorrector as bc
import pandas as pd
import sys

raw_codes = pd.read_table(sys.argv[1], names = ['barcode', 'count'])

bc.error_correct_file_or_df(
    fileOrDf = raw_codes,
    outfile = "barcodes_clustered.csv",
    bc_cols = ["barcode"],
    count_col = "count"
)


