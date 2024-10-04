#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys


path = sys.argv[1] 
out_type = sys.argv[2]

if out_type == "bartender":
    mid = sys.argv[3]

seqs = [x for x in SeqIO.parse(path, 'fasta')]


if out_type == "cutadapt":
    up = [str(x.seq) for x in seqs if x.name == "BARCODEUP"]
    dn = [str(x.seq) for x in seqs if x.name == "BARCODEDN"]
    bc = f'{up[0]}...{dn[0]}'
elif out_type == "bartender":
    up = [str(x.seq)[-5:] for x in seqs if x.name == "BARCODEUP"]
    dn = [str(x.seq)[:5] for x in seqs if x.name == "BARCODEDN"]
    bc = f'{up[0]}[20]{mid}[20]{dn[0]}'
else:
    raise ValueError("Please provide a type of either 'cutadapt' or 'bartender'")

#sys.stdout.write(bc)

SeqIO.write(SeqRecord(Seq(bc), id = "barcode", description = ""), f'{out_type}_bc.fasta', "fasta-2line")