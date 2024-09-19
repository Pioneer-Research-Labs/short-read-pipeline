#!/usr/bin/env python

from Bio import SeqIO
import sys

path = sys.argv[1] 
mid = sys.argv[2]
seqs = [x for x in SeqIO.parse(path, 'fasta')]
up = [str(x.seq)[-5:] for x in seqs if x.name == "BARCODEUP"]
dn = [str(x.seq)[:5] for x in seqs if x.name == "BARCODEDN"]
bc = f'{up[0]}[20]{mid}[20]{dn[0]}'
sys.stdout.write(bc)