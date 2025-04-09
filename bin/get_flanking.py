#!/usr/bin/env python

from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import re


def extract_flanks(path, feature_regex=r'BARCODE[0-9]{0,2}(UP|DN)'):
    """
    Parse the SnapGene file and extract the flanking sequences based on the provided regex
    :param path: Path to the SnapGene file
    :param feature_regex: Regex to identify the feature annotations in the SnapGene file
    :return: Written Genbank file with the flanking sequences
    """
    con = snapgene_file_to_dict(path)
    part_labels = [feature_regex]
    part_list = {}
    for x in con['features']:
        if any(re.match(reg, x['name']) for reg in part_labels):
            part_list[x['name']] = (x['start'], x['end'])

    flanking_seqs = []
    for x,y in part_list.items():
        seq = con['seq'][y[0]:y[1]]
        # We need location information later on to correctly use annotations in cutadapt
        f = FeatureLocation(y[0], y[1])
        feature = SeqFeature(f, type="misc_feature")
        record = SeqRecord(Seq(seq), id=x, description="", features=[feature], annotations={"molecule_type": "DNA"})
        flanking_seqs.append(record)

    return flanking_seqs



construct_file = sys.argv[1]
flanks = extract_flanks(construct_file)
SeqIO.write(flanks, "flanking.gb", "genbank")