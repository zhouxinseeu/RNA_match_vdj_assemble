import pysam
import pandas as pd
from collections import defaultdict


def get_barcode_reads(fqfile, sample):
    barcode_reads_dict = defaultdict(list)
    barcode_reads_count_dict = {}
    with pysam.FastxFile('HuG_TCR1225_2T_2.fq.gz') as fq:
        for entry in fq:
            attr = entry.name.split('_')
            barcode = attr[0]
            barcode_reads_dict[barcode].append(entry)
        for barcode in barcode_reads_dict.keys():
            barcode_reads_count_dict[barcode] = len(barcode_reads_dict[barcode])
    return barcode_reads_dict
    return barcode_reads_count_dict


def get_barcode_use(bclist, topn = None):
    with open(bclist, 'r') as bclist:
        if topn and topn != 'None':

        else:
            topn = int(topn)
        



