import pysam
import pandas as pd
from collections import defaultdict
import os


def barcodes_splitfq(fqfile, fastq_dir, bclist=None, topn=None):
    barcode_reads_dict = defaultdict(list)
    barcode_reads_count_dict = {}
    with pysam.FastxFile(fqfile) as fq:
        for entry in fq:
            attr = entry.name.split('_')
            barcode = attr[0]
            barcode_reads_dict[barcode].append(entry)
        for barcode in barcode_reads_dict.keys():
            barcode_reads_count_dict[barcode] = len(barcode_reads_dict[barcode])

    if bclist and bclist != 'None':
        with open(bclist, 'r') as bclist:
            barcode_annotation = list(bclist)
            barcodes = list(barcode_reads_count_dict.keys())
            barcodes_to_use = barcode_annotation.intersection(barcodes)
        if topn and topn != 'None':
            if int(topn) < len(barcodes_to_use):
                barcodes_to_use = sorted(barcodes_to_use.items(), key=lambda item: item[1], reverse=True)
                barcode_useful = barcodes_to_use[:int(topn)]
            else:
                barcode_useful = barcodes_to_use
        else:
            barcode_useful = barcodes_to_use
    else:
        barcode_reads_count_dict = sorted(barcode_reads_count_dict.items(), key=lambda item: item[1], reverse=True)
        if topn and topn != 'None':
            if int(topn) < len(barcode_reads_count_dict):
                barcode_useful = barcode_reads_count_dict[:int(topn)]
            else:
                barcode_useful = barcode_reads_count_dict
        else:
            barcode_useful = barcode_reads_count_dict

    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)
    for barcode in barcode_useful:
        fastq_file = f'{fastq_dir}/{barcode}.fq'
        with open(fastq_file, 'w') as f:
            for entry in barcode_reads_dict[barcode]:
                f.write(str(entry) + '\n')
