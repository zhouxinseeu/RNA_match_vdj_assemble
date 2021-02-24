import pysam
import pandas as pd
from collections import defaultdict
import os
from Bio.Seq import Seq


def get_fastq_to_assemble(fqfile, fastq_dir, reversed, bclist=None, topn=None):
    barcode_reads_dict = defaultdict(list)  # all barcodes from BCR fqfile paired with reads
    reads_count_dict = {}  # all barcodes and reads num for each barcode
    all_barcodes = []  # all barcodes
    with pysam.FastxFile(fqfile) as fq:
        for entry in fq:
            attr = entry.name.split('_')
            barcode = attr[0]
            all_barcodes.append(barcode)
            barcode_reads_dict[barcode].append(entry)
        for barcode in list(barcode_reads_dict.keys()):
            reads_count_dict[barcode] = len(barcode_reads_dict[barcode])

    if bclist and bclist != 'None':  # if bclist provided, compared with fqfile barcode and get intersection
        with open(bclist, 'r') as bclist:
            barcodes_from_matchfile = list(bclist)
            barcodes_for_match = []
            if reversed == 'True':  # whether reversed?
                for barcode in barcodes_from_matchfile:
                    barcode = barcode.strip('\n')
                    barcode = Seq(barcode)
                    barcode = barcode.reverse_complement()
                    barcodes_for_match.append(barcode)
            elif reversed == 'False':
                for barcode in barcodes_from_matchfile:
                    barcode = barcode.strip('\n')
                    barcodes_for_match.append(barcode)
            barcodes_to_use = [barcode for barcode in barcodes_for_match if barcode in all_barcodes]
            # barcodes in both RNA data and BCR data

            if topn and topn != 'None':  # if topn provided, copmared with length of barcode_to_use
                topn = int(topn)
                if topn < len(barcodes_to_use):
                    reads_count_dict_to_sort = {barcode: reads_count_dict[barcode] for barcode in barcodes_to_use}
                    reads_count_dict_to_sort = sorted(reads_count_dict_to_sort.items(), key=lambda item: item[1],
                                                      reverse=True)
                    barcode_useful = reads_count_dict_to_sort[:topn]
                    barcode_reads_useful = {i[0]: barcode_reads_dict[i[0]] for i in barcode_useful}
                    # matched barcodes and reads
                else:
                    barcode_reads_useful = {barcode: barcode_reads_dict[barcode] for barcode in barcodes_to_use}
            else:
                barcode_reads_useful = {barcode: barcode_reads_dict[barcode] for barcode in barcodes_to_use}
    else:  # if bclist not provided, compare all_barcode with topn
        if topn and topn != 'None':
            topn = int(topn)
            reads_count_dict = sorted(reads_count_dict.items(), key=lambda item: item[1], reverse=True)
            if topn < len(reads_count_dict):
                reads_count_dict = reads_count_dict[:topn]
                barcode_reads_useful = {i[0]: barcode_reads_dict[i[0]] for i in reads_count_dict}
            else:
                barcode_reads_useful = barcode_reads_dict
        else:
            barcode_reads_useful = barcode_reads_dict

    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)
    i = 1
    for barcode in list(barcode_reads_useful.keys()):
        fastq_file = f'{fastq_dir}/{i}.fq'
        with open(fastq_file, 'w') as f:
            for entry in barcode_reads_useful[barcode]:
                f.write(str(entry) + '\n')
        i += 1
