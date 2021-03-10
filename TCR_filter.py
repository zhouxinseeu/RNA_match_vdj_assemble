import pandas as pd
import os
import argparse
from Bio.Seq import Seq
from itertools import islice
from collections import defaultdict


def tpm_count(outdir):
    rec = pd.read_csv(f'{outdir}/tracer/filtered_TCRAB_summary/recom'
                      f'binants.txt', sep='\t')
    productive = rec[rec['productive'] == True]
    productive['TPM'] = ''
    indx = productive.index.to_list()
    for i in indx:
        cell_name = productive.at[i, 'cell_name']
        rec_id = productive.at[i, 'recombinant_id']
        with open(f'{outdir}/tracer/{cell_name}/expression_quantification/abundance.tsv') as tsvf:
            for line in tsvf:
                if rec_id in line:
                    line = line.rstrip()
                    line = line.split('\t')
                    tpm = float(line[4])
                    productive.loc[i, 'TPM'] = tpm

    productive.to_csv(f'{outdir}/productive.txt', sep='\t')


def tcr_filter(outdir):
    data = pd.read_csv(f'{outdir}/productive.txt', sep='\t', index_col=0)
    cell_name = set(list(data['cell_name']))
    filtered = pd.DataFrame()
    for name in cell_name:
        count_data = data[data['cell_name'] == name]
        tra = count_data[count_data['locus'] == 'A']
        trb = count_data[count_data['locus'] == 'B']
        if tra.empty is not True:
            tra = tra.sort_values(by='TPM', ascending=False)
            tra = tra.head(1)
            filtered = filtered.append(tra, ignore_index=True)
        if trb.empty is not True:
            trb = trb.sort_values(by='TPM', ascending=False)
            trb = trb.head(1)
            filtered = filtered.append(trb, ignore_index=True)
    filtered.to_csv(f'{outdir}/filtered.txt', sep='\t')

    count_a = len(filtered[filtered['locus'] == 'A'])
    count_b = len(filtered[filtered['locus'] == 'B'])
    paired_cell = pd.DataFrame(filtered['cell_name'].value_counts())
    unpaired_cell = paired_cell[paired_cell['cell_name'] == 1]
    paired_cell = paired_cell[paired_cell['cell_name'] == 2]
    paired_cell = paired_cell.index.to_list()
    string1 = f'productive TRA:{count_a}\nproductive TRB:{count_b}\npaired TRA and TRB:{len(paired_cell)}\n'
    with open(f'{outdir}/stat.txt', 'w') as f:
        f.write(string1)
    
    aaseqs = []
    for cell in paired_cell:
        temp = filtered[filtered['cell_name'] == cell]
        temp_loci = list(temp['locus'])
        temp_aaseq = list(temp['CDR3aa'])
        string = 'TR{}:{};TR{}:{}'.format(temp_loci[0], temp_aaseq[0], temp_loci[1], temp_aaseq[1])
        aaseqs.append(string)

    for cell in list(unpaired_cell.index):
        temp = filtered[filtered['cell_name'] == cell]
        temp_loci = list(temp['locus'])
        temp_aaseq = list(temp['CDR3aa'])
        string = 'TR{}:{}'.format(temp_loci[0], temp_aaseq[0])
        aaseqs.append(string)

    per_count_data = pd.DataFrame()
    per_count_data['CDR3aa'] = aaseqs
    per_count = pd.DataFrame(per_count_data['CDR3aa'].value_counts())
    per_count.columns = ["frequency"]
    proportation = []
    sum = per_count['frequency'].sum()
    for f in list(per_count['frequency']):
        p = f/sum
        proportation.append(p)
    per_count['proportation'] = proportation
    per_count = per_count.reset_index()
    per_count.rename(columns={'index': 'CDR3aa'}, inplace=True)
    per_count.to_csv(f'{outdir}/count.tsv', sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help='tracer output dir', required=True)
    args = parser.parse_args()
    tpm_count(args.outdir)
    tcr_filter(args.outdir)


if __name__ == "__main__":
    main()