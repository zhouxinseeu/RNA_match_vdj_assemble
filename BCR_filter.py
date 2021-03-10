import pandas as pd
import argparse
from Bio.Seq import Seq


def filter_by_tpm(bracer_outdir, match_dir):
    data = pd.read_csv(f'{bracer_outdir}/filtered_BCR_summary/changeodb.tab', sep='\t')
    data = data[data['FUNCTIONAL'] == True]
    cell_name = set(list(data['CELL']))
    filtered = pd.DataFrame()
    for name in cell_name:
        count_cell = data[data['CELL'] == name]
        count_h = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'H'])
        count_k = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'K'])
        count_l = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'L'])
        count_k_l = count_k.append(count_l)
        if count_h.empty is not True:
            count_h = count_h.sort_values(by='TPM', ascending=False)
            count_h = count_h.head(1)
            filtered = filtered.append(count_h, ignore_index=True)
        if count_k_l.empty is not True:
            count_k_l = count_k_l.sort_values(by='TPM', ascending=False)
            count_k_l = count_k_l.head(1)
            filtered = filtered.append(count_k_l, ignore_index=True)
    filtered.to_csv(f'{bracer_outdir}/filtered_BCR_summary/filtered_count.tsv', sep='\t')

    filtered_h = filtered[filtered['LOCUS'] == 'H']
    filtered_k = filtered[filtered['LOCUS'] == 'K']
    filtered_l = filtered[filtered['LOCUS'] == 'L']
    filtered_h_count = filtered_h.shape[0]
    filtered_k_count = filtered_k.shape[0]
    filtered_l_count = filtered_l.shape[0]
    stat_string_1 = "BCR_H reconstruction:\t{}\nBCR_K reconstruction:\t{}\nBCR_L reconstruction:\t{}\n".format(
        filtered_h_count, filtered_k_count, filtered_l_count)

    paired_cell = pd.DataFrame(filtered['CELL'].value_counts())
    unpaired_cell = paired_cell[paired_cell['CELL'] == 1]
    paired_cell = paired_cell[paired_cell['CELL'] == 2]
    paired_k = 0
    paired_l = 0

    clones = pd.DataFrame()
    cells = list(paired_cell.index)
    aaseqs = []

    for cell in cells:
        if 'K' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
            paired_k += 1
        elif 'L' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
            paired_l += 1
        tep = filtered[filtered['CELL'] == cell]
        tep_loci = list(tep['LOCUS'])
        cdr3 = list(tep['JUNCTION'])
        aaseq = []
        for seq in cdr3:
            seq = Seq(seq)
            seq = seq.translate()
            aaseq.append(seq)
        string = 'IG{}:{};IG{}:{}'.format(tep_loci[0], aaseq[0], tep_loci[1], aaseq[1])
        aaseqs.append(string)

    for cell in list(unpaired_cell.index):
        cells.append(cell)
        locus = list(filtered[filtered['CELL'] == cell]['LOCUS'])
        cdr3 = list(filtered[filtered['CELL'] == cell]['JUNCTION'])
        seq = Seq(cdr3[0])
        seq = seq.translate()
        string = 'IG{}:{}'.format(locus[0], seq)
        aaseqs.append(string)

    clones['CELLS'] = cells

    stat_string_2 = "Paired HK productive reconstruction:\t{}\nPaired HL productive reconstruction:\t{}\n".format(
        paired_k, paired_l
    )

    clones["aaseq"] = aaseqs
    clone_count = pd.DataFrame(clones['aaseq'].value_counts())
    clone_count.columns = ["frequency"]
    proportation = []
    sum = clone_count['frequency'].sum()
    for f in list(clone_count['frequency']):
        p = f/sum
        proportation.append(p)
    clone_count['proportation'] = proportation
    clone_count = clone_count.reset_index()
    clone_count.rename(columns={'index': 'aaseq'}, inplace=True)
    clone_count.to_csv(f'{bracer_outdir}/filtered_BCR_summary/clone_count.tsv', sep='\t')

    data_10X = pd.read_csv(f'{match_dir}/clonotypes.csv', sep=',')
    data_10X_top10 = data_10X.head(10)
    matched_clone = 0
    matched_aa = []
    for aa in list(data_10X_top10['cdr3s_aa']):
        if aa in aaseqs:
            matched_clone += 1
            matched_aa.append(aa)
    stat_string_3 = 'matched clones: {}\n'.format(matched_clone)

    with open(f'{bracer_outdir}/filtered_BCR_summary/stat.txt', 'w') as s:
        s.write(stat_string_1)
        s.write(stat_string_2)
        s.write(stat_string_3)
        for aa in matched_aa:
            aa = aa + '\n'
            s.write(aa)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bracer_outdir", help='bracer summarize output dir', required=True)
    parser.add_argument('--match_dir', help='10X outdir', required=True)
    args = parser.parse_args()
    filter_by_tpm(args.bracer_outdir, args.match_dir)


if __name__ == "__main__":
    main()
