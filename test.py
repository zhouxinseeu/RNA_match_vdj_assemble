import pandas as pd
import argparse


gencode = {
      'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
      'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
      'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
      'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
      'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
      'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
      'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
      'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
      'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
      'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
      'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
      'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
      'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
      'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
      'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
      'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}


def translate_frameshifted(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate


def filter_by_tpm(bracer_outdir):
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
    paired_cell = paired_cell[paired_cell['CELL'] == 2]
    paired_k = 0
    paired_l = 0
    for cell in list(paired_cell.index):
        if 'K' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
            paired_k += 1
        elif 'L' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
            paired_l += 1
    stat_string_2 = "Paired HK productive reconstruction:\t{}\nPaired HL productive reconstruction:\t{}\n".format(
        paired_k, paired_l
    )
    clones = pd.DataFrame(filtered['CLONE'].value_counts())
    clones.columns = ['clone_count']
    aaseq = []
    for clone in list(clones.index):
        sequence = str(filtered[filtered['CLONE'] == str(clone)].head(1)['CDR3'])
        aa = translate_frameshifted(sequence)
        aaseq.append(aa)

    clones['aaseq'] = aaseq
    clones.to_csv(f'{bracer_outdir}/filtered_BCR_summary/clone_count.tsv', sep='\t')
    with open(f'{bracer_outdir}/filtered_BCR_summary/stat.txt', 'a') as s:
        s.write(stat_string_1)
        s.write(stat_string_2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bracer_outdir", help='bracer summarize output dir', required=True)
    args = parser.parse_args()
    filter_by_tpm(args.bracer_outdir)


if __name__ == "__main__":
    main()

