import pysam
import argparse
import pandas as pd


def imgt_assembly_compare(imgtfile, assemblyfile):
    seq_cal = pd.DataFrame()
    cdr3_imgt = []
    cdr3_length_imgt = {}
    with pysam.FastxFile(imgtfile) as i:
        for record in i:
            name = record.name
            tags = name.split("|")
            length = len(record.sequence)
            receptor = tags[1]
            locus = tags[2]
            contig_name = tags[3]
            vj_pair = tags[4]
            cdr3_imgt.append(vj_pair)
            cdr3_length_imgt[vj_pair] = length

    cdr3_assembly = []
    cdr3_length_assembly = {}
    with pysam.FastxFile(assemblyfile) as a:
        for record in a:
            name = record.name
            tags = name.split("|")
            length = len(record.sequence)
            receptor = tags[1]
            locus = tags[2]
            contig_name = tags[3]
            vj_pair = tags[4]
            cdr3_assembly.append(vj_pair)
            cdr3_length_assembly[vj_pair] = length

    imgt_assembly = [seq for seq in cdr3_imgt if seq in cdr3_assembly]
    seq_cal['vjpair'] = imgt_assembly
    imgt_len = [cdr3_length_imgt[v] for v in imgt_assembly]
    assembly_len = [cdr3_length_assembly[v] for v in imgt_assembly]
    seq_cal['imgt_len'] = imgt_len
    seq_cal['assembly_len'] = assembly_len
    seq_cal.to_csv('len_count.tsv', sep='\t')
    same = 'same:' + str(len(imgt_assembly)) + "\n"
    imgt = 'imgt:' + str(len(cdr3_imgt)) + "\n"
    assembly = 'assembly:' + str(len(cdr3_assembly)) + "\n"
    with open('stat.txt', 'w') as f:
        f.write(same)
        f.write(imgt)
        f.write(assembly)
    return imgt_assembly


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--imgtfile', help='imgt results', required=True)
    parser.add_argument('--assemblyfile', help='assembly results', required=True)
    args = parser.parse_args()
    imgt_assembly_compare(args.imgtfile, args.assemblyfile)


if __name__ == "__main__":
    main()
