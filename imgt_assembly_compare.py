import pysam
import argparse


def imgt_assembly_compare(imgtfile, assemblyfile):
    cdr3_imgt = []
    with pysam.FastxFile(imgtfile) as i:
        for record in i:
            name = record.name
            tags = name.split("|")
            receptor = tags[1]
            locus = tags[2]
            contig_name = tags[3]
            vj_pair = tags[4]
            cdr3_imgt.append(vj_pair)

    cdr3_assembly = []
    with pysam.FastxFile(assemblyfile) as a:
        for record in a:
            name = record.name
            tags = name.split("|")
            receptor = tags[1]
            locus = tags[2]
            contig_name = tags[3]
            vj_pair = tags[4]
            cdr3_assembly.append(vj_pair)

    imgt_assembly = [seq for seq in cdr3_imgt if seq in cdr3_assembly]
    same = 'same:' + str(len(imgt_assembly))
    imgt = 'imgt:' + str(len(cdr3_imgt))
    assembly = 'assembly:' + str(len(cdr3_assembly))
    with open('stat.txt', 'w') as f:
        f.write(same)
        f.write(imgt)
        f.write(assembly)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--imgtfile', help='imgt results', required=True)
    parser.add_argument('--assemblyfile', help='assembly results', required=True)
    args = parser.parse_args()
    imgt_assembly_compare(args.imgtfile, args.assemblyfile)


if __name__ == "__main__":
    main()
