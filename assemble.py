import argparse
import os
from os import listdir
from os.path import isfile, join
from concurrent.futures import ProcessPoolExecutor
from celescope.tools.utils import genDict, format_number, log, read_barcode_file
import datetime


TRACER_PATH = '/SGRNJ/Public/Software/tracer/tracer'
CONF_PATH = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/tcr_fl/20201103/tracer_SGR.conf'
CONDA = 'vdjpuzzle1'
CONDA_SUB = 'celescope_tracer'
BRACER_PATH = '/SGRNJ03/randd/zhouxin/software/bracer/bracer'
BRACER_CONDA = 'bracer'
BRACER_CONF = '/SGRNJ03/randd/zhouxin/software/bracer/bracer.conf'


# 开始组装


@log
def bracer_summarise(outdir, sample):
    bracer_outdir = f'{outdir}/{sample}/bracer'
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{BRACER_PATH} summarise '
        f'-c {BRACER_CONF} ' 
        f'{bracer_outdir} '
        )
    bracer_summarise.logger.info(cmd)
    os.system(cmd)


def bracer(fq, outdir, species, sample):
    prefix = os.path.basename(fq).strip('.fq')
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{BRACER_PATH} assemble '
        f'--fragment_length 150 '
        f'--fragment_sd 5 '
        f'--single_end '
        f'--species {species} '
        f'-c {BRACER_CONF} '
        f'{prefix} '
        f'{outdir}/{sample}/bracer '
        f'{fq} '
    )
    os.system(cmd)


def tracer_summarise(outdir, sample):
    tracer_outdir = f'{outdir}/{sample}/tracer'
    cmd = (
        f'source activate {CONDA_SUB}; '
        f'{TRACER_PATH} summarise '
        f'-c {CONF_PATH} '
        f'{tracer_outdir} '
    )
    tracer_summarise.logger.info(cmd)
    os.system(cmd)


def tracer(fq, outdir, species, sample):
    prefix = os.path.basename(fq).strip('.fq')
    cmd = (
        f'source activate {CONDA}; '
        f'{TRACER_PATH} assemble '
        f'--fragment_length 150 '
        f'--fragment_sd 5 '
        f'--single_end '
        f'--species {species} '
        f'-c {CONF_PATH} '
        f'{fq} '
        f'{prefix} '
        f'{outdir}/{sample}/tracer '
    )
    os.system(cmd)


@log
def run_tracer(sample, outdir, fastq_dir, species, thread):
    fqs = [join(fastq_dir, f) for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    outdirs = [outdir] * len(fqs)
    species = [species] * len(fqs)
    samples = [sample] * len(fqs)
    if not os.path.exists(f'{outdir}/{sample}/tracer'):
        os.makedirs(f'{outdir}/{sample}/tracer')

    all_res = []
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(tracer, fqs, outdirs, species, samples):
            all_res.append(res)

    tracer_summarise(outdir, sample)


def run_bracer(sample, outdir, fastq_dir, species, thread):
    fqs = [join(fastq_dir, f) for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    outdirs = [outdir] * len(fqs)
    species = [species] * len(fqs)
    samples = [sample] * len(fqs)
    if not os.path.exists(f'{outdir}/{sample}/bracer'):
        os.makedirs(f'{outdir}/{sample}/bracer')

    all_res = []
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(bracer, fqs, outdirs, species, samples):
            all_res.append(res)

    bracer_summarise(outdir, sample)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help="assemble outdir")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument("--fastq_dir", help="fastq file dir")
    parser.add_argument('--mod', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], default="Hsap")
    parser.add_argument('--thread', help='thread')
    args = parser.parse_args()
    start_time = datetime.datetime.now()
    if args.mod == 'TCR':
        run_tracer(args.sample, args.outdir, args.fastq_dir, args.species, int(args.thread))
    elif args.mod == 'BCR':
        run_bracer(args.sample, args.outdir, args.fastq_dir, args.species, int(args.thread))
    end_time = datetime.datetime.now()
    time_report = 'assemble--start_at_{}--end_at_{}'.format(start_time, end_time)
    print(time_report)


if __name__ == "__main__":
    main()
