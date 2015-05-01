"""
calculate coverage across a list of regions
"""
import os

from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
from ichwrapper import cluster, arguments
import pandas as pd
from collections import Counter, defaultdict

import pybedtools
import vcf

from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, splitext_plus, tmpfile, safe_makedir

# from ecov.total import _calc_total_exome_coverage
# from ecov.bias import calculate_bias_over_multiple_regions
# from ecov.variants import calc_variants_stats


def calculate_bam(args):
    """
    samtools flagstat output
    """
    stats = defaultdict(list)
    samples = []
    for bam in args.bams:
        # sample = os.path.splitext(bam)[0].split("-")[0]
        sample = os.path.basename(bam).split("-")[0]
        out_sample = sample + ".flagstat"
        if not file_exists(out_sample):
            with file_transaction(out_sample) as tx_out:
                cmd = ("samtools flagstat {bam} > {tx_out}")
                do.run(cmd.format(**locals()), "bam stats for %s" % bam)
        with open(out_sample) as in_handle:
            for line in in_handle:
                if line.find("mapQ") == -1:
                    stats[line.strip().split(" + 0 ")[1].split("(")[0].strip()].append(line.strip().split(" + ")[0])
                # print stats[sample]
        samples.append(sample)
    with open(args.out, 'w') as out_handle:
        out_handle.write("\t".join(['measure'] + samples) + '\n')
        for feature in stats:
            out_handle.write("\t".join([feature] + stats[feature]) + "\n")

def select_regions(args):
    """
    select regions and create coverage plots
    """
    region = os.path.abspath(args.region)
    workdir = 'select'
    safe_makedir(workdir)
    out_file = os.path.join(workdir, args.out)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            with open(tx_out, 'w') as out_handle:
                print >> out_handle, "chrom\tstart\tend\tcu\tcm\tstrand\tgene\tsample"
                for in_vcf in args.files:
                    sample = splitext_plus(os.path.basename(in_vcf))[0].split("_")[0]
                    res = pybedtools.BedTool(in_vcf).intersect(b=region, wo=True)
                    # cmd = ("bedtools intersect -u -a {in_vcf} -b {region} > {tx_tmp}")
                    # do.run(cmd.format(**locals()), "selecting %s" % in_vcf)

                    for record in res:
                        gene = record[-2]
                        chrom, pos, info, header, frmt = record[0], int(record[1]), record[7], record[8], record[9]
                        cs = info.split(';')[0].split('=')[1]
                        frmt = dict(zip(header.split(":"), frmt.split(':')))
                        if int(frmt['CU']) > 10 and int(frmt['CM']) > 10:
                            print >> out_handle, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, pos, pos + 1, frmt['CU'], frmt['CM'], cs, gene, sample)

def example_cluster(args):
    resources = {'name': 'bias', 'mem': 1, 'cores': 1}
    # cluster.send_job(calculate_bias_over_multiple_regions, args.bams, args, resources)


if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser = arguments.myargs(parser)
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--out", required=True, help="output file.")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=1, help="Calculate bam stats", choices=['select'])
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    args = parser.parse_args()

    if os.path.exists(args.out):
        os.remove(args.out)

    if args.run == 'select':
        select_regions(args)
