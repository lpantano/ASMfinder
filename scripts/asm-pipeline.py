"""
calculate coverage across a list of regions
"""
import os
import os.path as op
import yaml

from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
from ichwrapper import cluster, arguments
# import pandas as pd
# from collections import Counter, defaultdict

import pybedtools
# import vcf

# from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, splitext_plus, tmpfile, safe_makedir
from bcbio.install import _get_data_dir

# from ecov.total import _calc_total_exome_coverage
# from ecov.bias import calculate_bias_over_multiple_regions
# from ecov.variants import calc_variants_stats
from asm.trimming import prepare
from asm.align import create_bam
from asm.bissnp import call_variations
from asm.report import create_report
from asm.select import get_het, is_good_cpg, cpg_het_pairs


def _update_algorithm(data, resources):
    """
    Update algorithm dict with new cores set
    """
    new_data = []
    for sample in data:
        sample[0]['config']['algorithm'] = resources
        new_data.append(sample)
    return new_data

def _prepare_samples(args):
    """
    create dict for each sample having all information
    """
    if args.galaxy:
        system_config = args.galaxy
    else:
        system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    config = yaml.load(open(system_config))
    config['algorithm'] = {}
    data = []
    for sample in args.files:
        dt = {}
        dt['name'] = splitext_plus(op.basename(sample))[0]
        dt['config'] = config
        dt['fastq'] = op.abspath(sample)
        data.append([dt])
    return data

def select_regions(args):
    """
    select regions and create coverage plots
    """
    assert args.files, "Need a set of fastq files"
    assert args.out, "Need --out"
    region = os.path.abspath(args.region)
    workdir = 'select'
    safe_makedir(workdir)
    out_file = os.path.join(workdir, splitext_plus(args.out)[0] + "_cpg.bed")
    out_snp_file = os.path.join(workdir, splitext_plus(args.out)[0] + '_snp.bed')
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            with open(tx_out, 'w') as out_handle:
                # print >> out_handle, "chrom\tstart\tend\tcu\tcm\tstrand\tgene\tsample"
                for in_vcf in args.files:
                    snp_file = in_vcf.replace("rawcpg", "rawsnp")
                    sample = splitext_plus(os.path.basename(in_vcf))[0].split("_")[0]
                    get_het(snp_file, region, sample, out_snp_file)
                    res = pybedtools.BedTool(in_vcf).intersect(b=region, wo=True)
                    # cmd = ("bedtools intersect -u -a {in_vcf} -b {region} > {tx_tmp}")
                    # do.run(cmd.format(**locals()), "selecting %s" % in_vcf)

                    for record in res:
                        gene = record[-2]
                        chrom, pos, info, header, frmt = record[0], int(record[1]), record[7], record[8], record[9]
                        cs = info.split(';')[0].split('=')[1]
                        frmt = dict(zip(header.split(":"), frmt.split(':')))
                        if is_good_cpg(frmt):
                            tag = "%s-%s-%s-%s" % (frmt['CU'], frmt['CM'], gene, sample)
                            print >> out_handle, "%s\t%s\t%s\t%s\t.\t%s" % (chrom, pos, pos + 1, tag, cs)


def detect_positions(data, args):
    assert args.reference, "Need --reference"
    assert args.index, "Need --index"
    assert args.files, "Need a set of fastq files"
    assert args.snp, "Need --snp"

    resources = {'name': 'trimming', 'mem': 4, 'cores': 1}
    data = _update_algorithm(data, resources)
    data = cluster.send_job(prepare, data, args, resources)

    resources = {'name': 'align', 'mem': 2, 'cores': 8}
    data = _update_algorithm(data, resources)
    data = cluster.send_job(create_bam, data, args, resources)

    resources = {'name': 'bissnp', 'mem': 3, 'cores': 8}
    data = _update_algorithm(data, resources)
    data = cluster.send_job(call_variations, data, args, resources)

    resources = {'name': 'report', 'mem': 2, 'cores': 5}
    data = _update_algorithm(data, resources)
    data = cluster.send_job(create_report, data, args, resources)


def _find_bam(bam_files, sample):
    """
    Find the most similar file name
    """
    score = 0
    candidate = None
    for fn in bam_files:
        sc = sum(a == b for a, b in zip(op.basename(sample), op.basename(fn)))
        if sc > score:
            score = sc
            candidate = fn
    return candidate


def link_sites(args):
    assert args.files, "Need a set of fastq files"
    workdir = 'link'
    workdir = op.abspath(safe_makedir(workdir))
    vcf_files = [fn for fn in args.files if fn.endswith('vcf')]
    bam_files = [fn for fn in args.files if fn.endswith('bam')]
    for in_vcf in vcf_files:
        snp_file = in_vcf.replace("rawcpg", "rawsnp")
        sample = splitext_plus(os.path.basename(in_vcf))[0].split("_")[0].replace(".rawcpg", "")
        bam_file = _find_bam(bam_files, sample)
        assert bam_file, "No bam file associated to vcf %s" % in_vcf
        out_file = op.join(workdir, sample + "_pairs.tsv")
        cpg_het_pairs(in_vcf, snp_file, bam_file, out_file, workdir)

if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser = arguments.myargs(parser)
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--index", help="genome index for bismark.")
    parser.add_argument("--is_rrbs", action="store_true", help="RRBS data.")
    parser.add_argument("--is_directional", action="store_true", help="is directional sequencing.")
    parser.add_argument("--snp", help="SNPdb database.")
    # parser.add_argument("--galaxy", help="bcbio galaxy resources.")

    parser.add_argument("--out", help="output file.")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=1, help="Calculate bam stats", choices=['select', 'detection', 'link'])
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    args = parser.parse_args()

    if args.run == 'select':
        select_regions(args)
    if args.run == "detection":
        data = _prepare_samples(args)
        detect_positions(data, args)
    if args.run == 'link':
        link_sites(args)
