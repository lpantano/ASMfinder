"""
Functions to spot hemy regions
"""
import os.path as op
from collections import defaultdict, Counter
import pybedtools
import pysam

from bcbio.distributed.transaction import file_transaction
from bcbio.utils import append_stem


def is_good_cpg(frmt, record):
    if record[6] != "PASS":
        return False
    if int(frmt['CU']) > 10 and int(frmt['CM']) > 10:
        return True

def _genotype(alleles):
    if alleles[0] == alleles[1]:
        return "homoz"
    else:
        return "heteroz"


def is_good_het(frmt, record):
    depth = sum(map(int, frmt['DP4'].split(','))[1:])
    if _genotype(frmt['GT'].split("/")) == "heteroz" and int(frmt['DP']) > 10 and depth > 10 and record[6] == "PASS":
        return True


def cpg_het_pairs(cpgvcf, snpvcf, out_file, workdir):
    """
    Detect het close to hemi-met sites
    """
    cpg_filter = op.join(workdir, op.basename(append_stem(cpgvcf, "_filtered")))
    snp_filter = op.join(workdir, op.basename(append_stem(snpvcf, "_filtered")))

    with open(cpg_filter, 'w') as out_handle:
        with open(cpgvcf) as in_handle:
            for line in in_handle:
                if line.startswith("#"):
                    continue
                record = line.strip().split("\t")
                # print record
                header, frmt = record[8], record[9]
                frmt = dict(zip(header.split(":"), frmt.split(':')))
                if is_good_cpg(frmt, record):
                    print >>out_handle, line

    with open(snp_filter, 'w') as out_handle:
        with open(snpvcf) as in_handle:
            for line in in_handle:
                if line.startswith("#"):
                    continue
                record = line.strip().split("\t")
                header, frmt = record[8], record[9]
                frmt = dict(zip(header.split(":"), frmt.split(':')))
                if is_good_het(frmt, record):
                    print >>out_handle, line
    res = pybedtools.BedTool(cpg_filter).window(snp_filter, w=75)
    with open(out_file, 'w') as out_handle:
        for record in res:
            if record[1] != record[11]:
                print >>out_handle, record


def _make_linkage(bam_file, chrom, cpg, snp):
    pairs = defaultdict(dict)
    pairs = _pairs_matrix(bam_file, [chrom, cpg-1, cpg], pairs, 'cpg')
    pairs = _pairs_matrix(bam_file, [chrom, snp-1, snp], pairs, 'snp')
    link = Counter()
    for pair in pairs:
        print pair
        print pairs[pair]
        if len(pairs[pair].keys()) > 1:
            link["/".join(pairs[pair].values())] += 1
    print link


def _pairs_matrix(bam_file, region, pileup, tag):
    """
    Get reads from the cpg region and pairs
    cpg nt with snp nt
    """
    c, s, e = region
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for pileupcolumn in samfile.pileup(c, s, e):
        if pileupcolumn.pos == s:
            # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
                    pileup[pileupread.alignment.query_name].update({tag: pileupread.alignment.query_sequence[pileupread.query_position]})

    return pileup


def get_het(in_vcf, region, sample, out_file):
    res = pybedtools.BedTool(in_vcf).intersect(b=region, wo=True)
    with file_transaction(out_file) as tx_out:
        with open(tx_out, 'w') as out_handle:
            # print >> out_handle, "chrom\tstart\tend\tgen\dp4\tstrand\tgene\tsample"

            for record in res:
                gene = record[-2]
                chrom, pos, info, header, frmt = record[0], int(record[1]), record[7], record[8], record[9]
                # cs = info.split(';')[0].split('=')[1]
                frmt = dict(zip(header.split(":"), frmt.split(':')))
                # if _genotype(frmt['GT'].split("/")) == "heteroz" and int(frmt['DP']) > 10 and int(frmt['DP4']) > 10 and record[6] == "PASS":
                if is_good_het(frmt, record):
                    tag = "%s-%s-%s-%s" % (frmt['GT'], frmt['DP'], gene, sample)
                    print >> out_handle, "%s\t%s\t%s\t%s\t.\t+" % (chrom, pos, pos + 1, tag )


