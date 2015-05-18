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
    alt_depth = sum(map(int, frmt['DP4'].split(','))[2:])
    ref_depth = sum(map(int, frmt['DP4'].split(','))[:2])
    if record[6] != "PASS":
        return False
    if int(ref_depth) > 3 and int(alt_depth) > 3:
        return True


def _genotype(alleles):
    if alleles[0] == alleles[1]:
        return "homoz"
    else:
        return "heteroz"


def is_good_het(frmt, record):
    depth = sum(map(int, frmt['DP4'].split(','))[2:])
    # if _genotype(frmt['GT'].split("/")) == "heteroz" and int(frmt['DP']) > 3 and depth > 3 and record[6] == "PASS":
    if _genotype(frmt['GT'].split("/")) == "heteroz" and int(frmt['DP']) > 3:
        return True


def _get_strand(record):
    return record[7].split(";")[0].split("=")[1]


def cpg_het_pairs(cpgvcf, snpvcf, bam_file, out_file, workdir):
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
                # if record[1] == "19889634":
                link, link_as, align = _make_linkage(bam_file, record[0], int(record[1]), int(record[11]), _get_strand(record), record[13])
                res = "%s %s %s %s %s/%s %s %s" % (record[0], record[1], record[3], record[11], record[13], record[14], link, link_as)
                if len(link) > 1:
                    print >>out_handle, res
                    # print >>out_handle, '\n'.join(align)


def _make_linkage(bam_file, chrom, cpg, snp, cpg_st, snp_ref):
    start, end = [cpg-1, snp] if cpg-1 < snp else [snp, cpg-1]
    pairs = _pairs_matrix(bam_file, [chrom, start, end], cpg-1, snp-1)
    link = Counter()
    link_as = Counter()
    align = []
    for pair in pairs:
        # print pair
        # print pairs[pair]
        # print strand
        if len(pairs[pair].keys()) == 1:
            continue
        nts = [pairs[pair]['snp'].split(":")[0], pairs[pair]['cpg'].split(":")[0]]
        align.append("-".join(nts) if cpg < snp else "-".join(nts[::-1]))
        info_snp = pairs[pair]['snp'].split(":")
        if info_snp[1] == cpg_st:
            # print pairs[pair]
            if pairs[pair]['cpg']:
                info_cpg = pairs[pair]['cpg'].split(":")
                if info_cpg[1] == cpg_st:
                    link["v%s/c%s(%s)" % (info_snp[0], info_cpg[0], cpg_st)] += 1
        else:
            link_as["v%s(%s)" % (info_snp[0], cpg_st)] += 1
    # print link
    return link, link_as, align


def _pairs_matrix(bam_file, region, cpg, snp):
    """
    Get reads from the cpg region and pairs
    cpg nt with snp nt
    """
    pileup = defaultdict(dict)
    c, s, e = region
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for pileupcolumn in samfile.pileup(c, s, e):
        if pileupcolumn.pos == cpg or pileupcolumn.pos == snp:
            # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
                    strand = "-" if pileupread.alignment.is_reverse else "+"
                    tag = "cpg" if pileupcolumn.pos == cpg else "snp"
                    nt = pileupread.alignment.query_sequence[pileupread.query_position]
                    nt = nt.lower() if strand == "-" else nt
                    pileup[pileupread.alignment.query_name].update({tag: nt + ":%s" % strand})

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


