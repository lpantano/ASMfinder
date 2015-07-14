import os.path as op

from collections import defaultdict, Counter

import pysam
from bcbio.utils import splitext_plus


def _complementary(nt):
    if nt.lower() == "a":
        return "T"
    if nt.lower() == "t":
        return "A"
    if nt.lower() == "c":
        return "G"
    if nt.lower() == "g":
        return "C"

def _code(nt):
    if nt == "T":
        return "C"
    return "Cm"

def _pairs_matrix(bam_file, region, strand):
    """
    Get reads from the cpg region and pairs
    cpg nt with snp nt
    """
    pileup = defaultdict(dict)
    c, s, e = region
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for pileupcolumn in samfile.pileup(c, s, e):
        # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
                st = "-" if pileupread.alignment.is_reverse else "+"
                if st == strand:
                    nt = pileupread.alignment.query_sequence[pileupread.query_position]
                    nt = nt.lower() if strand == "-" else nt
                    pileup[pileupread.alignment.query_name].update({pileupcolumn.pos: nt})

    return pileup

def plot(bams, params, prefix):
    chrom, cpg, snp, strand = params.split(":")
    cpg, snp = int(cpg), int(snp)
    s, e = min(cpg, snp), max(cpg, snp)
    s, e = s - 10, e + 10
    with open(prefix + ".txt", 'w') as out_handle:
        for bam_file in bams:
            link = defaultdict(Counter)
            sample = op.basename(splitext_plus(bam_file)[0])
            pairs = _pairs_matrix(bam_file, (chrom, s, e), strand)
            for read in pairs:
                r = ""
                for pos in range(int(s), int(e)):
                    if pos in pairs[read]:
                        r += pairs[read][int(pos)]
                    else:
                        r += "-"

                    if cpg in pairs[read] and snp in pairs[read]:
                        link[pairs[read][snp]][pairs[read][cpg]] += 1

                print >>out_handle, r

            for snp_a in link:
                total = sum(link[snp_a].values())
                meth = link[snp_a]["C"] if "C" in link[snp_a] else 0
                print "%s %s %s" % (snp_a, float(meth)/total, sample)

