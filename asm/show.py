import os.path as op

from collections import defaultdict

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


def plot(bams, params):
    chrom, s, e, strand = params.split(":")
    s, e = int(s), int(e)
    s, e = s - 10, e + 10
    for bam_file in bams:
        sample = op.basename(splitext_plus(bam_file)[0])
        pairs = _pairs_matrix(bam_file, (chrom, s, e), strand)
        for read in pairs:
            r = ""
            for pos in range(int(s), int(e)):
                if pos in pairs[read]:
                    r += pairs[read][int(pos)]
                else:
                    r += "-"
            print r
            # print "%s %s %s" % (read, pos, pairs[read][pos])
            #    nt2 = _code(nt2)
            #else:
            #    nt1 = _code(nt1)
            #print "%s %s %s %s" % (nt1, " ".join(spaces), nt2, sample)
