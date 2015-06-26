import os.path as op

from bcbio.utils import splitext_plus
from asm.select import _make_linkage


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

def plot(bams, params):
    chrom, cpg, snp, strand = params.split(":")
    snp, cpg = int(snp), int(cpg)
    for bam_file in bams:
        sample = op.basename(splitext_plus(bam_file)[0])
        pairs = _make_linkage(bam_file, chrom, cpg, snp, strand)
        for p in pairs[2]:
            nt1, nt2 = p.split("-")
            if strand == "-":
                nt1, nt2 = _complementary(nt1), _complementary(nt2)
            spaces = ["-"] * abs(cpg - snp)
            if snp < cpg:
                nt2 = _code(nt2)
            else:
                nt1 = _code(nt1)
            print "%s %s %s %s" % (nt1, " ".join(spaces), nt2, sample)
