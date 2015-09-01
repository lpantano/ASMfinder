import os.path as op

from collections import defaultdict, Counter

import pysam
from pybedtools import BedTool
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
                    # nt = nt.lower() if strand == "-" else nt
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

def raw(bams, chrom, cpg, snp, strand, als):
    cpg, snp = int(cpg), int(snp)
    s, e = min(cpg, snp), max(cpg, snp)
    s, e = s - 1, e + 1
    res = []
    for bam_file in bams:
        link = defaultdict(Counter)
        sample = op.basename(splitext_plus(bam_file)[0])
        pairs = _pairs_matrix(bam_file, (chrom, s, e), strand)
        for read in pairs:
            if cpg in pairs[read] and snp in pairs[read]:
                link[pairs[read][snp]][pairs[read][cpg]] += 1

        nom_meth, meth = ('T', 'C') if strand == "+" else ('A', 'G')
        for snp_a in link:
            total = sum(link[snp_a].values())
            if snp_a in als:
                res.append("%s %s %s %s %s %s %s %s %s" % (chrom, snp, snp_a, total, cpg, strand,link[snp_a][meth], link[snp_a][nom_meth], sample))
    return res

NT_METH = {'C': '+', 'G': '-'}

def region_selection(pair_file, bams, prefix, bed_file=None, min_samples=1):
    common = Counter()
    print min_samples
    for fn in pair_file:
        print fn
        with open(fn) as in_handle:
            for line in in_handle:
                cols = line.strip().split()
                if cols[2] in NT_METH:
                    strand = NT_METH[cols[2]]
                    common["%s:%s:%s:%s:%s" % (cols[0], cols[1], cols[3], strand, cols[4])] += 1
    string = ""
    for c in common:
        if common[c] > int(min_samples):
            cols = c.split(":")
            string += "%s\t%s\t%s\t%s\n" % (cols[0], cols[1], int(cols[1]) + 1, c)

    if bed_file:
        valid_pairs = BedTool(string, from_string=True).intersect(BedTool(bed_file))
    else:
        valid_pairs = BedTool(string, from_string=True)

    with open(prefix + ".tsv", 'w') as out_handle:
        print >>out_handle, "chrom snp_pos allele counts cpg_pos strand meth_C non_meth_C sample"
        for pair in valid_pairs:
            c = pair[3].strip().split(":")
            print >>out_handle, "\n".join(raw(bams, c[0], int(c[1]) - 1, int(c[2]) - 1, c[3], c[4].split("/")))

def region_selection_by_read(pair_file, cpg_vcf, snp_vcf, bams, prefix, bed_file, min_samples=1):
    common = Counter()
    for fn in pair_file:
        print fn
        with open(fn) as in_handle:
            for line in in_handle:
                cols = line.strip().split()
                if cols[2] in NT_METH:
                    strand = NT_METH[cols[2]]
                    common["%s:%s:%s:%s:%s" % (cols[0], cols[1], cols[3], strand, cols[4])] += 1
    string = ""
    for c in common:
        if common[c] > int(min_samples):
            cols = c.split(":")
            string += "%s\t%s\t%s\t%s\n" % (cols[0], cols[1], int(cols[1]) + 1, c)

    bed_valid_asm = BedTool(string, from_string=True)
    with open(prefix + ".tsv", 'w') as out_handle:
        print >>out_handle, "sample read chrom position strand call type asm"
        for line in BedTool(bed_file):
            single_line = BedTool("\t".join(line), from_string=True)
            this_region_asm = BedTool(bed_valid_asm).intersect(single_line)

            asm_cpg, asm_snp = set(), set()
            for pair in this_region_asm:
                c = pair[3].strip().split(":")
                asm_cpg.add(str(int(c[1]) - 1))
                asm_snp.add(str(int(c[2]) - 1))

            for bam in bams:
                this_cpg_vcf = bam.replace("_recal1.bam", ".rawcpg.vcf")
                this_snp_vcf = bam.replace("_recal1.bam", ".rawsnp.vcf")
                this_region_cpg = [[l[1], l[7].split(";")[0][3:]] for l in BedTool(this_cpg_vcf).intersect(single_line)]
                this_region_snp = [l[1] for l in BedTool(this_snp_vcf).intersect(single_line)]
                res = _get_reads(bam, [line[0], int(line[1]), int(line[2])],
                                 asm_snp, asm_cpg, this_region_cpg, this_region_snp, c[3])
                # print "\n".join(res)
                print >>out_handle, "\n".join(res)

NT_VALID = {'C': '+', 'G': '-', 'T': '+', 'A': '-'}

def _get_reads(bam_file, region, asm_snp, asm_cpg, cpg_info, snp, strand):
    """
    Get reads from the cpg region and pairs
    cpg nt with snp nt
    """
    name = op.basename(op.splitext(bam_file)[0])
    line = []
    c, s, e = region
    samfile = pysam.AlignmentFile(bam_file, "rb")
    cpg = dict(zip([item[0] for item in cpg_info], [item[1] for item in cpg_info]))

    for pileupcolumn in samfile.pileup(c, s, e):
        # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
                st = "-" if pileupread.alignment.is_reverse else "+"
                # if st == strand:
                nt = pileupread.alignment.query_sequence[pileupread.query_position]
                pos = str(pileupread.alignment.pos + pileupread.query_position)
                # nt = nt.lower() if strand == "-" else nt
                nt_type = []
                if nt == "N":
                    continue
                if pos in asm_snp:
                    nt_type = [name, pileupread.alignment.query_name, c, pos, st, nt, "SNP", "ASM"]
                elif pos in asm_cpg:
                    nt_type = [name, pileupread.alignment.query_name, c, pos, st, nt, "CpG", "ASM"]
                elif pos in cpg and NT_VALID[nt] == cpg[pos] and st == cpg[pos]:
                    nt_type = [name, pileupread.alignment.query_name, c, pos, st, nt, "CpG", "None"]
                    # print [pos, st, pileupread.alignment.is_reverse, nt, cpg[pos]]
                elif pos in snp:
                    nt_type = [name, pileupread.alignment.query_name, c, pos, st, nt, "SNP", "None"]
                if nt_type:
                    line.append(" ".join(nt_type))

    return line
