"""
Functions to spot hemy regions
"""
import os.path as op
import datetime
from collections import defaultdict, Counter
from tabulate import tabulate

from scipy import stats
import numpy as np

import pybedtools
import pysam
import vcf

from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import append_stem, file_exists, splitext_plus, safe_makedir
from bcbio.variation.vcfutils import bgzip_and_index


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

def _snp_veracity_both_strand(sense, anti):
    """
    Only if SNPs is detected in both strand with two alleles
    """
    gen_plus = sense.keys()
    gen_minus = anti.keys()
    allels1 = [g.split(":")[0].split("/")[0] for g in gen_plus]
    allels2 = [g.split(":")[0] for g in gen_minus]
    if len(allels1) == len(allels2):
        return True

def _read_pairs(gt):
    # print "read_pairs %s" % gt
    gt1 = gt.split(":")[0].split("/")[0]
    if gt.find("/") > -1:
        gt2 = gt.split(":")[0].split("/")[1]
        return (gt1, gt2)

def _get_total(gts, total):
    return [total[_read_pairs(gts[0][1])[0]], total[_read_pairs(gts[1][1])[0]]]

def _top_gt(gts):
    total = Counter()
    first = _read_pairs(gts[0][1])
    top = None
    for gt in gts:
        pair = _read_pairs(gt[1])
        if pair:
            if pair[0] != first[0] and pair[1] != first[1]:
                top = [gts[0], gt]
            total[pair[0]] += gt[0]
    if top:
        total = _get_total(top, total)
        return top, total
    return False, False

def _above_prop(x, s, p=0.8):
    pvals = []
    for p in [0.8, 0.9, 1.0]:
        pvals.append(stats.binom_test(x, s, p))
    return max(pvals) > 0.70

def _prop(gt):
    sense_sorted = sorted(zip(gt.values(), gt.keys()), reverse=True)
    top_2, total = _top_gt(sense_sorted)
    # print "top_2 %s totla %s" % (top_2, total)
    if top_2:
        gt2_prop = float(top_2[1][0]) / total[1]
        gt1_prop = float(top_2[0][0]) / total[0]
        table = np.array([[top_2[1][0], total[1] - top_2[1][0]], [total[0] - top_2[0][0], top_2[0][0]]])
        # print "table\n%s\ntotals %s %s" % (table, gt1_prop, gt2_prop)
        # print stats.fisher_exact(table)
        if stats.fisher_exact(table)[1] < 0.05 and _above_prop(top_2[0][0], total[0]) and _above_prop(top_2[1][0], total[1]):
            return True
    return False

def _valid_test(link, link_as):
    """
    Only if top2 associated nt are equally represented
    """
    # print "link %s %s" % (link, link_as)
    if len(link) > 1:
        sense_pval = _prop(link)
    else:
        sense_pval = False
    # if len(link_as) > 1:
    #     anti_pval = _prop(link_as)
    # else:
    #     anti_pval = True
    if sense_pval:
        return True
    return False

def _valid(link, link_as):
    """
    Only if one snp allele is associated with the Cu/Cm
    """
    if len(link) == 2:
        gen = link.keys()
        allels1 = gen[0].split(":")[0].split("/")
        allels2 = gen[1].split(":")[0].split("/")
        if allels1[0] != allels2[0] and allels1[1] != allels2[1] and _snp_veracity(link, link_as):
            return True

def _format(link):
    """
    Give nice format to dict with alleles and reads supporting
    """
    cell = ''
    for allele in link:
        cell += "%s=%s;" % (allele, link[allele])
    return cell

def _change_to_cpg(line, tag):
    return line.replace(tag, "CpG%s" % tag).strip()

def _change_to_snp(line, tag):
    return line.replace(tag, "SNP%s" % tag).strip()

def _create_vcf_header(vcf_file, out_handle):
    """
    Create header for final vcf
    """
    print >>out_handle, "##fileformat=VCFv4.1"
    print >>out_handle, "##fileData=%s" % datetime.date.today().strftime('%y%m%d')
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("##reference"):
                print >>out_handle, line.strip()
            if line.startswith("##contig"):
                print >>out_handle, line.strip()
            if line.startswith("#CHROM"):
                print >>out_handle, line.strip()
            if line.startswith("##BisSNP"):
                print >>out_handle, line.strip()
            if line.startswith("##FILTER"):
                print >>out_handle, line.strip()
            if line.startswith("##FORMAT=<ID=GT"):
                print >>out_handle, line.strip()
            if line.startswith("##INFO=<ID=DP"):
                print >>out_handle, line.strip()
            if line.startswith("##FORMAT=<ID=BRC6"):
                print >>out_handle, _change_to_cpg(line, 'BRC6')
                print >>out_handle, _change_to_snp(line, 'BRC6')
            if line.startswith("##FORMAT=<ID=CM"):
                print >>out_handle, _change_to_cpg(line, 'CM')
                print >>out_handle, _change_to_snp(line, 'CM')
            if line.startswith("##FORMAT=<ID=CU"):
                print >>out_handle, _change_to_cpg(line, 'CU')
                print >>out_handle, _change_to_snp(line, 'CU')
            if line.startswith("##FORMAT=<ID=CP"):
                print >>out_handle, _change_to_cpg(line, 'CP')
                print >>out_handle, _change_to_snp(line, 'CP')
            if line.startswith("##FORMAT=<ID=DP"):
                print >>out_handle, _change_to_cpg(line, 'DP')
                print >>out_handle, _change_to_snp(line, 'DP')
            if line.startswith("##INFO=<ID=CS"):
                print >>out_handle, line.strip()

def _get_info(info, tag):
    """
    get value from info vcf field
    """
    return next((value.split("=")[1] for value in info.split(";") if value.startswith(tag)), None)

def _get_format(header, frmt):
    """
    get format field in dict instance
    """
    frmt = dict(zip(header.split(":"), frmt.split(':')))
    return frmt

def _format_vcf_value(frmt1, frmt2, tag):
    return {_change_to_cpg(tag, tag): frmt1[tag],
            _change_to_snp(tag, tag): frmt2[tag]}

def _get_vcf_line(record):
    """
    create new vcf file with CpG and SNP information
    """
    frmt = {}
    cs = _get_info(record[7], "CS")
    ref = "%s%s" % ("C", record[13])
    alt = "%s%s" % ("C", record[14])
    qual = (float(record[5]) + float(record[15])) / 2
    filter = "LowQual"
    dp = int(_get_info(record[7], "DP")) + int(_get_info(record[17], "DP"))
    info = ";".join(["DP=%s" % dp, "CS=%s" % cs])
    cpg = _get_format(record[8], record[9])
    snp = _get_format(record[18], record[19])
    for value in ["BRC6", "CM", "CU", "CP", "DP"]:
        frmt.update(_format_vcf_value(cpg, snp, value))
    format = "GT:" + ":".join(frmt.keys())
    sample = snp["GT"] + ":" + ":".join(frmt.values())
    return record[0], record[11], ref, alt, qual, filter, info, format, sample

def _correct_vcf(vcf_file):
    """
    sort by genome/position, bgzip and index
    """
    vcf_sort = append_stem(vcf_file, "_sort") + ".gz"
    if not file_exists(vcf_sort):
        with file_transaction(vcf_sort) as tx_out:
            cmd = "cat {vcf_file} |vcf-sort | bgzip  > {tx_out}"
            do.run(cmd.format(**locals()), "sort %s" % vcf_file)
            do.run("tabix -f {0}".format(tx_out), "")
    return vcf_sort

def cpg_het_pairs(cpgvcf, snpvcf, bam_file, out_file, workdir):
    """
    Detect het close to hemi-met sites
    """
    out_vcf = splitext_plus(out_file)[0] + ".vcf"
    cpg_filter = op.join(workdir, op.basename(append_stem(cpgvcf, "_filtered")))
    snp_filter = op.join(workdir, op.basename(append_stem(snpvcf, "_filtered")))

    if not file_exists(cpg_filter):
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
    if not file_exists(snp_filter):
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

    if not file_exists(out_vcf):
        res = pybedtools.BedTool(cpg_filter).window(snp_filter, w=75)
        with open(out_file, 'w') as out_handle, open(out_vcf, 'w') as vcf_handle:
            _create_vcf_header(cpgvcf, vcf_handle)
            print >>out_handle, "chrom\tCpG_pos\tCpG_nt\tSNP_pos\tAlleles\tassociation_plus\tSNP_reads_minus"
            for record in res:
                if record[1] != record[11]:
                    # if record[1] == "19889634":
                    link, link_as, align = _make_linkage(bam_file, record[0], int(record[1]), int(record[11]), _get_strand(record))
                    res = "%s\t%s\t%s\t%s\t%s/%s\t%s\t%s" % (record[0], record[1], record[3], record[11], record[13], record[14], _format(link), _format(link_as))
                    chrom, pos, ref, alt, qual, filt, info, frmt, sample = _get_vcf_line(record)
                    # print res
                    if _valid_test(link, link_as):
                        filt = "PASS"
                        print >>out_handle, res
                        # print res
                        # print >>out_handle, '\n'.join(align)

                    vcf_res = "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{frmt}\t{sample}".format(**locals())
                    print >>vcf_handle, vcf_res
    return _correct_vcf(out_vcf)

def _complement(nt):
    if nt == 'a':
        return 't'
    elif nt == 't':
        return 'a'
    elif nt == 'c':
        return 'g'
    elif nt == 'g':
        return 'c'

def _model(pileup, snp, cpg_st):
    c_pos = v_pos = []
    for read in pileup:
        if len(pileup[read].keys()) == 1:
            continue
        info_snp = pileup[read]['snp'].split(":")
        info_cpg = pileup[read]['cpg'].split(":")
        if info_cpg[1] == cpg_st:
            if cpg_st == "+":
                c_pos.append(info_cpg[0].lower())
                v_pos.append(info_snp[0].lower())
            else:
                c_pos.append(_complement(info_cpg[0].lower()))
                v_pos.append(_complement(info_snp[0].lower()))
        else:
            if info_snp[1] == "+":
                v_pos.append(info_snp[0].lower())
            else:
                v_pos.append(_complement(info_snp[0].lower()))

def _make_linkage(bam_file, chrom, cpg, snp, cpg_st):
    start, end = [cpg-1, snp-1] if cpg-1 < snp-1 else [snp-1, cpg-1]
    pairs = _pairs_matrix(bam_file, [chrom, start, end], cpg-1, snp-1)
    link = Counter()
    link_as = Counter()
    align = []
    for pair in pairs:
        if len(pairs[pair].keys()) == 1:
            continue
        nts = [pairs[pair]['cpg'].split(":")[0], pairs[pair]['snp'].split(":")[0]]
        align.append("-".join(nts) if cpg < snp else "-".join(nts[::-1]))
        info_snp = pairs[pair]['snp'].split(":")
        # if info_snp[1] == cpg_st:
        # print pairs[pair]
        if pairs[pair]['cpg']:
            info_cpg = pairs[pair]['cpg'].split(":")
            if info_cpg[1] == info_snp[1] and info_cpg[1] == cpg_st:
                link["v%s/c%s:%s" % (info_snp[0], info_cpg[0], cpg_st)] += 1
        # else:
        #    link_as["v%s:%s" % (info_snp[0], info_snp[1])] += 1
    # print "LINK\n%s\n" % link
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

def post_processing(vcf_res, vcf_merged, out):
    """
    merge list of vcf files and get stats
    """
    if len(vcf_res) == 1:
        return vcf_res
    if not file_exists(vcf_merged):
        cmd = "bcftools merge {0} > {1}".format(" ".join(vcf_res), vcf_merged)
        do.run(cmd, "merge files")

    vcf_reader = vcf.Reader(open(vcf_merged, 'r'))
    samples = vcf_reader.samples
    num_call = Counter()
    num_call_sample = Counter()
    for record in vcf_reader:
        if not record.FILTER:
            num_call[record.num_called] += 1
            # print record.num_called

            for sample in samples:
                if record.genotype(sample)['GT'] != "./.":
                    # print record.genotype(sample)['GT']
                    num_call_sample[sample] += 1

    with open(out + "_shared_stat.tsv", 'w') as stat_handle:
        print >>stat_handle, tabulate([[k, v] for k, v in num_call.iteritems()], headers=["# samples", "# of SNPs"])
    with open(out + "_stat.tsv", 'w') as stat_handle:
        print >>stat_handle, tabulate([[k, v] for k, v in num_call_sample.iteritems()], headers=["samples", "# of SNPs"])

def detect_asm(data, args):
    vcf_res = []
    in_vcf = data['fastq']
    bam_file = data['bam']
    sample = splitext_plus(op.basename(in_vcf))[0].split(".raw")[0].replace(".rawcpg", "")
    workdir = op.join(args.out, sample)
    safe_makedir(workdir)
    snp_file = in_vcf.replace("rawcpg", "rawsnp")
    assert bam_file, "No bam file associated to vcf %s" % in_vcf
    out_file = op.join(workdir, sample + "_pairs.tsv")
    vcf_res = cpg_het_pairs(in_vcf, snp_file, bam_file, out_file, workdir)
    data['asm'] = vcf_res
    return data
