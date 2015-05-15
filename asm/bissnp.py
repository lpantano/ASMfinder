import os.path as op
# from collections import Counter
from bcbio.utils import splitext_plus, file_exists, safe_makedir, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.bam import index

from ichwrapper import log


def _count_covars(in_bam, sample, workdir, snp, reference, config):
    """
    countcovars from BisSNP tool
    """
    bissnp = do.find_cmd("bissnp")
    basename = sample
    num_cores = config['algorithm'].get('cores', 1)
    memory = config['algorithm'].get('memory', 4)
    jvm_opts = "-Xms750m -Xmx%sg" % memory
    cmd = ("{bissnp} {jvm_opts} -R {reference} -I {in_bam} "
           "-T BisulfiteCountCovariates "
           "-knownSites {snp} "
           "-cov ReadGroupCovariate "
           "-cov QualityScoreCovariate "
           "-cov CycleCovariate "
           "-recalFile {tx_out} "
           "-nt {num_cores} ")
    with chdir(workdir):
        out_count = op.join(workdir, sample + "_recal1.csv")
        if not file_exists(out_count):
            with file_transaction(out_count) as tx_out:
                log.logger.debug(cmd.format(**locals()))
                do.run(cmd.format(**locals()), "BisSNP countcovarts in %s" % in_bam)
    return out_count


def _recal_BQ_score(in_bam, sample, workdir, counts_file, reference, config):
    """
    recalibration from BisSNP tool
    """
    bissnp = do.find_cmd("bissnp")
    basename = sample
    num_cores = config['algorithm'].get('cores', 1)
    memory = config['algorithm'].get('memory', 4)
    jvm_opts = "-Xms750m -Xmx%sg" % memory
    cmd = ("{bissnp} {jvm_opts} -R {reference} -I {in_bam} "
           "-T BisulfiteTableRecalibration "
           "-recalFile {counts_file} "
           "-o {tx_out} "
           "-maxQ 60 ")
    with chdir(workdir):
        out_recal = op.join(workdir, sample + "_recal1.bam")
        if not file_exists(out_recal):
            with file_transaction(out_recal) as tx_out:
                log.logger.debug(cmd.format(**locals()))
                do.run(cmd.format(**locals()), "BisSNP writerecal in %s" % in_bam)
        index(out_recal, config)
    return out_recal


def _call_vcf(in_bam, sample, workdir, reference, config):
    """
    recalibration from BisSNP tool
    """
    bissnp = do.find_cmd("bissnp")
    basename = sample
    num_cores = config['algorithm'].get('cores', 1)
    memory = config['algorithm'].get('memory', 4)
    jvm_opts = "-Xms750m -Xmx%sg" % memory
    cmd = ("{bissnp} {jvm_opts} -R {reference} -I {in_bam} "
           "-T BisulfiteGenotyper "
           "-vfn1 {tx_out} "
           "-vfn2 {out_vfn2} "
           "-stand_call_conf 20 "
           "-stand_emit_conf 0 "
           "-mmq 30 "
           "-mbq 0 "
           "-nt {num_cores}")
    with chdir(workdir):
        out_vfn1 = op.join(workdir, sample + ".rawcpg.vcf")
        out_vfn2 = op.join(workdir, sample + ".rawsnp.vcf")
        if not file_exists(out_vfn1):
            with file_transaction(out_vfn1) as tx_out:
                log.logger.debug(cmd.format(**locals()))
                do.run(cmd.format(**locals()), "BisSNP writerecal in %s" % in_bam)
    return out_vfn1, out_vfn2


def call_variations(data, args):
    """
    Run BisSNP tool
    """
    safe_makedir("bissnp")
    sample = data['name']
    workdir = op.abspath(safe_makedir(op.join("bissnp", sample)))
    counts_file = _count_covars(data['final_bam'], sample, workdir, args.snp, args.reference, data['config'])
    recal_bam = _recal_BQ_score(data['final_bam'], sample, workdir, counts_file, args.reference, data['config'])
    cpg, snp = _call_vcf(recal_bam, sample, workdir, args.reference, data['config'])
    return data
