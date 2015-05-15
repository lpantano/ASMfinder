import os
import shutil
import os.path as op

# from collections import Counter
from bcbio.utils import splitext_plus, file_exists, safe_makedir, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import broad
from bcbio.bam import index

from ichwrapper import log

def _align(in_fastq, sample, workdir, genome_index, is_directional, reference, config):
    """
    align with bismark
    """
    bismark = do.find_cmd("bismark")
    num_cores = config['algorithm'].get('cores', 1)
    basename = sample
    if is_directional:
        is_directional = ""
    cmd = "{bismark} --bowtie2 -p {num_cores} -n 1 -o {tx_dir} --basename {sample} --unmapped {is_directional} {genome_index} {in_fastq}"
    out_dir = op.join(workdir, sample)
    out_bam = op.join(out_dir, sample + ".bam")

    with chdir(workdir):
        if not file_exists(out_bam):
            with tx_tmpdir() as tx_dir:
                cmd = cmd.format(**locals())
                log.logger.debug(cmd)
                do.run(cmd, "bismark in %s" % in_fastq)
                shutil.move(tx_dir, out_dir)

        broad_runner = broad.runner_from_config(config)
        # out_bam, _ = broad_runner.run_fn("picard_formatconverter", out_sam)
        names = {'rg': in_fastq, 'library': 'RRBS_LIB', 'pl': 'Illumina', 'pu': 'R1', 'sm': in_fastq, 'sample': sample}
        out_fix_bam = broad_runner.run_fn("picard_fix_rgs", out_bam, names)
        order_bam = splitext_plus(out_fix_bam)[0] + "_order.bam"
        broad_runner.run_fn("picard_reorder", out_fix_bam, reference, order_bam)
        index(order_bam, config)
    return order_bam


def create_bam(data, args):
    """
    aligner and conversion to BAM file
    """
    workdir = safe_makedir("align")
    sample = data['name']
    # workdir = op.join("align", sample)
    data['final_bam'] = _align(data['trimmed'], sample, op.abspath(workdir), args.index, args.is_directional, args.reference, data['config'])
    return data
