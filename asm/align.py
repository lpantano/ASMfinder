import os
import os.path as op
# from collections import Counter
from bcbio.utils import splitext_plus, file_exists, safe_makedir
from bcbio.provenance import do, find_cmd
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio import broad
from bcbio.bam import index
import shutil


def _align(in_fastq, sample, workdir, genome_index, is_directional, reference, config):
    """
    align with bismark
    """
    trim_galore = find_cmd("bismark")
    basename = sample
    cmd = "bismark -n 1 -o {tx_dir} --basename {sample} --unmapped {is_directional} {genome_index} {in_fastq}"
    out_bam = op.join(workdir, sample + ".bam")
    if not file_exists(out_bam):
        with tx_tmpdir() as tx_dir:
            do.run(cmd.format(**locals()), "bismark in %s" % in_fastq)
            shutil.move(tx_dir, workdir)

    broad_runner = broad.runner_from_config(config)
    # out_bam, _ = broad_runner.run_fn("picard_formatconverter", out_sam)
    names = {'rg': in_fastq, 'library': 'RRBS_LIB', 'pl': 'Illumina', 'pu': 'R1', 'sm': in_fastq}
    out_fix_bam = broad_runner.run_fn("picard_fix_rgs", out_bam, names)
    order_bam = splitext_plus(out_fix_bam)[0] + "_order.bam"
    broad_runner.run_fn("picard_reorder", out_fix_bam, reference, order_bam)
    index(order_bam, config)
    return order_bam


def aling(data, args):
    """
    aligner and conversion to BAM file
    """
    safe_makedir("align")
    sample = data['name']
    workdir = safe_makedir(sample)
    data['final_bam'] = _align(data['trimmed'], sample, op.abspath(workdir), args.index, args.is_directional, args.reference, data['config'])
    return data
