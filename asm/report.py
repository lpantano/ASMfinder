"""
Some reports for the analysis
"""

import shutil
import os.path as op

from bcbio.utils import splitext_plus, file_exists, safe_makedir, chdir
from bcbio.provenance import do
from bcbio.provenance.do import find_cmd
from bcbio.distributed.transaction import file_transaction, tx_tmpdir


def _run_meth_extractor(bam_in, sample, workdir, config):
    """
    Run bismark_methylation_extractor command
    """
    bismark = do.find_cmd("bismark_methylation_extractor")
    cores = config['algorithm'].get('cores', 1)
    memory = config['algorithm'].get('mem', 5)
    cmd = "{bismark}  --no_overlap --comprehensive --multicore {cores} --buffer_size {memory}G --bedGraph --counts --gzip {bam_in}"
    out_dir = op.join(workdir, sample)
    mbias_file = op.join(out_dir, op.basename(splitext_plus(bam_in)[0]) + '.M-bias.txt')
    if not file_exists(mbias_file):
        with tx_tmpdir() as tx_dir:
            with chdir(tx_dir):
                do.run(cmd.format(**locals()), "bismark_methylation_extractor  in %s" % bam_in)
                shutil.move(tx_dir, out_dir)
    assert op.exists(mbias_file), "mbias report doesn't exists:%s" % mbias_file
    return mbias_file


def _run_report(bam_in, sample, biasm_file, workdir, config):
    """
    Run bismark2report command
    """
    bismark = do.find_cmd("bismark2report")
    bam_report = op.join(op.dirname(bam_in), sample) + '_SE_report.txt'
    cmd = "{bismark} --alignment_report  {bam_report} -o {tx_out} --mbias_report {biasm_file}"
    out_dir = op.join(workdir, sample)
    out_file = op.join(out_dir, sample + '.html')
    with chdir(out_dir):
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out:
                do.run(cmd.format(**locals()), "bismarkr2report  in %s" % bam_in)

    return out_dir


def create_report(data, args):
    workdir = op.abspath(safe_makedir('report'))
    sample = data['name']
    config = data['config']

    biasm_file = _run_meth_extractor(data['final_bam'], sample, workdir, config)
    report = _run_report(data['final_bam'], sample, biasm_file, workdir, config)

    return report
