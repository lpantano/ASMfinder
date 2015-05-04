import os.path as op
from bcbio.utils import splitext_plus, file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.provenance.do import find_cmd
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
import shutil


def _trimming(in_fastq, out_dir, sample, is_rrbs, is_directional):
    """
    Trimming reads using trim_galore
    """
    trim_galore = find_cmd("trim_galore")
    if is_rrbs:
        is_rrbs = "--rrbs"
    if is_directional:
        is_directional = ""

    cmd = "{trim_galore} {is_rrbs} {is_directional} --length 30 --quality 30 {in_fastq} -o {tx_dir}"
    if not file_exists(out_dir):
        with tx_tmpdir() as tx_dir:
            do.run(cmd.format(**locals()), "trim_galore in %s" % in_fastq)
            shutil.move(tx_dir, out_dir)
    trimming = op.join(out_dir, sample + "_trimmed.fq")
    assert op.exists(trimming), "trimming file doesn't exists"
    return trimming


def prepare(data, args):
    """
    Run trimming and fastq for each fastq file
    """
    workdir = op.abspath(safe_makedir("prepare"))
    sample = data['name']
    data['trimmed'] = _trimming(data['fastq'], op.join(workdir, sample), data['name'], args.is_rrbs, args.is_directional)
    return data