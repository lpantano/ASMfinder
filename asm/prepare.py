import os.path as op
from bcbio.utils import splitext_plus, file_exists, safe_makedir
from bcbio.provenance import do, find_cmd
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
import shutil


def _trimming(in_fastq, out_dir):
    """
    Trimming reads using trim_galore
    """
    trim_galore = find_cmd("trim_galore")
    cmd = "{trim_galore} {is_rrbs} {is_directional} --length 30 --quality 30 {in_fastq} -o {tx_dir}"
    if not file_exists(out_dir):
        with tx_tmpdir() as tx_dir:
            do.run(cmd.format(**locals()), "trim_galore in %s" % in_fastq)
            shutil.move(tx_dir, out_dir)
    # get trimnmed file


def prepare(data, args):
    """
    Run trimming and fastq for each fastq file
    """
    safe_makedir("prepare")
    sample = data['name']
    workdir = safe_makedir(sample)
    data['trimmed'] = _trimming(data['fastq'], op.abspath(workdir))
    return data
