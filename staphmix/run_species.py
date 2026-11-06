# SPDX-License-Identifier: GPL-3.0-or-later
import os
import subprocess

from .configure import db_path, exe


def species_identify(query_genome, prefix, outdir_abpath, threads):
    """Compare the query with the S. aureus reference and report the ANI."""
    name = prefix + '.ani.txt'
    outfile = os.path.join(outdir_abpath, name)
    cmd = [exe['fastANI'],
           '-q', query_genome,
           '-r', db_path('species_ref'),
           '-o', outfile,
           '-t', str(threads)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError('fastANI failed:\n{0}'.format(proc.stderr.strip()))

    if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:
        return [False, '<75%']
    with open(outfile, 'r') as ani_file:
        for line in ani_file:
            part = line.strip().split('\t')
            ani = part[2]
            return [float(ani) >= 95, ani]
    return [False, '<75%']
