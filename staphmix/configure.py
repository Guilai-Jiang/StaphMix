# SPDX-License-Identifier: GPL-3.0-or-later
import os
import sys
from shutil import which

import numpy as np

# When StaphMix lives in a conda environment, its companion tools (blast,
# diamond, fastANI, amrfinder) are usually in the same environment. Make them
# discoverable even when the environment was not activated.
_env_bin = os.path.join(sys.prefix, 'bin')
if os.path.isdir(_env_bin) and _env_bin not in os.environ.get('PATH', '').split(os.pathsep):
    os.environ['PATH'] = _env_bin + os.pathsep + os.environ.get('PATH', '')

TOOLS = ('makeblastdb', 'blastn', 'diamond', 'fastANI', 'amrfinder')
exe = {tool: which(tool) for tool in TOOLS}

# File names expected inside the database directory. Callers refer to the keys,
# so a file can be renamed here without touching the rest of the code.
DB_FILES = dict(
    refseqs='MLSTdb0.1.fasta',
    core_genes='staph_aureus.cgMLST',
    repr='new_profile_md5.csv',
    hcc='new_HierCC.HierCC',
    species_ref='GCF_000025145.1_ASM2514v1_genomic.fna',
    vfdb='VFDB_setB_pro.fas',
    scc_type_db='whole_cassette_SCCmec_database_EXTENDED_20171117.fasta',
)

scheme_info = dict(min_iden=0.65, min_frag=0.6, max_iden=0.95)

db_folder = None
db_list = {}


def _db_candidates():
    """Places where a database directory is looked for, in priority order."""
    env = os.environ.get('STAPHMIX_DB')
    if env:
        yield os.path.expanduser(env)
    here = os.path.dirname(os.path.abspath(__file__))
    # conda package layout: $CONDA_PREFIX/share/staphmix/db
    yield os.path.join(sys.prefix, 'share', 'staphmix', 'db')
    # source checkout: <repo>/db (package sits next to it)
    yield os.path.join(os.path.dirname(here), 'db')
    yield os.path.join(here, 'db')


def is_db_dir(path):
    return bool(path) and os.path.isfile(os.path.join(path, DB_FILES['core_genes']))


def missing_db_files(folder):
    return [name for name in DB_FILES.values()
            if not os.path.isfile(os.path.join(folder, name))]


def default_db_folder():
    for path in _db_candidates():
        if is_db_dir(path):
            return os.path.abspath(path)
    return None


def set_db_folder(path):
    global db_folder, db_list
    folder = os.path.abspath(os.path.expanduser(path))
    if not is_db_dir(folder):
        raise ValueError(
            "Not a StaphMix database directory: {0}\n"
            "It must contain {1}.".format(folder, DB_FILES['core_genes']))
    db_folder = folder
    db_list = {key: os.path.join(folder, name) for key, name in DB_FILES.items()}
    return db_folder


def db_path(key):
    if not db_list:
        folder = default_db_folder()
        if folder is None:
            raise RuntimeError(
                "StaphMix database not found. Point --db or $STAPHMIX_DB at a "
                "directory containing {0}.".format(DB_FILES['core_genes']))
        set_db_folder(folder)
    return db_list[key]


def check_tools():
    missing = [tool for tool, path in exe.items() if not path]
    if missing:
        raise RuntimeError(
            "Missing required executable(s): {0}. Install them and make sure "
            "they are on PATH.".format(', '.join(missing)))

blosum62 = np.array(
    [4., -2., 0., -2., -1., -2., 0., -2., -1., 0., -1., -1., -1., -2., 0., -1., -1., -1., 1., 0., -4., 0.,
     -3., 0., -2., -1., 0., 0., 0., 0., 0., 0., -2., 4., -3., 4., 1., -3., -1., 0., -3., 0., 0., -4.,
     -3., 3., 0., -2., 0., -1., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., 0., -3.,
     9., -3., -4., -2., -3., -3., -1., 0., -3., -1., -1., -3., 0., -3., -3., -3., -1., -1., -4., -1., -2., -2.,
     -2., -3., 0., 0., 0., 0., 0., 0., -2., 4., -3., 6., 2., -3., -1., -1., -3., 0., -1., -4., -3., 1.,
     0., -1., 0., -2., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., -1., 1., -4., 2.,
     5., -3., -2., 0., -3., 0., 1., -3., -2., 0., 0., -1., 2., 0., 0., -1., -4., -2., -3., -1., -2., 4.,
     0., 0., 0., 0., 0., 0., -2., -3., -2., -3., -3., 6., -3., -1., 0., 0., -3., 0., 0., -3., 0., -4.,
     -3., -3., -2., -2., -4., -1., 1., -1., 3., -3., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -2., -3.,
     6., -2., -4., 0., -2., -4., -3., 0., 0., -2., -2., -2., 0., -2., -4., -3., -2., -1., -3., -2., 0., 0.,
     0., 0., 0., 0., -2., 0., -3., -1., 0., -1., -2., 8., -3., 0., -1., -3., -2., 1., 0., -2., 0., 0.,
     -1., -2., -4., -3., -2., -1., 2., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -3., 0., -4., -3.,
     4., 0., -3., 2., 1., -3., 0., -3., -3., -3., -2., -1., -4., 3., -3., -1., -1., -3., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., -3., -1., 1., -3., -2., -1., -3., 0.,
     5., -2., -1., 0., 0., -1., 1., 2., 0., -1., -4., -2., -3., -1., -2., 1., 0., 0., 0., 0., 0., 0.,
     -1., -4., -1., -4., -3., 0., -4., -3., 2., 0., -2., 4., 2., -3., 0., -3., -2., -2., -2., -1., -4., 1.,
     -2., -1., -1., -3., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -2., 0., -3., -2., 1., 0., -1., 2.,
     5., -2., 0., -2., 0., -1., -1., -1., -4., 1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., 3.,
     -3., 1., 0., -3., 0., 1., -3., 0., 0., -3., -2., 6., 0., -2., 0., 0., 1., 0., -4., -3., -4., -1.,
     -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., -2., -3., -1.,
     -1., -4., -2., -2., -3., 0., -1., -3., -2., -2., 0., 7., -1., -2., -1., -1., -4., -2., -4., -2., -3., -1.,
     0., 0., 0., 0., 0., 0., -1., 0., -3., 0., 2., -3., -2., 0., -3., 0., 1., -2., 0., 0., 0., -1.,
     5., 1., 0., -1., -4., -2., -2., -1., -1., 3., 0., 0., 0., 0., 0., 0., -1., -1., -3., -2., 0., -3.,
     -2., 0., -3., 0., 2., -2., -1., 0., 0., -2., 1., 5., -1., -1., -4., -3., -3., -1., -2., 0., 0., 0.,
     0., 0., 0., 0., 1., 0., -1., 0., 0., -2., 0., -1., -2., 0., 0., -2., -1., 1., 0., -1., 0., -1.,
     4., 1., -4., -2., -3., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., -1., -1., -1., -1., -2., -2., -2.,
     -1., 0., -1., -1., -1., 0., 0., -1., -1., -1., 1., 5., -4., 0., -2., 0., -2., -1., 0., 0., 0., 0.,
     0., 0., -4., -4., -4., -4., -4., -4., -4., -4., -4., 0., -4., -4., -4., -4., 0., -4., -4., -4., -4., -4.,
     1., -4., -4., -4., -4., -4., 0., 0., 0., 0., 0., 0., 0., -3., -1., -3., -2., -1., -3., -3., 3., 0.,
     -2., 1., 1., -3., 0., -2., -2., -3., -2., 0., -4., 4., -3., -1., -1., -2., 0., 0., 0., 0., 0., 0.,
     -3., -4., -2., -4., -3., 1., -2., -2., -3., 0., -3., -2., -1., -4., 0., -4., -2., -3., -3., -2., -4., -3.,
     11., -2., 2., -3., 0., 0., 0., 0., 0., 0., 0., -1., -2., -1., -1., -1., -1., -1., -1., 0., -1., -1.,
     -1., -1., 0., -2., -1., -1., 0., 0., -4., -1., -2., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., -3.,
     -2., -3., -2., 3., -3., 2., -1., 0., -2., -1., -1., -2., 0., -3., -1., -2., -2., -2., -4., -1., 2., -1.,
     7., -2., 0., 0., 0., 0., 0., 0., -1., 1., -3., 1., 4., -3., -2., 0., -3., 0., 1., -3., -1., 0.,
     0., -1., 3., 0., 0., -1., -4., -2., -3., -1., -2., 4., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])