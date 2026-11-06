# StaphMix

A one-stop analysis pipeline for rapid cgMLST, virulence and resistance
profiling, and SCCmec typing of *Staphylococcus aureus*.

## Principles

1. **Species confirmation** — fastANI against a *S. aureus* reference genome.
2. **cgMLST** — allele calling for a 1,119-gene core genome scheme, plus
   HierCC cluster assignment (HC6 and HC650).
3. **Virulence genes** — DIAMOND blastx against VFDB set B.
4. **Antimicrobial resistance** — AMRFinderPlus.
5. **SCCmec typing** — BLASTn against an extended whole-cassette SCCmec
   database.

## Requirements

Python 3.8+ and the following executables on `PATH`:

| Tool | Version |
|---|---|
| makeblastdb / blastn | >= 2.14 |
| diamond | >= 2.2 |
| fastANI | >= 1.33 |
| amrfinder (AMRFinderPlus) | >= 3.11 |

Python dependencies (installed automatically): numpy, pandas, numba,
biopython, click.

## Installation

### Conda (recommended)

The code and the reference databases are distributed as two packages, because
the database is about 1.5 GB:

```bash
conda install -c <your-channel> -c conda-forge -c bioconda staphmix staphmix-db
```

Without the database you can still install the code package and supply a
database directory yourself (see below).

### From source

```bash
git clone https://github.com/Guilai-Jiang/StaphMix.git
cd StaphMix
conda env create -n staphmix -f environment.yml   # or install the tools yourself
conda activate staphmix
pip install .
staphmix --help
```

## Database

StaphMix looks for its database in this order:

1. `--db PATH`
2. `$STAPHMIX_DB`
3. `$CONDA_PREFIX/share/staphmix/db` (installed `staphmix-db` package)
4. `./db` next to the source checkout

The directory must contain:

| File | Content |
|---|---|
| `staph_aureus.cgMLST` | list of the 1,119 core genes |
| `MLSTdb0.1.fasta` | reference allele sequences for every core gene |
| `new_profile_md5.csv` | md5 allele profiles of 31,410 reference genomes |
| `new_HierCC.HierCC` | HierCC cluster levels for the same genomes |
| `GCF_000025145.1_ASM2514v1_genomic.fna` | *S. aureus* reference for fastANI |
| `VFDB_setB_pro.fas` | VFDB set B protein database |
| `whole_cassette_SCCmec_database_EXTENDED_20171117.fasta` | SCCmec cassette database |

A parsed copy of the md5 profile (`new_profile_md5.csv.pkl`) is created on the
first cgMLST run to make later runs start in seconds instead of minutes.

### Building the database package

```bash
./scripts/build_db_archive.sh db                 # creates staphmix-db-1.0.tar.gz
# upload the tarball as a GitHub release asset, then fill url/sha256 in
# conda/recipe-db/meta.yaml and build:
conda build conda/recipe-db
anaconda upload -u <your-anaconda-user> \
  /path/to/conda-bld/noarch/staphmix-db-1.0-*.tar.bz2
```

Please check the redistribution terms of the bundled databases (VFDB,
SCCmec, cgMLST/HierCC) before publishing a database package.

## Usage

```bash
staphmix -q genome.fna.gz -s -c -m -v -a -p sample -o results -n 8
```

| Option | Meaning |
|---|---|
| `-q, --query` | input genome (fasta/fastq, optionally gzipped) |
| `-s, --species` | run the ANI check (stops if not *S. aureus*) |
| `-c, --cgmlst` | run cgMLST typing and HierCC assignment |
| `-m, --sccmec` | run SCCmec typing |
| `-v, --vf` | run virulence gene comparison |
| `-a, --amr` | run antimicrobial resistance comparison |
| `-p, --prefix` | prefix for the output files |
| `-o, --outdir` | output directory (created if needed) |
| `-n, --n_thread` | number of threads (default 8) |
| `--db` | database directory (see above) |

## Output

| File | Content |
|---|---|
| `*.summary.txt` | one row with all results: HC6, HC650, SCCmec type, virulence scores, MRSA/MSSA, resistance gene counts |
| `*.ani.txt` | fastANI result |
| `*.alleles.profile` | md5 allele profile for 1,119 genes (GrapeTree compatible) |
| `*.sccmec.txt` | SCCmec BLAST hits |
| `*.vf.txt` | virulence gene hits |
| `*.amr.txt` | resistance gene hits |

## Notes

- AMRFinderPlus does not bundle its own database. Install it with
  `amrfinder -u`, or point `AMRFINDER_DB` at a downloaded release. The conda
  package works normally once the environment is activated.
- The tool versions above are the minimum tested; newer releases of BLAST,
  DIAMOND, fastANI and AMRFinderPlus are compatible.
- On a cluster with a per-user thread limit, set `OMP_NUM_THREADS=1` and
  control parallelism with `-n`.

## Development

```bash
./scripts/build_conda.sh          # build the code package from this checkout
BUILD_DB=1 ./scripts/build_conda.sh   # also build the database package
```

`scripts/build_db_archive.sh` creates the database tarball used by
`staphmix-db`. See [docs/PUBLISHING.md](docs/PUBLISHING.md) for the release
process.

## License

GNU General Public License v3.0 or later; see `COPYING`.

## Contact

Bug reports and suggestions: [gljiang@stu.suda.edu.cn](mailto:gljiang@stu.suda.edu.cn)
