# StaphMix

[![build](https://github.com/Guilai-Jiang/StaphMix/actions/workflows/build-conda.yml/badge.svg)](https://github.com/Guilai-Jiang/StaphMix/actions/workflows/build-conda.yml)
[![license](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)
![python](https://img.shields.io/badge/python-3.8%2B-blue)

Genome-to-report profiling for *Staphylococcus aureus*: species check, cgMLST
and HierCC assignment, SCCmec type, virulence genes and acquired resistance
genes in a single pass.

| Step | Method |
|---|---|
| Species check | fastANI against a *S. aureus* reference (≥ 95 % ANI) |
| cgMLST / HierCC | 1,119 core genes, md5 allele profiles, 31,410-genome reference set |
| SCCmec | BLASTn against an extended whole-cassette database |
| Virulence | DIAMOND blastx against VFDB set B |
| Resistance | AMRFinderPlus |

Runtime is roughly 4–5 minutes for a 2.8 Mb genome with 8 threads.

## Requirements

Python ≥ 3.8 plus these executables on `PATH`:

| Tool | Version |
|---|---|
| `makeblastdb`, `blastn` | ≥ 2.14 |
| `diamond` | ≥ 2.2 |
| `fastANI` | ≥ 1.33 |
| `amrfinder` | ≥ 3.11 |

## Install

From a release (no Anaconda account needed; dependencies come from
conda-forge and bioconda):

```bash
conda install -c conda-forge -c bioconda \
  https://github.com/Guilai-Jiang/StaphMix/releases/download/v0.1.0/staphmix-0.1.0-py_0.conda
```

From source:

```bash
git clone https://github.com/Guilai-Jiang/StaphMix.git
cd StaphMix
conda env create -f environment.yml
conda activate staphmix
pip install .
```

## Database

The reference database is about 1.5 GB and ships as a release asset, not in the
repository:

```bash
curl -LO https://github.com/Guilai-Jiang/StaphMix/releases/download/db-v1.0/staphmix-db-1.0.tar.gz
mkdir -p ~/staphmix-db && tar xzf staphmix-db-1.0.tar.gz -C ~/staphmix-db
```

It is located in this order: `--db` → `$STAPHMIX_DB` →
`$CONDA_PREFIX/share/staphmix/db` → `./db` next to the source checkout.

| File | Content |
|---|---|
| `staph_aureus.cgMLST` | 1,119 core gene identifiers |
| `MLSTdb0.1.fasta` | reference allele for every core gene |
| `new_profile_md5.csv` | md5 allele profiles, 31,410 genomes |
| `new_HierCC.HierCC` | HierCC levels for the same genomes |
| `GCF_000025145.1_ASM2514v1_genomic.fna` | fastANI reference genome |
| `VFDB_setB_pro.fas` | VFDB set B proteins |
| `whole_cassette_SCCmec_database_EXTENDED_20171117.fasta` | SCCmec cassettes |

The first cgMLST run also writes `new_profile_md5.csv.pkl` (233 MB); later runs
start in ~2 s instead of ~20 s.

## Usage

```bash
staphmix -q genome.fna.gz -s -c -m -v -a -p sample -o results -n 8
```

| Option | Meaning |
|---|---|
| `-q, --query` | input genome: fasta/fastq, optionally gzipped |
| `-s, --species` | ANI check; stops if the genome is not *S. aureus* |
| `-c, --cgmlst` | cgMLST calling and HierCC assignment |
| `-m, --sccmec` | SCCmec typing |
| `-v, --vf` | virulence gene comparison |
| `-a, --amr` | resistance gene comparison |
| `-p, --prefix` | output file prefix |
| `-o, --outdir` | output directory |
| `-n, --n_thread` | threads (default 8) |
| `--db` | database directory |

## Output

| File | Content |
|---|---|
| `*.summary.txt` | one row: HC6, HC650, SCCmec type, virulence scores, MRSA/MSSA, resistance counts |
| `*.ani.txt` | fastANI result |
| `*.alleles.profile` | md5 allele profile, 1,119 loci (GrapeTree compatible) |
| `*.sccmec.txt` | SCCmec BLAST hits |
| `*.vf.txt` | virulence gene hits |
| `*.amr.txt` | resistance gene hits |

## Notes

- AMRFinderPlus does not ship its database. Run `amrfinder -u` once or set
  `AMRFINDER_DB`.
- On clusters with a per-user thread limit, export `OMP_NUM_THREADS=1` and
  control parallelism with `-n`.

## Citation

> Jiang G, Guo C, Hua Y, Wu C, Li S, Wang Y, Zhao J, Ji L, Cheng Y, Zhou Z,
> Wu X, Li H. **Large-scale genomic analysis places Chinese CC398 as a
> persistent human-associated MSSA lineage apart from the dominant global
> LA-MRSA clade.** *mSystems*. 2026;11(10).
> [doi:10.1128/msystems.00621-26](https://doi.org/10.1128/msystems.00621-26)

## License

GPL-3.0-or-later. See [LICENSE](LICENSE) for the grant and [COPYING](COPYING)
for the full GPLv3 text. The reference databases are distributed under their
own terms; see [docs/DATABASE-LICENSES.md](docs/DATABASE-LICENSES.md).

## Contact

[gljiang@alu.suda.edu.cn](mailto:gljiang@alu.suda.edu.cn)
