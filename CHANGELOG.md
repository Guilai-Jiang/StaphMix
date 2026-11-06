# Changelog

All notable changes to this project are documented here.

## [0.1.0] - 2026-11-06

### Added

- Species check (fastANI), cgMLST and HierCC assignment, SCCmec typing,
  VFDB virulence genes and AMRFinderPlus resistance genes.
- Database auto-discovery through `--db`, `$STAPHMIX_DB`, the conda share
  directory or `./db`.
- Pickle cache for the md5 profile table (about 20 s down to 2 s per run).
- Conda package built by CI and attached to each release; `staphmix-db`
  recipe for the reference data.
- Unit tests, CITATION.cff and licence notices.

### Fixed

- fastq input was converted to invalid FASTA (missing `>`).
- SCCmec output header did not match the data columns.
- The `alleles.profile` query column contained the output directory path.
- Dead branches in the SCCmec and VFDB overlap filters.
