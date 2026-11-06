# Publishing StaphMix

StaphMix is distributed as two conda packages:

- `staphmix` — the code (a few tens of KB), built from this repository.
- `staphmix-db` — the reference databases (~1.5 GB uncompressed), built from a
  release tarball because Git repositories should not carry files that large.

## 1. Push the code

```bash
git add .
git commit -m "Package StaphMix"
git push origin main
git tag v0.1.0 && git push origin v0.1.0
```

The GitHub Actions workflow in `.github/workflows/build-conda.yml` builds the
`staphmix` package for every pull request and publishes it to anaconda.org for
tags (requires the `ANACONDA_API_TOKEN` secret and the `ANACONDA_USER`
repository variable).

## 2. Create and publish the database tarball

```bash
./scripts/build_db_archive.sh db
# staphmix-db-1.0.tar.gz + sha256
```

Create a GitHub release tagged `db-v1.0` and attach the tarball. Copy the
printed sha256 into `conda/recipe-db/meta.yaml`. The default URL in that recipe
already points at:

```
https://github.com/Guilai-Jiang/StaphMix/releases/download/db-v1.0/staphmix-db-1.0.tar.gz
```

GitHub allows release assets up to 2 GB, so the archive fits. Zenodo, OSF or
an institutional server work equally well — only the URL matters.

## 3. Build and upload the packages

```bash
# code package: the script stages a copy without db/ (the recipe's local
# source path would otherwise copy the whole 1.5 GB database directory)
./scripts/build_conda.sh

conda build conda/recipe-db --output-folder dist   # database

anaconda upload -u <your-anaconda-user> dist/conda/noarch/staphmix-*.conda
anaconda upload -u <your-anaconda-user> dist/noarch/staphmix-db-*.conda
```

The database package is large; anaconda.org accepts it, but uploads can take a
while. If you prefer, keep `staphmix-db` on GitHub releases only and document
the `--db` option instead.

## 4. Install as a user

```bash
conda install -c <your-anaconda-user> -c conda-forge -c bioconda staphmix staphmix-db
staphmix --version
```

## Checklist

- [ ] `db/` and `db.zip` are not committed (see `.gitignore`)
- [ ] `scripts/build_db_archive.sh` only packs the seven files used at runtime
- [ ] the sha256 in `conda/recipe-db/meta.yaml` matches the uploaded tarball
- [ ] redistribution terms of VFDB, the SCCmec database and the cgMLST/HierCC
      data allow bundling, or the database is downloaded by the user instead
- [ ] `staphmix --help`, `staphmix --version` work in a clean environment
- [ ] a test genome runs end to end
