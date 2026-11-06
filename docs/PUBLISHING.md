# Publishing StaphMix

StaphMix can be distributed entirely through GitHub, without an Anaconda
account:

- `staphmix` — the code package (a few tens of KB), built by CI and attached to
  each release.
- `staphmix-db` — the reference databases (~1.5 GB uncompressed), distributed
  as a release tarball because Git repositories should not carry files that
  large. A conda package is optional.

## 1. Push the code and tag a release

```bash
git add .
git commit -m "Package StaphMix"
git push origin StaphMix

git tag v0.1.0
git push origin v0.1.0
```

The workflow in `.github/workflows/build-conda.yml` runs on every pull request
and, for tags, builds `staphmix` and attaches the `.conda` file to the GitHub
release. No secret needs to be configured: it uses the default
`GITHUB_TOKEN` with `contents: write`.

After the run finishes, the release contains:

```
staphmix-0.1.0-py_0.conda
```

Users install it directly by URL:

```bash
conda install -c conda-forge -c bioconda \
  https://github.com/Guilai-Jiang/StaphMix/releases/download/v0.1.0/staphmix-0.1.0-py_0.conda
```

## 2. Create and publish the database

```bash
./scripts/build_db_archive.sh db
# staphmix-db-1.0.tar.gz + sha256
```

Create a release tagged `db-v1.0` and attach the tarball. GitHub allows release
assets up to 2 GB, so the archive fits. Zenodo, OSF or an institutional server
work equally well.

Users then run:

```bash
curl -LO https://github.com/Guilai-Jiang/StaphMix/releases/download/db-v1.0/staphmix-db-1.0.tar.gz
mkdir -p ~/staphmix-db && tar xzf staphmix-db-1.0.tar.gz -C ~/staphmix-db
export STAPHMIX_DB=~/staphmix-db
```

### Optional: database as a conda package

Put the tarball URL and sha256 into `conda/recipe-db/meta.yaml`, then:

```bash
conda build conda/recipe-db --output-folder dist
```

Attach `dist/noarch/staphmix-db-1.0-0.conda` to the `db-v1.0` release as well.
It can be installed by URL like the code package, which puts the files in
`$CONDA_PREFIX/share/staphmix/db` automatically.

## 3. Optional: anaconda.org

An anaconda.org account is not required. If you want the shorter
`conda install -c <user> staphmix` form, upload the built packages:

```bash
anaconda upload -u <anaconda-user> dist/noarch/staphmix-*.conda
anaconda upload -u <anaconda-user> dist/noarch/staphmix-db-*.conda
```

## Checklist

- [ ] `LICENSE` (grant) and `COPYING` (GPLv3 text) are both present and the
      version in the metadata is `GPL-3.0-or-later`
- [ ] `db/` and `db.zip` are not committed (see `.gitignore`)
- [ ] `scripts/build_db_archive.sh` only packs the seven files used at runtime
- [ ] the database tarball is attached to a `db-*` release
- [ ] the tag release contains `staphmix-*.conda`
- [ ] redistribution terms of VFDB, the SCCmec database and the cgMLST/HierCC
      data allow bundling — otherwise distribute only the code and let users
      supply `--db`
- [ ] `staphmix --help`, `staphmix --version` work in a clean environment
- [ ] a test genome runs end to end
