# Database licensing

The StaphMix **code** is licensed under GPL-3.0-or-later (see `LICENSE` and
`COPYING`). That licence does not extend to the reference data, which are
separate works with their own terms. If you redistribute the database bundle
(for example as the `staphmix-db` package or a release asset), confirm the
terms below first.

| File | Origin | When publishing |
|---|---|---|
| `staph_aureus.cgMLST`, `MLSTdb0.1.fasta`, `new_profile_md5.csv`, `new_HierCC.HierCC` | Core-genome scheme, allele reference and HierCC profiles assembled for the CC398 study (Jiang et al., *mSystems* 2026). HierCC method: Zhou et al., *Genome Research* 2021. | Cite the study and the HierCC paper. Confirm with the corresponding author and the journal's data policy that redistribution is allowed. |
| `GCF_000025145.1_ASM2514v1_genomic.fna` | NCBI RefSeq assembly. | Public sequence data; attribution to NCBI/GenBank expected. |
| `VFDB_setB_pro.fas` | VFDB set B (Liu et al., *Nucleic Acids Research* 2022). | VFDB is provided for academic use; redistribution may require permission. Contact the VFDB maintainers before bundling. |
| `whole_cassette_SCCmec_database_EXTENDED_20171117.fasta` | Extended SCCmec cassette database used by SCCmecFinder (Kaya et al., *mSphere* 2018). | No explicit licence was found in the upstream repository at the time of writing. Obtain permission from the authors before redistributing. |
| AMRFinderPlus database | Downloaded separately by `amrfinder -u`. | Not bundled with StaphMix; governed by NCBI terms. |

## If redistribution is not permitted

Publish the code only. Users can build the database directory themselves and
point `--db` at it. The required file names are listed in the README; the code
does not care where the files came from as long as the names and formats match.

## Citation

- Jiang G, et al. Large-scale genomic analysis places Chinese CC398 as a
  persistent human-associated MSSA lineage apart from the dominant global
  LA-MRSA clade. *mSystems*. 2026;11(10).
  doi:10.1128/msystems.00621-26
- Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification
  scheme for bacterial virulence factors. *Nucleic Acids Res*.
  2022;50(D1):D912–D917.
- Kaya H, et al. SCCmecFinder, a web-based tool for typing of staphylococcal
  cassette chromosome mec in Staphylococcus aureus using whole-genome sequence
  data. *mSphere*. 2018;3(1):e00612-17.
- Zhou Z, Charlesworth J, Achtman M. HierCC: a multi-level clustering scheme
  for assignment and population genetics. *Genome Research*. 2021.
