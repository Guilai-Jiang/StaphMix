# StaphMix
A one-stop analysis pipeline for rapid cgMLST, virulence and resistance profiling, and SCCmec typing of *Staphylococcus aureus*.

## Required Tools

The following tools must be installed in advance:

- **makeblastdb** (Developer version 2.14.0+)
- **blastn** (Developer version 2.14.0+)
- **diamond** (Developer version 2.1.8.162)
- **fastANI** (Developer version 1.33)
- **amrfinder** (Developer version 3.12.8)

 Additionally, **Python 3.6+** is required.

## Usage

Run the following command to see the options:

```bash
python staphmix.py --help
Options:
-q, --query TEXT
    Path to the input genome file. Supports fasta(.gz)/fastq(.gz) file format. [required]
-s, --species
    Flag to run species identification.
-c, --cgmlst
    Flag to run cgMLST typing. The species must be Staphylococcus aureus.
-m, --sccmec
    Flag to run SCCmec typing.
-v, --vf
    Flag to run virulence genes comparison.
-a, --amr
    Flag to run antibiotic resistance genes comparison.
-p, --prefix TEXT
    Prefix for output files. [required]
-o, --outdir TEXT
    Output directory for results. [required]
-n, --n_thread INTEGER
    Number of threads to use.
--help
    Show this message and exit.
```

# StaphMix Analysis Example Output

## Current location:
`/titan/guilai/staph_aureus/staph_cgmlst/cc398_hc650-135/seq`

## Database directory:
`/titan/guilai/software/StaphMix/db`

## The working directory does not exist. Created directory:
`/titan/guilai/staph_aureus/staph_cgmlst/cc398_hc650-135/seq/test`

---

### 2025-02-14 06:24:24 Running ANI identification of *Staphylococcus aureus* species
- **Staphylococcus aureus** : True  
- Reference strain: **GCF000025145**  
- **ANI value** : 97.6805

---

### 2025-02-14 06:24:25 Running cgMLST typing and hierarchical clustering
- **Reference genome** : GCA_000637155
- **Allelic distance** : 151
- **Hierarchical cluster level** - HC650 : 135
- **Hierarchical cluster level** - HC6 : ND

---

### 2025-02-14 06:25:07 Running Staphylococcal chromosomal cassettes mec (SCCmec) typing
- **SCCmec type** : SCCmec_type_IV(2B)

---

### 2025-02-14 06:25:09 Running virulence genes comparison
- **Enterotoxin genes score**: 0
- **Blood infection genes score**: 0

---

### 2025-02-14 06:25:13 Running antibiotic resistance genes comparison
- **This strain is MRSA**
- **Antibiotic resistance genes score**: 1
- **Antibiotic resistance genes class**: 5
- **Antibiotic resistance genes number**: 10


