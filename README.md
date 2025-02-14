# StaphMix
A one-stop analysis pipeline for rapid cgMLST, virulence and resistance profiling, and SCCmec typing of *Staphylococcus aureus*.

# **Principle**

1. **ANI-based simple species identification**  
   Based on ANI (Average Nucleotide Identity), a simple method for species identification.

2. **cgMLST scheme based on 1,119 core genes**  
   Uses a cgMLST scheme developed from 1,119 core genes for bacterial typing.

3. **Virulence gene prediction based on VFDB database**  
   Predicts virulence genes by utilizing the VFDB (Virulence Factor Database).

4. **Antibiotic resistance gene prediction based on AMRFinder software**  
   Uses AMRFinder software to predict antibiotic resistance genes.

5. **SCCmec typing based on SCCmecFinder software with the whole_cassette_SCCmec_database_EXTENDED**  
   Predicts SCCmec types based on the SCCmecFinder software and the extended whole_cassette_SCCmec database.

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
