# StaphMix

         ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄         ▄  ▄▄       ▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄       ▄
        ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░▌     ▐░░▌▐░░░░░░░░░░░▌▐░▌     ▐░▌ 
        ▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌       ▐░▌▐░▌░▌   ▐░▐░▌ ▀▀▀▀█░█▀▀▀▀  ▐░▌   ▐░▌ 
        ▐░▌               ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌       ▐░▌▐░▌▐░▌ ▐░▌▐░▌     ▐░▌       ▐░▌ ▐░▌ 
        ▐░█▄▄▄▄▄▄▄▄▄      ▐░▌     ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌ ▐░▐░▌ ▐░▌     ▐░▌        ▐░▐░▌ 
        ▐░░░░░░░░░░░▌     ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌  ▐░▌  ▐░▌     ▐░▌         ▐░▌ 
         ▀▀▀▀▀▀▀▀▀█░▌     ▐░▌     ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░▌   ▀   ▐░▌     ▐░▌        ▐░▌░▌ 
                  ▐░▌     ▐░▌     ▐░▌       ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌     ▐░▌       ▐░▌ ▐░▌ 
         ▄▄▄▄▄▄▄▄▄█░▌     ▐░▌     ▐░▌       ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌ ▄▄▄▄█░█▄▄▄▄  ▐░▌   ▐░▌ 
        ▐░░░░░░░░░░░▌     ▐░▌     ▐░▌       ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌     ▐░▌ 
         ▀▀▀▀▀▀▀▀▀▀▀       ▀       ▀         ▀  ▀            ▀         ▀  ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀       ▀ 
A one-stop analysis pipeline for rapid cgMLST, virulence and resistance profiling, and SCCmec typing of *Staphylococcus aureus*.
(A pre-release. This release will be labeled as non-production ready)

# **Principle**

1. **ANI-based simple species identification**  

2. **cgMLST scheme based on 1,119 core genes**  

3. **Virulence gene prediction based on VFDB database**  

4. **Antibiotic resistance gene prediction based on AMRFinder software**  

5. **SCCmec typing based on SCCmecFinder software with the whole cassette SCCmec database EXTENDED**  

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
# **Results**

1. **alleles.profile**  
   cgMLST typing generates the allele profile, which can later be used with [Grapetree](https://achtman-lab.github.io/GrapeTree/MSTree_holder.html)
 software to construct phylogenetic trees or for visualization.

2. **amr.txt**  
   Antibiotic resistance gene comparison results.

3. **ani.txt**  
   Species similarity comparison results.

4. **sccmec.txt**  
   SCCmec type comparison results.

5. **vf.txt**  
   Virulence gene comparison results.

6. **summary.txt**  
   A summary of all results.

## Contact

If you find any bugs or have better suggestions, feel free to contact [gljiang@stu.suda.edu.cn](mailto:gljiang@stu.suda.edu.cn).
