### targeted AMOVA

A script is provided to estimate and store the among group component of AMOVA 
across genomic sites. The script takes genotype data stored in format `.geno`.
Genotye file is positional and provided without flag. No limit is set on number of 
genotype files but must be labelled for chromosome as such: `chr` followed by number.


**Other input files required:**

- **.bim** and **.fam** files (plink format) indicated by flags `--bim` and `--fam` respectively.

- **ref** and **admx** (optional) files: two column file {code; ID} positional, no headers. Labels in ref file used for 
AMOVA calculations if `--supervised` is passed.

- **Coordinates** file (optional): passed using `--aims` argument. Required columns: {ID; start; end; CHR}, header= True.


**Other features**
- `--mrg`: integer. Extends windows used by number of base pairs provided (devpd for known coordinates). 
- `--randN`: integer. Number of random windows to survey across chromosomes provided. 
- `--w`: integer. Size of random windows to extract. In Snp number.
- `--amova`: No value. Calculate and store amova across windows surveyed.


**Example command**


    qsub -N amovaEasy -V -q normal.q -l mem_free=30G -b y python -u AMOVA_PCA.py \
    /gs7k1/home/jgarcia/olde_home/FILTER_GENO/SplitGENO/GENO_chr08.geno \
    --fam /gs7k1/home/jgarcia/olde_home/FILTER_GENO/F3G/NG_001.fam \
    --bim /gs7k1/home/jgarcia/olde_home/FILTER_GENO/F3G/NG_001.bim \
    --ref /gs7k1/home/jgarcia/PO_Core/refs_CORE.txt \
    --admx /gs7k1/home/jgarcia/PO_Core/admx_CORE.txt \
    --amova \
    --supervised \
    -w 150 \
    --dr PCA \
    --random \
    -w 150 \
    --randN 20 \
    --id AMOVA_gitexample \
    --out AMOVA_gitexample

