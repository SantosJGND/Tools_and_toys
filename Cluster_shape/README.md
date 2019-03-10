### Allele frequency distribution

Estimate allele frequency from real data. 

Population genetic simulation usually relies on the Beta distribution as a model 
of allele frequencies of stationary populations. Realistically it is rarelly the 
case that the conditions required to achieve this state are met in nature. 
One alternative is to model frequencies as a function of prior knowledge of 
the history of the populations studied.

This repository proposes to extract the allele frequencies of local genomic windows from real data.

The usual problem with this approach is the degree of admixture in natural populations,
casting doubt on allele frequency estimates calculated along the genome, which might
inadvertently capture local genetic exchanges, biasing estimations.

We turn to mean shift, a density based clustering algorithm, to circumvent this issue.

The script `cluster_frequencies.py` was created to sample 
allele frequencies across genomic data sets in `.geno` format. Genotype files must 
contain the string `_chr` followed by the respective chromosome number. Accompanying sample and
marker information is required in the form `.fam` and `.bim` files respectively 
([plink](http://zzz.bwh.harvard.edu/plink/) format). A reference file is required, 
indicating samples to use and their respective population labels ([ex. ref file](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_test/refs_sim.txt)).
A file of non-reference accessions can also be provided using the flag `--admx` (ref file format). 
These will not be used for supervised analyses, but will be included in KDE estimates. 

The argument `--aims` takes a file of genomic regions to survey. The file must contain
columns **CHR**, **start**, **end** and **ID**, but the function that performs the reading
can be passed alternative column names as arguments. A margin, in bp, upstream and 
downstream of regions in file can be set using the argument `--mrg` (defaults to 1000).

If the `--aims` argument is not given windows will be sampled randomly across the genotype files provided. 
Window size is determined by `-w` (defaults to 150) in SNPs, the number of windows by the argument 
`--randN` (default 1000). 

If the argument `--amova` is passed, the among population variance component is calculated
across genomic windows. Defaults to calculations using euclidean distances in PCA feature space.
Alternatively, the user can specify using jaccard or hamming distances on haplotypes directly.
Metric option is set using the argument `--dist_metric`, possible arguments 'euclidean, 'hamming', 'jaccard'.
In theory takes any string afforded by the argument `metric` of the function 
[scipy.spatial.distance.pdist](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html).

If the argument `--supervised` is passed, then AMOVA, frequency and KDE estimates will 
be calculated at each window using reference labels. Otherwise, analyses are 
calculated on locally estimated mean shift clusters. 

An example job: [command](CLshape_command.txt)

Summary analysis and application:

>cluster frequency analysis: [data view notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Notebooks/cluster_freqs.ipynb)

>euclidean distances using observed frequency vectors: [data view notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Notebooks/Euclidian_to_fst.ipynb)

>example application - simulation and kde study: [simulation notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Notebooks/custom_manipulation.ipynb)

>input adaptation: [structure simulator](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Notebooks/Structure_simulator.ipynb)
