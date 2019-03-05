# Tools for Genetic and Genomic analysis

**Disclaimer**: Analyses and algorithms presented below deal with high quality haploid or phased data.

Genomics data sets that combine dense marker data with breadth of sampling are becoming availble for an increasing 
number of species. Be it for historical inference, to study recombination, selection or to attempt genomic prediction,
geneticists face in genomic data sets two main problems: (i) the more straighforward issue of size, i.e. computational 
load and speed, (ii) the more abstract problem of harnessing statistics descriptive of local genomic structure to produce 
informative outputs. 

We need to work with measurements that are light and that we understand, that are informative of patterns in fundamental
properties of the data. And we need to deploy these statistics and relate the distributions produced to that of genomic 
position and/or other layers of information.

This repository is structured to reflect these two steps. The first section, **Genetic Analysis**, is devoted to 
descriptions of variation in single data sets. Alongside a conversion tool, functions are introduced that combine powerful
python modules and dimensionality reduction for added speed in calculating correlation, using the fixation index *Fst*, and 
among population variance, using AMOVA. 

The second section, **Genomic Analysis**, focuses on conducting analyses across data sets with different characteristics.
In *Simulating Genomes: direct manipulation of known-truth*, the simulation procedure used to study the behaviour of the 
statistics introduced in section I (see [Generating Samples](https://nbviewer.jupyter.org/github/SantosJGND/Genetic-data-analysis/blob/master/Notebooks/1.%20Generating_haplotypes.ipynb) for details) is developed to allow for more complex 
patterns of data structure and admixture. The simulator is modular, to allow for easy break down of component functions and 
facilitate development by the user. Under **MS_target** are introduced functions and protocol to use mean shift, a 
density-based clustering algorithm, to characterise and visualize specific patterns, and to use this 
information in novel data sets.


## Genetic Analysis

- **Geno to Darwin** - *a script to extract genotypic data to DataMatrix format for dissimilarity analysis and representation using the software DARwin.* 

>plink and .geno format (combined). 

>[README](/DARwin_kde)

- **PCA and Fst** - *calculate Fst from euclidean distances in PCA space. Use for prediction*
 
>link to tutorial: [PCA and Fst](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/PCA_Fst_estimation/Euclidian_to_fst.ipynb)

>DataMatrix to Fst: [Notebook tool](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/PCA_Fst_estimation/DataMatrix_Fst.ipynb)

- **AMOVA calculation** - *calculate among group variance using different metrics*

>link to tutorial: [AMOVA_tutorial](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/AMOVA/AMOVA_tutorial.ipynb)

>Working Example [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/AMOVA/AMOVA_mRNA_supervised/AMOVA_output.ipynb)

>command line [README](AMOVA/)


## Genomic Analysis

- **Simulating Genomes: direct manipulation of known-truth** - *control pairwise Fst across data sets, allow for recombination and admixture*

>link to tutorial: [Genome simulator](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Simulate_genomes/Genomic%20structure%20Simulator.ipynb)

- **MeanShift for structure analysis** - *use meanshift to relate clusters across data sets. Adapted to the study of local genomic variation.*

>link to tutorial: [MStarget_tutorial](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/MS_target/MStarget_tutorial.ipynb)

- **Cluster allele frequencies** - use mean shift to control for admixture in calculating local allele frequencies.

>cluster frequency analysis: [data view notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/cluster_freqs.ipynb)

>euclidean distances using observed frequency vectors: [data view notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Euclidian_to_fst.ipynb)

>example application - simulation and kde study: [simulation notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/custom_manipulation.ipynb)

>input adaptation: [structure simulator](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Structure_simulator.ipynb)

> [README](/Cluster_shape)

- **From vcf** - analysis of global and local variation, supervised assignment, use for simulation.

>PCA, global and window based. Supervised classification: [vcf_analysis](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Simu_03-03-2019/vcf_analysis.ipynb)

>Frequency feed: [simulation basic](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Simu_03-03-2019/custom_manipulation.ipynb)