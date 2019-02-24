### Allele frequency distribution

Estimate allele frequency from real data. 

Population genetic simulation usually relies on the Beta distribution as a model 
of allele frequencies of stationary populations. Realistically it is rarelly the 
case that natural populations are sufficiently ideal that genetic structure is 
stationary. The alternative commonly resorted to is to model frequencies as a function
of prior knowledge on the history of the populations studied.

This repository proposes to model allele frequency distributions from real data.
The usual problem with this approach is the degree of admixture in natural populations,
casting doubt on allele frequency estimates calculated along the genome, which might
inadvertently capture local genetic exchanges, biasing estimations.

We turn to mean shift, a density based clustering algorithm, to circumvent this issue.

By calculating allele frequencies of local clusters, we ensure that the estimates 
are free from admixture. 

>example application - simulation: [simulation notebook](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/custom_manipulation.ipynb)
