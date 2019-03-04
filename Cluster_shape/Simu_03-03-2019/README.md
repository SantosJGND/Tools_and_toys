## .VCF analysis

This repository presents jupyter notebooks for the local analysis of .vcf files.

Our test example is data set of data simulated using the software ms (R version).

Different branches of analysis are proposed, based on the tools developped in the 
[Tools shop](https://github.com/SantosJGND/Tools_and_toys) repository.

### I. Structure and association analysis.

Analyse the global structure in the data set provided. Perform clustering to
derive reference population labels. Then, split the data into consecutive windows 
of a given size. Analyse genetic structure at a subset of local genomic windows. 

Finally, use supervised labels to conduct supervised classification across local
genomic windows. The classification procedure is based on the KDE of reference labels.
Allows for comparison threshold to identify conflicting assignments and outlier 
threshold to identify differentiation.

>analysis notebook: [vcf analysis](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Simu_03-03-2019/vcf_analysis.ipynb)

### Simulation input.

Simulate data sets using allele frequencies calculated from data stored in a .vcf file.
Begin by performing global structure and clustering to identify ref pops. Then,
subdivide the data into windows. At each window, calculate allele frequencies either 
using the supervised labels, calculated using the result of the global clustering, or
from Meanshfit clusters at each window. 

Perform an analysis of the allele frequencies extracted, and of the principal component
analysis projections of samples across local windows.

Finally, use the extracted allele frequencies to perform simulations. Sort pairs of allele 
frequency vectors by Fst and chose when to sample from each.

>analysis notebook: [custom manipulation](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/Cluster_shape/Simu_03-03-2019/custom_manipulation.ipynb)