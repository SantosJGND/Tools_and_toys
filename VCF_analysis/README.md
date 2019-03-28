### VCF Analysis

Jupyter notebooks to perform cluster analysis of genetic data in .vcf format. 

The tools presented in this repository represent a condensed form of the pipeline of analysis Galaxy_KDE 
(Santos *et al*. 2019 - *Awaiting review*, link to repository: [Galaxy KDE repo](https://github.com/SantosJGND/Galaxy_KDE_classifier)). 

> Galaxy_KDE: Pipeline developped to study the genetic diversity of local haplotypes across genomic data sets of homozygous or haploid data.


#### Supervised Layer

Jupyter notebooks are equiped to tackle the first layer of the Galaxy_KDE pipeline:

- Global clustering;
- KDE extraction across local genomic windows. 
- Classification.
- Summary Analysis. 

**Steps**:

>i. Perform global structure analysis using full data set. 

>ii. Sub-divide data into local windows.

>iii. Perform PCA and KDE of reference groups across windows. Extract diversity statistics. 

>iv. Display genetic differentiation and among group variance statistics across windows.

>v. Display patterns of assignment against local genomic position.


#### Reference substructure.

The notebook `Gy_pipeline` is further equiped with the second layer of the pipeline:

>i. KDE estimation of mean shift clusters across windows.

>ii. Structure analysis of p-value matrix. Clustering and profile analysis.

>iii. Classification of local haplotypes using cluster similarity. 

#### Curation

The notebook `kde_curate` uses the supervised classification of local haplotypes 
to capture local clusters and manually modify the labels of the recovered haplotypes.

Use for informed evolutionary inference of local patterns of structure.


- [Coalescence - simulation](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/VCF_analysis/Simu_17-03-2019/vcf_analysis.ipynb)

- [Structure manipulation - simulation](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/VCF_analysis/vcf_analysis/vcf_analysis.ipynb)

- [Extract Rice Data - Gy_pipeline](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/VCF_analysis/Extract/vcf_analysis.ipynb)

- [Extract Rice Data - kde_curate](https://nbviewer.jupyter.org/github/SantosJGND/Tools_and_toys/blob/master/VCF_analysis/Extract/kde_curate.ipynb)