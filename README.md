# mPackage ver 0.1.0

Temporal and functional analysis of human gut microbiome changes




### Introduction

This package will allow you to explore temporal and functional analysis of human gut-associated metagenomics species pangenome (MSP).



### Installing

This package can be installed through github repository

```
install.packages("devtools")
install_github("sysbiomelab/mPackage")
```
 
### Running example

```
#identifying functional clusters from mapping matrix
mspKoModules = makeFuncCluster(mspKoMat) 

#identifying highly covered MSPs by functional clusters
mspByModuleTab = getCoveredMspByModules(mspKoModules, mspKoMat) 

#identifying inflow/outflow probablities 
inflowStats = getInflowStats(metaTab, mgsMat)
outflowStats = getOutflowStats(metaTab, mgsMat)



```
### Citation

Saeed Shoaie, Sunjae Lee, Gholamreza Bidkhori et al., in preparation, 2020

### Authors

[Sunjae Lee](https://github.com/SunjaeLee)

[Gholamreza Bidkhori](https://scholar.google.com/citations?user=LrCrMp0AAAAJ&hl=en)

### Contacts
sunjae.lee@kcl.ac.uk
 
### Acknowledgments

*  [INRA](http://www.mgps.eu/index.php?id=accueil)
*  [SciLifeLab](https://www.scilifelab.se/)
*  [KCL](https://www.sysbiomelab.com/)
*  [KTH](https://www.kth.se/en)
*  [Karolinska Institute](https://ki.se/en)
*  [KAIST](https://www.kaist.ac.kr/en/)


