# microeco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128353984-942c7dad-ccc6-4e5b-8672-8325d3d576f8.png" width=150 align="right" ></a>

An R package for data mining in microbial community ecology

[![CRAN](https://www.r-pkg.org/badges/version/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
![](https://img.shields.io/badge/Release-0.6.5-orange.svg) ![](https://img.shields.io/badge/Test-0.6.6-red.svg)

## Background
In microbial community ecology, with the development of high-throughput sequencing techniques,
the increasing data amount and complexity make the community data analysis and management a challenge.
There has been a lot of R packages created for the microbiome profiling analysis.
However, it is still difficult to perform data mining fast and efficiently.
Therefore, we created R microeco package.

## Main Features
  + R6 Class to store and analyze data; fast, flexible and modularized
  + Taxonomic abundance analysis
  + Venn diagram
  + Alpha diversity
  + Beta diversity
  + Differential abundance analysis
  + Indicator species analysis
  + Environmental data analysis
  + Null model analysis
  + Network analysis
  + Functional analysis


## Install R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)

Put R in the computer env PATH, for example your_directory\R-4.0.0\bin\x64 

Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

## Install microeco

Install microeco package from CRAN directly.

```r
install.packages("microeco")
```

Or install the latest development version from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco
devtools::install_github("ChiLiubio/microeco")
```


## Tutorial
See the detailed package **tutorial** (https://chiliubio.github.io/microeco_tutorial/) and the help documentations (e.g. ?microtable).
Contructing the basic microtable object from other tools/platforms (e.g. QIIME, QIIME2, HUMAnN, Kraken2 and phyloseq) 
can be easily achieved with the package file2meco (https://github.com/ChiLiubio/file2meco).


## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li and Minjie Yao. 2021. microeco: an R package for data mining in microbial community ecology.
FEMS Microbiology Ecology, 97(2): fiaa255. https://doi.org/10.1093/femsec/fiaa255


## Contributing

We welcome any contribution, including but not limited to code, idea and [tutorial](https://chiliubio.github.io/microeco_tutorial/).
Please report errors and questions on github [Issues](https://github.com/ChiLiubio/microeco/issues).
Any contribution via [Pull requests](https://github.com/ChiLiubio/microeco/pulls) or Email(liuchi0426@126.com) will be appreciated.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CODE_OF_CONDUCT.md).


## References
  - Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272. DOI: 10.1126/science.aaf4507
  - Nguyen, N. H., Song, Z., Bates, S. T., Branco, S., Tedersoo, L., Menke, J., … Kennedy, P. G. (2016). 
    FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology, 20(1), 241–248.
  - Põlme, S., Abarenkov, K., Henrik Nilsson, R. et al. FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. Fungal Diversity 105, 1–16 (2020). DOI: 10.1007/s13225-020-00466-2
  - Aßhauer, K. P., Wemheuer, B., Daniel, R., & Meinicke, P. (2015). Tax4Fun: Predicting functional profiles from metagenomic 16S rRNA data. Bioinformatics, 31(17), 2882–2884.
  - Wemheuer, F., Taylor, J.A., Daniel, R. et al. Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences. Environmental Microbiome 15, 11 (2020). DOI: 10.1186/s40793-020-00358-7
  - Liu, C., Yao, M., Stegen, J. C., Rui, J., Li, J., & Li, X. (2017). Long-term nitrogen addition affects the phylogenetic turnover of soil microbial community responding to moisture pulse. Scientific Reports, 7(1), 17492.
  - Segata, N., Izard, J., Waldron, L., Gevers, D., Miropolsky, L., Garrett, W. S., & Huttenhower, C. (2011). Metagenomic biomarker discovery and explanation. Genome Biology, 12(6), R60.
  - Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255.
  - An, J., Liu, C., Wang, Q., Yao, M., Rui, J., Zhang, S., & Li, X. (2019). Soil bacterial community structure in Chinese wetlands. Geoderma, 337, 290–299.
  - Tackmann, J., Matias Rodrigues, J. F., & Mering, C. von. (2019). Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell Systems, 9(3), 286–296 e8.
  - White, J., Nagarajan, N., & Pop, M. (2009). Statistical methods for detecting differentially abundant features in clinical metagenomic samples. PLoS Computational Biology, 5(4), e1000352. 
  - Kurtz ZD, Muller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA. Sparse and compositionally robust inference of microbial ecological networks. PLoS Comput Biol 2015; 11: e1004226. 
  - McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8(4): e61217. 
  - Paulson, J., Stine, O., Bravo, H. et al. Differential abundance analysis for microbial marker-gene surveys. Nat Methods 10, 1200–1202 (2013). DOI: 10.1038/nmeth.2658
  - Deng Y, Jiang Y-H, Yang Y, He Z, Luo F, Zhou J. Molecular ecological network analyses. BMC bioinformatics 2012; 13: 113. 
  - Oksanen J, Blanchet FG, Friendly M, Kindt R, Legendre P, McGlinn D, et al. Vegan: Community ecology package. 2019. 
  - Picante: R tools for integrating phylogenies and ecology. Bioinformatics 2010; 26: 1463–1464.



