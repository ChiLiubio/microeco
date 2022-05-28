# microeco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128353984-942c7dad-ccc6-4e5b-8672-8325d3d576f8.png" width=150 align="right" ></a>

An R package for data mining in microbial community ecology

[![CRAN](https://www.r-pkg.org/badges/version/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
![](https://img.shields.io/badge/Release-0.10.0-orange.svg) ![](https://img.shields.io/badge/Test-0.10.1-red.svg)

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
  + Differential abundance test
  + Machine learning
  + Null model analysis
  + Network analysis
  + Environmental data analysis
  + Functional prediction


## Install R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)

Put R in the computer env PATH, for example your_directory\R-4.1.0\bin\x64 

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
See the detailed package **tutorial** (https://chiliubio.github.io/microeco_tutorial/).
The tutorial can also be downloaded to the computer to open (https://github.com/ChiLiubio/microeco_tutorial/releases).
Please use the class name to search the help documents (e.g., `?microtable`).
Creating the basic microtable object from other tools/platforms (e.g. QIIME, QIIME2, HUMAnN, Kraken2 and phyloseq) 
can be easily achieved with the package file2meco (https://github.com/ChiLiubio/file2meco).


## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li and Minjie Yao. 2021. _microeco_: an R package for data mining in microbial community ecology.
FEMS Microbiology Ecology, 97(2): fiaa255. https://doi.org/10.1093/femsec/fiaa255


## Contributing

We welcome any contribution, including but not limited to code, idea and [tutorial](https://chiliubio.github.io/microeco_tutorial/).
Please report errors and questions on github [Issues](https://github.com/ChiLiubio/microeco/issues).
Any contribution via [Pull requests](https://github.com/ChiLiubio/microeco/pulls) will be appreciated.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CODE_OF_CONDUCT.md).


## References
  - https://chiliubio.github.io/microeco_tutorial/references.html#references

