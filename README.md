# microeco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://user-images.githubusercontent.com/20815519/128353984-942c7dad-ccc6-4e5b-8672-8325d3d576f8.png" width=150 align="right" ></a>

An R package for data mining in microbial community ecology

[![CRAN](https://www.r-pkg.org/badges/version/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
![](https://img.shields.io/badge/Release-1.13.0-orange.svg) ![](https://img.shields.io/badge/Test-1.13.1-red.svg)

## Background
With the development of high-throughput sequencing techniques,
the increasing data amount and complexity make the microbiome omics data analysis and management a challenge.
Though there has been a lot of R packages in this filed, 
it is still difficult to perform data mining fast, efficiently and comprehensively.
Therefore, we created R microeco package (abbreviated and pronounced as **_[miːkəu]_**).

## Main Features
  + R6 Class to store and analyze data: flexible and modularized
  + Data normalization
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
If you do not already have R/RStudio installed, follow these steps:

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://posit.co/downloads/)

Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

## Install microeco

Install microeco package from CRAN.

```r
install.packages("microeco")
```

Or install the latest development version from Github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("ChiLiubio/microeco")
```


## Tutorial
Detailed online **tutorial** (https://chiliubio.github.io/microeco_tutorial/) is released along with the package.
The tutorial can also be downloaded to the computer to open (https://github.com/ChiLiubio/microeco_tutorial/releases).
When you are in an R session and want to have a look on those links, please run the command: `?microeco`.
Please use a class name to search its help document (e.g., `?microtable`) instead of the function therein.
Before creating a new issue in the [Issues](https://github.com/ChiLiubio/microeco/issues), 
please read the guideline (https://chiliubio.github.io/microeco_tutorial/notes.html#github-issues).
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

