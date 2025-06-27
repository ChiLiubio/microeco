# microeco <a href="https://chiliubio.github.io/microeco_tutorial/"><img src="https://private-user-images.githubusercontent.com/20815519/459707743-39462f01-05d7-4051-9df0-b0f9098ac588.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTA5ODgyNTIsIm5iZiI6MTc1MDk4Nzk1MiwicGF0aCI6Ii8yMDgxNTUxOS80NTk3MDc3NDMtMzk0NjJmMDEtMDVkNy00MDUxLTlkZjAtYjBmOTA5OGFjNTg4LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA2MjclMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNjI3VDAxMzIzMlomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTVkYTE1OTMxOGRjNzhhMTVlZDE5M2UwODQyZmUyOTI2N2FjOGIwZDA3OGEzZjk1ZmY3MGJhZGQ0OThhYmUwYjkmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.1s2Q6LYptVrT6wkBcCVZ7VRDfEUYvx_xzB6zNbqlPNE" width=150 align="right" ></a>

An R package for data mining in microbial community ecology

[![CRAN](https://www.r-pkg.org/badges/version/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
[![CRAN](https://cranlogs.r-pkg.org/badges/grand-total/microeco)](https://cran.r-project.org/web/packages/microeco/index.html)
![](https://img.shields.io/badge/Release-1.15.0-orange.svg) ![](https://img.shields.io/badge/Test-1.16.0-red.svg)

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
When the user is in an R session and want to have a look on those links, please run the command: `?microeco`.
Please use a class name to search its help document (e.g., `?microtable`) instead of the function therein.
Before creating a new issue in the [Issues](https://github.com/ChiLiubio/microeco/issues), 
please read the guideline (https://chiliubio.github.io/microeco_tutorial/notes.html#github-issues).
Creating the basic microtable object from other tools/platforms (e.g. QIIME2, HUMAnN, Kraken2 and phyloseq) 
can be easily achieved with the package file2meco (https://github.com/ChiLiubio/file2meco).
Additionally, to facilitate users' familiarity and utilization of the microeco package, 
we have also developed a systematic protocol based on multi-omics data (https://github.com/ChiLiubio/microeco_protocol_v1).


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

