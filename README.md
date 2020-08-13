# microeco
An R package for data analysis in microbial community ecology

## Background
In microbial community ecology, with the development of the high-throughput sequencing techniques,
the increasing data amount and complexity make the community data analysis and management a challenge.
There has been a lot of R packages created for the community data analysis in microbial ecology, such as phyloseq,
microbiomeSeq, ampvis2, mare and microbiome.
However, it is still difficult to perform data mining fast and efficiently.
Based on this, we created R package microeco.

## Main Features
  + R6 Class to store and analyze data; fast, flexible and modularized
  + Plotting the taxonomic abundance
  + Venn diagram
  + Alpha diversity
  + Beta diversity
  + Differential abundance analysis
  + Indicator species analysis
  + Environmental data analysis
  + Network analysis
  + Null model analysis
  + Functional analysis


## Installing R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)
3. With Windows, install also [Rtools](https://cran.r-project.org/bin/windows/Rtools/)  

Put R and Rtools in the computer env PATH, for example your_directory\R-3.6.3\bin\x64, your_directory\Rtools\bin and your_directory\Rtools\mingw_64\bin  
Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

## Install microeco
Directly install microeco online.
```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco
devtools::install_github("ChiLiubio/microeco")
```
If the installation of microeco is failed because of the bad internet, download the package first, then install it.
```r
devtools::install_local("microeco-master.zip")
```

## Use
See the detailed package tutorial (https://chiliubio.github.io/microeco/) and the help documentations.
If you want to run the codes in the tutorial completely, you need to install some additional packages, see the following **Notes** part.



## Notes

### packages important
To keep the start and use of microeco package simplified, 
the installation of microeco only depend on several packages, which are compulsory-installed and very useful in the data analysis.
These packages include R6, ape, vegan, rlang, data.table, magrittr, dplyr, tibble, reshape2, scales, grid, ggplot2, RColorBrewer, Rcpp, RcppArmadillo and RcppEigen.
So the question is that you may encounter an error when using a class or function that invoke an additional package like this:

```r
library(microeco)
data(sample_info_16S)
data(otu_table_16S)
data(taxonomy_table_16S)
data(phylo_tree_16S)
dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
dataset$tidy_dataset()
dataset$cal_betadiv(unifrac = TRUE)
```


```html
Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...
```


<br>
The reason is that calculating unifrac distance require ‘GUniFrac’ package.

The solutions:

1. install the package when encounter such an error. Indeed, it's very easy to install in Rstudio. Just try it.

2. install the packages in advance. We recommend this solution if you are interest in many methods of the microeco package. We first show some packages that are necessary in some functions.


<div id="content-wrapper">
  <div class="inner clearfix">
    <section id="main-content">
<table>
<colgroup>
<col width="19%"></col>
<col width="38%"></col>
<col width="42%"></col>
</colgroup>
<thead>
<tr class="header">
<th align="center">Package</th>
<th align="center">where</th>
<th align="center">description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">GUniFrac</td>
<td align="center">cal_betadiv</td>
<td align="center">UniFrac distance matrix</td>
</tr>
<tr class="even">
<td align="center">picante</td>
<td align="center">cal_alphadiv</td>
<td align="center">Faith’s phylogenetic alpha diversity</td>
</tr>
<tr class="odd">
<td align="center">agricolae</td>
<td align="center">cal_diff(method = anova)</td>
<td align="center">multiple comparisons</td>
</tr>
<tr class="even">
<td align="center">ggpubr</td>
<td align="center">plot_alpha</td>
<td align="center">some plotting functions</td>
</tr>
<tr class="odd">
<td align="center">ggdendro</td>
<td align="center">plot_clustering</td>
<td align="center">plotting clustering dendrogram</td>
</tr>
<tr class="even">
<td align="center">MASS</td>
<td align="center">trans_diff$new(method = “lefse”,…)</td>
<td align="center">linear discriminant analysis</td>
</tr>
<tr class="odd">
<td align="center">randomForest</td>
<td align="center">trans_diff$new(method = “rf”,…)</td>
<td align="center">random forest analysis</td>
</tr>
<tr class="even">
<td align="center">ggrepel</td>
<td align="center">trans_rda</td>
<td align="center">reduce the text overlap in the plot</td>
</tr>
<tr class="odd">
<td align="center">pheatmap</td>
<td align="center">plot_corr(pheatmap = TRUE)</td>
<td align="center">correlation heatmap with clustering dendrogram</td>
</tr>
<tr class="even">
<td align="center">igraph</td>
<td align="center">trans_network class</td>
<td align="center">network related operations</td>
</tr>
<tr class="odd">
<td align="center">rgexf</td>
<td align="center">save_network</td>
<td align="center">save network with gexf style</td>
</tr>
<tr class="even">
<td align="center">VGAM</td>
<td align="center">trans_network class</td>
<td align="center">generate Dirichlet random variates in SparCC</td>
</tr>
<tr class="odd">
<td align="center">RJSONIO</td>
<td align="center">trans_func</td>
<td align="center">the dependency of biom package</td>
</tr>
<tr class="even">
<td align="center">ggalluvial</td>
<td align="center">plot_bar(use_alluvium = TRUE)</td>
<td align="center">alluvial plot</td>
</tr>
</tbody>
</table>
    </section>
  </div>
</div>


Then, if you want to install these packages or some of them, you can do like this:

```r
# If a package is not installed, it will be installed from CRAN.
# First select the packages of interest
packages <- c("GUniFrac", "picante", "agricolae", "ggpubr", "ggdendro", "MASS", "randomForest", 
	"ggrepel", "pheatmap", "igraph", "rgexf", "VGAM", "RJSONIO", "ggalluvial")
# Now check or install
lapply(packages, function(x) {
	if(!require(x, character.only = TRUE)) {
		install.packages(x, dependencies = TRUE)
	}})
```

#### WGCNA
In the correlation-based network, when the species number is large,
the correlation algorithm in WGCNA is very fast compared to the cor function in R base.
WGCNA depends on several packages in Bioconductor, including GO.db and impute.
So if you want to install WGCNA, first install  GO.db (https://bioconductor.org/packages/release/data/annotation/html/GO.db.html)
and impute (http://www.bioconductor.org/packages/release/bioc/html/impute.html) with the following code.

```r
# install GO.db and impute
# First check and install BiocManager package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# install GO.db and impute
BiocManager::install("GO.db")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
# then install WGCNA.
install.packages("WGCNA", dependencies = TRUE)
```

#### ggtree
Plotting the cladogram from the LEfSe result requires the ggtree package in bioconductor (https://bioconductor.org/packages/release/bioc/html/ggtree.html).
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
```

#### chorddiag

The R package chorddiag is used for the chord plot in the network analysis and can be installed from Github https://github.com/mattflor/chorddiag

#### Tax4Fun
Tax4Fun is an R package used for the prediction of functional potential of microbial communities.

1. install Tax4Fun package
```r
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
```
2. download SILVA123 reference data from http://tax4fun.gobics.de/
　unzip SILVA123.zip , move it to a place you can remember


#### python
Predicting the functional potential on the biogeochemical cycles require python 2.7 and several packages.

1. download python 2.7 from https://www.python.org/downloads/release
2. With windows, put python in the computer env PATH manually, 
　such as your_directory_path\python and your_directory_path\python\Scripts
3. Open terminal or cmd or Powershell , run

```python
pip install numpy
pip install argparse
```

If the installation is too slow and failed, use -i select the appropriate mirror, for example, in China, you can use:

```python
pip install numpy -i https://pypi.douban.com/simple/
pip install argparse -i https://pypi.douban.com/simple/
```

#### FlashWeave
FlashWeave is a julia package used for network analysis.
It predicts ecological interactions among microbes from large-scale compositional abundance data (i.e. OTU tables constructed from sequencing data) 
through statistical co-occurrence.

1. download and install julia from https://julialang.org/downloads/
2. Put julia in the computer env PATH, such as  your_directory_path\Julia\bin
3. Open terminal or cmd or Powershell, input julia, install FlashWeave following the operation in https://github.com/meringlab/FlashWeave.jl  

#### Gephi
Gephi is used to open saved network file, i.e. network.gexf in the [tutorial](https://chiliubio.github.io/microeco/).
You can download Gephi and learn how to use it from https://gephi.org/users/download/

## plotting
Most of the plotting in the package rely on the ggplot2 package system.
We provide some parameters to change the corresponding plot.
If you want to change the output plot, you can also assign the output a name and use the ggplot2-style grammer to modify it as you need.
Of course, you can also directly modify the function or class to reload them.

## read the raw files
In this part, we show how to construct the object of microtable class using the raw otu file from QIIME.

```r
# use the raw data files stored inside the package
otu_file_path <- system.file("extdata", "otu_table_raw.txt", package="microeco")
# the example sample table is csv style
sample_file_path <- system.file("extdata", "sample_info.csv", package="microeco")
# phylogenetic tree
phylo_file_path <- system.file("extdata", "rep_phylo.tre", package="microeco")
# load microeco and qiimer, if qiimer is not installed, see Tax4Fun part to install qiimer package
library(microeco)
library(qiimer)
# read and parse otu_table_raw.txt; this file does not have the first commented line, so we use commented = FALSE
otu_raw_table <- read_qiime_otu_table(otu_file_path, commented=FALSE)
# obtain the otu table data.frame
otu_table_1 <- as.data.frame(otu_raw_table[[3]])
colnames(otu_table_1) <- unlist(otu_raw_table[[1]])
# obtain the taxonomic table  data.frame
taxonomy_table_1 <- as.data.frame(split_assignments(unlist(otu_raw_table[[4]])))
# read sample metadata table, data.frame
sample_info <- read.csv(sample_file_path, row.names = 1, stringsAsFactors = FALSE)
# read the phylogenetic tree
phylo_tree <- read.tree(phylo_file_path)
# check whether the tree is rooted, if unrooted, transform to rooted
if(!is.rooted(phylo_tree)){
	phylo_tree <- multi2di(phylo_tree)
}
# make the taxonomic table clean, this is very important
taxonomy_table_1 %<>% tidy_taxonomy
# create a microtable object
dataset <- microtable$new(sample_table = sample_info, otu_table = otu_table_1, tax_table = taxonomy_table_1, phylo_tree = phylo_tree)
# for other operations, see the tutorial (https://chiliubio.github.io/microeco/) and the help documentations
# the class documentation include the function links, see the microtable class, input:
?microtable
```

## QQ
If the user has problems or suggestions, feel free to join the QQ group for discussions.  
QQ group: 207510995


## Acknowledgement
  - [R6](https://github.com/r-lib/R6), The
    main class system in this package.
  - [lefse python
    script](https://bitbucket.org/biobakery/biobakery/wiki/lefse), The
    main lefse codes are translated from **lefse python script**.
  - [phyloseq](https://github.com/joey711/phyloseq), the idea of data
    structures of microtable class in microeco comes from
    `phyloseq-class` in package **phyloseq**.
  - [microbiomeSeq](https://github.com/umerijaz/microbiomeSeq), 
    the method that calculates the roles of nodes within- and among- modules connectivity is 
    modified from the package **microbiomeSeq**.
  - [SpiecEasi](https://github.com/zdk123/SpiecEasi), 
    the method that calculates SparCC is
    modified from the package **SpiecEasi**.
  - [microbiomeMarker](https://github.com/yiluheihei/microbiomeMarker), 
    the method that plots the LEfSe cladogram is
    modified from the package **microbiomeMarker**.

## References
  - If the user use cal_spe_func function for the functional identification of prokaryotes or cal_FAPROTAX function in trans_func class, 
    please also cite the original FAPROTAX database paper: 
    Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272. DOI: 10.1126/science.aaf4507.
  - If the user use cal_spe_func function for the functional identification of fungi,
    please also cite the original FUNGuild database paper: 
    Nguyen, N. H., Song, Z., Bates, S. T., Branco, S., Tedersoo, L., Menke, J., … Kennedy, P. G. (2016). 
    FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology, 20(1), 241–248.
  - If the user use trans_nullmodel class, please also cite the original paper: 
    Liu, C., Yao, M., Stegen, J. C., Rui, J., Li, J., & Li, X. (2017). Long-term nitrogen addition affects the phylogenetic turnover of soil microbial community responding to moisture pulse. Scientific Reports, 7(1), 17492. Journal Article.
  - If the user use LefSe method in trans_diff class, please also cite the original paper:
    Segata, N., Izard, J., Waldron, L., Gevers, D., Miropolsky, L., Garrett, W. S., & Huttenhower, C. (2011). Metagenomic biomarker discovery and explanation. Genome Biology, 12(6), R60.
  - An, J., Liu, C., Wang, Q., Yao, M., Rui, J., Zhang, S., & Li, X. (2019). Soil bacterial community structure in Chinese wetlands. Geoderma, 337, 290–299.
  - Aßhauer, K. P., Wemheuer, B., Daniel, R., & Meinicke, P. (2015). Tax4Fun: Predicting functional profiles from metagenomic 16S rRNA data. Bioinformatics, 31(17), 2882–2884.
  - Tackmann, J., Matias Rodrigues, J. F., & Mering, C. von. (2019). Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell Systems, 9(3), 286–296 e8.
  - White, J., Nagarajan, N., & Pop, M. (2009). Statistical methods for detecting differentially abundant features in clinical metagenomic samples. PLoS Computational Biology, 5(4), e1000352. 


















