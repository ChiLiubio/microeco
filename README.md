# microeco
An R package for data mining in microbial community ecology

## Background
In microbial community ecology, with the development of the high-throughput sequencing techniques,
the increasing data amount and complexity make the community data analysis and management a challenge.
There has been a lot of R packages created for the microbiome profiling analysis, such as phyloseq,
microbiomeSeq, ampvis2, mare and microbiome.
However, it is still difficult to perform data mining fast and efficiently.
Based on this, we created an R package microeco.

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


## Installing R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)

Put R and Rtools in the computer env PATH, for example your_directory\R-4.0.0\bin\x64 .
Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

## Install microeco
Directly install microeco online from CRAN.

```r
install.packages("microeco")
```

Install microeco from github (beta version).
We suggest this method because there are some minor improvements compared to the version in CRAN.

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

## Tutorial
See the detailed package tutorial (https://chiliubio.github.io/microeco/) and the help documentations.
If the tutorial website can not be opened because of bad internet connection, you can download the microeco-master.zip and open the index.html to see the tutorial 
using your browser directly.
If you want to run the codes in the tutorial and README completely, you need to install some additional packages, see the following **Notes** part.
It is notable that there are some features that exist in other microbiome analysis software or platforms and currently not implemented in microeco package.
In the microeco package, we provide the functions for the conversions between microtable object and phyloseq object.
Beside the data conversion for the basic object, we suggest that users save the intermediate files in the operations of other objects in microeco package and
apply those files to other tools according to the format requirement.
The main files stored in the objects of microeco package is the commonly used data.frame format.
So the intermediate and result files are easily saved and checked for the use of other tools in microbial ecology.

## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: An R package for data mining in microbial community ecology, FEMS Microbiology Ecology, 
2020;, fiaa255, https://doi.org/10.1093/femsec/fiaa255


## Notes

### packages important
To keep the start and use of microeco package simplified, 
the installation of microeco depend on some packages, which are compulsory-installed from CRAN and useful in the data analysis.
These packages include
R6, stats, ape, vegan, rlang, data.table, magrittr, dplyr, tibble, reshape2, scales, grid, ggplot2, VGAM, MASS, RColorBrewer.
So the question is that you may encounter an error when using a class or function that invoke an additional package like this:

```r
library(microeco)
data(dataset)
t1 <- trans_network$new(dataset = dataset, cal_cor = NA, taxa_level = "OTU", filter_thres = 0.0005)
t1$cal_network(network_method = "SpiecEasi")
```


```html
Error in t1$cal_network(network_method = "SpiecEasi"): igraph package not installed ...
```


<br>
The reason is that network construction require igraph package. We donot put the igraph and some other packages (in bioconductor or github) on the "Imports" part of microeco package.

The solutions:

1. install the package when encounter such an error. Indeed, it's very easy to install the packages of CRAN in Rstudio. Just try it.

2. install the packages in advance. We recommend this solution if you are interest in most of the methods in the microeco package.

We show several packages that are published in CRAN and not installed automatically.



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
<tr class="even">
<td align="center">reshape2</td>
<td align="center">microtable class</td>
<td align="center">data transformation</td>
</tr>
<tr class="odd">
<td align="center">MASS</td>
<td align="center">trans_diff$new(method = "lefse",…)</td>
<td align="center">linear discriminant analysis</td>
</tr>
<tr class="even">
<td align="center">GUniFrac</td>
<td align="center">cal_betadiv()</td>
<td align="center">UniFrac distance matrix</td>
</tr>
<tr class="odd">
<td align="center">ggpubr</td>
<td align="center">plot_alpha()</td>
<td align="center">some plotting functions</td>
</tr>
<tr class="even">
<td align="center">randomForest</td>
<td align="center">trans_diff$new(method = "rf",…)</td>
<td align="center">random forest analysis</td>
</tr>
<tr class="odd">
<td align="center">ggdendro</td>
<td align="center">plot_clustering()</td>
<td align="center">plotting clustering dendrogram</td>
</tr>
<tr class="even">
<td align="center">ggrepel</td>
<td align="center">trans_rda class</td>
<td align="center">reduce the text overlap in the plot</td>
</tr>
<tr class="odd">
<td align="center">agricolae</td>
<td align="center">cal_diff(method = anova)</td>
<td align="center">multiple comparisons</td>
</tr>
<tr class="even">
<td align="center">gridExtra</td>
<td align="center">trans_diff class</td>
<td align="center">merge plots</td>
</tr>
<tr class="odd">
<td align="center">picante</td>
<td align="center">cal_alphadiv()</td>
<td align="center">Faith’s phylogenetic alpha diversity</td>
</tr>
<tr class="even">
<td align="center">pheatmap</td>
<td align="center">plot_corr(pheatmap = TRUE)</td>
<td align="center">correlation heatmap with clustering dendrogram</td>
</tr>
<tr class="odd">
<td align="center">tidytree</td>
<td align="center">trans_diff class</td>
<td align="center">plot the taxonomic tree</td>
</tr>
<tr class="even">
<td align="center">mice</td>
<td align="center">trans_env class</td>
<td align="center">Insert missing value in env data</td>
</tr>
<tr class="odd">
<td align="center">phyloseq</td>
<td align="center">meco2phyloseq()</td>
<td align="center">convert between microtable and phyloseq</td>
</tr>
<tr class="even">
<td align="center">qiime2R</td>
<td align="center">qiimed2meco()</td>
<td align="center">QIIME2 files to microtable object</td>
</tr>
<tr class="odd">
<td align="center">igraph</td>
<td align="center">trans_network class</td>
<td align="center">network related operations</td>
</tr>
<tr class="even">
<td align="center">rgexf</td>
<td align="center">save_network</td>
<td align="center">save network with gexf style</td>
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
packages <- c("reshape2", "GUniFrac", "MASS", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "gridExtra", "picante", "pheatmap", "igraph", "rgexf", "RJSONIO", "ggalluvial")
# Now check or install
lapply(packages, function(x) {
	if(!require(x, character.only = TRUE)) {
		install.packages(x, dependencies = TRUE)
	}})
```

There are also some packages that may be useful in some parameters or functions. These packages may be R packages published in github or bioconductor,
or packages written by other languages.


#### ggtree
Plotting the cladogram from the LEfSe result requires the ggtree package in bioconductor (https://bioconductor.org/packages/release/bioc/html/ggtree.html).
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
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


#### SpiecEasi

The R package SpiecEasi can be used for the network construction using SPIEC-EASI (SParse InversE Covariance Estimation for Ecological Association Inference) approach.
The package can be installed from Github https://github.com/zdk123/SpiecEasi

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


#### FlashWeave
FlashWeave is a julia package used for network analysis.
It predicts ecological interactions among microbes from large-scale compositional abundance data (i.e. OTU tables constructed from sequencing data) 
through statistical co-occurrence.

1. download and install julia from https://julialang.org/downloads/
2. Put julia in the computer env PATH, such as  your_directory_path\Julia\bin
3. Open terminal or cmd or Powershell, input julia, install FlashWeave following the operation in https://github.com/meringlab/FlashWeave.jl  

#### Gephi
Gephi is an excellent network visualization tool and used to open the saved network file, i.e. network.gexf in the [tutorial](https://chiliubio.github.io/microeco/).
You can download Gephi and learn how to use it from https://gephi.org/users/download/

## Plotting
Most of the plotting in the package rely on the ggplot2 package system.
We provide some parameters to change the corresponding plot.
If you want to change the output plot, you can also assign the output a name and use the ggplot2-style grammer to modify it as you need.
Each data table used for plotting is stored in the object and can be downloaded for the personalized analysis and plotting.
Of course, you can also directly modify the function or class to reload them.

## Read QIIME files
In this part, we first show how to construct the object of microtable class using the raw OTU file from QIIME.

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

## Read QIIME2 files
We provide a function qiimed2meco() to convert QIIME2 file to microtable object directly.
If you want to run the codes, first download the QIIME2 data files from https://docs.qiime2.org/2020.8/tutorials/pd-mice/


```r
# Please first install qiime2R from github, see https://github.com/jbisanz/qiime2R
libraray(qiime2R)
library(magrittr)
qiimed2meco <- function(ASV_data, sample_data, taxonomy_data, phylo_tree = NULL){
	# Read ASV data
	ASV <- as.data.frame(read_qza(ASV_data)$data)
	#  Read metadata
	metadata <- read_q2metadata(sample_data)
	rownames(metadata) <- as.character(metadata[, 1])
	# Read taxonomy table
	taxa_table <- read_qza(taxonomy_data)
	taxa_table <- parse_taxonomy(taxa_table$data)
	# Make the taxonomic table clean, this is very important.
	taxa_table %<>% tidy_taxonomy
	# Read phylo tree
	if(!is.null(phylo_tree)){
		phylo_tree <- read_qza(phylo_tree)$data
	}
	dataset <- microtable$new(sample_table = metadata, tax_table = taxa_table, otu_table = ASV, phylo_tree = phylo_tree)
	dataset
}
# first download QIIME2 files from https://docs.qiime2.org/2020.8/tutorials/pd-mice/
library(microeco)
meco_dataset <- qiimed2meco(ASV_data = "dada2_table.qza", sample_data = "sample-metadata.tsv", taxonomy_data = "taxonomy.qza", phylo_tree = "tree.qza")
meco_dataset
```


## Conversion between microtable and phyloseq
We provide two functions meco2phyloseq() and phyloseq2meco() for the conversion between microtable object and phyloseq object (phyloseq package).

```r
# Please first install phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
```

```r
# from microtable to phyloseq object
library(microeco)
library(phyloseq)
data("dataset")
physeq <- meco2phyloseq(dataset)
physeq
```

```r
# from phyloseq to microtable object
library(phyloseq)
library(microeco)
data("GlobalPatterns")
meco_dataset <- phyloseq2meco(GlobalPatterns)
meco_dataset
```




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
  - Kurtz ZD, Muller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA. Sparse and compositionally robust inference of microbial ecological networks. PLoS Comput Biol 2015; 11: e1004226. 

















