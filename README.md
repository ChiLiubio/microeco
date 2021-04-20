# microeco
An R package for data mining in microbial community ecology

![](https://img.shields.io/badge/Release-Ver0.3.3-blue.svg) ![](https://img.shields.io/badge/Test-Ver0.4.0-red.svg)

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

Put R in the computer env PATH, for example your_directory\R-4.0.0\bin\x64 

Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

## Install microeco

Install the latest microeco from github.
This is the **best** way as there are some minor improvements here compared to the version in CRAN.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
# then install microeco
devtools::install_github("ChiLiubio/microeco")
```

If the installation from github is failed because of the bad internet, download the package first, then install it locally.

```r
devtools::install_local("microeco-master.zip")
```

It is also good to install the microeco from CRAN directly. There may be some update delay in the CRAN version.

```r
install.packages("microeco")
```

## Tutorial
See the detailed package tutorial (https://chiliubio.github.io/microeco/) and the help documentations.
If the tutorial website can not be opened because of bad internet connection, you can download the microeco-master.zip and open the index.html to see the tutorial 
using your browser directly.
If you want to run the codes in the tutorial and README completely, you need to install some additional packages, see the following **Notes** part.
Beside the demonstration in the tutorial, users can also save the intermediate files in each object and
apply those files to other tools according to the format requirement.
Main files stored in the objects of microeco package is the commonly used data.frame format.
So the intermediate and result files are easily saved and used for other tools in microbial ecology.
Contructing the basic microtable object from other tools/platforms (e.g. QIIME, QIIME2 and phyloseq) 
can be easily achieved using the package file2meco (https://github.com/ChiLiubio/file2meco).
There are some important features and approaches that exist in other microbiome analysis software/platforms and currently not implemented in microeco package.
More approaches for the file conversion will be provided in the package file2meco.


## Citation
Chi Liu, Yaoming Cui, Xiangzhen Li, Minjie Yao, microeco: an R package for data mining in microbial community ecology, 
FEMS Microbiology Ecology, Volume 97, Issue 2, February 2021, fiaa255, https://doi.org/10.1093/femsec/fiaa255


## Notes

### packages important
To keep the start and use of microeco package simplified, 
the installation of microeco depend on some packages, which are compulsory-installed from CRAN and useful in the data analysis.
These packages include
R6, stats, ape, vegan, rlang, data.table, magrittr, dplyr, tibble, scales, grid, ggplot2, RColorBrewer.
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

1. install the package when encounter such an error. Actually, it's very easy to install the packages of CRAN in Rstudio. Just try it.

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
packages <- c("reshape2", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "gridExtra", "picante", "pheatmap", "igraph", "rgexf", "ggalluvial")
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

#### SpiecEasi

The R package SpiecEasi can be used for the network construction using SPIEC-EASI (SParse InversE Covariance Estimation for Ecological Association Inference) approach.
The package can be installed from Github https://github.com/zdk123/SpiecEasi

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

#### chorddiag

The R package chorddiag is used for the chord plot in the network analysis and can be installed from Github https://github.com/mattflor/chorddiag

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


#### Tax4Fun
Tax4Fun is an R package used for the prediction of functional potential of prokaryotic communities.

1. install Tax4Fun package
```r
install.packages("RJSONIO")
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
```
2. download SILVA123 reference data from http://tax4fun.gobics.de/
　unzip SILVA123.zip , move it to a place you can remember


#### Tax4Fun2
Tax4Fun2 is another R package for the the prediction of functional profiles and functional gene redundancies of prokaryotic communities.
It has higher accuracies than PICRUSt and Tax4Fun. The Tax4Fun2 approach implemented in microeco is a little different from the original package.
Using Tax4Fun2 approach require the representative fasta file.
The user do not need to install Tax4Fun2 R package.
The only thing need to do is to download the blast tool and Ref99NR/Ref100NR.
Downlaod blast tools from "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+" ; e.g. ncbi-blast-\*\*\*\*-x64-win64.tar.gz  for windows system.
Downlaod Ref99NR.zip from "https://cloudstor.aarnet.edu.au/plus/s/DkoZIyZpMNbrzSw/download"  or Ref100NR.zip from "https://cloudstor.aarnet.edu.au/plus/s/jIByczak9ZAFUB4/download" .
Uncompress all the folders. The final folders should be like these structures:

blast tools:  
　|-- ncbi-blast-2.11.0+  
　　|---- bin  
　　　|------ blastn.exe  
　　　|------ makeblastdb.exe  
　　　|------ ......  

Ref99NR/Ref100NR:  
　|-- Tax4Fun2_ReferenceData_v2  
　　|---- Ref99NR  
　　　|------ otu000001.tbl.gz  
　　　|------ ......  
　　　|------ Ref99NR.fasta  
　　　|------ Ref99NR.tre  

The path "ncbi-blast-2.11.0+/bin" and "Tax4Fun2_ReferenceData_v2" will be required in the trans_func$cal_tax4fun2() function.

```r
# seqinr should be installed for reading and writing fasta file
install.packages("seqinr", dependencies = TRUE)
# Now we show how to read the fasta file
# see https://github.com/ChiLiubio/file2meco if you do not have installed file2meco
rep_fasta_path <- system.file("extdata", "rep.fna", package="file2meco")
rep_fasta <- seqinr::read.fasta(rep_fasta_path)
# then see the help document of microtable class about the rep_fasta in microtable$new().
```

## Plotting
Most of the plotting in the package rely on the ggplot2 package system.
We provide some parameters to change the corresponding plot.
If you want to change the output plot, you can also assign the output a name and use the ggplot2-style grammer to modify it as you need.
Each data table used for plotting is stored in the object and can be downloaded for the personalized analysis and plotting.
Of course, you can also directly modify the function or class to reload them.

## Files from other tools to microtable object
Previous descriptions on how to construct microtable object from QIIME, QIIME2 and phyloseq have been moved to the package file2meco (https://github.com/ChiLiubio/file2meco)
The package file2meco is designed to transform files from other tools/platforms into microtable object.


## sample_table in microtable
The rownames of sample_table are used for selecting samples/groups in all the related operations in the package.
Before you create microtable object, make sure that the rownames of sample_table are the sample names.


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
  - [microbiomeSeq](https://github.com/umerijaz/microbiomeSeq)
  - [microbiomeMarker](https://github.com/yiluheihei/microbiomeMarker)













