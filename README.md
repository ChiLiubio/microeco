# microeco
An R package for ecological analysis of microbial communities

### Installing R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)
3. With Windows, install also [Rtools](https://cran.r-project.org/bin/windows/Rtools/)  

Put R and Rtools in the computer env PATH: your_directory\R-3.6.3\bin\x64, your_directory\Rtools\bin and your_directory\Rtools\mingw_64\bin  
Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

### Install microeco
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

### Use
See the detailed package [tutorial](https://chiliubio.github.io/microeco/) and the help documentations.
If you want to run the codes in the tutorial completely, you need to install some additional packages, see the following Notes part.


### QQ
If the user has problems or suggestions, feel free to join the QQ group for discussions.  
QQ group: 277434916

### Notes

#### Packages important
To keep the start and use of the package simplified, 
the installation of microeco package only depend on several packages, which are compulsory-installed and very useful in the data analysis.
These packages include R6, ape, vegan, rlang, data.table, magrittr, dplyr, tibble, reshape2, scales, grid, ggplot2, RColorBrewer, Rcpp, RcppArmadillo and RcppEigen.
So the question is that you may encounter an error when using a class or function that invoke an additional package like this:

```r
library(microeco)
data(sample_info)
data(otu_table)
data(taxonomy_table)
data(phylo_tree)
dataset <- microtable$new(sample_table = sample_info, otu_table = otu_table, tax_table = taxonomy_table, phylo_tree = phylo_tree)
dataset$tidy_dataset()
dataset$cal_betadiv(unifrac = TRUE)
```


```html
Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...
```


<br>

The solutions:

1. install the package when encounter such an error. Indeed, it's very easy to install in Rstudio. Just try it.

2. install the packages in advance. We recommend this solution if you are interest at many methods of the microeco package. We first show some packages that are necessary in some functions.


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
<td align="center">UniFrac beta diversity matrix</td>
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
<td align="center">WGCNA</td>
<td align="center">trans_network$new(cal_cor = “WGCNA”,…)</td>
<td align="center">invoke WGCNA package to calcuate correlations</td>
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
<td align="center">VGAM</td>
<td align="center">trans_corr class</td>
<td align="center">Generates Dirichlet random variates in SparCC</td>
</tr>
<tr class="even">
<td align="center">RJSONIO</td>
<td align="center">trans_func</td>
<td align="center">the dependency of biom package</td>
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
	"ggrepel", "pheatmap", "WGCNA", "igraph", "rgexf", "VGAM", "RJSONIO")
# Now check or install
lapply(packages,
	function(x) {
		if (!require(x, character.only = TRUE)) {
			install.packages(x, dependencies = TRUE)
		}
	}
)
```

Besides, WGCNA also depends on the GO.db package,
which can be installed from Bioconductor (https://bioconductor.org/packages/release/data/annotation/html/GO.db.html).

##### chorddiag
The R package chorddiag can be installed from Github https://github.com/mattflor/chorddiag

##### Tax4Fun 
1. install Tax4Fun package
```r
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
```
2. download SILVA123 reference data from http://tax4fun.gobics.de/
　unzip SILVA123.zip , move it to a place you can remember


##### python
Install python 

1. download python 2.7 from https://www.python.org/downloads/release
2. With windows, put python in the computer env PATH manually, 
　such as your_directory_path\python and your_directory_path\python\Scripts
3. Open terminal or cmd or Powershell , run

```python
pip install numpy
pip install argparse
```

If the installation is too slow to be failed, use -i select the appropriate mirror, for example, in China, you can use:

```python
pip install numpy -i https://pypi.douban.com/simple/
pip install argparse -i https://pypi.douban.com/simple/
```

##### FlashWeave

1. download and install julia from https://julialang.org/downloads/
2. Put julia in the computer env PATH, such as  your_directory_path\Julia\bin
3. Open terminal or cmd or Powershell, input julia, install FlashWeave following the operation in https://github.com/meringlab/FlashWeave.jl  
	or 
```julia
using Pkg
Pkg.add("FlashWeave")
```

### Acknowledgement
































