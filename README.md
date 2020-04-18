# microeco
An R package for ecological analysis of microbial communities

### Installing R/RStudio
If you do not already have R/RStudio installed, do as follows.

1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://rstudio.com/)
3. With Windows, install also [Rtools](https://cran.r-project.org/bin/windows/Rtools/)  

Put R and Rtools in the computer env path.  
Open RStudio...Tools...Global Options...Packages, select the appropriate mirror in Primary CRAN repository.

### Install microeco
Directly install microeco online.
```r
# require devtools package
devtools::install_github("ChiLiubio/microeco")
```
If failed because of the bad internet, download the package first, then install it.
```r
devtools::install_local("microeco-master.zip")
```

### Use
See the detailed package [tutorial](https://chiliubio.github.io/microeco/) and the help documentations.
If you want to run the codes in the tutorial completely, you need to install some packages, see the following Notes part.


### QQ
If the user has problems or suggestions, feel free to join the QQ group for discussions.  
QQ group: 277434916

### Notes

#### Packages important
To keep the start and use of the package simplified, 
the installation of microeco package only depend on several packages, which are compulsory-installed and very useful in the data analysis.
So the question is that you may encounter an error when using a class or function like this:

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

<div style="background-color: #FFFF00; color:red">
<strong>Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...</strong>
</div>

<span style="background-color: #FFFF00; color:red">
<strong>Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...</strong>
</span>

The solutions:

1. install the package when encounter such error. It's not troublesome.

2. install the packages in advance. We first show some packages that are necessary in some functions.


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
</tbody>
</table>
    </section>
  </div>
</div>


Then, if you want to install these packages or some of them, you can do like this:

```r
# If a package is not installed, it will be installed from CRAN.
# First select the packages of interest
packages <- c("GUniFrac", "picante", "agricolae", "ggpubr", "ggdendro", "MASS", "randomForest", "ggrepel", "pheatmap", "WGCNA", "igraph", "rgexf")
# Now check or install
package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
  }
)
```

Besides, some dependency packages of WGCNA are stored in Bioconductor. The tax4fun package can be downloaded from http://tax4fun.gobics.de/
and the correponding SILVA123 ref data is also required.















