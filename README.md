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

<div id="content-wrapper">
  <div class="inner clearfix">
    <section id="main-content">
<p><span style="color:red"> 这个非常重要，一定不能搞错。 </span></p>
    </section>
  </div>
</div>

```html
<p>
<span style="color:red">
<strong>Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...</strong>
</span>
</p>
```

```html
<span style="background-color: #FFFF00; color:red">
<strong>Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...</strong>
</span>

```

<div id="content-wrapper">
  <div class="inner clearfix">
    <section id="main-content">
<p><img src="data:image/png;base64,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" width="600px" style="display: block; margin: auto;" /></p>
    </section>
  </div>
</div>

<br>

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















