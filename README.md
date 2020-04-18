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

<font color='red'> Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ... </font>

The solutions:

1. install the package when encounter such error. It's not troublesome.

2. install the packages in advance. We first show some packages that are necessary in some functions.

-----------------------------------------------------------------------------
   Package                 where                       description           
-------------- ----------------------------- --------------------------------
   GUniFrac             cal_betadiv           UniFrac beta diversity matrix  

   picante             cal_alphadiv             Faith's phylogenetic alpha   
                                                        diversity            

  agricolae      cal_diff(method = anova)          multiple comparisons      

    ggpubr              plot_alpha               some plotting functions     

   ggdendro           plot_clustering         plotting clustering dendrogram 

     MASS         trans_diff$new(method =      linear discriminant analysis  
                       "lefse",...)                                          

 randomForest     trans_diff$new(method =         random forest analysis     
                         "rf",...)                                           

   ggrepel               trans_rda            reduce the text overlap in the 
                                                           plot              

   pheatmap     plot_corr(pheatmap = TRUE)       correlation heatmap with    
                                                  clustering dendrogram      

    WGCNA       trans_network$new(cal_cor =      invoke WGCNA package to     
                       "WGCNA",...)               calcuate correlations      

    igraph          trans_network class         network related operations   

    rgexf              save_network            save network with gexf style  
-----------------------------------------------------------------------------






















