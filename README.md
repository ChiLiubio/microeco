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

## Some packages important
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

Error in loadNamespace(name) : there is no package called ‘GUniFrac’ ...

The solutions:

1. install the package when encounter such error. It's not troublesome.

2. install the packages in advance. We first show some packages that are necessary in some functions.

----------------------------------------------------
    Package        function            description  
-------------- ------------- ----------------------
 **GUniFrac**   cal_betadiv       GUniFrac beta diversity matrix

 **OTU_236**    k__Bacteria      p__Chloroflexi       

 **OTU_399**    k__Bacteria    p__Proteobacteria   

 **OTU_1556**   k__Bacteria     p__Acidobacteria    

  **OTU_32**    k__Archaea      p__Miscellaneous       
                              Crenarchaeotic Group                         
----------------------------------------------------























