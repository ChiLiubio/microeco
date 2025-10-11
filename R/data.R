#' The dataset structured with microtable class for the demonstration of examples
#'
#' The dataset arose from 16S rRNA gene amplicon sequencing of wetland soils in China <doi:10.1016/j.geoderma.2018.09.035>. 
#' In \code{dataset$sample_table}, the 'Group' column means Chinese inland wetlands (IW), coastal wetland (CW) and Tibet plateau wetlands (TW).
#' The column 'Type' denotes the sampling region: northeastern region (NE), northwest region (NW), North China area (NC), 
#' middle-lower reaches of the Yangtze River (YML), southern coastal area (SC), upper reaches of the Yangtze River (YU) and Qinghai-Tibet Plateau (QTP). 
#' The column 'Saline' represents the saline soils and non-saline soils.
#'
#' \itemize{
#'   \item sample_table: sample information table
#'   \item otu_table: species-sample abundance table
#'   \item tax_table: taxonomic table
#'   \item phylo_tree: phylogenetic tree
#'   \item taxa_abund: taxa abundance list with several tables for Phylum...Genus
#'   \item alpha_diversity: alpha diversity table
#'   \item beta_diversity: list with several beta diversity distance matrix
#' }
#'
#' @docType data
#' @keywords datasets, R6
#' @name dataset
#' @usage data(dataset)
#' @format R6 class object
NULL

#' The environmental factors for the 16S example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name env_data_16S
#' @usage data(env_data_16S)
#' @format data.frame
NULL


#' The OTU table of the 16S example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name otu_table_16S
#' @usage data(otu_table_16S)
#' @format data.frame
NULL


#' The OTU table of the ITS example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name otu_table_ITS
#' @usage data(otu_table_ITS)
#' @format data.frame
NULL


#' The phylogenetic tree of 16S example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name phylo_tree_16S
#' @usage data(phylo_tree_16S)
#' @format data.frame
NULL


#' Customized FAPROTAX trait database
#'
#'
#' @docType data
#' @keywords datasets
#' @name prok_func_FAPROTAX
#' @usage data(prok_func_FAPROTAX)
#' @format list
NULL


#' The NJC19 database
#'
#'
#' @docType data
#' @keywords datasets
#' @name prok_func_NJC19_list
#' @usage data(prok_func_NJC19_list)
#' @format list
NULL


#' The FUNGuild database for fungi trait prediction
#'
#'
#' @docType data
#' @keywords datasets
#' @name fungi_func_FUNGuild
#' @usage data(fungi_func_FUNGuild)
#' @format data.frame
NULL


#' The FungalTraits database for fungi trait prediction
#'
#'
#' @docType data
#' @keywords datasets
#' @name fungi_func_FungalTraits
#' @usage data(fungi_func_FungalTraits)
#' @format data.frame
NULL


#' The KEGG data files used in the trans_func class
#'
#'
#' @docType data
#' @keywords datasets
#' @name Tax4Fun2_KEGG
#' @usage data(Tax4Fun2_KEGG)
#' @format list
NULL


#' The sample information of 16S example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name sample_info_16S
#' @usage data(sample_info_16S)
#' @format data.frame
NULL


#' The sample information of ITS example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name sample_info_ITS
#' @usage data(sample_info_ITS)
#' @format data.frame
NULL


#' The taxonomic information of 16S example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name taxonomy_table_16S
#' @usage data(taxonomy_table_16S)
#' @format data.frame
NULL


#' The taxonomic information of ITS example data
#'
#'
#' @docType data
#' @keywords datasets
#' @name taxonomy_table_ITS
#' @usage data(taxonomy_table_ITS)
#' @format data.frame
NULL
