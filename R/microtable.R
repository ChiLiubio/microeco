#' @title
#' Create microtable object to store and manage all the basic files.
#'
#' @description
#' This class is a wrapper for a series of operations on the original files and the basic manipulations,
#' including the microtable object creation, data reduction, data rarefaction based on Paul et al. (2013) <doi:10.1371/journal.pone.0061217>, taxa abundance calculation, 
#' alpha and beta diversity calculation based on the An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and 
#' Lozupone et al. (2005) <doi:10.1128/AEM.71.12.8228–8235.2005> and other basic operations.
#'
#' @export
microtable <- R6Class(classname = "microtable",
	public = list(
		#' @param otu_table data.frame; necessary; The feature abundance table, rows are features, e.g. species, cols are samples.
		#' @param sample_table data.frame; default NULL; The sample information table, rows are samples, cols are sample metadata; 
		#' 	If not provided, the function can generate a table automatically according to the sample names in otu_table.
		#' @param tax_table data.frame; default NULL; The taxonomic information table, rows are species, cols are taxonomic classes.
		#' @param phylo_tree phylo; default NULL; The phylogenetic tree; use read.tree function in ape package for input.
		#' @param rep_fasta list; default NULL; The representative sequences; use read.fasta function in seqinr package for input.
		#' @return an object of class "microtable" with the following components:
		#' \describe{
		#'   \item{\code{sample_table}}{The sample information table.}
		#'   \item{\code{otu_table}}{The OTU table.}
		#'   \item{\code{tax_table}}{The taxonomic table.}
		#'   \item{\code{phylo_tree}}{The phylogenetic tree.}
		#'   \item{\code{rep_fasta}}{The representative sequence.}
		#'   \item{\code{taxa_abund}}{default NULL; use cal_abund function to calculate.}
		#'   \item{\code{alpha_diversity}}{default NULL; use cal_alphadiv function to calculate.}
		#'   \item{\code{beta_diversity}}{default NULL; use cal_betadiv function to calculate.}
		#' }
		#' @format microtable.
		#' @examples
		#' data(otu_table_16S)
		#' data(taxonomy_table_16S)
		#' data(sample_info_16S)
		#' data(phylo_tree_16S)
		#' dataset <- microtable$new(otu_table = otu_table_16S)
		#' dataset <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, 
		#'   tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
		#' # trim the files in the dataset
		#' dataset$tidy_dataset()
		initialize = function(otu_table, sample_table = NULL, tax_table = NULL, phylo_tree = NULL, rep_fasta = NULL)
			{
			if(!all(sapply(otu_table, is.numeric))){
				stop("Some columns in otu_table are not numeric vector! Please check the otu_table!")
			}else{
				otu_table <- private$check_abund_table(otu_table)
				self$otu_table <- otu_table
			}
			if(is.null(sample_table)){
				message("No sample_table provided, automatically use colnames in otu_table to create one ...")
				self$sample_table <- data.frame(SampleID = colnames(otu_table), Group = colnames(otu_table)) %>% 
					`row.names<-`(.$SampleID)
			}else{
				self$sample_table <- sample_table
			}
			# check whether phylogenetic tree is rooted
			if(!is.null(phylo_tree)){
				if(!ape::is.rooted(phylo_tree)){
					phylo_tree <- ape::multi2di(phylo_tree)
				}
			}
			self$tax_table <- tax_table
			self$phylo_tree <- phylo_tree
			self$rep_fasta <- rep_fasta
			self$taxa_abund <- NULL
			self$alpha_diversity <- NULL
			self$beta_diversity <- NULL
		},
		#' @description
		#' Print the microtable object.
		print = function(){
			cat("microtable class:\n")
			cat(paste("sample_table have", nrow(self$sample_table), "rows and", ncol(self$sample_table), "columns\n"))
			cat(paste("otu_table have", nrow(self$otu_table), "rows and", ncol(self$otu_table), "columns\n"))
			if(!is.null(self$tax_table)) cat(paste("tax_table have", nrow(self$tax_table), "rows and", ncol(self$tax_table), "columns\n"))
			if(!is.null(self$phylo_tree)) cat(paste("phylo_tree have", length(self$phylo_tree$tip.label), "tips\n"))
			if(!is.null(self$rep_fasta)) cat(paste("rep_fasta have", length(self$rep_fasta), "sequences\n"))
			if(!is.null(self$taxa_abund)) cat(paste("Taxa abundance: calculated for", paste0(names(self$taxa_abund), collapse = ","), "\n"))
			if(!is.null(self$alpha_diversity)) cat(paste("Alpha diversity: calculated for", paste0(colnames(self$alpha_diversity), collapse = ","), "\n"))
			if(!is.null(self$beta_diversity)) cat(paste("Beta diversity: calculated for", paste0(names(self$beta_diversity), collapse = ","), "\n"))
			invisible(self)
		},
		#' @description
		#' Filter the taxa considered as pollution from tax_table.
		#' This operation will remove any line of the tax_table containing any the word in taxa parameter regardless of word case.
		#'
		#' @param taxa default: c("mitochondria", "chloroplast"); filter mitochondria and chloroplast, or others as needed.
		#' @return None
		#' @examples 
		#' dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
		filter_pollution = function(taxa = c("mitochondria", "chloroplast")){
			if(is.null(self$tax_table)){
				stop("The tax_table in the microtable object is NULL ! Please check it!")
			}else{
				tax_table_use <- self$tax_table
			}
			if(length(taxa) > 1){
				taxa <- paste0(taxa, collapse = "|")
			}
			tax_table_use %<>% base::subset(unlist(lapply(data.frame(t(.)), function(x) !any(grepl(taxa, x, ignore.case = TRUE)))))
			filter_num <- nrow(self$tax_table) - nrow(tax_table_use)
			message(paste("Total", filter_num, "taxa are removed from tax_table ..."))
			self$tax_table <- tax_table_use
		},
		#' @description
		#' Rarefy communities to make all samples have same species number. See also \code{\link{rrarefy}} for the alternative method.
		#'
		#' @param sample.size default:NULL; species number, If not provided, use minimum number of all samples.
		#' @param rngseed random seed; default: 123.
		#' @param replace default: TRUE; See \code{\link{sample}} for the random sampling.
		#' @return None; rarefied dataset.
		#' @examples
		#' \donttest{
		#' dataset$rarefy_samples(sample.size = min(dataset$sample_sums()), replace = TRUE)
		#' }
		rarefy_samples = function(sample.size = NULL, rngseed = 123, replace = TRUE){
			set.seed(rngseed)
			self$tidy_dataset()
			if(is.null(sample.size)){
				sample.size <- min(self$sample_sums())
				message("Use the minimum number across samples: ", sample.size)
			}
			if (length(sample.size) > 1) {
				stop("`sample.size` had more than one value !")
			}
			if (sample.size <= 0) {
				stop("sample.size less than or equal to zero. ", "Need positive sample size to work !")
			}
			if (max(self$sample_sums()) < sample.size){
				stop("sample.size is larger than the maximum of sample sums, pleasure check input sample.size !")
			}
			if (min(self$sample_sums()) < sample.size) {
				rmsamples <- self$sample_names()[self$sample_sums() < sample.size]
				message(length(rmsamples), " samples removed, ", "because they contained fewer reads than `sample.size`.")
				self$sample_table <- base::subset(self$sample_table, ! self$sample_names() %in% rmsamples)
				self$tidy_dataset()
			}
			newotu <- self$otu_table
			newotu <- as.data.frame(apply(newotu, 2, private$rarefaction_subsample, sample.size = sample.size, replace = replace))
			rownames(newotu) <- self$taxa_names()
			self$otu_table <- newotu
			# remove OTUs with 0 sequence
			rmtaxa <- self$taxa_names()[self$taxa_sums() == 0]
			if (length(rmtaxa) > 0) {
				message(length(rmtaxa), " OTUs were removed because they are no longer present in any sample after random subsampling ...")
				self$tax_table <- base::subset(self$tax_table, ! self$taxa_names() %in% rmtaxa)
				self$tidy_dataset()
			}
		},
		#' @description
		#' Tidy the object of microtable Class.
		#' Trim files in the object to make taxa and samples consistent across all files in the object. So the results are intersections.
		#'
		#' @param main_data TRUE or FALSE, if TRUE, only basic files in microtable object is trimmed, otherwise, all files, 
		#' 	  including taxa_abund, alpha_diversity and beta_diversity, are all trimed.
		#' @return None, Object of microtable itself cleaned up. 
		#' @examples
		#' dataset$tidy_dataset(main_data = TRUE)
		tidy_dataset = function(main_data = TRUE){
			# check whether there is 0 abundance in otu_table
			self$otu_table <- private$check_abund_table(self$otu_table)
			
			sample_names <- intersect(rownames(self$sample_table), colnames(self$otu_table))
			if(length(sample_names) == 0){
				stop("No same sample name found between rownames of sample_table and colnames of otu_table! Please check whether the rownames of sample_table are sample names!")
			}
			# keep the sample order same with raw sample table
			sample_names <- rownames(self$sample_table) %>% .[. %in% sample_names]
			self$sample_table %<>% .[sample_names, , drop = FALSE]
			self$otu_table %<>% .[ , sample_names, drop = FALSE]
			# trim taxa
			self$otu_table %<>% {.[apply(., 1, sum) > 0, , drop = FALSE]}
			taxa_list <- list(rownames(self$otu_table), rownames(self$tax_table), self$phylo_tree$tip.label) %>% 
				.[!unlist(lapply(., is.null))]
			taxa_names <- Reduce(intersect, taxa_list)
			self$otu_table %<>% .[taxa_names, , drop = FALSE]
			if(!is.null(self$tax_table)){
				self$tax_table %<>% .[taxa_names, , drop = FALSE]
			}
			if(!is.null(self$phylo_tree)){
				self$phylo_tree %<>% ape::drop.tip(., base::setdiff(.$tip.label, taxa_names))
			}
			if(!is.null(self$rep_fasta)){
				self$rep_fasta %<>% .[taxa_names]
			}
			# other files will also be changed if main_data FALSE
			if(main_data == F){
				if(!is.null(self$taxa_abund)){
					self$taxa_abund %<>% lapply(., function(x) x[, sample_names, drop = FALSE])
				}
				if(!is.null(self$alpha_diversity)){
					self$alpha_diversity %<>% .[sample_names, , drop = FALSE]
				}
				if(!is.null(self$beta_diversity)){
					self$beta_diversity %<>% lapply(., function(x) x[sample_names, sample_names, drop = FALSE])
				}
			}
		},
		#' @description
		#' Add the rownames of tax_table as the last column of tax_table. 
		#' This is especially useful when the rownames of tax_table are required as a taxonomic level for the following taxa_abund calculation and biomarker idenfification.
		#'
		#' @param use_name default "OTU"; The column name used in the tax_table.
		#' @return new tax_table stored in object.
		#' @examples
		#' \donttest{
		#' dataset$add_rownames2taxonomy()
		#' }
		add_rownames2taxonomy = function(use_name = "OTU"){
			if(is.null(self$tax_table)){
				stop("The tax_table in the microtable object is NULL ! However it is necessary!")
			}else{
				tax_table_use <- self$tax_table
			}
			tax_table_use <- data.frame(tax_table_use, rownames(tax_table_use), check.names = FALSE, stringsAsFactors = FALSE)
			if(use_name %in% colnames(tax_table_use)){
				stop("The input use_name: ", use_name, " has been used in the raw tax_table! Please check it!")
			}
			colnames(tax_table_use)[ncol(tax_table_use)] <- use_name
			self$tax_table <- tax_table_use
		},
		#' @description
		#' Calculate the taxonomic abundance at each taxonomic rank.
		#'
		#' @param select_cols default NULL; numeric vector or character vector of colnames of tax_table; used to select columns to merge and calculate abundances.
		#'   This is very useful if there are commented columns or some columns with multiple structure that cannot be used directly.
		#' @param rel default TRUE; if TRUE, relative abundance is used; if FALSE, absolute abundance will be summed.
		#' @param split_group default FALSE; if TRUE, split the rows to multiple rows according to one or more columns in tax_table. Very useful when multiple mapping info exist.
		#' @param split_by default "&&"; Separator delimiting collapsed values; only useful when split_group == TRUE; see sep in separate_rows function.
		#' @param split_column default NULL; character vector or list; only useful when split_group == TRUE; character vector: 
		#'     fixed column or columns used for the splitting in tax_table in each abundance calculation; 
		#'     list: containing more character vectors to assign the column names to each calculation, such as list(c("Phylum"), c("Phylum", "Class")).
		#' @return taxa_abund in object.
		#' @examples
		#' \donttest{
		#' dataset$cal_abund()
		#' }
		cal_abund = function(
			select_cols = NULL, 
			rel = TRUE, 
			split_group = FALSE, 
			split_by = "&&", 
			split_column = NULL
			){
			taxa_abund = list()
			if(is.null(self$tax_table)){
				stop("No tax_table found! Please check your data!")
			}
			# check data corresponding
			if(nrow(self$tax_table) != nrow(self$otu_table)){
				message("The row number of tax_table is not equal to that of otu_table ...")
				message("Automatically applying tidy_dataset() function to trim the data ...")
				self$tidy_dataset()
				print(self)
			}
			# check whether no row in tax_table
			if(nrow(self$tax_table) == 0){
				stop("0 rows in tax_table! Please check your data!")
			}
			if(is.null(select_cols)){
				select_cols <- 1:ncol(self$tax_table)
			}else{
				if(!is.numeric(select_cols)){
					if(any(! select_cols %in% colnames(self$tax_table))){
						stop("Part of the input names are not in the tax_table!")
					}else{
						select_cols <- match(select_cols, colnames(self$tax_table))
					}
				}
			}
			if(split_group){
				if(is.null(split_column)){
					stop("Spliting rows by one or more columns require split_column parameter! Please set split_column and try again!")
				}
			}
			for(i in seq_along(select_cols)) {
				taxrank <- colnames(self$tax_table)[select_cols[i]]
				# assign the columns used for the splitting.
				if(!is.null(split_column)){
					if(is.list(split_column)){
						use_split_column <- split_column[[i]]
					}else{
						use_split_column <- split_column
					}
				}
				taxa_abund[[taxrank]] <- private$transform_data_proportion(
											self, 
											columns = select_cols[1:i], 
											rel = rel, 
											split_group = split_group, 
											split_by = split_by, 
											split_column = use_split_column)
			}
			self$taxa_abund <- taxa_abund
			message('The result is stored in object$taxa_abund ...')
		},
		#' @description
		#' Save taxonomic abundance to the computer local place.
		#'
		#' @param dirpath default "taxa_abund"; directory name to save the taxonomic abundance files.
		#' @examples
		#' \dontrun{
		#' dataset$save_abund(dirpath = "taxa_abund")
		#' }
		save_abund = function(dirpath = "taxa_abund"){
			if(!dir.exists(dirpath)){
				dir.create(dirpath)
			}
			for(i in names(self$taxa_abund)){
				write.csv(self$taxa_abund[[i]], file = paste0(dirpath, "/", i, "_abund.csv"), row.names = TRUE)
			}
		},
		#' @description
		#' Sum the species number for each sample.
		#'
		#' @return species number of samples.
		#' @examples
		#' \donttest{
		#' dataset$sample_sums()
		#' }
		sample_sums = function(){
			colSums(self$otu_table)
		},
		#' @description
		#' Sum the species number for each taxa.
		#'
		#' @return species number of taxa.
		#' @examples
		#' \donttest{
		#' dataset$taxa_sums()
		#' }
		taxa_sums = function(){
			rowSums(self$otu_table)
		},
		#' @description
		#' Show sample names.
		#'
		#' @return sample names.
		#' @examples
		#' \donttest{
		#' dataset$sample_names()
		#' }
		sample_names = function(){
			rownames(self$sample_table)
		},
		#' @description
		#' Show taxa names of tax_table.
		#'
		#' @return taxa names.
		#' @examples
		#' \donttest{
		#' dataset$taxa_names()
		#' }
		taxa_names = function(){
			rownames(self$tax_table)
		},
		#' @description
		#' Rename the taxa, including the rownames of otu_table, rownames of tax_table, tip labels of phylogenetic tree and representative sequences.
		#'
		#' @param newname_prefix default "ASV_"; the prefix of new names; new names will be newname_prefix + numbers according to the rowname order of otu_table.
		#' @return renamed dataset.
		#' @examples
		#' \donttest{
		#' dataset$rename_taxa()
		#' }
		rename_taxa = function(newname_prefix = "ASV_"){
			self$tidy_dataset()
			# extract old names for futher matching
			old_names <- rownames(self$otu_table)
			new_names <- paste0(newname_prefix, seq_len(nrow(self$otu_table)))
			rownames(self$otu_table) <- new_names
			if(!is.null(self$tax_table)){
				rownames(self$tax_table) <- new_names
			}
			if(!is.null(self$phylo_tree)){
				self$phylo_tree$tip.label[match(old_names, self$phylo_tree$tip.label)] <- new_names
			}
			if(!is.null(self$rep_fasta)){
				names(self$rep_fasta)[match(old_names, names(self$rep_fasta))] <- new_names
			}
		},
		#' @description
		#' Merge samples according to specific group to generate a new microtable.
		#'
		#' @param use_group the group column in sample_table.
		#' @return a new merged microtable object.
		#' @examples
		#' \donttest{
		#' dataset$merge_samples(use_group = "Group")
		#' }
		merge_samples = function(use_group){
			otu_table <- self$otu_table
			sample_table <- self$sample_table
			if(!is.null(self$tax_table)){
				tax_table <- self$tax_table
			}else{
				tax_table <- NULL
			}
			if(!is.null(self$phylo_tree)){
				phylo_tree <- self$phylo_tree
			}else{
				phylo_tree <- NULL
			}
			if(!is.null(self$rep_fasta)){
				rep_fasta <- self$rep_fasta
			}else{
				rep_fasta <- NULL
			}
			otu_table_new <- rowsum(t(otu_table), as.factor(as.character(sample_table[, use_group]))) %>% t %>% as.data.frame
			sample_table_new <- data.frame(SampleID = unique(as.character(sample_table[, use_group]))) %>% `row.names<-`(.[,1])
			# return a new microtable object
			microtable$new(sample_table = sample_table_new, otu_table = otu_table_new, tax_table = tax_table, phylo_tree = phylo_tree, rep_fasta = rep_fasta)
		},
		#' @description
		#' Merge taxa according to specific taxonomic rank to generate a new microtable.
		#'
		#' @param taxa the specific rank in tax_table.
		#' @return a new merged microtable object.
		#' @examples
		#' \donttest{
		#' dataset$merge_taxa(taxa = "Genus")
		#' }
		merge_taxa = function(taxa = "Genus"){
			# Agglomerate all OTUs by given taxonomic level
			ranknumber <- which(colnames(self$tax_table) %in% taxa)
			sampleinfo <- self$sample_table
			abund <- self$otu_table
			tax <- self$tax_table
			tax <- tax[, 1:ranknumber, drop=FALSE]
			# concatenate taxonomy in case of duplicated taxa names
			merged_taxonomy <- apply(tax, 1, paste, collapse="|")
			abund1 <- cbind.data.frame(Display = merged_taxonomy, abund) %>% 
				reshape2::melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
			abund1 <- data.table(abund1)[, sum_abund:=sum(Abundance), by=list(Display, Sample)] %>% 
				.[, c("Abundance"):=NULL] %>% 
				setkey(Display, Sample) %>% 
				unique() %>% 
				as.data.frame()
			# use dcast to generate table
			new_abund <- as.data.frame(data.table::dcast(data.table(abund1), Display~Sample, value.var= list("sum_abund"))) %>% 
				`row.names<-`(.[,1]) %>% 
				.[,-1, drop = FALSE]
			new_abund <- new_abund[order(apply(new_abund, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]
			# choose OTU names with highest abundance to replace the long taxonomic information in names
			name1 <- cbind.data.frame(otuname = rownames(tax), Display = merged_taxonomy, abundance = apply(abund[rownames(tax), ], 1, sum))
			name1 <- data.table(name1)[, max_abund:=max(abundance), by = Display]
			name1 <- name1[max_abund == abundance] %>% 
				.[, c("abundance", "max_abund"):=NULL] %>% 
				setkey(Display) %>% 
				unique() %>% 
				as.data.frame()
			name1 <- name1[!duplicated(name1$Display), ] %>% 
				`row.names<-`(.$Display)
			rownames(new_abund) <- name1[rownames(new_abund), "otuname"]
			new_tax <- tax[rownames(new_abund), , drop = FALSE]
			microtable$new(sample_table = sampleinfo, otu_table = new_abund, tax_table = new_tax)
		},
		#' @description
		#' Calculate alpha diversity in microtable object.
		#'
		#' @param measures default NULL; one or more indexes from "Observed", "Coverage", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "PD"; 
		#'   If null, use all those measures.
		#' @param PD TRUE or FALSE, whether phylogenetic tree should be calculated, default FALSE.
		#' @return alpha_diversity stored in object.
		#' @examples
		#' \donttest{
		#' dataset$cal_alphadiv(measures = NULL, PD = FALSE)
		#' class(dataset$alpha_diversity)
		#' }
		cal_alphadiv = function(measures = NULL, PD = FALSE){
			if (!any(self$otu_table == 1)){
				warning("The data you have provided does not have\n", 
					"any singletons. This is highly suspicious. \n", 
					"Results of richness estimates are probably unreliable, or wrong.")
			}
			OTU <- as.data.frame(t(self$otu_table), check.names = FALSE)
			renamevec    <-     c("Observed", "Coverage", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
			names(renamevec) <- c("S.obs", "coverage", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher")
			if (is.null(measures)){
				measures <- as.character(renamevec)
			}
			if (any(measures %in% names(renamevec))){
				measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% measures]
			}
			if (!any(measures %in% renamevec)){
				stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
			}
			outlist <- vector("list")
			estimRmeas <- c("Chao1", "Observed", "ACE")
			if (any(estimRmeas %in% measures)){
				outlist <- c(outlist, list(t(data.frame(vegan::estimateR(OTU), check.names = FALSE))))
			}
			if ("Shannon" %in% measures){
				outlist <- c(outlist, list(shannon = vegan::diversity(OTU, index = "shannon")))
			}
			if ("Simpson" %in% measures){
				outlist <- c(outlist, list(simpson = vegan::diversity(OTU, index = "simpson")))
			}
			if ("InvSimpson" %in% measures){
				outlist <- c(outlist, list(invsimpson = vegan::diversity(OTU, index = "invsimpson")))
			}
			if ("Fisher" %in% measures) {
				fisher = tryCatch(vegan::fisher.alpha(OTU, se = TRUE), warning = function(w) {
					suppressWarnings(vegan::fisher.alpha(OTU, se = TRUE)[, c("alpha", "se")])
				})
				if (!is.null(dim(fisher))) {
					colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
					outlist <- c(outlist, list(fisher))
				}
				else {
					outlist <- c(outlist, Fisher = list(fisher))
				}
			}
			outlist <- c(outlist, list(coverage = private$goods(OTU)))
			if(PD == T){
				if(is.null(self$phylo_tree)){
					stop("Please provide phylogenetic tree for PD calculation!")
				}else{
					outlist <- c(outlist, list(PD = picante::pd(OTU, self$phylo_tree)[,"PD", drop=TRUE]))
				}
			}
			res <- do.call("cbind", outlist)
			namechange <- base::intersect(colnames(res), names(renamevec))
			colnames(res)[colnames(res) %in% namechange] <- renamevec[namechange]
			self$alpha_diversity <- as.data.frame(res)
			message('The result is stored in object$alpha_diversity ...')
		},
		#' @description
		#' Save alpha diversity table to the computer.
		#'
		#' @param dirpath default "alpha_diversity"; directory name to save the alpha_diversity.csv file.
		save_alphadiv = function(dirpath = "alpha_diversity"){
			if(!dir.exists(dirpath)){
				dir.create(dirpath)
				# stop("The directory is not found, please first create it!")
			}
			write.csv(self$alpha_diversity, file = paste0(dirpath, "/", "alpha_diversity.csv"), row.names = TRUE)
		},
		#' @description
		#' Calculate beta diversity in microtable object, including Bray-Curtis, Jaccard, and UniFrac.
		#' See An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and Lozupone et al. (2005) <doi:10.1128/AEM.71.12.8228–8235.2005>.
		#'
		#' @param method default NULL; a character vector with one or more elements; If default, "bray" and "jaccard" will be used; 
		#'   see \code{\link{vegdist}} function and method parameter in vegan package. 
		#' @param unifrac default FALSE; TRUE or FALSE, whether unifrac index should be calculated.
		#' @param binary default FALSE; TRUE is used for jaccard and unweighted unifrac; optional for other indexes.
		#' @param ... parameters passed to \code{\link{vegdist}} function.
		#' @return beta_diversity stored in object.
		#' @examples
		#' \donttest{
		#' dataset$cal_betadiv(unifrac = FALSE)
		#' class(dataset$beta_diversity)
		#' }
		cal_betadiv = function(method = NULL, unifrac = FALSE, binary = FALSE, ...){
			res <- list()
			eco_table <- t(self$otu_table)
			sample_table <- self$sample_table
			if(is.null(method)){
				method <- c("bray", "jaccard")
			}
			for(i in method){
				if(i == "jaccard"){
					binary_use <- TRUE
				}else{
					binary_use <- binary
				}
				res[[i]] <- as.matrix(vegan::vegdist(eco_table, method = i, binary = binary_use, ...))
			}
			
			if(unifrac == T){
				if(is.null(self$phylo_tree)){
					stop("No phylogenetic tree provided, please change the parameter unifrac to FALSE")
				}
				phylo_tree <- self$phylo_tree
				# require GUniFrac package; do not consider too much about alpha parameter
				unifrac1 <- GUniFrac::GUniFrac(eco_table, phylo_tree, alpha = c(0, 0.5, 1))
				unifrac2 <- unifrac1$unifracs
				wei_unifrac <- unifrac2[,, "d_1"]
				res$wei_unifrac <- wei_unifrac
				unwei_unifrac <- unifrac2[,, "d_UW"]
				res$unwei_unifrac <- unwei_unifrac
			}
			self$beta_diversity <- res
			message('The result is stored in object$beta_diversity ...')
		},
		#' @description
		#' Save beta diversity matrix to the computer.
		#'
		#' @param dirpath default "beta_diversity"; directory name to save the beta diversity matrix files.
		save_betadiv = function(dirpath = "beta_diversity"){
			if(!dir.exists(dirpath)){
				dir.create(dirpath)
				# stop("The directory is not found, please first create it!")
			}
			for(i in names(self$beta_diversity)){
				write.csv(self$beta_diversity[[i]], file = paste0(dirpath, "/", i, ".csv"), row.names = TRUE)
			}
		}
		),
	private = list(
		# check whether there is OTU or sample with 0 abundance
		# input and return are both otu_table
		check_abund_table = function(otu_table){
			if(any(apply(otu_table, 1, sum) == 0)){
				remove_num <- sum(apply(otu_table, 1, sum) == 0)
				message(remove_num, " taxa are removed from the otu_table, as the abundance is 0 ...")
				otu_table %<>% .[apply(., 1, sum) > 0, , drop = FALSE]
			}
			if(any(apply(otu_table, 2, sum) == 0)){
				remove_num <- sum(apply(otu_table, 2, sum) == 0)
				message(remove_num, " samples are removed from the otu_table, as the abundance is 0 ...")
				otu_table %<>% .[, apply(., 2, sum) > 0, drop = FALSE]
			}
			if(ncol(otu_table) == 0){
				stop("No sample have abundance! Please check you data!")
			}
			if(nrow(otu_table) == 0){
				stop("No taxon have abundance! Please check you data!")
			}
			otu_table
		},
		# taxa abundance calculation function
		transform_data_proportion = function(
			input, 
			columns, 
			rel, 
			split_group = FALSE, 
			split_by = "&&", 
			split_column
			){
			sampleinfo <- input$sample_table
			abund <- input$otu_table
			tax <- input$tax_table
			tax <- tax[, columns, drop=FALSE]
			# split rows to multiple rows if multiple correspondence
			if(split_group){
				merge_abund <- cbind.data.frame(tax, abund)
				split_merge_abund <- tidyr::separate_rows(merge_abund, all_of(split_column), sep = split_by)
				new_tax <- split_merge_abund[, columns, drop = FALSE]
				new_abund <- split_merge_abund[, (ncol(tax) + 1):(ncol(merge_abund)), drop = FALSE]
				abund1 <- cbind.data.frame(Display = apply(new_tax, 1, paste, collapse="|"), new_abund)
			}else{
				abund1 <- cbind.data.frame(Display = apply(tax, 1, paste, collapse="|"), abund)
			}
			# first convert table to long format
			# then sum abundance by sample and taxonomy
			abund1 <- abund1 %>% 
				data.table() %>% 
				data.table::melt(id.vars = "Display", value.name= "Abundance", variable.name = "Sample") %>%
				.[, sum_abund:=sum(Abundance), by=list(Display, Sample)] %>% 
				.[, c("Abundance"):=NULL] %>% 
				setkey(Display, Sample) %>% 
				unique()
			if(rel == T){
				abund1 <- abund1[, res_abund:=sum_abund/sum(sum_abund), by=list(Sample)] %>% 
					.[, c("sum_abund"):=NULL]
			}else{
				colnames(abund1)[colnames(abund1) == "sum_abund"] <- "res_abund"
			}
			# dcast the table
			abund2 <- as.data.frame(data.table::dcast(abund1, Display~Sample, value.var= list("res_abund"))) %>%
				`row.names<-`(.[,1]) %>% 
				.[,-1, drop = FALSE]
			abund2 <- abund2[order(apply(abund2, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]
			abund2
		},
		rarefaction_subsample = function(x, sample.size, replace=FALSE){
			# Adapted from the rarefy_even_depth() in phyloseq package, see Paul et al. (2013) <doi:10.1371/journal.pone.0061217>.
			# All rights reserved.
			# Create replacement species vector
			rarvec <- numeric(length(x))
			# Perform the sub-sampling. Suppress warnings due to old R compat issue.
			if(sum(x) <= 0){
				# Protect against, and quickly return an empty vector, 
				# if x is already an empty count vector
				return(rarvec)
			}
			if(replace){
				suppressWarnings(subsample <- sample(1:length(x), sample.size, replace = TRUE, prob=x))
			} else {
				# resample without replacement
				obsvec <- apply(data.frame(OTUi=1:length(x), times=x), 1, function(x){
					rep_len(x["OTUi"], x["times"])
				})
				obsvec <- unlist(obsvec, use.names=FALSE)
				suppressWarnings(subsample <- sample(obsvec, sample.size, replace = FALSE))
			}
			sstab <- table(subsample)
			# Assign the tabulated random subsample values to the species vector
			rarvec[as(names(sstab), "integer")] <- sstab
			return(rarvec)
		},
		# define coverage function
		goods = function(com){
			no.seqs <- rowSums(com)
			sing <- com==1
			no.sing <- apply(sing, 1, sum)
			goods <- 1-no.sing/no.seqs
			return(goods)
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
