#' @title
#' Create \code{microtable} object to store and manage all the basic files.
#'
#' @description
#' This class is a wrapper for a series of operations on the input files and basic manipulations,
#' including microtable object creation, data trimming, data filtering, rarefaction based on Paul et al. (2013) <doi:10.1371/journal.pone.0061217>, taxonomic abundance calculation, 
#' alpha and beta diversity calculation based on the An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and 
#' Lozupone et al. (2005) <doi:10.1128/AEM.71.12.8228-8235.2005> and other basic operations.\cr
#' \cr
#' Online tutorial: \href{https://chiliubio.github.io/microeco_tutorial/}{https://chiliubio.github.io/microeco_tutorial/} \cr
#' Download tutorial: \href{https://github.com/ChiLiubio/microeco_tutorial/releases}{https://github.com/ChiLiubio/microeco_tutorial/releases}
#' 
#' @export
microtable <- R6Class(classname = "microtable",
	public = list(
		#' @param otu_table data.frame; The feature abundance table; rownames are features (e.g. OTUs/ASVs/species/genes); column names are samples.
		#' @param sample_table data.frame; default NULL; The sample information table; rownames are samples; columns are sample metadata; 
		#' 	 If not provided, the function can generate a table automatically according to the sample names in otu_table.
		#' @param tax_table data.frame; default NULL; The taxonomic information table; rownames are features; column names are taxonomic classes.
		#' @param phylo_tree phylo; default NULL; The phylogenetic tree that must be read with the \code{read.tree} function of ape package.
		#' @param rep_fasta \code{DNAStringSet}, \code{list} or \code{DNAbin} format; default NULL; The sequences.
		#'   The sequences should be read with the \code{readDNAStringSet} function in \code{Biostrings} package (DNAStringSet class), 
		#'   \code{read.fasta} function in \code{seqinr} package (list class),
		#'   or \code{read.FASTA} function in \code{ape} package (DNAbin class).
		#' @param auto_tidy default FALSE; Whether tidy the files in the \code{microtable} object automatically.
		#'   If TRUE, the function can invoke the \code{tidy_dataset} function.
		#' @return an object of class \code{microtable} with the following components:
		#' \describe{
		#'   \item{\code{sample_table}}{The sample information table.}
		#'   \item{\code{otu_table}}{The feature table.}
		#'   \item{\code{tax_table}}{The taxonomic table.}
		#'   \item{\code{phylo_tree}}{The phylogenetic tree.}
		#'   \item{\code{rep_fasta}}{The sequence.}
		#'   \item{\code{taxa_abund}}{default NULL; use \code{cal_abund} function to calculate.}
		#'   \item{\code{alpha_diversity}}{default NULL; use \code{cal_alphadiv} function to calculate.}
		#'   \item{\code{beta_diversity}}{default NULL; use \code{cal_betadiv} function to calculate.}
		#' }
		#' @format microtable.
		#' @examples
		#' data(otu_table_16S)
		#' data(taxonomy_table_16S)
		#' data(sample_info_16S)
		#' data(phylo_tree_16S)
		#' m1 <- microtable$new(otu_table = otu_table_16S)
		#' m1 <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, 
		#'   tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)
		#' # trim the files in the dataset
		#' m1$tidy_dataset()
		initialize = function(otu_table, sample_table = NULL, tax_table = NULL, phylo_tree = NULL, rep_fasta = NULL, auto_tidy = FALSE)
			{
			if(missing(otu_table)){
				stop("otu_table must be provided!")
			}
			if(!inherits(otu_table, "data.frame")){
				stop("The input otu_table must be data.frame format!")
			}
			otu_table <- private$check_abund_table(otu_table)
			self$otu_table <- otu_table
			if(is.null(sample_table)){
				message("No sample_table provided, automatically use colnames in otu_table to create one ...")
				self$sample_table <- data.frame(SampleID = colnames(otu_table), Group = colnames(otu_table)) %>% 
					`row.names<-`(.$SampleID)
			}else{
				if(!inherits(sample_table, "data.frame")){
					stop("Input sample_table must be data.frame format!")
				}
				if(inherits(sample_table, "tbl_df")){
					stop("Input sample_table is of tbl_df class! It may be created using tibble package! Please convert it to traditional data.frame class!")
				}
				self$sample_table <- sample_table
			}
			if(!is.null(tax_table)){
				if(!inherits(tax_table, "data.frame")){
					stop("The input tax_table must be data.frame format!")
				}
			}
			if(!is.null(phylo_tree)){
				if(!inherits(phylo_tree, "phylo")){
					stop("The input phylo_tree must be phylo format! Please use read.tree function in ape package to read the phylogenetic tree!")
				}
				if(!ape::is.rooted(phylo_tree)){
					phylo_tree <- ape::multi2di(phylo_tree)
				}
			}
			if(!is.null(rep_fasta)){
				if(!(inherits(rep_fasta, "list") | inherits(rep_fasta, "DNAbin") | inherits(rep_fasta, "DNAStringSet"))){
					stop("Unknown fasta format! Must be one of DNAStringSet (from readDNAStringSet function of Biostrings package), ", 
						"list (from read.fasta function of seqinr package), and DNAbin (from read.FASTA function of ape package)!")
				}
			}
			self$tax_table <- tax_table
			self$phylo_tree <- phylo_tree
			self$rep_fasta <- rep_fasta
			self$taxa_abund <- NULL
			self$alpha_diversity <- NULL
			self$beta_diversity <- NULL
			self$auto_tidy <- auto_tidy
			if(self$auto_tidy) self$tidy_dataset()
		},
		#' @description
		#' Filter the features considered pollution in \code{microtable$tax_table}.
		#' This operation will remove any line of the \code{microtable$tax_table} containing any the word in taxa parameter regardless of word case.
		#'
		#' @param taxa default \code{c("mitochondria", "chloroplast")}; filter mitochondria and chloroplast, or others as needed.
		#' @return None
		#' @examples 
		#' m1$filter_pollution(taxa = c("mitochondria", "chloroplast"))
		filter_pollution = function(taxa = c("mitochondria", "chloroplast")){
			if(is.null(self$tax_table)){
				stop("The tax_table in the microtable object is NULL ! Please check it!")
			}
			tax_table_use <- self$tax_table
			if(length(taxa) > 1){
				taxa <- paste0(taxa, collapse = "|")
			}
			tax_table_use %<>% base::subset(unlist(lapply(data.frame(t(.)), function(x) !any(grepl(taxa, x, ignore.case = TRUE)))))
			filter_num <- nrow(self$tax_table) - nrow(tax_table_use)
			message("Total ", filter_num, " features are removed from tax_table ...")
			self$tax_table <- tax_table_use
			if(self$auto_tidy) self$tidy_dataset()
		},
		#' @description
		#' Filter the feature with low abundance and/or low occurrence frequency.
		#'
		#' @param rel_abund default 0; the relative abundance threshold, such as 0.0001.
		#' @param freq default 1; the occurrence frequency threshold. 
		#' 	 For example, the number 2 represents filtering the feature that occurs less than 2 times.
		#' 	 A number smaller than 1 is also allowable. 
		#' 	 For instance, the number 0.1 represents filtering the feature that occurs in less than 10\% samples.
		#' @param include_lowest default TRUE; whether include the feature with the threshold.
		#' @return None
		#' @examples
		#' \donttest{
		#' d1 <- clone(m1)
		#' d1$filter_taxa(rel_abund = 0.0001, freq = 0.2)
		#' }
		filter_taxa = function(rel_abund = 0, freq = 1, include_lowest = TRUE){
			raw_otu_table <- self$otu_table
			if(rel_abund != 0){
				if(rel_abund >= 1){
					stop("rel_abund must be smaller than 1!")
				}
				taxa_raw_abund <- self$taxa_sums()
				taxa_rel_abund <- taxa_raw_abund/sum(taxa_raw_abund)
				if(include_lowest){
					abund_names <- taxa_rel_abund[taxa_rel_abund < rel_abund] %>% names
				}else{
					abund_names <- taxa_rel_abund[taxa_rel_abund <= rel_abund] %>% names
				}
				if(length(abund_names) == nrow(raw_otu_table)){
					stop("No feature remained! Please check the rel_abund parameter!")
				}
				message(length(abund_names), " features filtered based on the abundance ...")
			}else{
				abund_names <- c()
			}
			if(freq != 0){
				if(freq < 1){
					message("freq smaller than 1; first convert it to an integer ...")
					freq <- round(ncol(raw_otu_table) * freq)
					message("Use converted freq integer: ", freq, " for the following filtering ...")
				}
				taxa_occur_num <- apply(raw_otu_table, 1, function(x){sum(x != 0)})
				if(include_lowest){
					freq_names <- taxa_occur_num[taxa_occur_num < freq] %>% names
				}else{
					freq_names <- taxa_occur_num[taxa_occur_num <= freq] %>% names
				}
				if(length(freq_names) == nrow(raw_otu_table)){
					stop("No feature remained! Please check the freq parameter!")
				}
				message(length(freq_names), " features filtered based on the occurrence...")
			}else{
				freq_names <- c()
			}
			filter_names <- c(abund_names, freq_names) %>% unique
			if(length(filter_names) == 0){
				new_table <- raw_otu_table
			}else{
				if(length(filter_names) == nrow(raw_otu_table)){
					stop("All features are filtered! Please adjust the parameters")
				}
				new_table <- raw_otu_table[! rownames(raw_otu_table) %in% filter_names, , drop = FALSE]
			}
			self$otu_table <- new_table
			self$tidy_dataset()
		},
		#' @description
		#' Rarefy communities to make all samples have same count number.
		#'
		#' @param method default c("rarefy", "SRS")[1]; "rarefy" represents the classical resampling like \code{rrarefy} function of \code{vegan} package.
		#'    "SRS" is scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <DOI:10.7717/peerj.9593>.
		#' @param sample.size default NULL; libray size. If not provided, use the minimum number across all samples. 
		#'    For "SRS" method, this parameter is passed to \code{Cmin} parameter of \code{SRS} function of SRS package.
		#' @param ... parameters pass to \code{norm} function of \code{\link{trans_norm}} class.
		#' @return None; rarefied dataset.
		#' @examples
		#' \donttest{
		#' m1$rarefy_samples(sample.size = min(m1$sample_sums()))
		#' }
		rarefy_samples = function(method = c("rarefy", "SRS")[1], sample.size, ...){
			self$tidy_dataset()
			if(method == "rarefying"){
				method <- "rarefy"
			}
			method <- match.arg(method, c("rarefy", "SRS"))
			tmp <- suppressMessages(trans_norm$new(self))
			tmp_new <- tmp$norm(method = method, sample.size = sample.size, ...)
			self$otu_table <- tmp_new$otu_table
			suppressMessages(self$tidy_dataset())
		},
		#' @description
		#' Trim all the data in the \code{microtable} object to make taxa and samples consistent. The results are intersections across data.
		#'
		#' @param main_data default FALSE; if TRUE, only basic data in \code{microtable} object is trimmed. Otherwise, all data, 
		#' 	  including \code{taxa_abund}, \code{alpha_diversity} and \code{beta_diversity}, are all trimed.
		#' @return None. The data in the object are tidied up. 
		#' 	  If \code{tax_table} is in object, its row names are totally same with the row names of \code{otu_table}.
		#' @examples
		#' m1$tidy_dataset(main_data = TRUE)
		tidy_dataset = function(main_data = FALSE){
			self <- private$tidy_samples(self)
			self$otu_table %<>% {.[apply(., 1, sum) > 0, , drop = FALSE]}
			taxa_list <- list(rownames(self$otu_table), rownames(self$tax_table), self$phylo_tree$tip.label) %>% 
				.[!unlist(lapply(., is.null))]
			taxa_names <- Reduce(intersect, taxa_list)
			if(length(taxa_names) == 0){
				if(is.null(self$phylo_tree)){
					stop("No same feature name found between rownames of otu_table and rownames of tax_table! Please check rownames of those tables!")
				}else{
					stop("No same feature name found among otu_table, tax_table and phylo_tree! Please check feature names in those objects!")
				}
			}
			self$otu_table %<>% .[taxa_names, , drop = FALSE]
			if(!is.null(self$tax_table)){
				self$tax_table %<>% .[taxa_names, , drop = FALSE]
			}
			if(!is.null(self$phylo_tree)){
				self$phylo_tree %<>% ape::drop.tip(., base::setdiff(.$tip.label, taxa_names))
			}
			if(!is.null(self$rep_fasta)){
				invisible(self$rep_fasta[1])
				fasta_names <- names(self$rep_fasta)
				if(is.null(fasta_names)){
					stop("The name of rep_fasta is NULL! Please provide a correct fasta file!")
				}
				if(!all(taxa_names %in% fasta_names)){
					stop("Some feature names are not found in the names of rep_fasta! Please provide a complete fasta file or manually check the names!")
				}
				self$rep_fasta %<>% .[taxa_names]
			}
			# check again whether a sample has 0 abundance after feature filtering
			self <- private$tidy_samples(self)
			if(!main_data){
				sample_names <- rownames(self$sample_table)
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
		#' Add the rownames of \code{microtable$tax_table} as its last column. 
		#' This is especially useful when the rownames of \code{microtable$tax_table} are required as a taxonomic level 
		#' 	 for the taxonomic abundance calculation and biomarker idenfification.
		#'
		#' @param use_name default "OTU"; The column name used in the \code{tax_table}.
		#' @return NULL, a new tax_table stored in the object.
		#' @examples
		#' \donttest{
		#' m1$add_rownames2taxonomy()
		#' }
		add_rownames2taxonomy = function(use_name = "OTU"){
			self$tidy_dataset()
			if(is.null(self$tax_table)){
				message("The tax_table in the microtable object is NULL! Create one ...")
				tax_table_use <- data.frame(rownames(self$otu_table), stringsAsFactors = FALSE)
				rownames(tax_table_use) <- rownames(self$otu_table)
			}else{
				tax_table_use <- self$tax_table
				tax_table_use <- data.frame(tax_table_use, rownames(tax_table_use), check.names = FALSE, stringsAsFactors = FALSE)
				if(use_name %in% colnames(tax_table_use)){
					stop("The input use_name: ", use_name, " has been used in the raw tax_table! Please check it!")
				}
			}
			colnames(tax_table_use)[ncol(tax_table_use)] <- use_name
			message("Use ", use_name, " as the name of new column in tax_table ...")
			self$tax_table <- tax_table_use
		},
		#' @description
		#' Sum the species number for each sample.
		#'
		#' @return species number of samples.
		#' @examples
		#' \donttest{
		#' m1$sample_sums()
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
		#' m1$taxa_sums()
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
		#' m1$sample_names()
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
		#' m1$taxa_names()
		#' }
		taxa_names = function(){
			rownames(self$otu_table)
		},
		#' @description
		#' Rename the features, including the rownames of \code{otu_table}, rownames of \code{tax_table}, tip labels of \code{phylo_tree} and \code{rep_fasta}.
		#'
		#' @param newname_prefix default "ASV_"; the prefix of new names; new names will be newname_prefix + numbers according to the rownames order of \code{otu_table}.
		#' @return None; renamed dataset.
		#' @examples
		#' \donttest{
		#' m1$rename_taxa()
		#' }
		rename_taxa = function(newname_prefix = "ASV_"){
			self$tidy_dataset()
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
		#' Merge samples according to specific group to generate a new \code{microtable}.
		#'
		#' @param group a column name in \code{sample_table} of \code{microtable} object.
		#' @return a new merged microtable object.
		#' @examples
		#' \donttest{
		#' m1$merge_samples("Group")
		#' }
		merge_samples = function(group){
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
			if(! group %in% colnames(sample_table)){
				stop("Provided parameter group must be a column name of sample_table !")
			}
			otu_table_new <- rowsum(t(otu_table), as.factor(as.character(sample_table[, group]))) %>% t %>% as.data.frame
			sample_table_new <- data.frame(SampleID = unique(as.character(sample_table[, group]))) %>% `row.names<-`(.[,1])
			microtable$new(
				sample_table = sample_table_new, 
				otu_table = otu_table_new, 
				tax_table = tax_table, 
				phylo_tree = phylo_tree, 
				rep_fasta = rep_fasta,
				auto_tidy = self$auto_tidy
			)
		},
		#' @description
		#' Merge taxa according to specific taxonomic rank to generate a new \code{microtable}.
		#'
		#' @param taxa default "Genus"; the specific rank in \code{tax_table}.
		#' @return a new merged \code{microtable} object.
		#' @examples
		#' \donttest{
		#' m1$merge_taxa(taxa = "Genus")
		#' }
		merge_taxa = function(taxa = "Genus"){
			check_tax_level(taxa, self)
			ranknumber <- which(colnames(self$tax_table) %in% taxa)
			sampleinfo <- self$sample_table
			abund <- self$otu_table
			tax <- self$tax_table
			tax <- tax[, 1:ranknumber, drop=FALSE]
			# concatenate taxonomy in case of duplicated taxa names
			merged_taxonomy <- apply(tax, 1, paste, collapse = "|")
			abund1 <- cbind.data.frame(Display = merged_taxonomy, abund) %>% 
				reshape2::melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
			abund1 <- data.table(abund1)[, sum_abund:=sum(Abundance), by=list(Display, Sample)] %>% 
				.[, c("Abundance"):=NULL] %>% 
				setkey(Display, Sample) %>% 
				unique() %>% 
				as.data.frame()
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
			microtable$new(sample_table = sampleinfo, otu_table = new_abund, tax_table = new_tax, auto_tidy = self$auto_tidy)
		},
		#' @description
		#' Save each basic data in microtable object as local file.
		#'
		#' @param dirpath default "basic_files"; directory to save the tables, phylogenetic tree and sequences in microtable object. It will be created if not found.
		#' @param sep default ","; the field separator string, used to save tables. Same with \code{sep} parameter in \code{\link{write.table}} function.
		#' 	  default \code{','} correspond to the file name suffix 'csv'. The option \code{'\t'} correspond to the file name suffix 'tsv'. For other options, suffix are all 'txt'.
		#' @param ... parameters passed to \code{\link{write.table}}.
		#' @examples
		#' \dontrun{
		#' m1$save_table()
		#' }
		save_table = function(dirpath = "basic_files", sep = ",", ...){
			if(!dir.exists(dirpath)){
				dir.create(dirpath)
			}
			suffix <- switch(sep, ',' = "csv", '\t' = "tsv", "txt")
			tmp <- self$otu_table
			tmp <- data.frame(ID = rownames(tmp), tmp)
			save_path <- paste0(dirpath, "/feature_table.", suffix)
			write.table(tmp, file = save_path, row.names = FALSE, sep = sep, ...)
			message('Save feature abundance to ', save_path, ' ...')
			tmp <- self$sample_table
			tmp <- data.frame(ID = rownames(tmp), tmp)
			save_path <- paste0(dirpath, "/sample_table.", suffix)
			write.table(tmp, file = save_path, row.names = FALSE, sep = sep, ...)
			message('Save metadata to ', save_path, ' ...')
			if(!is.null(self$tax_table)){
				tmp <- self$tax_table
				tmp <- data.frame(ID = rownames(tmp), tmp)
				save_path <- paste0(dirpath, "/tax_table.", suffix)
				write.table(tmp, file = save_path, row.names = FALSE, sep = sep, ...)
				message('Save taxonomic information to ', save_path, ' ...')
			}
			if(!is.null(self$phylo_tree)){
				tmp <- self$phylo_tree
				save_path <- file.path(dirpath, "phylo_tree.tre")
				ape::write.tree(tmp, file = save_path)
				message('Save phylogenetic tree to ', save_path, ' ...')
			}
			if(!is.null(self$rep_fasta)){
				tmp <- self$rep_fasta
				save_path <- file.path(dirpath, "rep_fasta.fasta")
				if(inherits(tmp, "list")){
					seqinr::write.fasta(tmp, names = names(tmp), file.out = save_path)
				}else{
					if(inherits(tmp, "DNAStringSet")){
						Biostrings::writeXStringSet(x = tmp, filepath = save_path)
					}else{
						if(inherits(tmp, "DNAbin")){
							ape::write.FASTA(tmp, file = save_path)
						}else{
							stop("Unknown fasta format! Must be one of DNAStringSet (from readDNAStringSet function of Biostrings package), ", 
								"list (from read.fasta function of seqinr package), and DNAbin (from read.FASTA function of ape package)!")
						}
					}
				}
				message('Save sequences to ', save_path, ' ...')
			}
		},
		#' @description
		#' Calculate the taxonomic abundance at each taxonomic level or selected levels.
		#'
		#' @param select_cols default NULL; numeric vector (column sequences) or character vector (column names of \code{microtable$tax_table}); 
		#'   applied to select columns to calculate abundances according to ordered hierarchical levels.
		#'   This parameter is very useful when only part of the columns are needed to calculate abundances.
		#' @param rel default TRUE; if TRUE, relative abundance is used; if FALSE, absolute abundance (i.e. raw values) will be summed.
		#' @param merge_by default "|"; the symbol to merge and concatenate taxonomic names of different levels.
		#' @param split_group default FALSE; if TRUE, split the rows to multiple rows according to one or more columns in \code{tax_table} 
		#'   when there is multiple mapping information.
		#' @param split_by default "&"; Separator delimiting collapsed values; only available when \code{split_group = TRUE}.
		#' @param split_column default NULL; one column name used for the splitting in tax_table for each abundance calculation; 
		#'   only available when \code{split_group = TRUE}. If not provided, the function will split each column that containing the \code{split_by} character.
		#' @param split_special_char default "&&"; special character that will be used forcibly to split multiple mapping information in \code{tax_table} by default
		#'   no matter \code{split_group} setting.
		#' @return \code{taxa_abund} list in object.
		#' @examples
		#' \donttest{
		#' m1$cal_abund()
		#' }
		cal_abund = function(
			select_cols = NULL, 
			rel = TRUE, 
			merge_by = "|",
			split_group = FALSE, 
			split_by = "&", 
			split_column = NULL,
			split_special_char = "&&"
			){
			taxa_abund <- list()
			if(is.null(self$tax_table)){
				stop("No tax_table found! Please check your data!")
			}
			if(nrow(self$tax_table) != nrow(self$otu_table)){
				message("The row number of tax_table is not equal to that of otu_table ...")
				message("Automatically applying tidy_dataset() function to trim the data ...")
				self$tidy_dataset()
				print(self)
			}
			if(nrow(self$sample_table) != ncol(self$otu_table)){
				message("The sample numbers of sample_table is not equal to that of otu_table ...")
				message("Automatically applying tidy_dataset() function to trim the data ...")
				self$tidy_dataset()
				print(self)
			}
			if(nrow(self$tax_table) == 0){
				stop("0 row in tax_table! Please check your data!")
			}
			if(is.null(select_cols)){
				select_cols <- seq_len(ncol(self$tax_table))
			}else{
				if(!is.numeric(select_cols)){
					if(any(! select_cols %in% colnames(self$tax_table))){
						stop("Part of input names of select_cols are not in the tax_table!")
					}else{
						select_cols <- match(select_cols, colnames(self$tax_table))
					}
				}
			}
			# built-in characters, such as those in MetaCyc mapping file
			if(any(grepl(split_special_char, self$tax_table[, select_cols], fixed = TRUE))){
				split_group <- TRUE
				split_by <- split_special_char
			}
			for(i in seq_along(select_cols)){
				tmp_mt <- clone(self)
				sel <- select_cols[1:i]
				taxrank <- colnames(self$tax_table)[select_cols[i]]
				
				test_doubleand <- lapply(1:i, function(x){any(grepl(split_special_char, self$tax_table[, select_cols[x]], fixed = TRUE))}) %>% unlist
				if(any(test_doubleand)){
					if(sum(test_doubleand) > 1){
						use_split_column <- taxrank
						sel <- select_cols[i]
					}else{
						use_split_column <- colnames(self$tax_table)[select_cols[1:i]]
					}
				}else{
					if(split_group){
						# assign the columns used for the splitting
						if(is.null(split_column)){
							use_split_column <- colnames(self$tax_table)[select_cols[1:i]]
						}else{
							use_split_column <- split_column
						}
					}
				}
				
				tmp_mt$tax_table %<>% .[, sel, drop = FALSE]
				taxa_abund[[taxrank]] <- private$transform_data_proportion(
											input = tmp_mt, 
											rel = rel, 
											merge_by = merge_by,
											split_group = split_group, 
											split_by = split_by, 
											split_column = use_split_column)
			}
			self$taxa_abund <- taxa_abund
			message('The result is stored in object$taxa_abund ...')
		},
		#' @description
		#' Save taxonomic abundance as local file.
		#'
		#' @param dirpath default "taxa_abund"; directory to save the taxonomic abundance files. It will be created if not found.
		#' @param merge_all default FALSE; Whether merge all tables into one. The merged file format is generally called 'mpa' style.
		#' @param rm_un default FALSE; Whether remove unclassified taxa in which the name ends with '__' generally.
		#' @param rm_pattern default "__$"; The pattern searched through the merged taxonomic names. See also \code{pattern} parameter in \code{\link{grepl}} function. 
		#' 	  Only available when \code{rm_un = TRUE}. The default "__$" means removing the names end with '__'.
		#' @param sep default ","; the field separator string. Same with \code{sep} parameter in \code{\link{write.table}} function.
		#' 	  default \code{','} correspond to the file name suffix 'csv'. The option \code{'\t'} correspond to the file name suffix 'tsv'. For other options, suffix are all 'txt'.
		#' @param ... parameters passed to \code{\link{write.table}}.
		#' @examples
		#' \dontrun{
		#' m1$save_abund(dirpath = "taxa_abund")
		#' m1$save_abund(merge_all = TRUE, rm_un = TRUE, sep = "\t")
		#' }
		save_abund = function(dirpath = "taxa_abund", merge_all = FALSE, rm_un = FALSE, rm_pattern = "__$", sep = ",", ...){
			if(!dir.exists(dirpath)){
				dir.create(dirpath)
			}
			suffix <- switch(sep, ',' = "csv", '\t' = "tsv", "txt")
			if(merge_all){
				res <- data.frame()
				for(i in names(self$taxa_abund)){
					res %<>% rbind(., self$taxa_abund[[i]])
				}
				res <- data.frame(Taxa = rownames(res), res)
				if(rm_un){
					res %<>% .[!grepl(rm_pattern, .$Taxa), ]
				}
				save_path <- paste0(dirpath, "/mpa_abund.", suffix)
				write.table(res, file = save_path, row.names = FALSE, sep = sep, ...)
				message('Save abundance to ', save_path, ' ...')
			}else{
				for(i in names(self$taxa_abund)){
					tmp <- self$taxa_abund[[i]]
					if(rm_un){
						tmp %<>% .[!grepl(rm_pattern, rownames(.)), ]
					}
					tmp <- data.frame(Taxa = rownames(tmp), tmp)
					save_path <- paste0(dirpath, "/", i, "_abund.", suffix)
					write.table(tmp, file = save_path, row.names = FALSE, sep = sep, ...)
				}
			}
		},
		#' @description
		#' Calculate alpha diversity.
		#'
		#' @param measures default NULL; one or more indexes in \code{c("Observed", "Coverage", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Pielou")}; 
		#'   The default NULL represents that all the measures are calculated. 'Shannon', 'Simpson' and 'InvSimpson' are calculated based on \code{vegan::diversity} function;
		#'   'Chao1' and 'ACE' depend on the function \code{vegan::estimateR}.
		#'   'Fisher' index relies on the function \code{vegan::fisher.alpha}.
		#'   "Observed" means the observed species number in a community, i.e. richness.
		#'   "Coverage" represents good's coverage. It is defined:
		#' 	   	     \deqn{Coverage = 1 - \frac{f1}{n}} 
		#'    where \emph{n} is the total abundance of a sample, and \emph{f1} is the number of singleton (species with abundance 1) in the sample.
		#'   "Pielou" denotes the Pielou evenness index. It is defined:
		#' 	   	     \deqn{J = \frac{H'}{\ln(S)}}
		#'    where \emph{H'} is Shannon index, and \emph{S} is the species number.
		#' @param PD default FALSE; whether Faith's phylogenetic diversity is calculated. The calculation depends on the function \code{picante::pd}.
		#'   Note that the phylogenetic tree (\code{phylo_tree} object in the data) is required for PD.
		#' @return alpha_diversity stored in the object. The se.chao1 and se.ACE are the standard erros of Chao1 and ACE, respectively.
		#' @examples
		#' \donttest{
		#' m1$cal_alphadiv(measures = NULL, PD = FALSE)
		#' class(m1$alpha_diversity)
		#' }
		cal_alphadiv = function(measures = NULL, PD = FALSE){
			# modified based on the alpha diversity analysis of phyloseq package
			OTU <- as.data.frame(t(self$otu_table), check.names = FALSE)
			renamevec    <-     c("Observed", "Coverage", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Pielou")
			names(renamevec) <- c("S.obs", "coverage", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher", "pielou")
			if(is.null(measures)){
				use_measures <- as.character(renamevec)
			}else{
				use_measures <- measures
			}
			if(any(use_measures %in% names(renamevec))){
				use_measures[use_measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% use_measures]
			}
			if(!any(use_measures %in% renamevec)){
				stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
			}
			outlist <- vector("list")
			estimRmeas <- c("Chao1", "Observed", "ACE")
			if(any(estimRmeas %in% use_measures)){
				est <- tryCatch(vegan::estimateR(OTU), error = function(e){c("Skip 'Chao1', 'ACE' and 'Observed' ...")})
				if(is.numeric(est)){
					outlist <- c(outlist, list(t(data.frame(est, check.names = FALSE))))
				}else{
					message(est)
				}
			}
			if("Shannon" %in% use_measures){
				outlist <- c(outlist, list(shannon = vegan::diversity(OTU, index = "shannon")))
			}
			if("Simpson" %in% use_measures){
				outlist <- c(outlist, list(simpson = vegan::diversity(OTU, index = "simpson")))
			}
			if("InvSimpson" %in% use_measures){
				outlist <- c(outlist, list(invsimpson = vegan::diversity(OTU, index = "invsimpson")))
			}
			if("Fisher" %in% use_measures){
				fisher <- tryCatch(vegan::fisher.alpha(OTU, se = TRUE), 
					warning = function(w){suppressWarnings(vegan::fisher.alpha(OTU, se = TRUE)[, c("alpha", "se")])},
					error = function(e){c("Skip the index Fisher because of an error ...")}
					)
				if(!is.null(dim(fisher))) {
					colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
					outlist <- c(outlist, list(fisher))
				}else{
					if(is.numeric(fisher)){
						outlist <- c(outlist, Fisher = list(fisher))
					}else{
						message(fisher)
					}
				}
			}
			if("Pielou" %in% use_measures){
				outlist <- c(outlist, list(pielou = vegan::diversity(OTU, index = "shannon")/log(vegan::specnumber(OTU))))
			}
			if("Coverage" %in% use_measures){
				outlist <- c(outlist, list(coverage = private$goods(OTU)))
			}
			if(PD){
				if(is.null(self$phylo_tree)){
					stop("Please provide phylogenetic tree for PD calculation!")
				}else{
					outlist <- c(outlist, list(PD = picante::pd(OTU, self$phylo_tree)[, "PD", drop=TRUE]))
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
			}
			write.csv(self$alpha_diversity, file = paste0(dirpath, "/", "alpha_diversity.csv"), row.names = TRUE)
		},
		#' @description
		#' Calculate beta diversity dissimilarity matrix, such as Bray-Curtis, Jaccard, and UniFrac.
		#' See An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and Lozupone et al. (2005) <doi:10.1128/AEM.71.12.8228–8235.2005>.
		#'
		#' @param method default NULL; a character vector with one or more elements; \code{c("bray", "jaccard")} is used when \code{method = NULL}; 
		#'   See the \code{method} parameter in \code{vegdist} function for more available options, such as 'aitchison' and 'robust.aitchison'. 
		#' @param unifrac default FALSE; whether UniFrac indexes (weighted and unweighted) are calculated. Phylogenetic tree is necessary when \code{unifrac = TRUE}.
		#' @param binary default FALSE; Whether convert abundance to binary data (presence/absence) when \code{method} is not "jaccard". 
		#'   TRUE is used for "jaccard" automatically.
		#' @param ... parameters passed to \code{vegdist} function of vegan package.
		#' @return beta_diversity list stored in the object.
		#' @examples
		#' \donttest{
		#' m1$cal_betadiv(unifrac = FALSE)
		#' class(m1$beta_diversity)
		#' }
		cal_betadiv = function(method = NULL, unifrac = FALSE, binary = FALSE, ...){
			res <- list()
			eco_table <- t(self$otu_table)
			sample_table <- self$sample_table
			if(is.null(method)){
				method <- c("bray", "jaccard")
			}
			vegdist_methods <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", "morisita", "horn", "mountford", 
				"jaccard", "raup", "binomial", "chao", "altGower", "cao", "mahalanobis", "clark", "chisq", "chord", "hellinger", 
				"aitchison", "robust.aitchison")
			for(i in method){
				i <- match.arg(i, vegdist_methods)
				if(i == "jaccard"){
					binary_use <- TRUE
				}else{
					binary_use <- binary
				}
				if(i == "aitchison"){
					if(any(eco_table == 0)){
						eco_table <- eco_table + 1
					}
				}
				res[[i]] <- as.matrix(vegan::vegdist(eco_table, method = i, binary = binary_use, ...))
			}
			if(unifrac){
				if(is.null(self$phylo_tree)){
					stop("No phylogenetic tree provided, please change the parameter unifrac to FALSE")
				}
				phylo_tree <- self$phylo_tree
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
			}
			for(i in names(self$beta_diversity)){
				write.csv(self$beta_diversity[[i]], file = paste0(dirpath, "/", i, ".csv"), row.names = TRUE)
			}
		},
		#' @description
		#' Print the microtable object.
		print = function(){
			cat("microtable-class object:\n")
			cat(paste("sample_table have", nrow(self$sample_table), "rows and", ncol(self$sample_table), "columns\n"))
			cat(paste("otu_table have", nrow(self$otu_table), "rows and", ncol(self$otu_table), "columns\n"))
			if(!is.null(self$tax_table)) cat(paste("tax_table have", nrow(self$tax_table), "rows and", ncol(self$tax_table), "columns\n"))
			if(!is.null(self$phylo_tree)) cat(paste("phylo_tree have", length(self$phylo_tree$tip.label), "tips\n"))
			if(!is.null(self$rep_fasta)) cat(paste("rep_fasta have", length(self$rep_fasta), "sequences\n"))
			if(!is.null(self$taxa_abund)) cat(paste("Taxa abundance: calculated for", paste0(names(self$taxa_abund), collapse = ","), "\n"))
			if(!is.null(self$alpha_diversity)) cat(paste("Alpha diversity: calculated for", paste0(colnames(self$alpha_diversity), collapse = ","), "\n"))
			if(!is.null(self$beta_diversity)) cat(paste("Beta diversity: calculated for", paste0(names(self$beta_diversity), collapse = ","), "\n"))
			invisible(self)
		}
		),
	private = list(
		# check and remove OTU or sample with 0 abundance
		check_abund_table = function(otu_table){
			if(!all(sapply(otu_table, is.numeric))){
				stop("Some columns in otu_table are not numeric class! Please check the input data!")
			}
			if(any(apply(otu_table, 1, sum) == 0)){
				remove_num <- sum(apply(otu_table, 1, sum) == 0)
				message(remove_num, " taxa with 0 abundance are removed from the otu_table ...")
				otu_table %<>% .[apply(., 1, sum) > 0, , drop = FALSE]
			}
			if(any(apply(otu_table, 2, sum) == 0)){
				remove_num <- sum(apply(otu_table, 2, sum) == 0)
				message(remove_num, " samples with 0 abundance are removed from the otu_table ...")
				otu_table %<>% .[, apply(., 2, sum) > 0, drop = FALSE]
			}
			if(ncol(otu_table) == 0){
				stop("No available sample! Please check the data!")
			}
			if(nrow(otu_table) == 0){
				stop("No available taxon! Please check the data!")
			}
			otu_table
		},
		tidy_samples = function(microtable_obj){
			microtable_obj$otu_table <- private$check_abund_table(microtable_obj$otu_table)
			sample_names <- intersect(rownames(microtable_obj$sample_table), colnames(microtable_obj$otu_table))
			if(length(sample_names) == 0){
				stop("No same sample name found between rownames of sample_table and colnames of otu_table! ",
					"Please first check whether the rownames of sample_table are sample names! Then check through the sample names of each table!")
			}
			# keep the sample order same with original sample table
			sample_names <- rownames(microtable_obj$sample_table) %>% .[. %in% sample_names]
			microtable_obj$sample_table %<>% .[sample_names, , drop = FALSE]
			microtable_obj$otu_table %<>% .[ , sample_names, drop = FALSE]
			microtable_obj
		},
		transform_data_proportion = function(
			input,
			rel,
			merge_by,
			split_group,
			split_by,
			split_column
			){
			sampleinfo <- input$sample_table
			abund <- input$otu_table
			tax <- input$tax_table
			# split rows to multiple rows for multiple mapping
			if(split_group){
				merge_abund <- cbind.data.frame(tax, abund)
				split_merge_abund <- tidyr::separate_longer_delim(merge_abund, cols = all_of(split_column), delim = split_by)
				new_tax <- split_merge_abund[, 1:ncol(tax), drop = FALSE]
				new_abund <- split_merge_abund[, (ncol(tax) + 1):(ncol(split_merge_abund)), drop = FALSE]
				abund1 <- cbind.data.frame(Display = apply(new_tax, 1, paste, collapse = merge_by), new_abund)
			}else{
				abund1 <- cbind.data.frame(Display = apply(tax, 1, paste, collapse = merge_by), abund)
			}
			# first convert table to long format
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
			abund2 <- as.data.frame(data.table::dcast(abund1, Display ~ Sample, value.var = list("res_abund"))) %>%
				`row.names<-`(.[,1]) %>% 
				.[,-1, drop = FALSE]
			abund2 <- abund2[order(apply(abund2, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]
			abund2
		},
		goods = function(com){
			no.seqs <- rowSums(com)
			sing <- com==1
			no.sing <- apply(sing, 1, sum)
			goods <- 1 - no.sing/no.seqs
			goods
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
