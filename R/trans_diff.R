#' @title 
#' Create trans_diff object for the differential analysis on the taxonomic abundance.
#'
#' @description
#' This class is a wrapper for a series of differential abundance test and indicator analysis methods, including 
#'  LEfSe based on the Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>,
#'  random forest <doi:10.1016/j.geoderma.2018.09.035>, metastat based on White et al. (2009) <doi:10.1371/journal.pcbi.1000352>,
#'  the method in R package metagenomeSeq Paulson et al. (2013) <doi:10.1038/nmeth.2658>, non-parametric Kruskal-Wallis Rank Sum Test,
#'  Dunn's Kruskal-Wallis Multiple Comparisons based on the FSA package, Wilcoxon Rank Sum and Signed Rank Tests, t test and anova.
#' 
#' Authors: Chi Liu, Yang Cao, Chenhao Li
#' 
#' @export
trans_diff <- R6Class(classname = "trans_diff",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param method default "lefse"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'lefse'}}{LEfSe method based on Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>}
		#'     \item{\strong{'rf'}}{random forest and non-parametric test method based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>}
		#'     \item{\strong{'metastat'}}{Metastat method for all paired groups based on White et al. (2009) <doi:10.1371/journal.pcbi.1000352>}
		#'     \item{\strong{'metagenomeSeq'}}{zero-inflated log-normal model-based differential test method from metagenomeSeq package.}
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons when group number > 2; see dunnTest function in FSA package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'anova'}}{Duncan's multiple range test for anova}
		#'     \item{\strong{'ANCOMBC'}}{Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC) <doi:10.1038/s41467-020-17041-7>}
		#'   }
		#' @param group default NULL; sample group used for the comparision; a colname of microtable$sample_table.
		#' @param taxa_level default "all"; 'all' represents using abundance data at all taxonomic ranks; 
		#' 	  For testing at a specific rank, provide taxonomic rank name, such as "Genus".
		#' 	  If the provided taxonomic name is neither 'all' nor a colname in tax_table of dataset, the function will use the features in otu_table automatically.
		#' @param filter_thres default 0; the relative abundance threshold, such as 0.0005; only useful when method != "metastat", "metagenomeSeq" or "ANCOMBC". 
		#' @param alpha default 0.05; differential significance threshold for method = "lefse" or "rf"; used to select taxa with significance across groups.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of p.adjust function for other available options; 
		#'    NULL mean disuse the p value adjustment; So when p_adjust_method = NULL, P.adj is same with P.unadj.
		#' @param lefse_subgroup default NULL; sample sub group used for sub-comparision in lefse; Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>.
		#' @param lefse_min_subsam default 10; sample numbers required in the subgroup test.
		#' @param lefse_norm default 1000000; scale value in lefse.
		#' @param nresam default 0.6667; sample number ratio used in each bootstrap for method = "lefse" or "rf".
		#' @param boots default 30; bootstrap test number for method = "lefse" or "rf".
		#' @param rf_ntree default 1000; see ntree in randomForest function of randomForest package when method = "rf".
		#' @param group_choose_paired default NULL; a vector used for selecting the required groups for paired testing, only used for method = "metastat" or "metagenomeSeq".
		#' @param metagenomeSeq_count default 1; Filter features to have at least 'counts' counts.; see the count parameter in MRcoefs function of metagenomeSeq package.
		#' @param ANCOMBC_formula default NULL; same with the formula parameter in ANCOMBC::ancombc function; If NULL; assign the group parameter to it automatically.
		#' @param ... parameters passed to cal_diff function of trans_alpha class when method is one of "KW", "KW_dunn", "wilcox", "t.test" and "anova";
		#' 	 passed to ANCOMBC::ancombc function when method is "ANCOMBC" (except formula and global parameters; please see ANCOMBC_formula parameter).
		#' @return res_diff and res_abund.\cr
		#'   \strong{res_abund} includes mean abudance of each taxa (Mean), standard deviation (SD), standard error (SE) and sample number (N) in the group (Group).\cr
		#'   \strong{res_diff} is the detailed differential test result, containing:\cr
		#'     \strong{"Comparison"}: The groups for the comparision, maybe all groups or paired groups. If this column is not found, all groups used;\cr
		#'     \strong{"Group"}: Which group has the maximum median or mean value across the test groups; 
		#'        For non-parametric methods, median value; For t.test, mean value;\cr
		#'     \strong{"Taxa"}: which taxa is used in this comparision;\cr
		#'     \strong{"Method"}: Test method used in the analysis depending on the method input;\cr
		#'     \strong{"LDA" or "MeanDecreaseGini"}: LDA: linear discriminant score in LEfSe; MeanDecreaseGini: mean decreasing gini index in random forest;\cr
		#'     \strong{"P.unadj" and "P.adj"}: raw p value; P.adj: adjusted p value;\cr
		#'     \strong{"qvalue"}: qvalue for metastat analysis.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group")
		#' t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group")
		#' t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Group", taxa_level = "Genus")
		#' t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group")
		#' }
		initialize = function(
			dataset = NULL,
			method = c("lefse", "rf", "metastat", "metagenomeSeq", "KW", "KW_dunn", "wilcox", "t.test", "anova", "ANCOMBC")[1],
			group = NULL,
			taxa_level = "all",
			filter_thres = 0,
			alpha = 0.05,
			p_adjust_method = "fdr",
			lefse_subgroup = NULL,
			lefse_min_subsam = 10,
			lefse_norm = 1000000,
			nresam = 0.6667,
			boots = 30,
			rf_ntree = 1000,
			group_choose_paired = NULL,
			metagenomeSeq_count = 1,
			ANCOMBC_formula = NULL,
			...
			){
			if(is.null(dataset)){
				stop("No dataset provided!")
			}
			# in case of dataset changes
			tmp_dataset <- clone(dataset)
			sampleinfo <- tmp_dataset$sample_table
			if(is.null(group)){
				stop("The group parameter need to be provided for differential test among groups!")
			}else{
				if(length(group) > 1){
					stop("Please provide only one colname of sample_table for group parameter!")
				}else{
					if(! group %in% colnames(sampleinfo)){
						stop("Please provide a correct colname of sample_table for group parameter!")
					}
				}
			}
			method <- match.arg(method, c("lefse", "rf", "metastat", "metagenomeSeq", "KW", "KW_dunn", "wilcox", "t.test", "anova", "ANCOMBC"))
			
			if(is.factor(sampleinfo[, group])){
				self$group_order <- levels(sampleinfo[, group])
				sampleinfo[, group] %<>% as.character
			}
			
			################################
			# generate abudance table
			if(is.null(tmp_dataset$taxa_abund)){
				message("No taxa_abund found. First calculate it with cal_abund function ...")
				tmp_dataset$cal_abund()
			}
			# make sure the taxa_level can be found from tax_table
			if(grepl("all", taxa_level, ignore.case = TRUE)){
				abund_table <- do.call(rbind, unname(tmp_dataset$taxa_abund))
			}else{
				# make sure taxa_level can be extracted from taxa_abund
				if(! taxa_level %in% names(tmp_dataset$taxa_abund)){
					# recalculate taxa_abund with rownames as features in otu_table
					message("Provided taxa_level: ", taxa_level, " not in tax_table of dataset; use features in otu_table ...")
					tmp_dataset$add_rownames2taxonomy(use_name = taxa_level)
					suppressMessages(tmp_dataset$cal_abund(rel = TRUE))
				}
				abund_table <- tmp_dataset$taxa_abund[[taxa_level]]
			}
			if(filter_thres > 0){
				if(filter_thres >= 1){
					stop("Parameter filter_thres represents relative abudance. It should be smaller than 1 !")
				}else{
					abund_table %<>% .[apply(., 1, mean) > filter_thres, ]
				}
			}
			if(grepl("lefse", method, ignore.case = TRUE)){
				abund_table %<>% {. * lefse_norm}
				self$lefse_norm <- lefse_norm
			}
			abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]}
			################################
			
			if(method %in% c("KW", "KW_dunn", "wilcox", "t.test", "anova")){
				tem_data <- clone(tmp_dataset)
				# use test method in trans_alpha
				tem_data$alpha_diversity <- as.data.frame(t(abund_table))
				tem_data1 <- suppressMessages(trans_alpha$new(dataset = tem_data, group = group))
				suppressMessages(tem_data1$cal_diff(method = method, p_adjust_method = p_adjust_method, ...))
				output <- tem_data1$res_diff
				if(method != "anova"){
					colnames(output)[colnames(output) == "Measure"] <- "Taxa"
				}else{
					output <- rownames_to_column(output, var = "Group")
					output <- reshape2::melt(output, id.vars = "Group", variable.name = "Taxa", value.name = "Significance")
				}
			}
			if(method %in% c("lefse", "rf")){
				# differential test
				group_vec <- sampleinfo[, group] %>% as.factor
				comparisions <- paste0(levels(group_vec), collapse = " - ")

				message("Start Kruskal-Wallis rank sum test for ", group, " ...")
				res_class <- suppressWarnings(lapply(seq_len(nrow(abund_table)), function(x) private$test_mark(abund_table[x, ], group_vec, method = "kruskal.test")))
				
				pvalue_raw <- unlist(lapply(res_class, function(x) x$p_value))
				names(pvalue_raw) <- rownames(abund_table)
				pvalue_raw[is.nan(pvalue_raw)] <- 1
				message(sum(pvalue_raw < alpha), " taxa found significant ...")
				if(is.null(p_adjust_method)){
					pvalue <- pvalue_raw
				}else{
					pvalue <- p.adjust(pvalue_raw, method = p_adjust_method)
				}
				# select significant taxa
				sel_taxa <- pvalue < alpha
				message("After P value adjustment, ", sum(sel_taxa), " taxa found significant ...")
				if(sum(sel_taxa) == 0){
					stop("No significant taxa found! Stop running!")
				}

				abund_table_sub <- abund_table[sel_taxa, ]
				pvalue_sub <- pvalue[sel_taxa]
				class_taxa_median_sub <- lapply(res_class, function(x) x$med) %>% do.call(cbind, .) %>% .[, sel_taxa]
			}
			if(method == "rf"){
				# change the name in case of additional problem from the taxonomic names
				nametable <- cbind.data.frame(name = rownames(abund_table_sub), repl = paste0("t", 1:nrow(abund_table_sub)), stringsAsFactors = FALSE)
				rownames(nametable) <- nametable$repl
				predictors <- t(abund_table_sub)
				colnames(predictors) <- nametable$repl
				res <- NULL
				if(all(table(as.character(sampleinfo[, group])) < 4)){
					if(nresam < 1){
						message("The sample number in all your group less than 4, automatically set nresam = 1 for random forest analysis !")
						nresam <- 1
					}
				}
				for(num in seq_len(boots)){
					# resampling
					sample_names_resample <- rownames(predictors)[base::sample(1:nrow(predictors), size = ceiling(nrow(predictors) * nresam))]
					predictors_sub <- predictors[sample_names_resample, ]
					sampleinfo_resample <- sampleinfo[sample_names_resample, , drop = FALSE]
					# make sure the groups and samples numbers right
					if(length(unique(sampleinfo_resample[, group])) != length(unique(sampleinfo[, group])) | min(table(sampleinfo_resample[, group])) < 2){
						next
					}
					rf_data <- data.frame(response = as.factor(sampleinfo_resample[, group]), predictors_sub, stringsAsFactors = FALSE)
					tem_classify <- randomForest::randomForest(response~., data = rf_data, ntree = rf_ntree)
					# use importance to evaluate
					imp <- randomForest::importance(tem_classify)
					colnames(imp)[1] <- num
					if(is.null(res)){
						res <- imp
					}else{
						res <- cbind(res, imp[rownames(res), , drop = FALSE])
					}
				}
				res <- apply(res, 1, mean) %>% as.data.frame
				use_method <- ifelse(length(levels(group_vec)) > 2, "Kruskal-Wallis rank sum test & Random Forest", "Wilcoxon rank sum test & Random Forest")
				Taxa_name <- nametable[rownames(res), "name"]
				res <- data.frame(
					Comparison = comparisions,
					Taxa = Taxa_name,
					Group = apply(class_taxa_median_sub, 2, function(x) rownames(class_taxa_median_sub)[which.max(x)])[Taxa_name], 
					Method = use_method,
					MeanDecreaseGini = res[, 1], 
					stringsAsFactors = FALSE)
				output <- dplyr::arrange(res, dplyr::desc(MeanDecreaseGini))
				rownames(output) <- output$Taxa
				output$P.unadj <- pvalue_raw[as.character(output$Taxa)]
				output$P.adj <- pvalue_sub[as.character(output$Taxa)]
				output$Significance <- cut(output$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
			}
			if(method == "lefse"){
				all_class_pairs <- combn(unique(as.character(group_vec)), 2)
				# check the difference among subgroups
				if(!is.null(lefse_subgroup)){
					message("Start lefse subgroup biomarkers check for ", lefse_subgroup, " ...")
					all_sub_number <- as.data.table(sampleinfo)[, .N, by = c(group, lefse_subgroup)] %>% as.data.frame %>% .$N
					if(all(all_sub_number < lefse_min_subsam)){
						stop("All sample numbers for subgroups < ", lefse_min_subsam, "! Please consider using small lefse_min_subsam parameter!")
					}
					remove_list_total <- list()
					# for each group paires
					for(i in seq_len(ncol(all_class_pairs))){
						y1 <- all_class_pairs[, i]
						y1_sub_pairs <- expand.grid(
							unique(sampleinfo[sampleinfo[, group] == y1[1], lefse_subgroup]), 
							unique(sampleinfo[sampleinfo[, group] == y1[2], lefse_subgroup]), 
							stringsAsFactors = FALSE
							) %>% t
						y1_sub_pairs <- y1_sub_pairs[, unlist(lapply(1:ncol(y1_sub_pairs), function(x){
							ifelse(any(c(sum(sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[1, x]) < lefse_min_subsam, 
								sum(sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[2, x]) < lefse_min_subsam)), FALSE, TRUE)
						})), drop = FALSE]
						if(ncol(y1_sub_pairs) == 0) next
						res_sub_total <- list()
						# check each subgroup pairs under fixed group pair condition
						for(j in 1:ncol(y1_sub_pairs)){
							y2 <- y1_sub_pairs[, j]
							abund_table_sub_y2 <- abund_table_sub[, c(
								rownames(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], ]), 
								rownames(sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], ])
								)]
							group_vec_sub2 <- c(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], group], 
								sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], group])
							res_sub <- lapply(seq_len(nrow(abund_table_sub_y2)), function(x) private$test_mark(abund_table_sub_y2[x,], group_vec_sub2))
							res_sub_total[[j]] <- res_sub
						}
						raw_median <- class_taxa_median_sub[y1, ] %>% 
							{.[1, ] > .[2, ]} %>% 
							as.vector
						check_median_sub <- sapply(res_sub_total, function(x) unlist(lapply(x, function(y) {y$med[y1, 1] %>% {.[1] > .[2]}}))) %>% as.data.frame
						check_median_sub[] <- lapply(check_median_sub, function(x) x == raw_median)
						check_p_sub <- sapply(res_sub_total, function(x) unlist(lapply(x, function(y) y$p_value))) %>% as.data.frame
						remove_list <- unlist(lapply(seq_len(nrow(check_median_sub)), function(x){
							if(all(unlist(check_median_sub[x, ]))){
								FALSE
							}else{
								if(any(check_p_sub[x, !unlist(check_median_sub[x, ])] < alpha)){
									TRUE
								}else{
									FALSE
								}
							}
						}))
						remove_list_total[[i]] <- remove_list
					}
					if(!identical(remove_list_total, list())){
						remove_list_total %<>% do.call(cbind, .) %>% apply(., 1, any)
						message("Remove ", sum(remove_list_total), " biomarkers after subgroup check ...")
						abund_table_sub %<>% .[!remove_list_total, ]
						if(nrow(abund_table_sub) == 0){
							stop("No biomarkers remained after subgroup check! stop running!")
						}
						pvalue_sub %<>% .[!remove_list_total]
						class_taxa_median_sub %<>% .[, !remove_list_total]
					}
				}
				res_lda <- list()
				# bootstrap default 30 times
				for(num in seq_len(boots)){
					res_lda_pair <- list()
					# resampling samples
					sample_names_resample <- colnames(abund_table_sub)[base::sample(1:ncol(abund_table_sub), size = ceiling(ncol(abund_table_sub) * nresam))]
					abund_table_sub_resample <- abund_table_sub[, sample_names_resample]
					sampleinfo_resample <- sampleinfo[sample_names_resample, , drop = FALSE]
					# make sure the groups and samples numbers available					
					if(sum(table(as.character(sampleinfo_resample[, group])) > 1) < 2){
						res_lda[[num]] <- NA
						next
					}
					# cycle all paired groups
					for(i in seq_len(ncol(all_class_pairs))){
						sel_samples <- sampleinfo_resample[, group] %in% all_class_pairs[, i]
						# make sure the groups and samples numbers available
						if(length(table(as.character(sampleinfo_resample[sel_samples, group]))) < 2){
							res_lda_pair[[i]] <- NA
							next
						}
						if(min(table(as.character(sampleinfo_resample[sel_samples, group]))) < 2){
							res_lda_pair[[i]] <- NA
							next
						}
						group_vec_lda <- sampleinfo_resample[sel_samples, group] %>% as.character %>% as.factor
						abund_table_sub_lda <- abund_table_sub_resample[, sel_samples]
						abund_table_sub_lda %<>% .[apply(., 1, sd) > 1.0e-10, ]
						if(is.null(lefse_subgroup)){
							abund1 <- cbind.data.frame(t(abund_table_sub_lda), Group = group_vec_lda)
						}else{
							subgroup_vec <- sampleinfo_resample[sel_samples, lefse_subgroup] %>% as.character %>% as.factor
							# consider subgroup as a independent variable
							abund1 <- cbind.data.frame(t(abund_table_sub_lda), Group = group_vec_lda, lefse_subgroup = subgroup_vec)
						}
						# LDA analysis
						check_res <- tryCatch(mod1 <- MASS::lda(Group ~ ., abund1, tol = 1.0e-10), error = function(e) { skip_to_next <- TRUE})
						if(rlang::is_true(check_res)) {
							res_lda_pair[[i]] <- NA
							next
						}else{
							# calculate effect size
							w <- mod1$scaling[,1]
							w_unit <- w/sqrt(sum(w^2))
							w_unit %<>% {.[!grepl("lefse_subgroup", names(.))]}
							ss <- abund1[, !colnames(abund1) %in% c("Group", "lefse_subgroup")]
							xy_matrix <- as.matrix(ss)
							LD <- xy_matrix %*% w_unit
							effect_size <- tapply(LD, group_vec_lda, mean) %>% as.vector %>% {.[1] - .[2]} %>% abs
							coeff <- w_unit * effect_size %>% abs
							coeff[is.nan(coeff)] <- 0
							names(coeff) %<>% gsub("^`|`$", "", .)
							rres <- mod1$means %>% as.data.frame
							colnames(rres) %<>% gsub("^`|`$", "", .)
							rres <- rres[, rownames(abund_table_sub_lda)]
							rres1 <- apply(rres, 2, function(x) abs(x[1] - x[2]))
							res_lda_pair[[i]] <- (rres1 + coeff[names(rres1)]) *0.5
						}
					}
					res_lda[[num]] <- res_lda_pair
				}
				# obtain the final lda value
				res <- sapply(rownames(abund_table_sub), function(k){
					unlist(lapply(seq_len(ncol(all_class_pairs)), function(p){
						unlist(lapply(res_lda, function(x){ x[[p]][k]})) %>% .[!is.na(.)] %>% mean
					})) %>% 
					.[!is.na(.)] %>% 
					.[!is.nan(.)] %>% 
					max
				})
				res <- sapply(res, function(x) {log10(1 + abs(x)) * ifelse(x > 0, 1, -1)})
				output <- cbind.data.frame(Group = apply(class_taxa_median_sub, 2, function(x) rownames(class_taxa_median_sub)[which.max(x)]), 
					LDA = res,
					P.unadj = pvalue_raw[names(pvalue_sub)],
					P.adj = pvalue_sub)
				output %<>% .[order(.$LDA, decreasing = TRUE), ]
				output <- cbind.data.frame(Comparison = comparisions, Taxa = rownames(output), Method = "LEfSe", output)
				message("Minimum LDA score: ", range(output$LDA)[1], " maximum LDA score: ", range(output$LDA)[2])
				output$Significance <- cut(output$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
			}
			#######################################
			# output taxonomic abundance mean and sd for the final res_abund and enrich group finding in metagenomeSeq or ANCOMBC
			if(grepl("lefse", method, ignore.case = TRUE)){
				res_abund <- reshape2::melt(rownames_to_column(abund_table_sub/lefse_norm, "Taxa"), id.vars = "Taxa")
			}else{
				if(grepl("rf", method, ignore.case = TRUE)){
					res_abund <- reshape2::melt(rownames_to_column(abund_table_sub, "Taxa"), id.vars = "Taxa")
				}else{
					res_abund <- reshape2::melt(rownames_to_column(abund_table, "Taxa"), id.vars = "Taxa")
				}
			}
			colnames(res_abund) <- c("Taxa", "Sample", "Abund")
			res_abund <- suppressWarnings(dplyr::left_join(res_abund, rownames_to_column(sampleinfo), by = c("Sample" = "rowname")))
			res_abund <- microeco:::summarySE_inter(res_abund, measurevar = "Abund", groupvars = c("Taxa", group))
			colnames(res_abund)[colnames(res_abund) == group] <- "Group"
			#######################################
			if(method %in% c("metastat", "metagenomeSeq", "ANCOMBC")){
				if(taxa_level == "all"){
					stop("Please provide the taxa_level instead of 'all', such as 'Genus' !")
				}
				if(is.null(group_choose_paired)){
					all_name <- combn(unique(as.character(sampleinfo[, group])), 2)
				}else{
					all_name <- combn(unique(group_choose_paired), 2)
				}
				output <- data.frame()			
			}
			if(method == "metastat"){
				# transform data
				ranknumber <- which(colnames(tmp_dataset$tax_table) %in% taxa_level)
				abund <- tmp_dataset$otu_table
				tax <- tmp_dataset$tax_table[, 1:ranknumber, drop=FALSE]
				merged_taxonomy <- apply(tax, 1, paste, collapse="|")
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

				message("Total ", ncol(all_name), " paired group for calculation ...")
				for(i in 1:ncol(all_name)) {
					message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ..."))
					use_data <- new_abund[ , unlist(lapply(as.character(all_name[,i]), function(x) which(as.character(sampleinfo[, group]) %in% x)))]
					use_data %<>% .[!grepl("__$", rownames(.)), ]
					use_data <- use_data[apply(use_data, 1, sum) != 0, ]
					g <- sum(as.character(sampleinfo[, group]) == as.character(all_name[1, i])) + 1
					# calculate metastat
					res <- private$calculate_metastat(inputdata = use_data, g = g)
					add_name <- paste0(as.character(all_name[, i]), collapse = " - ") %>% rep(., nrow(res))
					res <- cbind.data.frame(compare = add_name, res)
					output <- rbind.data.frame(output, res)
				}
				output %<>% dropallfactors(unfac2num = TRUE)
				colnames(output)[1:2] <- c("Comparison", "Taxa")
				max_group <- lapply(seq_len(nrow(output)), function(x){
					group_select <- unlist(strsplit(output[x, "Comparison"], split = " - "))
					ifelse(output[x, "mean(group1)"] > output[x, "mean(group2)"], group_select[1], group_select[2])
				}) %>% unlist
				output$Significance <- cut(output$qvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
				output <- data.frame(output, Group = max_group)
			}
			if(method == "metagenomeSeq"){
				if(!require("metagenomeSeq")){
					stop("metagenomeSeq package not installed !")
				}
				message("Total ", ncol(all_name), " paired group for calculation ...")
				for(i in 1:ncol(all_name)) {
					message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ...\n"))
					use_dataset <- clone(tmp_dataset)
					use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[,i]), ]
					use_dataset$tidy_dataset()
					suppressMessages(use_dataset$cal_abund(rel = FALSE))
					newdata <- microtable$new(otu_table = use_dataset$taxa_abund[[taxa_level]], 
						sample_table = use_dataset$sample_table)
					newdata$tidy_dataset()
					obj <- newMRexperiment(
						newdata$otu_table, 
						phenoData= AnnotatedDataFrame(newdata$sample_table)
#						featureData = AnnotatedDataFrame(use_dataset$tax_table)
						)
					## Normalization and Statistical testing
					obj_1 <- cumNorm(obj)
					pd <- pData(obj)
					colnames(pd)[which(colnames(pd) == group)] <- "Group"
					# construct linear model
					mod <- model.matrix(~1 + Group, data = pd)
					objres1 <- fitFeatureModel(obj_1, mod)
					# extract the result
					tb <- data.frame(logFC = objres1@fitZeroLogNormal$logFC, se = objres1@fitZeroLogNormal$se)
					p <- objres1@pvalues
					if(is.null(p_adjust_method)){
						padj <- p
					}else{
						if (p_adjust_method == "ihw-ubiquity" | p_adjust_method == "ihw-abundance") {
							padj <- MRihw(objres1, p, p_adjust_method, 0.1)
						} else {
							padj <- p.adjust(p, method = p_adjust_method)
						}
					}
					srt <- order(p, decreasing = FALSE)
					valid <- 1:length(padj)
					if (metagenomeSeq_count > 0) {
						np <- rowSums(objres1@counts)
						valid <- intersect(valid, which(np >= metagenomeSeq_count))
					}
					srt <- srt[which(srt %in% valid)]
					res <- cbind(tb[, 1:2], p)
					res <- cbind(res, padj)
					res <- as.data.frame(res[srt, ])
					colnames(res) <- c(colnames(tb)[1:2], "P.unadj", "P.adj")
					res <- cbind.data.frame(feature = rownames(res), res)
					rownames(res) <- NULL
					add_name <- paste0(as.character(all_name[, i]), collapse = " - ") %>% rep(., nrow(res))
					res <- cbind.data.frame(compare = add_name, res)
					output <- rbind.data.frame(output, res)
				}
			}
			if(method == "ANCOMBC"){
				if(!require("ANCOMBC")){
					stop("ANCOMBC package is not installed !")
				}
				if(!require("file2meco")){
					stop("Please install file2meco package! The function meco2phyloseq is required!")
				}
				if(is.null(ANCOMBC_formula)){
					ANCOMBC_formula <- group
				}
				if(ncol(all_name) > 1){
					message("First run the global test ...")
					use_dataset <- clone(tmp_dataset)
					suppressMessages(use_dataset$cal_abund(rel = FALSE))
					newdata <- microtable$new(otu_table = use_dataset$taxa_abund[[taxa_level]], sample_table = use_dataset$sample_table)
					newdata$tidy_dataset()
					newdata <- file2meco::meco2phyloseq(newdata)
					res_raw <- ancombc(phyloseq = newdata, group = group, formula = ANCOMBC_formula, global = TRUE, ...)
					res <- res_raw$res_global
					colnames(res) <- c("W", "P.unadj", "P.adj", "diff_abn")
					res <- cbind.data.frame(feature = rownames(res), res)
					group_vec <- use_dataset$sample_table[, group] %>% as.character %>% unique
					comparisions <- paste0(group_vec, collapse = " - ")
					res <- cbind.data.frame(compare = comparisions, res)
					output <- rbind.data.frame(output, res)
				}
				message("Total ", ncol(all_name), " paired group for test ...")
				for(i in 1:ncol(all_name)) {
					message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ...\n"))
					use_dataset <- clone(tmp_dataset)
					use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[,i]), ]
					use_dataset$tidy_dataset()
					suppressMessages(use_dataset$cal_abund(rel = FALSE))
					newdata <- microtable$new(otu_table = use_dataset$taxa_abund[[taxa_level]], sample_table = use_dataset$sample_table)
					newdata$tidy_dataset()
					newdata <- file2meco::meco2phyloseq(newdata)
					res_raw <- ancombc(phyloseq = newdata, group = group, formula = ANCOMBC_formula, ...)
					res_raw2 <- res_raw$res
					res <- data.frame(feature = rownames(res_raw2$W), W = res_raw2$W[, 1], P.unadj = res_raw2$p_val[, 1], 
						P.adj = res_raw2$q_val[, 1], diff_abn = res_raw2$diff_abn[, 1])
					rownames(res) <- NULL
					add_name <- paste0(as.character(all_name[, i]), collapse = " - ") %>% rep(., nrow(res))
					res <- cbind.data.frame(compare = add_name, res)
					output <- rbind.data.frame(output, res)
				}
			}
			if(method %in% c("metagenomeSeq", "ANCOMBC")){
				output %<>% dropallfactors(unfac2num = TRUE)
				colnames(output)[1:2] <- c("Comparison", "Taxa")
				output$Significance <- cut(output$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
				# filter the unknown taxa in output
				output %<>% .[.$Taxa %in% res_abund$Taxa, ]
				# get enriched group
				output$Group <- lapply(seq_along(output$Taxa), function(x){
					select_group_split <- strsplit(output[x, "Comparison"], split = " - ") %>% unlist
					res_abund[res_abund$Taxa == output[x, "Taxa"] & res_abund$Group %in% select_group_split, ] %>%
					{.[which.max(.$Mean), "Group"]}
				}) %>% unlist
			}
			self$res_diff <- output
			message(method , " analysis result is stored in object$res_diff ...")
			self$res_abund <- res_abund
			message('Taxa abundance data is stored in object$res_abund ...')
			# save abund_table for the cladogram
			self$abund_table <- abund_table
			self$method <- method
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Plotting the abundance of differential taxa.
		#'
		#' @param use_number default 1:20; numeric vector; the taxa numbers (1:n) used in the plot; 
		#'   If the n is larger than the number of total significant taxa, automatically use all the taxa.		
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors palette.
		#' @param select_group default NULL; this is used to select the paired groups. 
		#'   This parameter is especially useful when the comparision methods is applied to paired groups;
		#'   The input select_group must be one of object$res_diff$Comparison.
		#' @param select_taxa default NULL; character vector to provide taxa names. 
		#' 	 The taxa names should be same with the names shown in the plot, not the 'Taxa' column names in object$res_diff$Taxa.
		#' @param simplify_names default TRUE; whether use the simplified taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param group_order default NULL; a vector to order groups, i.e. reorder the legend and colors in plot; 
		#' 	  If NULL, the function can first check whether the group column of sample_table is factor. If yes, use the levels in it.
		#' 	  If provided, overlook the levels in the group of sample_table.
		#' @param barwidth default 0.9; the bar width in plot.
		#' @param use_se default TRUE; whether use SE in plot, if FALSE, use SD.
		#' @param add_sig default FALSE; whether add the significance label to the plot.
		#' @param add_sig_label default "Significance"; select a colname of object$res_diff for the label text, such as 'P.adj' or 'Significance'.
		#' @param add_sig_label_color default "black"; the color for the label text when add_sig = TRUE.
		#' @param add_sig_tip_length default 0.01; the tip length for the added line when add_sig = TRUE.
		#' @param y_start default 1.01; the y axis position from which to add the label; the default 1.01 means 1.01 * Value;
		#'   For method != "anova", all the start positions are same, i.e. Value = max(Mean+SD or Mean+SE); 
		#'   For method = "anova"; the stat position is calculated for each point, i.e. Value = Mean+SD or Mean+SE.
		#' @param y_increase default 0.05; the increasing y axia space to add label for paired groups; the default 0.05 means 0.05 * y_start * Value; 
		#' 	  In addition, this parameter is also used to label the letters of anova result with the fixed (1 + y_increase) * y_start * Value.
		#' @param text_y_size default 10; the size for the y axis text.
		#' @param coord_flip default TRUE; whether flip cartesian coordinates so that horizontal becomes vertical, and vertical, horizontal.
		#' @param ... parameters passed to ggsignif::stat_signif when add_sig = TRUE.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Genus")
		#' t1$plot_diff_abund(use_number = 1:10)
		#' t1$plot_diff_abund(use_number = 1:10, add_sig = TRUE)
		#' t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group")
		#' t1$plot_diff_abund(use_number = 1:20)
		#' t1$plot_diff_abund(use_number = 1:20, add_sig = TRUE)
		#' t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group")
		#' t1$plot_diff_abund(use_number = 1:20)
		#' t1$plot_diff_abund(use_number = 1:20, add_sig = TRUE)
		#' }
		plot_diff_abund = function(
			use_number = 1:20,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			select_group = NULL,
			select_taxa = NULL,
			simplify_names = TRUE,
			keep_prefix = TRUE,
			group_order = NULL,
			barwidth = 0.9,
			use_se = TRUE,
			add_sig = FALSE,
			add_sig_label = "Significance",
			add_sig_label_color = "black",
			add_sig_tip_length = 0.01,
			y_start = 1.01,
			y_increase = 0.05,
			text_y_size = 10,
			coord_flip = TRUE,
			...
			){
			abund_data <- self$res_abund
			method <- self$method
			diff_data <- self$res_diff

			# sort according to different columns
			if(method == "metastat"){
				message('Reorder taxa according to qvalue in res_diff from low to high ...')
				diff_data %<>% .[order(.$qvalue, decreasing = FALSE), ]
				# diff_data %<>% .[.$qvalue < 0.05, ]
			}else{
				# lefse and rf are ordered
				if(! method %in% c("lefse", "rf", "anova")){
					message('Reorder taxa according to P.adj in res_diff from low to high ...')
					diff_data %<>% .[order(.$P.adj, decreasing = FALSE), ]
					# diff_data %<>% .[.$P.adj < 0.05, ]
				}
			}

			if(!is.null(select_group)){
				if(length(select_group) > 1){
					stop("The select_group parameter should only have one element! Please check the input!")
				}
				if(! select_group %in% diff_data$Comparison){
					stop("The select_group parameter must be one of elements of object$res_diff$Comparison!")
				}
				diff_data %<>% .[.$Comparison %in% select_group, ]
				select_group_split <- strsplit(select_group, split = " - ") %>% unlist
				abund_data %<>% .[.$Group %in% select_group_split, ]
			}
			if(nrow(diff_data) == 0){
				stop("No significant taxa can be used to plot the abudance!")
			}
			if(simplify_names == T){
				diff_data$Taxa %<>% gsub(".*\\|", "", .)
				abund_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				diff_data$Taxa %<>% gsub(".__", "", .)
				abund_data$Taxa %<>% gsub(".__", "", .)
			}
			if(is.null(select_taxa)){
				if(length(use_number) > length(unique(as.character(diff_data$Taxa)))){
					message("The length of use_number is larger than taxa number in object$diff_data. Use all taxa ...")
					use_number <- 1:length(unique(as.character(diff_data$Taxa)))
				}
				diff_data %<>% .[.$Taxa %in% unique(as.character(diff_data$Taxa))[use_number], ]
				diff_data$Taxa %<>% factor(., levels = rev(unique(as.character(.))))
			}else{
				diff_data %<>% .[.$Taxa %in% select_taxa, ]
				if(nrow(diff_data) == 0){
					stop("No significant taxa can be used to plot the abudance!")
				}
				diff_data$Taxa %<>% factor(., levels = rev(select_taxa))
			}
			abund_data %<>% .[.$Taxa %in% levels(diff_data$Taxa), ]
			abund_data$Taxa %<>% factor(., levels = levels(diff_data$Taxa))
			if(is.null(group_order)){
				if((!is.null(self$group_order)) & (length(unique(abund_data$Group)) == length(self$group_order))){
					abund_data$Group %<>% factor(., levels = rev(self$group_order))
				}else{
					abund_data$Group %<>% as.character %>% as.factor
				}
			}else{
				abund_data$Group %<>% factor(., levels = rev(group_order))
			}
			if(length(color_values) < length(levels(abund_data$Group))){
				stop("Please provide color_values parameter with more colors!")
			}else{
				color_values %<>% .[1:length(levels(abund_data$Group))] %>% rev
			}
			if(!coord_flip){
				abund_data$Group %<>% factor(., levels = rev(levels(.)))
				abund_data$Taxa %<>% factor(., levels = rev(levels(.)))
				diff_data$Taxa %<>% factor(., levels = rev(levels(.)))
				color_values %<>% rev
			}
			
			# get labels info
			if(add_sig){
				# assign labels by factor orders
				x_axis_order <- levels(abund_data$Group)
				if(! add_sig_label %in% colnames(diff_data)){
					stop("add_sig_label parameter must be one of colnames of object$res_diff!")
				}
				if(is.factor(diff_data[, add_sig_label])){
					diff_data[, add_sig_label] %<>% as.character
				}else{
					if(is.numeric(diff_data[, add_sig_label])){
						diff_data[, add_sig_label] %<>% round(., 4)
					}
				}
				if(use_se){
					y_start_use <- max((abund_data$Mean + abund_data$SE)) * y_start
				}else{
					y_start_use <- max((abund_data$Mean + abund_data$SD)) * y_start
				}
				all_taxa <- levels(abund_data$Taxa)

				if(length(levels(abund_data$Group)) > 2 & method %in% c("lefse", "rf", "KW")){
					add_letter_text <- diff_data[match(all_taxa, diff_data$Taxa), add_sig_label]
					textdf <- data.frame(
						x = all_taxa, 
						y = y_start_use, 
						add = add_letter_text, 
						stringsAsFactors = FALSE
						)
				}else{
					if(method != "anova"){
						if(any(grepl("\\s-\\s", x_axis_order))){
							stop("The group names have ' - ' characters, which can impede the group recognition and mapping in the plot! Please rename groups and rerun!")
						}
						annotations <- c()
						x_min <- c()
						x_max <- c()
						y_position <- c()

						start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
						increase_bar_mid <- barwidth/length(x_axis_order)

						for(j in all_taxa){
							select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Taxa == j, ]
							for(i in seq_len(nrow(select_use_diff_data))){
								# first determine the bar range
								mid_num <- match(j, all_taxa) - 1
								annotations %<>% c(., select_use_diff_data[i, add_sig_label])
								x_min %<>% c(., mid_num + 
									(start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\1", select_use_diff_data[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
								x_max %<>% c(., mid_num + 
									(start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\2", select_use_diff_data[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
								y_position %<>% c(., y_start_use * (1 + i * y_increase))
							}
						}
					}else{
						x_mid <- c()
						annotations <- c()
						y_position <- c()

						start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
						increase_bar_mid <- barwidth/length(x_axis_order)

						for(j in all_taxa){
							select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Taxa == j, ]
							for(i in seq_len(nrow(select_use_diff_data))){
								# first determine the bar range
								mid_num <- match(j, all_taxa) - 1
								annotations %<>% c(., select_use_diff_data[i, add_sig_label])
								x_mid %<>% c(., mid_num + (start_bar_mid + (match(select_use_diff_data[i, "Group"], x_axis_order) - 1) * increase_bar_mid))
								abund_data_select <- abund_data[abund_data$Group == select_use_diff_data[i, "Group"] & abund_data$Taxa == j, ]
								if(use_se){
									y_position %<>% c(., y_start * (abund_data_select$Mean + abund_data_select$SE) + y_increase * max(abund_data$Mean))						
								}else{
									y_position %<>% c(., y_start * (abund_data_select$Mean + abund_data_select$SD) + y_increase * max(abund_data$Mean))
								}
							}
						}
						textdf <- data.frame(
							x = x_mid, 
							y = y_position, 
							add = annotations, 
							stringsAsFactors = FALSE
							)
					}
				}
			}
			
			p <- ggplot(abund_data, aes(x = Taxa, y = Mean, color = Group, fill = Group)) +
				theme_bw() +
				geom_bar(stat="identity", position = position_dodge(), width = barwidth)
			if(use_se == T){
				p <- p + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.45, position=position_dodge(barwidth), color = "black")
			}else{
				p <- p + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.45, position=position_dodge(barwidth), color = "black")
			}
			if(add_sig){
				if(length(levels(abund_data$Group)) > 2 & method %in% c("lefse", "rf", "KW")){
					p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
				}else{
					if(method != "anova"){
						p <- p + ggsignif::geom_signif(
							annotations = annotations,
							y_position = y_position, 
							xmin = x_min, 
							xmax = x_max,
							color = add_sig_label_color,
							tip_length = add_sig_tip_length,
							...
							)
					}else{
						p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
					}
				}
			}
			p <- p +
				scale_color_manual(values = color_values) +
				scale_fill_manual(values = color_values) +
				ylab("Relative abundance") +
				theme(legend.position = "right") +
				theme(panel.border = element_blank(), panel.background=element_rect(fill="white")) +
				theme(axis.title = element_text(size = 17))
				
			if(coord_flip){
				p <- p + coord_flip() + guides(fill = guide_legend(reverse = TRUE, ncol = 1), color = "none") +
					theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) +
					theme(axis.title.y = element_blank(), axis.text.y = element_text(size = text_y_size, color = "black"))
			}else{
				p <- p + guides(fill = guide_legend(reverse = FALSE, ncol=1), color = "none") +
					theme(axis.text.x = element_text(angle = 45, colour = "black", vjust = 1, hjust = 1, size = text_y_size)) +
					theme(axis.title.x = element_blank(), axis.text.x = element_text(size = text_y_size, color = "black")) +
					theme(plot.margin = unit(c(.1, .1, .1, 1), "cm"))
			}
			p
		},
		#' @description
		#' Bar plot for indicator index, such as LDA score and P value.
		#'
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors palette for different groups.
		#' @param use_number default 1:10; numeric vector; the taxa numbers used in the plot, i.e. 1:n.
		#' @param threshold default NULL; threshold value of indicators for selecting taxa, such as 3 for LDA score of LEfSe.
		#' @param select_group default NULL; this is used to select the paired group when multiple comparisions are generated;
		#'   The input select_group must be one of object$res_diff$Comparison.
		#' @param simplify_names default TRUE; whether use the simplified taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param group_order default NULL; a vector to order the legend and colors in plot; 
		#' 	  If NULL, the function can first check whether the group column of sample_table is factor. If yes, use the levels in it.
		#' 	  If provided, this parameter can overwrite the levels in the group of sample_table.
		#' @param axis_text_y default 12; the size for the y axis text.
		#' @param plot_vertical default TRUE; whether use vertical bar plot or horizontal.
		#' @param ... parameters pass to \code{\link{geom_bar}}
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_diff_bar(use_number = 1:20)
		#' }
		plot_diff_bar = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			use_number = 1:10,
			threshold = NULL,
			select_group = NULL,
			simplify_names = TRUE,
			keep_prefix = TRUE,
			group_order = NULL,
			axis_text_y = 12,
			plot_vertical = TRUE,
			...
			){
			use_data <- self$res_diff
			method <- self$method
			if(method == "lefse"){
				colnames(use_data)[colnames(use_data) == "LDA"] <- "Value"
				ylab_title <- "LDA score"
			}else{
				if(method == "rf"){
					colnames(use_data)[colnames(use_data) == "MeanDecreaseGini"] <- "Value"
					ylab_title <- "MeanDecreaseGini"
				}else{
					if(method == "metastat"){
						use_data %<>% .[.$qvalue < 0.05, ]
						use_data$Value <- 1 - use_data$qvalue
						ylab_title <- "1 - qvalue"
					}else{
						if(method != "anova"){
							use_data %<>% .[.$P.adj < 0.05, ]
							use_data$Value <- 1 - use_data$P.adj
							ylab_title <- "1 - P.adjust"
						}else{
							stop("This function can not be used to ", method," currently!")
						}
					}
				}
			}
			if(!method %in% c("lefse", "rf")){
				if(length(unique(use_data$Comparison)) > 1){
					# make sure the Group not replicated for multiple comparisions
					if(is.null(select_group)){
						message('Multiple comparisions found. But select_group parameter not provided. Select the first group pair to show ...')
						select_group <- unique(use_data$Comparison)[1]
					}else{
						if(length(select_group) > 1){
							stop("The select_group parameter should only have one element! Please check the input!")
						}
						if(! select_group %in% use_data$Comparison){
							stop("The select_group parameter must be one of elements of object$res_diff$Comparison!")
						}
					}
					use_data %<>% .[.$Comparison %in% select_group, ]
				}
			}
			if(nrow(use_data) == 0){
				stop("No significant taxa can be used to show!")
			}
			if(simplify_names == T){
				use_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				use_data$Taxa %<>% gsub(".__", "", .)
			}
			# make sure no duplication
			use_data %<>% .[!duplicated(.$Taxa), ]
			if(nrow(use_data) == 0){
				stop("No available data can be used to show!")
			}
			if(is.null(threshold)){
				sel_num <- use_number
			}else{
				sel_num <- sum(use_data$Value > threshold)
				if(sel_num == 0){
					stop("Too large threshold provided, no data selected!")
				}
				sel_num <- 1:sel_num
			}
			if(length(sel_num) > nrow(use_data)){
				sel_num <- 1:nrow(use_data)
			}
			use_data %<>% .[sel_num, ]
			if(is.null(group_order)){
				if((!is.null(self$group_order)) & (length(unique(use_data$Group)) == length(self$group_order))){
					use_data$Group %<>% factor(., levels = self$group_order)
				}else{
					use_data$Group %<>% as.character %>% as.factor
				}
			}else{
				use_data$Group %<>% factor(., levels = group_order)
			}
			if(length(color_values) < length(levels(use_data$Group))){
				stop("Please provide color_values parameter with more colors!")
			}
			
			# rearrange orders
			if(length(levels(use_data$Group)) == 2){
				use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$Group), function(x){
					if(x == levels(use_data$Group)[1]){
						use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$Value, decreasing = TRUE), "Taxa"]
					}else{
						use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$Value, decreasing = FALSE), "Taxa"]
					}
				})))))
				use_data[use_data$Group == levels(use_data$Group)[2], "Value"] %<>% {. * -1}
			}else{
				use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$Group), function(x){
					use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$Value, decreasing = TRUE), "Taxa"]
				})))))
			}
			self$plot_diff_bar_taxa <- levels(use_data$Taxa) %>% rev
			
			p <- ggplot(use_data, aes(x = Taxa, y = Value, color = Group, fill = Group, group = Group)) +
				geom_bar(stat="identity", position = position_dodge(), ...) +
				theme_bw() +
				scale_color_manual(values = color_values) +
				scale_fill_manual(values = color_values) +
				ylab(ylab_title) +
				xlab("") +
				theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = axis_text_y, color = "black")) +
				theme(axis.text.x = element_text(size = 10)) +
				theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()) +
				theme(panel.border = element_blank()) +
				theme(axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square"))
				
			if(plot_vertical == T){
				p <- p + coord_flip()
			}
			p
		},
		#' @description
		#' Plot the cladogram using taxa with significant difference.
		#'
		#' @param color default RColorBrewer::brewer.pal(8, "Dark2"); color palette used in the plot.
		#' @param group_order default NULL; a vector to order the legend in plot; 
		#' 	  If NULL, the function can first check whether the group column of sample_table is factor. If yes, use the levels in it.
		#' 	  If provided, this parameter can overwrite the levels in the group of sample_table. 
		#' 	  If the number of provided group_order is less than the number of groups in res_diff$Group, the function will select the groups of group_order automatically.
		#' @param use_taxa_num default 200; integer; The taxa number used in the background tree plot; select the taxa according to the mean abundance .
		#' @param filter_taxa default NULL; The mean relative abundance used to filter the taxa with low abundance.
		#' @param use_feature_num default NULL; integer; The feature number used in the plot; 
		#'	  select the features according to the LDA score (method = "lefse") or MeanDecreaseGini (method = "rf") from high to low.
		#' @param clade_label_level default 4; the taxonomic level for marking the label with letters, root is the largest.
		#' @param select_show_labels default NULL; character vector; The features to show in the plot with full label names, not the letters.
		#' @param only_select_show default FALSE; whether only use the the select features in the parameter select_show_labels.
		#' @param sep default "|"; the seperate character in the taxonomic information.
		#' @param branch_size default 0.2; numberic, size of branch.
		#' @param alpha default 0.2; shading of the color.
		#' @param clade_label_size default 2; basic size for the clade label; please also see clade_label_size_add and clade_label_size_log
		#' @param clade_label_size_add default 5; added basic size for the clade label; see the formula in clade_label_size_log parameter.
		#' @param clade_label_size_log default exp(1); the base of log function for added size of the clade label; the size formula: 
		#'   clade_label_size + log(clade_label_level + clade_label_size_add, base = clade_label_size_log); 
		#'   so use clade_label_size_log, clade_label_size_add and clade_label_size
		#'   can totally control the label size for different taxonomic levels.
		#' @param node_size_scale default 1; scale for the node size.
		#' @param node_size_offset default 1; offset for the node size.
		#' @param annotation_shape default 22; shape used in the annotation legend.
		#' @param annotation_shape_size default 5; size used in the annotation legend.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, select_show_labels = NULL)
		#' }
		plot_diff_cladogram = function(
			color = RColorBrewer::brewer.pal(8, "Dark2"),
			group_order = NULL,
			use_taxa_num = 200,
			filter_taxa = NULL,
			use_feature_num = NULL,
			clade_label_level = 4,
			select_show_labels = NULL,
			only_select_show = FALSE,
			sep = "|",
			branch_size = 0.2,
			alpha = 0.2,
			clade_label_size = 2,
			clade_label_size_add = 5,
			clade_label_size_log = exp(1),
			node_size_scale = 1,
			node_size_offset = 1,
			annotation_shape = 22,
			annotation_shape_size = 5
			){
			# developed based on microbiomeMarker 
			abund_table <- self$abund_table
			marker_table <- self$res_diff %>% dropallfactors
			method <- self$method
			
			if(! method %in% c("lefse", "rf")){
				stop("This function currently can only be used for method = 'lefse' or 'rf' !")
			}
			if(self$taxa_level != "all"){
				stop("This function is only useful when taxa_level = 'all' !")
			}
			if(!is.null(use_feature_num)){
				marker_table %<>% .[1:use_feature_num, ]
			}
			if(only_select_show == T){
				marker_table %<>% .[.$Taxa %in% select_show_labels, ]
			}
			# color legend order settings
			if(is.null(group_order)){
				if(! is.null(self$group_order)){
					color_groups <- self$group_order
				}else{
					color_groups <- marker_table$Group %>% as.character %>% as.factor %>% levels
				}
			}else{
				color_groups <- group_order
			}
			
			# filter redundant groups
			color_groups %<>% .[. %in% unique(marker_table$Group)]
			marker_table %<>% .[.$Group %in% color_groups, ]
			# get the color palette
			if(length(color) < length(unique(marker_table$Group))){
				stop("Please provide enough color palette! There are ", length(unique(marker_table$Group)), 
					" groups, but only ", length(color), " colors provideed in color parameter!")
			}else{
				color <- color[1:length(unique(marker_table$Group))]
			}

			# filter the taxa with unidentified classification or with space, in case of the unexpected error in the following operations
			abund_table %<>% {.[!grepl("\\|.__\\|", rownames(.)), ]} %>%
				{.[!grepl("\\s", rownames(.)), ]} %>%
				# also filter uncleared classification to make it in line with the lefse above
				{.[!grepl("Incertae_sedis|unculture", rownames(.), ignore.case = TRUE), ]}

			if(!is.null(use_taxa_num)){
				abund_table %<>% .[names(sort(apply(., 1, mean), decreasing = TRUE)[1:use_taxa_num]), ]
			}
			if(!is.null(filter_taxa)){
				abund_table %<>% .[apply(., 1, mean) > (self$lefse_norm * filter_taxa), ]
			}
			abund_table %<>% .[sort(rownames(.)), ]

			tree_table <- data.frame(taxa = row.names(abund_table), abd = rowMeans(abund_table), stringsAsFactors = FALSE) %>%
				dplyr::mutate(taxa =  paste("r__Root", .data$taxa, sep = sep), abd = .data$abd/max(.data$abd)*100)
			taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
			nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)

			# check whether some nodes duplicated from bad classification
			if(any(duplicated(nodes))){
				del <- nodes %>% .[duplicated(.)] %>% unique
				for(i in del){
					tree_table %<>% .[!grepl(paste0("\\|", i, "($|\\|)"), .$taxa), ]
				}
				taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
				nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
			}

			# add root node
			nodes %<>% c("r__Root", .)
			# levels used for extend of clade label
			label_levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
				factor(levels = rev(unlist(lapply(taxa_split, function(x) gsub("(.)__.*", "\\1", x))) %>% .[!duplicated(.)]))

			# root must be a parent node
			nodes_parent <- purrr::map_chr(taxa_split, ~ .x[length(.x) - 1]) %>% c("root", .)

			## add index for nodes
			is_tip <- !nodes %in% nodes_parent
			index <- vector("integer", length(is_tip))
			index[is_tip] <- 1:sum(is_tip)
			index[!is_tip] <- (sum(is_tip)+1):length(is_tip)

			edges <- cbind(parent = index[match(nodes_parent, nodes)], child = index)
			edges <- edges[!is.na(edges[, 1]), ]
			# not label the tips
			node_label <- nodes[!is_tip]
			phylo <- structure(list(
				edge = edges, 
				node.label = node_label, 
				tip.label = nodes[is_tip], 
				edge.length = rep(1, nrow(edges)), 
				Nnode = length(node_label)
				), class = "phylo")
			mapping <- data.frame(
				node = index, 
				abd = c(100, tree_table$abd),
				node_label = nodes, 
				stringsAsFactors = FALSE)
			mapping$node_class <- label_levels
			tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
			tree <- ggtree::ggtree(tree, size = 0.2, layout = 'circular')
			
			annotation <- private$generate_cladogram_annotation(marker_table, tree = tree, color = color, color_groups = color_groups)
			# backgroup hilight
			annotation_info <- dplyr::left_join(annotation, tree$data, by = c("node" = "label")) %>%
				dplyr::mutate(label = .data$node, id = .data$node.y, level = as.numeric(.data$node_class))
			hilight_para <- dplyr::transmute(
				annotation_info,
				node = .data$id,
				fill = .data$color,
				alpha = alpha,
				extend = private$get_offset(.data$level)
			)
			hilights_g <- purrr::pmap(hilight_para, ggtree::geom_hilight)
			tree <- purrr::reduce(hilights_g, `+`, .init = tree)

			# hilight legend
			hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group, .data$color)
			hilights_df$x <- 0
			hilights_df$y <- 1
			# resort the table used for the legend color and text
			hilights_df %<>% `row.names<-`(.$enrich_group) %>% .[color_groups, ]
			# make sure the right order in legend
			hilights_df$enrich_group %<>% factor(., levels = color_groups)

			# add legend
			tree <- tree + 
				geom_rect(aes_(xmin = ~x, xmax = ~x, ymax = ~y, ymin = ~y, fill = ~enrich_group), data = hilights_df, inherit.aes = FALSE) +
				guides(fill = guide_legend(title = NULL, order = 1, override.aes = list(fill = hilights_df$color)))

			# set nodes color and size
			nodes_colors <- rep("white", nrow(tree$data))
			nodes_colors[annotation_info$id] <- annotation_info$color
			node_size <- node_size_scale*log(tree$data$abd) + node_size_offset
			tree$data$node_size <- node_size
			tree <- tree + ggtree::geom_point2(aes(size = I(node_size)), fill = nodes_colors, shape = 21)

			## add clade labels
			clade_label <- dplyr::transmute(
				annotation_info,
				node = .data$id,
				offset = private$get_offset(.data$level)-0.4,
				offset.text = 0,
				angle = purrr::map_dbl(.data$id, private$get_angle, tree = tree),
				label = .data$label,
				fontsize = clade_label_size + log(.data$level + clade_label_size_add, base = clade_label_size_log),
				barsize = 0,
				extend = 0.2,
				hjust = 0.5,
				level = .data$level
			) %>% dplyr::arrange(desc(.data$level))

			clade_label$offset.text <- unlist(lapply(seq_len(nrow(clade_label)), function(x){
				if(clade_label$angle[x] < 180){ 0.2 }else{ 0 }}))
			clade_label$angle <- unlist(lapply(clade_label$angle, function(x){
				if(x < 180){ x - 90 }else{ x + 90 }}))
			clade_label_new <- clade_label

			# add letters label to replace long taxonomic label
			if(is.null(select_show_labels)){
				# outer circle --> larger level; label smaller levels, i.e. finer taxonomy
				ind <- clade_label$level < clade_label_level	
			}else{
				ind <- ! clade_label$label %in% select_show_labels
			}
			ind_num <- sum(ind)

			if(ind_num > 0){
				if(ind_num < 27){
					use_letters <- letters
				}else{
					if(ind_num < 326){
						use_letters <- apply(combn(letters, 2), 2, function(x){paste0(x, collapse = "")})
					}else{
						stop("Too much features to be labelled with letters, consider to use use_feature_num parameter to reduce the number!")
					}
				}
				clade_label_new$label_legend <- clade_label_new$label_show <- clade_label_new$label_raw <- clade_label_new$label
				clade_label_new$label_show[ind] <- use_letters[1:ind_num]
				clade_label_new$label_legend[ind] <- paste0(clade_label_new$label_show[ind], ": ", clade_label_new$label[ind])
				clade_label_new$label <- clade_label_new$label_show
				# delete redundant columns to avoid warnings
				clade_label <- clade_label_new %>% .[, which(! colnames(.) %in% c("label_raw", "label_show", "label_legend", "level"))]
			}

			clade_label_g <- purrr::pmap(clade_label, ggtree::geom_cladelabel)
			tree <- purrr::reduce(clade_label_g, `+`, .init = tree)

			# if letters are used, add guide labels
			if(ind_num > 0){
				guide_label <- clade_label_new[ind, ] %>%
					dplyr::mutate(color = annotation_info$color[match(.data$label_raw, annotation_info$label)])
				tree <- tree + 
					geom_point(data = guide_label, inherit.aes = FALSE, aes_(x = 0, y = 0, shape = ~label_legend), size = 0, stroke = 0) +
						scale_shape_manual(values = rep(annotation_shape, nrow(guide_label))) +
						guides(shape = guide_legend(override.aes = list(
							size = annotation_shape_size, shape = annotation_shape, fill = guide_label$color)))
			}
			tree <- tree + theme(legend.position = "right", legend.title = element_blank())
			tree

		},
		#' @description
		#' Print the trans_alpha object.
		print = function() {
			cat("trans_diff object:\n")
			cat(paste("res_diff have", ncol(self$res_diff), "columns: ", paste0(colnames(self$res_diff), collapse = ", "), "\n"))
			if(!is.null(self$res_abund)) cat(paste("res_abund have", ncol(self$res_abund), "columns: ", paste0(colnames(self$res_abund), collapse = ", "), "\n"))
			invisible(self)
		}
		),
	private = list(
		# group test for lefse or rf
		test_mark = function(dataframe, group, min_num_nonpara = 1, method = NULL){
			d1 <- as.data.frame(t(dataframe))
			taxaname <- colnames(d1)[1]
			d1$Group <- group
			colnames(d1)[1] <- "Value"
			formu <- reformulate("Group", "Value")
			if(any(table(as.character(group))) < min_num_nonpara){
				list(p_value = NA, med = NA)
			}else{
				if(! is.null(method)){
					method <- match.arg(method, c("wilcox.test", "kruskal.test"))
					if(method == "wilcox.test"){
						res1 <- wilcox.test(formula = formu, data = d1)
					}else{
						res1 <- kruskal.test(formula = formu, data = d1)
					}
				}else{
					if(length(unique(as.character(group))) == 2){
						res1 <- wilcox.test(formula = formu, data = d1)
					}else{
						res1 <- kruskal.test(formula = formu, data = d1)
					}
				}
				if(is.nan(res1$p.value)){
					res1$p.value <- 1
				}
				med <- tapply(d1[,1], group, median) %>% as.data.frame
				colnames(med) <- taxaname
				list(p_value = res1$p.value, med = med)
			}
		},
		# generate the cladogram annotation table
		generate_cladogram_annotation = function(marker_table, tree, color, color_groups, sep = "|") {
			use_marker_table <- marker_table
			feature <- use_marker_table$Taxa
			label <- strsplit(feature, split = sep, fixed = TRUE) %>% 
				purrr::map_chr(utils::tail, n =1)
			plot_color <- use_marker_table$Group %>% 
				as.character
			for(i in seq_along(color_groups)){
				plot_color[plot_color == color_groups[i]] <- color[i]
			}
			annotation <- data.frame(
				node = label,
				color = plot_color,
				enrich_group = use_marker_table$Group,
				stringsAsFactors = FALSE
			)
			# filter the feature with bad classification
			annotation %<>% .[label %in% tree$data$label, ]
			annotation
		},
		get_angle = function(tree, node){
			if (length(node) != 1) {
				stop("The length of `node` must be 1")
			}
			tree_data <- tree$data
			sp <- tidytree::offspring(tree_data, node)$node
			sp2 <- c(sp, node)
			sp.df <- tree_data[match(sp2, tree_data$node),]
			mean(range(sp.df$angle))
		},
		get_offset = function(x) {(x*0.2+0.2)^2},
		# metastat input
		calculate_metastat = function(inputdata, g, pflag = FALSE, threshold = NULL, B = NULL){
			trans_data <- private$load_frequency_matrix(input = inputdata)
			res <- private$detect_differentially_abundant_features(jobj = trans_data, g = g, pflag = pflag, threshold = threshold, B = B)
			res
		},
		#*****************************************************************************************************
		#Modified from raw metastat code
		#load up the frequency matrix from a file
		#*****************************************************************************************************
		# Note sep
		load_frequency_matrix = function(input){
			# load names 
			subjects <- array(0,dim=c(ncol(input)));
			for(i in 1:length(subjects)) {
				subjects[i] <- as.character(colnames(input)[i]);
			}
			# load taxa
			taxa <- array(0,dim=c(nrow(input)));
			for(i in 1:length(taxa)) {
				taxa[i] <- as.character(rownames(input)[i]);
			}
			# load remaining counts
			matrix <- array(0, dim=c(length(taxa),length(subjects)));
			for(i in 1:length(taxa)){
				for(j in 1:length(subjects)){ 
					matrix[i,j] <- as.numeric(input[i,j]);
				}
			}
			jobj <- list(matrix=matrix, taxa=taxa)
			return(jobj)
		},
		#  Modified from metastat raw codes
		# http://metastats.cbcb.umd.edu/detect_DA_features.r
		#*****************************************************************************************************
		#  Author: james robert white, whitej@umd.edu, Center for Bioinformatics and Computational Biology.
		#  University of Maryland - College Park, MD 20740
		#
		#  This software is designed to identify differentially abundant features between two groups
		#  Input is a matrix of frequency data. Several thresholding options are available.
		#  See documentation for details.
		#*****************************************************************************************************
		#*****************************************************************************************************
		#  detect_differentially_abundant_features:
		#  the major function - inputs an R object "jobj" containing a list of feature names and the 
		#  corresponding frequency matrix, the argument g is the first column of the second group. 
		#  
		#  -> set the pflag to be TRUE or FALSE to threshold by p or q values, respectively
		#  -> threshold is the significance level to reject hypotheses by.
		#  -> B is the number of bootstrapping permutations to use in estimating the null t-stat distribution.
		#*****************************************************************************************************
		#*****************************************************************************************************
		detect_differentially_abundant_features = function(jobj, g, pflag = NULL, threshold = NULL, B = NULL){
			#**********************************************************************************
			# ************************ INITIALIZE COMMAND-LINE ********************************
			# ************************        PARAMETERS       ********************************
			#**********************************************************************************
			qflag = FALSE;
			if (is.null(B)){
				B = 1000;
			}
			if (is.null(threshold)){
				threshold = 0.05;
			}
			if (is.null(pflag)){
				pflag = TRUE;
				qflag = FALSE;
			}
			if (pflag == TRUE){
				qflag = FALSE;
			}
			if (pflag == FALSE){
				qflag = TRUE;
			}
			#********************************************************************************
			# ************************ INITIALIZE PARAMETERS ********************************
			#********************************************************************************
			Fmatrix <- jobj$matrix;                   # the feature abundance matrix
			taxa <- jobj$taxa;                        # the taxa/(feature) labels of the TAM
			nrows = nrow(Fmatrix);                   
			ncols = ncol(Fmatrix);
			Pmatrix <- array(0, dim=c(nrows,ncols));  # the relative proportion matrix
			C1 <- array(0, dim=c(nrows,3));           # statistic profiles for class1 and class 2
			C2 <- array(0, dim=c(nrows,3));           # mean[1], variance[2], standard error[3]   
			T_statistics <- array(0, dim=c(nrows,1)); # a place to store the true t-statistics 
			pvalues <- array(0, dim=c(nrows,1));      # place to store pvalues
			qvalues <- array(0, dim=c(nrows,1));      # stores qvalues
			#*************************************
			#  convert to proportions
			#  generate Pmatrix
			#*************************************
			totals <- array(0, dim=c(ncol(Fmatrix)));
			for (i in 1:ncol(Fmatrix)) {
				# sum the ith column 
				totals[i] = sum(Fmatrix[,i]);
			}
			for (i in 1:ncols) {   # for each subject
				for (j in 1:nrows) { # for each row
					Pmatrix[j,i] = Fmatrix[j,i]/totals[i];
				}
			}
			#********************************************************************************
			# ************************** STATISTICAL TESTING ********************************
			#********************************************************************************
			if (ncols == 2){  # then we have a two sample comparison
				#************************************************************
				#  generate p values using chisquared or fisher's exact test
				#************************************************************
				for (i in 1:nrows){           # for each feature
					f11 = sum(Fmatrix[i,1]);
					f12 = sum(Fmatrix[i,2]);
					f21 = totals[1] - f11;
					f22 = totals[2] - f12;
					C1[i,1] = f11/totals[1];                       # proportion estimate
					C1[i,2] = (C1[i,1]*(1-C1[i,1]))/(totals[1]-1); # sample variance
					C1[i,3] = sqrt(C1[i,2]);                       # sample standard error
					C2[i,1] = f12/totals[2];
					C2[i,2] = (C2[i,1]*(1-C2[i,1]))/(totals[2]-1);
					C2[i,3] = sqrt(C2[i,2]);
					
					#  f11  f12
					#  f21  f22  <- contigency table format
					contingencytable <- array(0, dim=c(2,2));
					contingencytable[1,1] = f11;
					contingencytable[1,2] = f12;
					contingencytable[2,1] = f21;
					contingencytable[2,2] = f22;
					if (f11 > 20 && f22 > 20){
						csqt <- stats::chisq.test(contingencytable);
						pvalues[i] = csqt$p.value;
					}else{
						ft <- stats::fisher.test(contingencytable, workspace = 8e6, alternative = "two.sided", conf.int = FALSE);
						pvalues[i] = ft$p.value;
					}
				}
				#*************************************
				#  calculate q values from p values
				#*************************************
				qvalues <- private$calc_qvalues(pvalues);
			}else{ # we have multiple subjects per population
				#*************************************
				#  generate statistics mean, var, stderr    
				#*************************************
				for (i in 1:nrows){ # for each taxa
					# find the mean of each group
					C1[i,1] = mean(Pmatrix[i, 1:g-1]);  
					C1[i,2] = stats::var(Pmatrix[i, 1:g-1]); # variance
					C1[i,3] = C1[i,2]/(g-1);    # std err^2 (will change to std err at end)

					C2[i,1] = mean(Pmatrix[i, g:ncols]);  
					C2[i,2] = stats::var(Pmatrix[i, g:ncols]);  # variance
					C2[i,3] = C2[i,2]/(ncols-g+1); # std err^2 (will change to std err at end)
				}
				#*************************************
				#  two sample t-statistics
				#*************************************
				for (i in 1:nrows){                   # for each taxa
					xbar_diff = C1[i,1] - C2[i,1]; 
					denom = sqrt(C1[i,3] + C2[i,3]);
					T_statistics[i] = xbar_diff/denom;  # calculate two sample t-statistic
				}
				#*************************************
				# generate initial permuted p-values
				#*************************************
				pvalues <- private$permuted_pvalues(Pmatrix, T_statistics, B, g, Fmatrix);
				#*************************************
				#  generate p values for sparse data 
				#  using fisher's exact test
				#*************************************
				for (i in 1:nrows){                   # for each taxa
					if (sum(Fmatrix[i,1:(g-1)]) < (g-1) && sum(Fmatrix[i,g:ncols]) < (ncols-g+1)){
						# then this is a candidate for fisher's exact test
						f11 = sum(Fmatrix[i,1:(g-1)]);
						f12 = sum(Fmatrix[i,g:ncols]);
						f21 = sum(totals[1:(g-1)]) - f11;
						f22 = sum(totals[g:ncols]) - f12;
						#  f11  f12
						#  f21  f22  <- contigency table format
						contingencytable <- array(0, dim=c(2,2));
						contingencytable[1,1] = f11;
						contingencytable[1,2] = f12;
						contingencytable[2,1] = f21;
						contingencytable[2,2] = f22;
						ft <- fisher.test(contingencytable, workspace = 8e6, alternative = "two.sided", conf.int = FALSE);
						pvalues[i] = ft$p.value; 
					}  
				}

				#*************************************
				#  calculate q values from p values
				#*************************************
				qvalues <- private$calc_qvalues(pvalues);
				#*************************************
				#  convert stderr^2 to std error
				#*************************************
				for (i in 1:nrows){
					C1[i,3] = sqrt(C1[i,3]);
					C2[i,3] = sqrt(C2[i,3]);
				}
			}
			#*************************************
			#  threshold sigvalues and print
			#*************************************
			sigvalues <- array(0, dim=c(nrows,1));
			if (pflag == TRUE){  # which are you thresholding by?
				sigvalues <- pvalues;
			}else{
				sigvalues <- qvalues;
			}
			s = sum(sigvalues <= threshold);
			Differential_matrix <- array(0, dim=c(s,9));
			dex = 1;
			for (i in 1:nrows){
				if (sigvalues[i] <= threshold){
					Differential_matrix[dex,1]   = jobj$taxa[i];
					Differential_matrix[dex,2:4] = C1[i,];
					Differential_matrix[dex,5:7] = C2[i,];
					Differential_matrix[dex,8]   = pvalues[i];  
					Differential_matrix[dex,9]   = qvalues[i];
					dex = dex+1;  
				}
			}
			# show(Differential_matrix);
			Total_matrix <- array(0, dim=c(nrows,9));
			for (i in 1:nrows){
				Total_matrix[i,1]   = jobj$taxa[i];
				Total_matrix[i,2:4] = C1[i,];
				Total_matrix[i,5:7] = C2[i,];
				Total_matrix[i,8]   = pvalues[i];
				Total_matrix[i,9]   = qvalues[i];
			}
			colnames(Total_matrix) <- c("taxa", "mean(group1)", "variance(group1)", "std.err(group1)", "mean(group2)", "variance(group2)", "std.err(group2)", "pvalue", "qvalue")
			Total_matrix <- Total_matrix[order(Total_matrix[, "pvalue"], decreasing = FALSE), ]
			Total_matrix
		},
		#*****************************************************************************************************
		# takes a matrix, a permutation vector, and a group division g.
		# returns a set of ts based on the permutation.
		#*****************************************************************************************************
		permute_and_calc_ts = function(Imatrix, y, g){
			nr = nrow(Imatrix);
			nc = ncol(Imatrix);
			# first permute the rows in the matrix
			Pmatrix <- Imatrix[,y[1:length(y)]];
			Ts <- private$calc_twosample_ts(Pmatrix, g, nr, nc);
			return (Ts)
		},
		#*****************************************************************************************************
		#  function to calculate qvalues.
		#  takes an unordered set of pvalues corresponding the rows of the matrix
		#*****************************************************************************************************
		calc_qvalues = function(pvalues){
			nrows = length(pvalues);
			# create lambda vector
			lambdas <- seq(0,0.95,0.01);
			pi0_hat <- array(0, dim=c(length(lambdas)));
			# calculate pi0_hat
			for (l in 1:length(lambdas)){ # for each lambda value
				count = 0;
				for (i in 1:nrows){ # for each p-value in order
					if (pvalues[i] > lambdas[l]){
						count = count + 1; 	
					}
					pi0_hat[l] = count/(nrows*(1-lambdas[l]));
				}
			}
			f <- unclass(stats::smooth.spline(lambdas,pi0_hat,df=3));
			f_spline <- f$y;
			pi0 = f_spline[length(lambdas)];   # this is the essential pi0_hat value
			# order p-values
			ordered_ps <- order(pvalues);
			pvalues <- pvalues;
			qvalues <- array(0, dim=c(nrows));
			ordered_qs <- array(0, dim=c(nrows));
			ordered_qs[nrows] <- min(pvalues[ordered_ps[nrows]]*pi0, 1);
			for(i in (nrows-1):1) {
				p = pvalues[ordered_ps[i]];
				new = p*nrows*pi0/i;

				ordered_qs[i] <- min(new,ordered_qs[i+1],1);
			}
			# re-distribute calculated qvalues to appropriate rows
			for (i in 1:nrows){
				qvalues[ordered_ps[i]] = ordered_qs[i];
			}
			################################
			# plotting pi_hat vs. lambda
			################################
			# plot(lambdas,pi0_hat,xlab=expression(lambda),ylab=expression(hat(pi)[0](lambda)),type="p");
			# lines(f);
			return (qvalues);
		},
		##################################################################################
		# metastat code from White et al. (2009) <doi:10.1371/journal.pcbi.1000352>.
		#************************************************************************
		# ************************** SUBROUTINES ********************************
		#*****************************************************************************************************
		#  calc two sample two statistics
		#  g is the first column in the matrix representing the second condition
		#*****************************************************************************************************
		calc_twosample_ts = function(Pmatrix, g, nrows, ncols)
		{
			C1 <- array(0, dim=c(nrows,3));  # statistic profiles
			C2 <- array(0, dim=c(nrows,3)); 
			Ts <- array(0, dim=c(nrows,1));
			if (nrows == 1){
				C1[1,1] = mean(Pmatrix[1:g-1]);
				C1[1,2] = stats::var(Pmatrix[1:g-1]); # variance
				C1[1,3] = C1[1,2]/(g-1);    # std err^2

				C2[1,1] = mean(Pmatrix[g:ncols]);
				C2[1,2] = stats::var(Pmatrix[g:ncols]);  # variance
				C2[1,3] = C2[1,2]/(ncols-g+1); # std err^2
			}else{
				# generate statistic profiles for both groups
				# mean, var, stderr
				for (i in 1:nrows){ # for each taxa
					# find the mean of each group
					C1[i,1] = mean(Pmatrix[i, 1:g-1]);  
					C1[i,2] = stats::var(Pmatrix[i, 1:g-1]); # variance
					C1[i,3] = C1[i,2]/(g-1);    # std err^2

					C2[i,1] = mean(Pmatrix[i, g:ncols]);  
					C2[i,2] = stats::var(Pmatrix[i, g:ncols]);  # variance
					C2[i,3] = C2[i,2]/(ncols-g+1); # std err^2
				}
			}
			# permutation based t-statistics
			for (i in 1:nrows){ # for each taxa
				xbar_diff = C1[i,1] - C2[i,1]; 
				denom = sqrt(C1[i,3] + C2[i,3]);
				Ts[i] = xbar_diff/denom;  # calculate two sample t-statistic 
			}
			return (Ts)
		},
		#*****************************************************************************************************
		#  function to calculate permuted pvalues from Storey and Tibshirani(2003)
		#  B is the number of permutation cycles
		#  g is the first column in the matrix of the second condition 
		#*****************************************************************************************************
		permuted_pvalues = function(Imatrix, tstats, B, g, Fmatrix)
		{
			# B is the number of permutations were going to use!
			# g is the first column of the second sample
			# matrix stores tstats for each taxa(row) for each permuted trial(column)
			M = nrow(Imatrix);
			ps <- array(0, dim=c(M)); # to store the pvalues
			if (is.null(M) || M == 0){
				return (ps);
			}
			permuted_ttests <- array(0, dim=c(M, B));
			ncols = ncol(Fmatrix);
			# calculate null version of tstats using B permutations.
			for (j in 1:B){
				trial_ts <- private$permute_and_calc_ts(Imatrix, sample(1:ncol(Imatrix)), g);
				permuted_ttests[,j] <- abs(trial_ts); 
			}
			# calculate each pvalue using the null ts
			if ((g-1) < 8 || (ncols-g+1) < 8){
				# then pool the t's together!
				# count how many high freq taxa there are
				hfc = 0;
				for (i in 1:M){                   # for each taxa
					if (sum(Fmatrix[i,1:(g-1)]) >= (g-1) || sum(Fmatrix[i,g:ncols]) >= (ncols-g+1)){
						hfc = hfc + 1;
					}
				}
				# the array pooling just the frequently observed ts  
				cleanedpermuted_ttests <- array(0, dim=c(hfc,B));
				hfc = 1;
				for (i in 1:M){
					if (sum(Fmatrix[i,1:(g-1)]) >= (g-1) || sum(Fmatrix[i,g:ncols]) >= (ncols-g+1)){
						cleanedpermuted_ttests[hfc,] = permuted_ttests[i,];
						hfc = hfc + 1;
					}
				}
				#now for each taxa
				for (i in 1:M){  
					ps[i] = (1/(B*hfc))*sum(cleanedpermuted_ttests > abs(tstats[i]));
				}
			}else{
				for (i in 1:M){
					ps[i] = (1/(B+1))*(sum(permuted_ttests[i,] > abs(tstats[i]))+1);
				}
			}
			return(ps)
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
