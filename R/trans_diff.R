#' @title 
#' Create trans_diff object for the differential analysis on the taxonomic abundance.
#'
#' @description
#' This class is a wrapper for a series of differential abundance test and indicator analysis methods, including non-parametric test, 
#' LEfSe based on the Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>, random forest, metastat based on White et al. (2009) <doi:10.1371/journal.pcbi.1000352> and
#' the method in R package metagenomeSeq Paulson et al. (2013) <doi:10.1038/nmeth.2658>.
#'
#' @export
trans_diff <- R6Class(classname = "trans_diff",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param method default "lefse"; "lefse", "rf", "metastat" or "mseq". "lefse": Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>; 
		#' 	  "rf" represents random forest; metastat: White et al. (2009) <doi:10.1371/journal.pcbi.1000352>; "mseq" represents the method in metagenomeSeq package.
		#' @param group default NULL; sample group used for main comparision.
		#' @param lefse_subgroup default NULL; sample sub group used for sub-comparision in lefse; Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>.
		#' @param alpha default .05; significance threshold.
		#' @param lefse_min_subsam default 10; sample numbers required in the subgroup test.
		#' @param lefse_norm default 1000000; scale value in lefse.
		#' @param nresam default .6667; sample number ratio used in each bootstrap or LEfSe or random forest.
		#' @param boots default 30; bootstrap test number for lefse or rf.
		#' @param rf_taxa_level default "all"; use all taxonomic rank data, if want to test a specific rank, provide taxonomic rank name, such as "Genus".
		#' @param rf_ntree default 1000; see ntree in randomForest function of randomForest package.
		#' @param metastat_taxa_level default "Genus"; taxonomic rank level used in metastat test; White et al. (2009) <doi:10.1371/journal.pcbi.1000352>.
		#' @param group_choose_paired default NULL; a vector used for selecting the required groups for paired testing, only used for metastat or mseq.
		#' @param mseq_adjustMethod default "fdr"; Method to adjust p-values by. Default is "fdr". 
		#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
		#' @param mseq_count default 1; Filter features to have at least 'counts' counts.; see the count parameter in MRcoefs function of metagenomeSeq package.
		#' @return res_rf, res_lefse, res_abund, res_metastat, or res_mseq in trans_diff object, depending on the method.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group")
		#' }
		initialize = function(
			dataset = NULL,
			method = c("lefse", "rf", "metastat", "mseq")[1],
			group = NULL,
			lefse_subgroup = NULL,
			alpha = 0.05,
			lefse_min_subsam = 10,
			lefse_norm = 1000000,
			nresam = 0.6667,
			boots = 30,
			rf_taxa_level = "all",
			rf_ntree = 1000,
			metastat_taxa_level = "Genus",
			group_choose_paired = NULL,
			mseq_adjustMethod = "fdr",
			mseq_count = 1
			){
			if(is.null(dataset)){
				stop("No dataset provided!")
			}
			if(is.null(dataset$taxa_abund)){
				stop("Please first calculate taxa_abund! see cal_abund function in microtable class!")
			}
			sampleinfo <- dataset$sample_table
			sampleinfo[, group] %<>% as.character
#			self$method <- method
			if(grepl("lefse|rf", method, ignore.case = TRUE)){
				if(grepl("lefse", method, ignore.case = TRUE)){
					abund_table <- do.call(rbind, unname(lapply(dataset$taxa_abund, function(x) x * lefse_norm)))
					self$lefse_norm <- lefse_norm
				}else{
					if(grepl("all", rf_taxa_level, ignore.case = TRUE)){
						abund_table <- do.call(rbind, unname(dataset$taxa_abund))
					}else{
						abund_table <- dataset$taxa_abund[[rf_taxa_level]]
					}
				}
				abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.)), ]}
				# differential test
				group_vec <- sampleinfo[, group] %>% as.factor
				message("Start differential test for ", group, " ...")
				res_class <- lapply(seq_len(nrow(abund_table)), function(x) private$test_mark(abund_table[x,], group_vec))
				pvalue <- unlist(lapply(res_class, function(x) x$p_value))
				pvalue[is.nan(pvalue)] <- 1
				# select significant taxa
				sel_taxa <- pvalue < alpha
				message("Total ", sum(sel_taxa), " biomarkers found ...")
				if(sum(sel_taxa) == 0){
					stop("No significant biomarkers found! stop running!")
				}
				# save abund_table in self for the cladogram
				self$abund_table <- abund_table
				abund_table_sub <- abund_table[sel_taxa, ]
				pvalue_sub <- pvalue[sel_taxa]
			}
			if(grepl("rf", method, ignore.case = TRUE)){
				names(pvalue_sub) <- rownames(abund_table_sub)
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
				res <- data.frame(Taxa = rownames(res), MeanDecreaseGini = res[, 1], stringsAsFactors = FALSE)
				res$Taxa <- nametable[res$Taxa, "name"]
				imp_sort <- dplyr::arrange(res, dplyr::desc(MeanDecreaseGini))
				rownames(imp_sort) <- imp_sort$Taxa
				imp_sort$pvalue <- pvalue_sub[as.character(imp_sort$Taxa)]
				self$res_rf <- imp_sort
				message('The result is stored in object$res_rf !')
			}
			if(grepl("lefse", method, ignore.case = TRUE)){
				class_taxa_median_sub <- lapply(res_class, function(x) x$med) %>% do.call(cbind, .) %>% .[, sel_taxa]
				all_class_pairs <- combn(unique(as.character(group_vec)), 2)
				# check the difference among subgroups
				if(!is.null(lefse_subgroup)){
					message("Start lefse subgroup biomarkers check for ", lefse_subgroup, " ...")
					all_sub_number <- as.data.table(sampleinfo)[, .N, by = c(group, lefse_subgroup)] %>% as.data.frame %>% .$N
					if(all(all_sub_number < lefse_min_subsam)){
						stop("All sample numbers for subgroups < ", lefse_min_subsam, "! Please consider use small lefse_min_subsam parameter!")
					}
					remove_list_total <- list()
					# for each group paires
					for(i in 1:ncol(all_class_pairs)){
						y1 <- all_class_pairs[, i]
						y1_sub_pairs <- expand.grid(unique(sampleinfo[sampleinfo[, group] == y1[1], lefse_subgroup]), 
							unique(sampleinfo[sampleinfo[, group] == y1[2], lefse_subgroup]), stringsAsFactors = FALSE) %>% t
						y1_sub_pairs <- y1_sub_pairs[, unlist(lapply(1:ncol(y1_sub_pairs), function(x){
							ifelse(any(c(sum(sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[1, x]) < lefse_min_subsam, 
								sum(sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[2, x]) < lefse_min_subsam)), FALSE, TRUE)
						})), drop = FALSE]
						if(ncol(y1_sub_pairs) == 0) next
						res_sub_total <- list()
						# check each subgroup pairs under fixed group pair condition
						for(j in 1:ncol(y1_sub_pairs)){
							y2 <- y1_sub_pairs[, j]
							abund_table_sub_y2 <- abund_table_sub[, c(rownames(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], ]), 
								rownames(sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], ]))]
							group_vec_sub2 <- c(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], group], 
								sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], group])
							res_sub <- lapply(seq_len(nrow(abund_table_sub_y2)), function(x) private$test_mark(abund_table_sub_y2[x,], group_vec_sub2))
							res_sub_total[[j]] <- res_sub
						}
						raw_median <- class_taxa_median_sub[y1, ] %>% {.[1, ] > .[2, ]} %>% as.vector
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
					})) %>% .[!is.na(.)] %>% .[!is.nan(.)] %>% max
				})
				res <- sapply(res, function(x) {log10(1 + abs(x)) * ifelse(x > 0, 1, -1)})
				res1 <- cbind.data.frame(Group = apply(class_taxa_median_sub, 2, function(x) rownames(class_taxa_median_sub)[which.max(x)]), 
					pvalue = pvalue_sub, LDA = res)
				res1 %<>% .[order(.$LDA, decreasing = TRUE), ]
				res1 <- cbind.data.frame(Taxa = rownames(res1), res1)
				message("Finished, minimum LDA score: ", range(res1$LDA)[1], " maximum LDA score: ", range(res1$LDA)[2])
				self$res_lefse <- res1
				message('The lefse result is stored in object$res_lefse !')
			}
			if(grepl("lefse|rf", method, ignore.case = TRUE)){
				if(grepl("lefse", method, ignore.case = TRUE)){
					res_abund <- reshape2::melt(rownames_to_column(abund_table_sub/lefse_norm, "Taxa"), id.vars = "Taxa")
				}else{
					res_abund <- reshape2::melt(rownames_to_column(abund_table_sub, "Taxa"), id.vars = "Taxa")				
				}
				colnames(res_abund) <- c("Taxa", "Sample", "Abund")
				res_abund <- suppressWarnings(dplyr::left_join(res_abund, rownames_to_column(sampleinfo), by = c("Sample" = "rowname")))
				res_abund <- microeco:::summarySE_inter(res_abund, measurevar = "Abund", groupvars = c("Taxa", group))
				colnames(res_abund)[colnames(res_abund) == group] <- "Group"
				self$res_abund <- res_abund
				message('The abundance is stored in object$res_abund !')
			}
			if(grepl("metastat|mseq", method, ignore.case = TRUE)){
				if(is.null(group_choose_paired)){
					all_name <- combn(unique(as.character(sampleinfo[, group])), 2)
				}else{
					all_name <- combn(unique(group_choose_paired), 2)
				}
				output <- data.frame()			
			}
			if(grepl("metastat", method, ignore.case = TRUE)){
				self$metastat_taxa_level <- metastat_taxa_level
				# transform data
				ranknumber <- which(colnames(dataset$tax_table) %in% metastat_taxa_level)
				abund <- dataset$otu_table
				tax <- dataset$tax_table[, 1:ranknumber, drop=FALSE]
				merged_taxonomy <- apply(tax, 1, paste, collapse="|")
				abund1 <- cbind.data.frame(Display = merged_taxonomy, abund) %>% 
					reshape2::melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
				abund1 <- data.table(abund1)[, sum_abund:=sum(Abundance), by=list(Display, Sample)] %>% 
					.[, c("Abundance"):=NULL] %>% setkey(Display, Sample) %>% unique() %>% as.data.frame()
				new_abund <- as.data.frame(data.table::dcast(data.table(abund1), Display~Sample, value.var= list("sum_abund"))) %>% 
					`row.names<-`(.[,1]) %>% .[,-1, drop = FALSE]
				new_abund <- new_abund[order(apply(new_abund, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]

				message("Total ", ncol(all_name), " paired group for calculation ...")
				for(i in 1:ncol(all_name)) {
					message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " vs "), " ...\n"))
					use_data <- new_abund[ , unlist(lapply(as.character(all_name[,i]), function(x) which(as.character(sampleinfo[, group]) %in% x)))]
					use_data %<>% .[!grepl("__$", rownames(.)), ]
					use_data <- use_data[apply(use_data, 1, sum) != 0, ]

					g <- sum(as.character(sampleinfo[, group]) == as.character(all_name[1, i])) + 1
					# calculate metastat
					res <- calculate_metastat(inputdata = use_data, g = g)
					add_name <- paste0(as.character(all_name[, i]), collapse = " vs ") %>% rep(., nrow(res))
					res <- cbind.data.frame(compare = add_name, res)
					output <- rbind.data.frame(output, res)
				}
				output %<>% dropallfactors(unfac2num = TRUE)
				self$res_metastat <- output
				message('The metastat result is stored in object$res_metastat !')
				self$res_metastat_group_matrix <- all_name
				message('The metastat group information is stored in object$res_metastat_group_matrix !')
			}
			if(grepl("mseq", method, ignore.case = TRUE)){
				if(!require(metagenomeSeq)){
					stop("metagenomeSeq package not installed")
				}
				message("Total ", ncol(all_name), " paired group for calculation ...")
				for(i in 1:ncol(all_name)) {
					message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " vs "), " ...\n"))
					use_dataset <- clone(dataset)
					use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[,i]), ]
					use_dataset$tidy_dataset()
					obj <- newMRexperiment(use_dataset$otu_table, phenoData= AnnotatedDataFrame(use_dataset$sample_table), featureData = AnnotatedDataFrame(use_dataset$tax_table))
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
					if (mseq_adjustMethod == "ihw-ubiquity" | mseq_adjustMethod == "ihw-abundance") {
						padj = MRihw(objres1, p, mseq_adjustMethod, 0.1)
					} else {
						padj = p.adjust(p, method = mseq_adjustMethod)
					}
					srt <- order(p, decreasing = FALSE)
					valid <- 1:length(padj)
					if (mseq_count > 0) {
						np = rowSums(objres1@counts)
						valid = intersect(valid, which(np >= mseq_count))
					}
					srt <- srt[which(srt %in% valid)]
					res <- cbind(tb[, 1:2], p)
					res <- cbind(res, padj)
					res <- as.data.frame(res[srt, ])
					colnames(res) <- c(colnames(tb)[1:2], "pvalues", "adjPvalues")
					res <- cbind.data.frame(feature = rownames(res), res)
					rownames(res) <- NULL
					add_name <- paste0(as.character(all_name[, i]), collapse = " vs ") %>% rep(., nrow(res))
					res <- cbind.data.frame(compare = add_name, res)
					output <- rbind.data.frame(output, res)
				}
				output %<>% dropallfactors(unfac2num = TRUE)
				self$res_mseq <- output
				message('The result is stored in object$res_mseq !')
				self$res_mseq_group_matrix <- all_name
				message('The group information is stored in object$res_mseq_group_matrix !')
			}
		},
		#' @description
		#' Plotting the abundance of differential taxa.
		#'
		#' @param method default NULL; "rf" or "lefse"; automatically check the method in the result.
		#' @param only_abund_plot default TRUE; if true, return only abundance plot; if false, return both indicator plot and abundance plot
		#' @param use_number default 1:10; vector, the taxa numbers used in the plot, 1:n.
		#' @param color_values colors for presentation.
		#' @param plot1_bar_color default "grey30"; the color for the plot 1.
		#' @param plot2_sig_color default "red"; the color for the significance in plot 2.
		#' @param plot2_sig_size default 1.5; the size for the significance in plot 2.
		#' @param axis_text_y default 12; the size for the y axis text.
		#' @param simplify_names default TRUE; whether use the simplified taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param group_order default NULL; a vector to order the legend in plot.
		#' @param plot2_barwidth default .9; the bar width in plot 2.
		#' @param add_significance default TRUE; whether add the significance asterisk; only available when only_abund_plot FALSE.
		#' @param use_se default TRUE; whether use SE in plot 2, if FALSE, use SD.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_diff_abund(use_number = 1:10)
		#' }
		plot_diff_abund = function(
			method = NULL,
			only_abund_plot = TRUE,
			use_number = 1:10,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			plot1_bar_color = "grey50",
			plot2_sig_color = "red",
			plot2_sig_size = 1.2,
			axis_text_y = 10,
			simplify_names = TRUE,
			keep_prefix = TRUE,
			group_order = NULL,
			plot2_barwidth = .9,
			add_significance = TRUE,
			use_se = TRUE
			){
			data2 <- self$res_abund
			if(is.null(method)){
				if(!is.null(self$res_lefse) & !is.null(self$res_rf)){
					method <- "lefse"
				}else{
					if(!is.null(self$res_lefse)){
						method <- "lefse"
					}else{
						if(!is.null(self$res_rf)){
							method <- "rf"
						}else{
							stop("No lefse or randomForest result found!")
						}
					}
				}
			}
			if(grepl("lefse", method, ignore.case = TRUE)){
				data1 <- self$res_lefse
				colnames(data1)[colnames(data1) == "LDA"] <- "Value"
				p1_xtile <- "LDA score"
			}else{
				if(grepl("rf", method, ignore.case = TRUE)){
					data1 <- self$res_rf
					colnames(data1)[colnames(data1) == "MeanDecreaseGini"] <- "Value"
					p1_xtile <- "MeanDecreaseGini"
				}else{
					stop("Provided method is not found, choose lefse or rf!")
				}
			}
			if(simplify_names == T){
				data1$Taxa %<>% gsub(".*\\|", "", .)
				data2$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				data1$Taxa %<>% gsub(".__", "", .)
				data2$Taxa %<>% gsub(".__", "", .)
			}
			if(length(use_number) > nrow(data1)){
				use_number <- 1:nrow(data1)
			}
			data1 %<>% .[use_number, ]
			data1$Taxa %<>% factor(., levels = rev(unique(as.character(.))))
			data2 %<>% .[.$Taxa %in% levels(data1$Taxa), ]
			data2$Taxa %<>% factor(., levels = levels(data1$Taxa))
			if(is.null(group_order)){
				data2$Group %<>% as.character %>% as.factor
			}else{
				data2$Group %<>% factor(., levels = rev(group_order))
			}

			p1 <- ggplot(data1, aes(x = Taxa, y = Value)) +
				theme_bw() +
				geom_bar(stat = "identity", fill = plot1_bar_color) +
				coord_flip() +
				xlab("") +
				ylab(p1_xtile) +
				theme(panel.border = element_blank(), panel.background=element_rect(fill="white")) +
				theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) + #, panel.grid.minor.x = element_blank())
				theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = axis_text_y, color = "black")) +
				theme(plot.margin = unit(c(.1, 0, .1, 0), "cm"))

			p2 <- ggplot(data2, aes(x=Taxa, y=Mean, color = Group, fill = Group, group = Group)) +
				geom_bar(stat="identity", position = position_dodge(), width = plot2_barwidth)
			if(use_se == T){
				p2 <- p2 + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.45, position=position_dodge(plot2_barwidth), color = "black")
			}else{
				p2 <- p2 + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.45, position=position_dodge(plot2_barwidth), color = "black")
			}
			p2 <- p2 + theme_bw() +
				coord_flip() +
				scale_color_manual(values=color_values) +
				scale_fill_manual(values=color_values) +
				ylab("Relative abundance") +
				theme(legend.position = "right") +
				theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.border = element_blank(), 
					panel.background=element_rect(fill="white")) +
				theme(axis.title = element_text(size = 17)) +
				guides(fill=guide_legend(reverse=TRUE, ncol=1), color = FALSE)
			
			if(only_abund_plot == T){
				p2 <- p2 + theme(axis.title.y=element_blank(), axis.text.y = element_text(size = axis_text_y, color = "black")) + 
					theme(plot.margin = unit(c(.1, 0, .1, 0), "cm"))
				return(p2)
			}else{
				if(add_significance == T){
					Significance <- rev(as.character(cut(data1$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))))
					p2 <- p2 + scale_x_discrete(labels=Significance) +
						theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), 
						axis.text.y = element_text(color = plot2_sig_color, size = rel(plot2_sig_size))) +
						theme(plot.margin = unit(c(.1, 0, .1, .8), "cm"))
				}else{
					p2 <- p2 + theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
						theme(plot.margin = unit(c(.1, 0, .1, 0), "cm"))
				}
				return(list(p1 = p1, p2 = p2))
			}
		},
		#' @description
		#' Bar plot for LDA score.
		#'
		#' @param use_number default 1:10; vector, the taxa numbers used in the plot, 1:n.
		#' @param color_values colors for presentation.
		#' @param LDA_score default NULL; numeric value as the threshold, such as 2, limited with use_number.
		#' @param simplify_names default TRUE; whether use the simplified taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param group_order default NULL; a vector to order the legend in plot.
		#' @param axis_text_y default 12; the size for the y axis text.
		#' @param plot_vertical default TRUE; whether use vertical bar plot or horizontal.
		#' @param ... parameters pass to \code{\link{geom_bar}}
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_lefse_bar(LDA_score = 4)
		#' }
		plot_lefse_bar = function(
			use_number = 1:10,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			LDA_score = NULL,
			simplify_names = TRUE,
			keep_prefix = TRUE,
			group_order = NULL,
			axis_text_y = 12,
			plot_vertical = TRUE,
			...
			){
			use_data <- self$res_lefse
			if(simplify_names == T){
				use_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				use_data$Taxa %<>% gsub(".__", "", .)
			}
			if(is.null(LDA_score)){
				sel_num <- use_number
			}else{
				sel_num <- sum(use_data$LDA > LDA_score)
				if(sel_num == 0){
					stop("Too large LDA_score provided, no data selected!")
				}
				sel_num <- 1:sel_num
			}
			if(length(sel_num) > nrow(use_data)){
				sel_num <- 1:nrow(use_data)
			}
			use_data %<>% .[sel_num, ]
			if(is.null(group_order)){
				use_data$Group %<>% as.character %>% as.factor
			}else{
				use_data$Group %<>% factor(., levels = group_order)
			}
			# rearrange orders
			if(length(levels(use_data$Group)) == 2){
				use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$Group), function(x){
					if(x == levels(use_data$Group)[1]){
						use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$LDA, decreasing = TRUE), "Taxa"]
					}else{
						use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$LDA, decreasing = FALSE), "Taxa"]
					}
				})))))
				use_data[use_data$Group == levels(use_data$Group)[2], "LDA"] %<>% {. * -1}
			}else{
				use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$Group), function(x){
					use_data[as.character(use_data$Group) %in% x, ] %>% .[order(.$LDA, decreasing = TRUE), "Taxa"]
				})))))
			}
			
			p <- ggplot(use_data, aes(x = Taxa, y = LDA, color = Group, fill = Group, group = Group)) +
				geom_bar(stat="identity", position = position_dodge(), ...) +
				theme_bw() +
				scale_color_manual(values=color_values) +
				scale_fill_manual(values=color_values) +
				ylab("LDA score") +
				xlab("") +
				theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = axis_text_y, color = "black")) +
				theme(axis.text.x = element_text(size = 12)) +
				theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()) +
				theme(panel.border = element_blank()) +
				theme(axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square"))
				
			if(plot_vertical == T){
				p <- p + coord_flip()
			}
			p
		},
		#' @description
		#' Plot the cladogram for LEfSe result similar with the python version. Codes are modified from microbiomeMarker 
		#'
		#' @param color default RColorBrewer::brewer.pal(8, "Dark2"); color used in the plot.
		#' @param use_taxa_num default 200; integer; The taxa number used in the background tree plot; select the taxa according to the mean abundance 
		#' @param filter_taxa default NULL; The mean relative abundance used to filter the taxa with low abundance
		#' @param use_feature_num default NULL; integer; The feature number used in the plot; select the features according to the LDA score
		#' @param group_order default NULL; a vector to order the legend in plot.
		#' @param clade_label_level default 4; the taxonomic level for marking the label with letters, root is the largest
		#' @param select_show_labels default NULL; character vector; The features to show in the plot with full label names, not the letters
		#' @param only_select_show default FALSE; whether only use the the select features in the parameter select_show_labels
		#' @param sep default "|"; the seperate character in the taxonomic information
		#' @param branch_size default 0.2; numberic, size of branch
		#' @param alpha default 0.2; shading of the color
		#' @param clade_label_size default 0.7; size for the clade label
		#' @param node_size_scale default 1; scale for the node size
		#' @param node_size_offset default 1; offset for the node size
		#' @param annotation_shape default 22; shape used in the annotation legend
		#' @param annotation_shape_size default 5; size used in the annotation legend
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_lefse_cladogram(use_taxa_num = 100, use_feature_num = 30, select_show_labels = NULL)
		#' }
		plot_lefse_cladogram = function(
			color = RColorBrewer::brewer.pal(8, "Dark2"),
			use_taxa_num = 200,
			filter_taxa = NULL,
			use_feature_num = NULL,
			group_order = NULL,
			clade_label_level = 4,
			select_show_labels = NULL,
			only_select_show = FALSE,
			sep = "|",
			branch_size = 0.2,
			alpha = 0.2,
			clade_label_size = 0.7,
			node_size_scale = 1,
			node_size_offset = 1,
			annotation_shape = 22,
			annotation_shape_size = 5
			){
			abund_table <- self$abund_table
			marker_table <- self$res_lefse %>% dropallfactors

			if(!is.null(use_feature_num)){
				marker_table %<>% .[1:use_feature_num, ]
			}
			if(only_select_show == T){
				marker_table %<>% .[.$Taxa %in% select_show_labels, ]
			}
			color <- color[1:length(unique(marker_table$Group))]

			# filter the taxa with unidentified classification or with space, in case of the unexpected error in the following operations
			abund_table %<>% {.[!grepl("\\|.__\\|", rownames(.)), ]}
			abund_table %<>% {.[!grepl("\\s", rownames(.)), ]}

			if(!is.null(use_taxa_num)){
				abund_table %<>% .[names(sort(apply(., 1, mean), decreasing = TRUE)[1:use_taxa_num]), ]
			}
			if(!is.null(filter_taxa)){
				abund_table %<>% .[apply(., 1, mean) > (self$lefse_norm * filter_taxa), ]
			}
			abund_table %<>% .[sort(rownames(abund_table)), ]

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
			nodes <- c("r__Root", nodes)
			# levels used for extend of clade label
			label_levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
				factor(levels = rev(unlist(lapply(taxa_split, function(x) gsub("(.)__.*", "\\1", x))) %>% .[!duplicated(.)]))

			nodes_parent <- purrr::map_chr(taxa_split, ~ .x[length(.x) - 1])
			# root must be a parent node
			nodes_parent <- c("root", nodes_parent)

			## add index for nodes
			is_tip <- !nodes %in% nodes_parent
			index <- vector("integer", length(is_tip))
			index[is_tip] <- 1:sum(is_tip)
			index[!is_tip] <- (sum(is_tip)+1):length(is_tip)

			edges <- cbind(parent = index[match(nodes_parent, nodes)], child = index)
			edges <- edges[!is.na(edges[, 1]), ]
			# not label the tips
			node_label <- nodes[!is_tip]
			phylo <- structure(list(edge = edges, node.label = node_label, tip.label = nodes[is_tip], edge.length = rep(1, nrow(edges)), Nnode = length(node_label)),
				class = "phylo")
			mapping <- data.frame(node = index, abd = c(100, tree_table$abd), node_label = nodes, stringsAsFactors = FALSE)
			mapping$node_class <- label_levels
			tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
			tree <- ggtree::ggtree(tree, size = 0.2, layout = 'circular')
			# color legend order settings
			if(is.null(group_order)){
				color_groups <- marker_table$Group %>% as.character %>% as.factor %>% levels
			}else{
				color_groups <- group_order
			}
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
			rownames(hilights_df) <- hilights_df$enrich_group
			hilights_df <- hilights_df[color_groups, ]
			# make sure the right order in legend
			hilights_df$enrich_group %<>% factor(., levels = color_groups)

			# add legend
			tree <- tree + geom_rect(aes_(xmin = ~x, xmax = ~x, ymax = ~y, ymin = ~y, fill = ~enrich_group), data = hilights_df, inherit.aes = FALSE) +
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
				fontsize = clade_label_size + log(.data$level + 20),
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
			}

			clade_label_g <- purrr::pmap(clade_label_new, ggtree::geom_cladelabel)
			p <- purrr::reduce(clade_label_g, `+`, .init = tree)

			# if letters are used, add guide labels
			if(ind_num > 0){
				guide_label <- clade_label_new[ind, ] %>%
					dplyr::mutate(color = annotation_info$color[match(.data$label_raw, annotation_info$label)])
				p <- p + geom_point(data = guide_label, inherit.aes = FALSE, aes_(x = 0, y = 0, shape = ~label_legend), size = 0, stroke = 0) +
					scale_shape_manual(values = rep(annotation_shape, nrow(guide_label))) +
					guides(shape = guide_legend(override.aes = list(size = annotation_shape_size, shape = annotation_shape, fill = guide_label$color)))
			}
			p <- p + theme(legend.position = "right", legend.title = element_blank())
			p
		},
		#' @description
		#' Bar plot for metastat.
		#'
		#' @param use_number default 1:10; vector, the taxa numbers used in the plot, 1:n.
		#' @param color_values colors for presentation.
		#' @param qvalue default .05; numeric value as the threshold of q value.
		#' @param choose_group default 1; which column in res_metastat_group_matrix will be used.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Group")
		#' t1$plot_metastat(use_number = 1:10, qvalue = 0.05, choose_group = 1)
		#' }
		plot_metastat = function(use_number = 1:10, color_values = RColorBrewer::brewer.pal(8, "Dark2"), qvalue = 0.05, choose_group = 1
			){
			use_data <- self$res_metastat
			group_char <- self$res_metastat_group_matrix[, choose_group]
			use_group <- use_data[,1] %>% unique %>% .[1]
			use_data %<>% .[.$qvalue < qvalue & .[,1] == use_group, ]
			use_data[,2] %<>% gsub(paste0(".*", tolower(substr(self$metastat_taxa_level, 1, 1)), "__(.*)$"), "\\1", .)
			use_data %<>% .[.[,2] != "", ]
			if(nrow(use_data) > length(use_number)){
				use_data %<>% .[use_number, ]
			}
			plot_data <- data.frame(taxa = rep(use_data[,2], 2), Mean = c(use_data[,3], use_data[,6]), SE = c(use_data[,5], use_data[,8]), 
				Group = rep(group_char, times = 1, each = nrow(use_data)))
			plot_data$taxa %<>% factor(., levels = unique(.))
			p <- ggplot(plot_data, aes(x=taxa, y=Mean, color = Group, fill = Group, group = Group)) +
				geom_bar(stat="identity", position = position_dodge()) +
				geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.45, position=position_dodge(.9), color = "black") +
				theme_classic() +
				scale_color_manual(values=color_values) +
				scale_fill_manual(values=color_values) +
				ylab("Relative abundance") +
				theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 9), legend.position = "top") +
				theme(axis.title = element_text(size = 15)) +
				xlab(self$metastat_taxa_level)
			p
		},
		#' @description
		#' Print the trans_diff object.
		print = function() {
			cat("trans_diff class:\n")
			if(!is.null(self$res_rf)) cat("Randomeforest has been calculated \n")
			if(!is.null(self$res_lefse)) cat("Lefse has been calculated \n")
			if(!is.null(self$res_metastat)) cat("Metastat has been calculated \n")
			invisible(self)
		}
		),
	private = list(
		# group test in lefse or rf
		test_mark = function(dataset, group, nonpara = TRUE, para = "anova", min_num_nonpara = 1){
			d1 <- as.data.frame(t(dataset))
			num_vals <- as.numeric(d1[,1])
			#shapiro test function  - check normality
			if(length(unique((num_vals[!is.na(num_vals)]))) > 3 
			   & length(unique((num_vals[!is.na(num_vals)]))) < 5000 ){
			  d1.shapiro <- shapiro.test(num_vals)$p.value
			  if (d1.shapiro > 0.05){
			    nonpara = F
			  }
			}else{
			  message("It was not possible to verify the normality of the taxa ", colnames(d1)[1], " !")
			} 

			if(nonpara == T){
				if(any(table(as.character(group))) < min_num_nonpara){
					list(p_value = NA, med = NA)
				}else{
					if(length(unique(as.character(group))) == 2){
						res1 <- wilcox.test(d1[,1], g=group)
					}else{
						res1 <- kruskal.test(d1[,1], g=group)
					}
					if(is.nan(res1$p.value)){
						res1$p.value <- 1
					}
					med <- tapply(d1[,1], group, median) %>% as.data.frame
					colnames(med) <- colnames(d1)
					list(p_value = res1$p.value, med = med)		
				}
			}else{
				if(para == "anova"){
				  colnames(d1)[1] <- "Abundance" #dataframe of 1 col
					d2 <- cbind.data.frame(d1, Group = group)
					res1 <- aov(Abundance ~ Group, d2)
					pvalue <- as.numeric(unlist(summary(res1))[9])
					list(p_value = pvalue) #default return is list
				}
			}
		},
		# generate the cladogram annotation table
		generate_cladogram_annotation = function(marker_table, tree, color, color_groups, sep = "|") {
			use_marker_table <- marker_table
			feature <- use_marker_table$Taxa
			label <- strsplit(feature, split = sep, fixed = TRUE) %>% purrr::map_chr(utils::tail, n =1)
			plot_color <- use_marker_table$Group %>% as.character
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
		get_offset = function(x) {(x*0.2+0.2)^2}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)

