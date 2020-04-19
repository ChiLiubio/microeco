#' Create trans_diff object for the difference tests on the taxonomic abundance and plotting.
#'
#' This class is a wrapper for a series of differential abundance test and indicator analysis methods.
#' The functions in this class include \code{\link{plot_diff_abund}}, \code{\link{plot_lefse_bar}}, \code{\link{plot_metastat}}
#'
#'
#' @param dataset the object of \code{\link{microtable}} Class.
#' @param method default "lefse"; one of c("lefse", "rf", "metastat").
#' @param group default NULL; sample group used for main comparision.
#' @param lefse_subgroup default NULL; sample sub group used in lefse.
#' @param alpha default .05; significance threshold.
#' @param lefse_min_subsam default 10; sample numbers required in the subgroup test.
#' @param lefse_norm default 1000000; scale value.
#' @param lefse_nresam default .6667; sample number ratio used in each bootstrap.
#' @param boots default 30; bootstrap test number for lefse or rf.
#' @param rf_taxa_level default "all"; use all taxonomic rank data, if want to test a specific rank, provide taxonomic rank name, such as "Genus".
#' @param rf_ntree default 1000; see ntree in \code{\link{randomForest}}.
#' @param metastat_taxa_level default "Genus"; taxonomic rank level used in metastat test.
#' @param metastat_group_choose default NULL; a vector used for selecting the required groups for testing.
#' @return res_rf res_lefse res_abund or res_metastat in trans_diff object.
#' @examples 
#' t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group")
#' @export
trans_diff <- R6Class(classname = "trans_diff",
	public = list(
		initialize = function(dataset = NULL, method = c("lefse", "rf", "metastat")[1],
			group = NULL, alpha = 0.05,
			lefse_subgroup = NULL, lefse_min_subsam = 10, lefse_norm = 1000000, lefse_nresam = 0.6667, boots = 30,
			rf_taxa_level = "all", rf_ntree = 1000,
			metastat_taxa_level = "Genus", metastat_group_choose = NULL
			){
			if(is.null(dataset)){
				stop("No dataset provided!")
			}
			if(is.null(dataset$taxa_abund)){
				stop("Please first calculate taxa_abund use cal_abund function!")
			}
			sampleinfo <- dataset$sample_table
			if(grepl("lefse|rf", method, ignore.case = TRUE)){
				if(grepl("lefse", method, ignore.case = TRUE)){
					abund_table <- do.call(rbind, unname(lapply(dataset$taxa_abund, function(x) x * lefse_norm)))
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
				res_class <- lapply(seq_len(nrow(abund_table)), function(x) private$test_mark(abund_table[x,], group_vec))
				pvalue <- unlist(lapply(res_class, function(x) x$p_value))
				pvalue[is.nan(pvalue)] <- 1
				# select significant taxa
				sel_taxa <- pvalue < alpha
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
				for(num in seq_len(boots)){
					# resampling
					sample_names_resample <- rownames(predictors)[base::sample(1:nrow(predictors), size = ceiling(nrow(predictors) * lefse_nresam))]
					predictors_sub <- predictors[sample_names_resample, ]
					sampleinfo_resample <- sampleinfo[sample_names_resample, ]
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
			}
			if(grepl("lefse", method, ignore.case = TRUE)){
				class_taxa_median_sub <- lapply(res_class, function(x) x$med) %>% do.call(cbind, .) %>% .[, sel_taxa]
				all_class_pairs <- combn(unique(as.character(group_vec)), 2)
				# check the difference among subgroups
				if(!is.null(lefse_subgroup)){
					remove_list_total <- list()
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
					remove_list_total %<>% do.call(cbind, .) %>% apply(., 1, any)
					abund_table_sub %<>% .[!remove_list_total, ]
					pvalue_sub %<>% .[!remove_list_total]
					class_taxa_median_sub %<>% .[, !remove_list_total]
				}
				res_lda <- list()
				# bootstrap default 30 times
				for(num in seq_len(boots)){
					res_lda_pair <- list()
					# resampling samples
					sample_names_resample <- colnames(abund_table_sub)[base::sample(1:ncol(abund_table_sub), size = ceiling(ncol(abund_table_sub) * lefse_nresam))]
					abund_table_sub_resample <- abund_table_sub[, sample_names_resample]
					sampleinfo_resample <- sampleinfo[sample_names_resample, ]
					# make sure the groups and samples numbers right
					if(length(unique(sampleinfo_resample[, group])) != length(unique(sampleinfo[, group])) | min(table(sampleinfo_resample[, group])) < 2){
						next
					}
					# cycle all paired groups
					for(i in seq_len(ncol(all_class_pairs))){
						sel_samples <- sampleinfo_resample[, group] %in% all_class_pairs[, i]
						abund_table_sub_lda <- abund_table_sub_resample[, sel_samples]
						abund_table_sub_lda %<>% .[apply(., 1, sd) > 1.0e-10, ]
						group_vec_lda <- sampleinfo_resample[sel_samples, group] %>% as.character %>% as.factor
						if(is.null(lefse_subgroup)){
							abund1 <- cbind.data.frame(t(abund_table_sub_lda), Group = group_vec_lda)
						}else{
							subgroup_vec <- sampleinfo_resample[sel_samples, lefse_subgroup] %>% as.character %>% as.factor
							# consider subgroup as a independent variable
							abund1 <- cbind.data.frame(t(abund_table_sub_lda), Group = group_vec_lda, lefse_subgroup = subgroup_vec)
						}
						mod1 <- MASS::lda(Group ~ ., abund1, tol = 1.0e-10)
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
					res_lda[[num]] <- res_lda_pair
				}
				# same with the python lefse
				res <- sapply(rownames(abund_table_sub), function(k){
					unlist(lapply(seq_len(ncol(all_class_pairs)), function(p){
						unlist(lapply(res_lda, function(x){ x[[p]][k]})) %>% .[!is.na(.)] %>% mean
					})) %>% .[!is.nan(.)] %>% max
				})
				res <- sapply(res, function(x) {log10(1 + abs(x)) * ifelse(x > 0, 1, -1)})
				res1 <- cbind.data.frame(Group = apply(class_taxa_median_sub, 2, function(x) rownames(class_taxa_median_sub)[which.max(x)]), pvalue = pvalue_sub, LDA = res)
				res1 %<>% .[order(.$LDA, decreasing = TRUE), ]
				res1 <- cbind.data.frame(Taxa = rownames(res1), res1)
				self$res_lefse <- res1
			}
			if(grepl("lefse|rf", method, ignore.case = TRUE)){
				if(grepl("lefse", method, ignore.case = TRUE)){
					tran_abund <- reshape2::melt(rownames_to_column(abund_table_sub/lefse_norm, "Taxa"), id.vars = "Taxa")
				}else{
					tran_abund <- reshape2::melt(rownames_to_column(abund_table_sub, "Taxa"), id.vars = "Taxa")				
				}
				colnames(tran_abund) <- c("Taxa", "Sample", "Abund")
				tran_abund <- suppressWarnings(left_join(tran_abund, rownames_to_column(sampleinfo), by = c("Sample" = "rowname")))
				tran_abund <- summarySE(tran_abund, measurevar = "Abund", groupvars = c("Taxa", group))
				self$res_abund <- tran_abund
			}
			
			if(grepl("metastat", method, ignore.case = TRUE)){
				self$metastat_taxa_level <- metastat_taxa_level
				# transform data
				ranknumber <- which(colnames(dataset$tax_table) %in% metastat_taxa_level)
				abund <- dataset$otu_table
				tax <- dataset$tax_table[, 1:ranknumber, drop=FALSE]
				merged_taxonomy <- apply(tax, 1, paste, collapse="|")
				abund1 <- cbind.data.frame(Display = merged_taxonomy, abund) %>% reshape2::melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
				abund1 <- data.table(abund1)[, sum_abund:=sum(Abundance), by=list(Display, Sample)] %>% .[, c("Abundance"):=NULL] %>% setkey(Display, Sample) %>% unique() %>% as.data.frame()
				new_abund <- as.data.frame(data.table::dcast(data.table(abund1), Display~Sample, value.var= list("sum_abund"))) %>% `row.names<-`(.[,1]) %>% .[,-1, drop = FALSE]
				new_abund <- new_abund[order(apply(new_abund, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]
				if(is.null(metastat_group_choose)){
					all_name <- combn(unique(as.character(sampleinfo[, group])), 2)
				}else{
					all_name <- combn(unique(metastat_group_choose), 2)
				}
				output <- data.frame()
				for(i in 1:ncol(all_name)) {
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
				self$res_metastat <- output
				self$res_metastat_group_matrix <- all_name				
			}
		},
		plot_diff_abund = function(method = NULL, only_abund_plot = TRUE, use_number = 1:10, color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			plot1_bar_color = "grey50", plot2_sig_color = "red", plot2_sig_size = 1.2,
			axis_text_y = 10, 
			simplify_names = TRUE, keep_prefix = TRUE, group_order = NULL, 
			plot2_barwidth = .9, add_significance = TRUE, use_se = TRUE
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
				theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.border = element_blank(), panel.background=element_rect(fill="white")) + #, panel.grid.minor.x = element_blank())
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
						theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(color = plot2_sig_color, size = rel(plot2_sig_size))) +
						theme(plot.margin = unit(c(.1, 0, .1, .8), "cm"))
				}else{
					p2 <- p2 + theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
						theme(plot.margin = unit(c(.1, 0, .1, 0), "cm"))
				}
				return(list(p1 = p1, p2 = p2))
			}
		},
		plot_lefse_bar = function(use_number = 1:10, color_values = RColorBrewer::brewer.pal(8, "Dark2"), LDA_score = NULL,
			simplify_names = TRUE, keep_prefix = TRUE, group_order = NULL, axis_text_y = 12, plot_vertical = TRUE, ...
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
		plot_metastat = function(use_number = 1:10, qvalue = 0.05, choose_group = 1, color_values = RColorBrewer::brewer.pal(8, "Dark2")
			){
			use_data <- self$res_metastat %>% dropallfactors
			group_char <- self$res_metastat_group_matrix[, choose_group]
			use_group <- use_data[,1] %>% unique %>% .[1]
			use_data %<>% .[.$qvalue < qvalue & .[,1] == use_group, ]
			use_data[,2] %<>% gsub(paste0(".*", tolower(substr(self$metastat_taxa_level, 1, 1)), "__(.*)$"), "\\1", .)
			use_data %<>% .[.[,2] != "", ]
			if(nrow(use_data) > length(use_number)){
				use_data %<>% .[use_number, ]
			}
			plot_data <- data.frame(taxa = rep(use_data[,2], 2), Mean = c(use_data[,3], use_data[,6]), SE = c(use_data[,5], use_data[,8]), Group = rep(group_char, times = 1, each = nrow(use_data)))
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
		print = function(...) {
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
			if(nonpara == T){
				if(any(table(as.character(group))) < min_num_nonpara){
					list(p_value = NA, med = NA)
				}else{
					if(length(unique(as.character(group))) == 2){
						res1 <- wilcox.test(d1[,1], g=group)
					}else{
						res1 <- kruskal.test(d1[,1], g=group)
					}
					med <- tapply(d1[,1], group, median) %>% as.data.frame
					colnames(med) <- colnames(d1)
					list(p_value = res1$p.value, med = med)		
				}
			}else{
				if(para == "anova"){
					rownames(d1)[1] <- "Abundance"
					d2 <- cbind.data.frame(d1, Group = group)
					res1 <- aov(Abundance ~ Group, d2)
					pvalue <- as.numeric(unlist(summary(res1))[9])
					pvalue
				}
			}
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)


#' Plotting the abundance of differencial taxa.
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
#' @return a list with two ggplot.
#' @examples
#' t1$plot_diff_abund(use_number = 1:10)

plot_diff_abund <- function(method = NULL, only_abund_plot = TRUE, use_number = 1:10, color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			plot1_bar_color = "grey50", plot2_sig_color = "red", plot2_sig_size = 1.2,
			axis_text_y = 10, 
			simplify_names = TRUE, keep_prefix = TRUE, group_order = NULL, 
			plot2_barwidth = .9, add_significance = TRUE, use_se = TRUE){
	dataset$plot_diff_abund()
}



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
#' t1$plot_lefse_bar(LDA_score = 2)

plot_lefse_bar <- function(use_number = 1:10, color_values = RColorBrewer::brewer.pal(8, "Dark2"), LDA_score = NULL,
			simplify_names = TRUE, keep_prefix = TRUE, group_order = NULL, axis_text_y = 12, plot_vertical = TRUE, ...){
	dataset$plot_lefse_bar()
}

#' Bar plot metastat.
#'
#' @param use_number default 1:10; vector, the taxa numbers used in the plot, 1:n.
#' @param color_values colors for presentation.
#' @param qvalue default .05; numeric value as the threshold of q value.
#' @param choose_group default 1; which column in res_metastat_group_matrix will be used.
#' @return ggplot.
#' @examples
#' t1$plot_metastat(use_number = 1:10, qvalue = 0.05, choose_group = 1)
plot_metastat <- function(use_number = 1:10, qvalue = 0.05, choose_group = 1, color_values = RColorBrewer::brewer.pal(8, "Dark2")){
	dataset$plot_metastat()
}









