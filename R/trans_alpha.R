#' @title Create \code{trans_alpha} object for alpha diversity statistics and plotting.
#'
#' @description
#' This class is a wrapper for a series of alpha diversity related analysis, including the statistics and plotting based on 
#' An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and Paul et al. (2013) <doi:10.1371/journal.pone.0061217>.
#'
#' @export
trans_alpha <- R6Class(classname = "trans_alpha",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param group default NULL; the sample column used for the statistics; If provided, can return \code{data_stat}.
		#' @param by_group default NULL; perform the differential test among groups within each by_group.
		#' @param order_x default NULL; a \code{sample_table} column name or a vector containg sample names; if provided, order samples by using \code{factor}.
		#' @return \code{data_alpha} and \code{data_stat} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Group")
		#' }
		initialize = function(dataset = NULL, group = NULL, by_group = NULL, order_x = NULL) {
			if(is.null(dataset)){
				message("Parameter dataset not provided. Please run the functions with your other customized data!")
				self$data_alpha <- NULL
			}else{
				if(is.null(dataset$alpha_diversity)){
					message("The alpha_diversity in dataset not found! Calculate it automatically ...")
					dataset$cal_alphadiv()
				}
				data_alpha <- dataset$alpha_diversity %>% 
					cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
					.[, !grepl("^se", colnames(.))] %>%
					# to long format
					reshape2::melt(id.vars = "Sample") %>%
					`colnames<-`(c("Sample", "Measure", "Value")) %>%
					dplyr::left_join(., rownames_to_column(dataset$sample_table), by = c("Sample" = "rowname"))
				if(!is.null(order_x)){
					if(length(order_x == 1)){
						data_alpha$Sample %<>% factor(., levels = unique(dataset$sample_table[, order_x]))
					} else {
						data_alpha$Sample %<>% factor(., levels = order_x)
					}
				}
				self$data_alpha <- data_alpha
				message('The transformed diversity data is stored in object$data_alpha ...')
			}
			if(! is.null(by_group)){
				if(! by_group %in% colnames(data_alpha)){
					stop("Provided by_group: ", by_group, " is not found in dataset$sample_table!")
				}
			}
			if(! is.null(group)){
				if(is.null(dataset)){
					stop("Parameter dataset not provided, but group is provided!")
				}
				if(! group %in% colnames(data_alpha)){
					stop("Provided group: ", group, " is not found in dataset$sample_table!")
				}
				self$data_stat <- microeco:::summarySE_inter(data_alpha, measurevar = "Value", groupvars = c(group, "Measure"))
				message('The group statistics are stored in object$data_stat ...')
			}else{
				self$data_stat <- NULL
			}
			self$group <- group
			self$by_group <- by_group
		},
		#' @description
		#' Test the difference of alpha diversity across groups.
		#'
		#' @param method default "KW"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons, see \code{dunnTest} function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'anova'}}{Duncan's multiple range test for anova}
		#'   }
		#' @param measure default NULL; a vector; if null, all indexes will be calculated; see names of \code{microtable$alpha_diversity}, 
		#' 	 e.g. Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage, PD.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of \code{p.adjust} function for available options; 
		#'    NULL can disuse the p value adjustment.
		#' @param anova_set default NULL; specified group set for anova, such as \code{'block + N*P*K'}, see \code{\link{aov}}.
		#' @param ... parameters passed to \code{kruskal.test} or \code{wilcox.test} function (\code{method = "KW"}) or \code{dunnTest} function of \code{FSA} package 
		#'   (\code{method = "KW_dunn"}) or \code{agricolae::duncan.test} (\code{method = "anova"}).
		#' @return \code{res_diff} in object. A \code{data.frame} generally. A list for anova when anova_set is assigned.
		#'   In the data frame, 'Group' column means that the group has the maximum median or mean value across the test groups;
		#'   For non-parametric methods, maximum median value; For t.test, maximum mean value.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(method = "KW")
		#' t1$cal_diff(method = "KW_dunn")
		#' t1$cal_diff(method = "anova")
		#' }
		cal_diff = function(method = c("KW", "KW_dunn", "wilcox", "t.test", "anova")[1], measure = NULL, p_adjust_method = "fdr", anova_set = NULL, ...){
			group <- self$group
			# 'by_group' for test inside each by_group instead of all groups in 'group'
			by_group <- self$by_group
			data_alpha <- self$data_alpha
			if(is.null(measure)){
				measure <- unique(as.character(data_alpha$Measure))
			}else{
				if(! all(measure %in% as.character(data_alpha$Measure))){
					stop("One or more measures input not in the object$data_alpha! Please check the input!")
				}
			}
			method <- match.arg(method, c("KW", "KW_dunn", "wilcox", "t.test", "anova"))
			if(is.null(by_group)){
				if(method %in% c("wilcox", "t.test") & length(unique(as.character(data_alpha[, group]))) > 10){
					stop("There are too many groups to do paired comparisons! please use method = 'KW' or 'KW_dunn' or 'anova' !")
				}
			}
			if(!is.null(anova_set)){
				method <- "anova"
			}
			if(!is.null(by_group)){
				all_bygroups <- as.character(data_alpha[, by_group])
				unique_bygroups <- unique(all_bygroups)
			}
			if(method %in% c("KW", "wilcox", "t.test")){
				res_list <- list()
				res_list$comnames <- c()
				res_list$p_value <- c()
				res_list$measure_use <- c()
				res_list$test_method <- c()
				res_list$max_group <- c()
				res_list$group_by <- c()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, c(group, "Value")]
						res_list <- private$kwwitt_test(method = method, input_table = div_table, group = group, res_list = res_list, measure = k, ...)
					}else{
						for(each_group in unique_bygroups){
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, c(by_group, group, "Value")]
							if(length(unique(as.character(div_table[, group]))) < 2){
								next
							}
							res_list <- private$kwwitt_test(method = method, input_table = div_table, group = group, res_list = res_list, 
								group_by = each_group, measure = k, ...)
						}
					}
				}
				if(is.null(p_adjust_method)){
					p_value_adjust <- res_list$p_value
				}else{
					p_value_adjust <- p.adjust(res_list$p_value, method = p_adjust_method)
				}
				if(is.null(by_group)){
					compare_result <- data.frame(res_list$comnames, res_list$measure_use, res_list$test_method, res_list$max_group, res_list$p_value, p_value_adjust)
					colnames(compare_result) <- c("Comparison", "Measure", "Test_method", "Group", "P.unadj", "P.adj")
				}else{
					compare_result <- data.frame(res_list$comnames, res_list$group_by, res_list$measure_use, res_list$test_method, 
						res_list$max_group, res_list$p_value, p_value_adjust)
					colnames(compare_result) <- c("Comparison", "by_group", "Measure", "Test_method", "Group", "P.unadj", "P.adj")
				}
				compare_result$Significance <- cut(compare_result$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
			}
			if(method == "KW_dunn"){
				if(is.null(by_group)){
					if(length(unique(data_alpha[, group])) == 2){
						stop("There are only 2 groups. Please select other method instead of KW_dunn !")
					}
				}
				compare_result <- data.frame()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, c(group, "Value")]
						res_table <- private$kdunn_test(input_table = div_table, group = group, measure = k, ...)
						compare_result %<>% rbind(., res_table)
					}else{
						for(each_group in unique_bygroups){
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, c(by_group, group, "Value")]
							if(length(unique(as.character(div_table[, group]))) < 3){
								next
							}
							res_table <- private$kdunn_test(input_table = div_table, group = group, group_by = each_group, measure = k, ...)
							compare_result %<>% rbind(., res_table)
						}
					}
				}
				compare_result$Significance <- cut(compare_result$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "ns"))
			}
			if(method == "anova"){
				if(is.null(anova_set)){
					compare_result <- NULL
				}else{
					compare_result <- list()
				}
				if(is.null(anova_set)){
					if(is.null(by_group)){
						for(k in measure){
							div_table <- data_alpha[data_alpha$Measure == k, ]
							compare_result <- private$anova_test(input_table = div_table, group = group, measure = k, compare_result = compare_result, ...)
						}
						colnames(compare_result)[1] <- "Group"
					}else{
						for(each_group in unique_bygroups){
							test <- data_alpha[all_bygroups == each_group, ]
							if(length(unique(as.character(test[, group]))) < 2){
								message("Skip the by_group: ", each_group, " as groups number < 2!")
								next
							}
							tmp_measure <- NULL
							for(k in measure){
								# regenerate table for test
								div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, ]
								tmp_measure <- private$anova_test(input_table = div_table, group = group, measure = k, compare_result = tmp_measure, ...)
							}
							if(!is.null(tmp_measure)){
								tmp_measure <- cbind.data.frame(by_group = each_group, tmp_measure)
								compare_result %<>% rbind.data.frame(., tmp_measure)
							}
						}
						colnames(compare_result)[2] <- "Group"
					}
				}else{
					for(k in measure){
						model <- aov(reformulate(anova_set, "Value"), div_table)
						compare_result[[k]] <- summary(model)
					}
				}
			}
			self$res_diff <- compare_result
			self$cal_diff_method <- method
			message('The result is stored in object$res_diff ...')
		},
		#' @description
		#' Plotting the alpha diversity.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color pallete for groups.
		#' @param measure default Shannon; alpha diversity measurement; see names of alpha_diversity of dataset, 
		#'   e.g., Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage, PD.
		#' @param group default NULL; group name used for the plot.
		#' @param add_sig default TRUE; wheter add significance label using the result of cal_diff function, i.e. \code{object$res_diff};
		#'   This is manily designed to add post hoc test of anova or Dunn's Kruskal-Wallis Multiple Comparisons to make the label adding easy.
		#' @param add_sig_label default "Significance"; select a colname of \code{object$res_diff} for the label text, such as 'P.adj' or 'Significance'.
		#' @param add_sig_text_size default 3.88; the size of text in added label.
		#' @param use_boxplot default TRUE; TRUE: boxplot; FALSE: mean-se plot.
		#' @param boxplot_color default TRUE; TRUE: use color_values, FALSE: use "black".
		#' @param boxplot_add default "jitter"; points type, see the add parameter in \code{ggpubr::ggboxplot}.
		#' @param order_x_mean default FALSE; whether order x axis by the means of groups from large to small.
		#' @param y_start default 1.01; the y axis value from which to add the label; the default 1.01 means \code{1.01 * max(values)}.
		#' @param y_increase default 0.05; the increasing y axia space to add label; the default 0.05 means \code{0.05 * y_start}; 
		#' 	  this parameter is also used to label the letters of anova result with the fixed \code{(1 + y_increase) * y_start space}.
		#' @param xtext_angle default NULL; number (e.g. 30) used to make x axis text generate angle.
		#' @param xtext_size default 15; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @param ... parameters pass to \code{ggpubr::ggboxplot} function.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_alpha(measure = "Shannon", group = "Group")
		#' }
		plot_alpha = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			measure = "Shannon",
			group = NULL,
			add_sig = TRUE,
			add_sig_label = "Significance",
			add_sig_text_size = 3.88,
			use_boxplot = TRUE,
			boxplot_color = TRUE,
			boxplot_add = "jitter",
			order_x_mean = FALSE,
			y_start = 1.01,
			y_increase = 0.05,
			xtext_angle = NULL,
			xtext_size = 15,
			ytitle_size = 17,
			...
			){
			cal_diff_method <- self$cal_diff_method
			if(is.null(group)){
				if(is.null(self$group)){
					stop("No group provided !")
				}else{
					group <- self$group
				}
			}
			if(add_sig){
				if(is.null(self$res_diff)){
					message("No object$res_diff found. Only plot the data. To use difference test result, ",
						"please first run cal_diff function to get the significance result!")
					add_sig <- FALSE
				}
			}
			if(! measure %in% self$data_alpha$Measure){
				stop("Please provide a correct measure parameter !")
			}else{
				use_data <- self$data_alpha[self$data_alpha$Measure == measure, ]
			}
			if(order_x_mean){
				mean_orders <- names(sort(tapply(use_data$Value, use_data[, group], mean), decreasing = TRUE))
				use_data[, group] %<>% factor(., levels = mean_orders)
			}else{
				if(!is.factor(use_data[, group])){
					use_data[, group] %<>% as.factor
				}
			}
			if(use_boxplot){
				if(boxplot_color){
					color_use <- group
				}else{
					color_use <- "black"
				}
				p <- ggpubr::ggboxplot(
					use_data, 
					x = group, 
					y= "Value", 
					color = color_use, 
					palette = color_values, 
					add = boxplot_add, 
					outlier.colour = "white", 
					...
					)
			}else{
				p <- ggplot(use_data, aes_string(x = group, y = "Value")) + 
					theme_minimal() +
					stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2) +
					stat_summary(fun.y=mean, geom="point", size = rel(3)) + 
					theme(
						axis.title = element_text(face = "bold",size = rel(1.8)),
						axis.line.x = element_line(colour="black"),
						axis.line.y = element_line(colour="black"),
						axis.ticks = element_line(),
						panel.grid.major = element_line(colour="#f0f0f0"),
						panel.grid.minor = element_blank(),
						plot.margin=unit(c(10,5,5,5),"mm")
						)
			}
			if(add_sig & group == self$group){
				x_axis_order <- levels(use_data[, group])
				if(cal_diff_method == "anova"){
					if(inherits(self$res_diff, "data.frame")){
						add_letter_text <- self$res_diff[x_axis_order, measure]
						group_position <- tapply(use_data$Value, use_data[, group], function(x) {res <- max(x); ifelse(is.na(res), x, res)}) %>% 
							{. + max(.) * y_increase}
						textdata <- data.frame(
							x = x_axis_order, 
							y = group_position[x_axis_order], 
							add = add_letter_text, 
							stringsAsFactors = FALSE
							)
						p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size)
					}
				}else{
					if(!(cal_diff_method == "KW" & length(unique(use_data[, group])) > 2)){
						use_diff_data <- self$res_diff %>% .[.$Measure == measure, ]
						if(! add_sig_label %in% colnames(use_diff_data)){
							stop("Please provide a correct add_sig_label parameter! Must be a colname of object$res_diff !")
						}else{
							annotations <- use_diff_data[, add_sig_label]
						}
						if(is.numeric(annotations)){
							annotations %<>% round(., 4)
						}
						y_start <- max(use_data$Value) * y_start
						y_position <- c()
						x_min <- c()
						x_max <- c()
						for(i in seq_len(nrow(use_diff_data))){
							x_min %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\1", use_diff_data[i, "Comparison"]), x_axis_order))
							x_max %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\2", use_diff_data[i, "Comparison"]), x_axis_order))
							y_position %<>% c(., y_start * (1 + i * y_increase))
						}
						p <- p + ggsignif::geom_signif(
							annotations = annotations,
							y_position = y_position, 
							xmin = x_min, 
							xmax = x_max,
							textsize = add_sig_text_size
						)
					}
				}
			}
			p <- p + ylab(measure) + xlab("") + theme(legend.position="none")
			p <- p + theme(
					axis.text.x = element_text(colour = "black", size = xtext_size),
					axis.title.y= element_text(size = ytitle_size),
					axis.text.y = element_text(size = rel(1.1)),
					axis.title.x = element_blank()
					)
			if(!is.null(xtext_angle)){
				p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", vjust = 1, hjust = 1, size = xtext_size))
			}
			p
		},
		#' @description
		#' Print the trans_alpha object.
		print = function() {
			cat("trans_alpha object:\n")
			cat(paste("data_alpha have", ncol(self$data_alpha), "columns: ", paste0(colnames(self$data_alpha), collapse = ", "), "\n"))
			cat(paste("data_alpha$Measure: ", paste0(unique(as.character(self$data_alpha$Measure)), collapse = ", "), "\n"))
			if(!is.null(self$data_stat)) cat(paste("data_stat have", ncol(self$data_stat), "columns: ", paste0(colnames(self$data_stat), collapse = ", "), "\n"))
			invisible(self)
		}
		),
	private = list(
		kwwitt_test = function(method = NULL, input_table = NULL, group = NULL, res_list = NULL, group_by = NULL, measure = NULL, ...){
			groupvec <- as.character(input_table[, group])
			use_comp_group_num <- unique(c(2, length(unique(groupvec))))
			for(i in use_comp_group_num){
				all_name <- combn(unique(groupvec), i)
				use_method <- switch(method, KW = "Kruskal-Wallis Rank Sum Test", wilcox = "Wilcoxon Rank Sum Test", t.test = "t.test")
				for(j in 1:ncol(all_name)){
					table_compare <- input_table[groupvec %in% as.character(all_name[, j]), ]
					table_compare[, group] %<>% factor(., levels = unique(as.character(.)))
					formu <- reformulate(group, "Value")
					if(i == 2){
						if(method == "t.test"){
							res1 <- t.test(formu, data = table_compare, ...)
							max_group_select <- tapply(table_compare$Value, table_compare[, group], mean) %>% {.[which.max(.)]} %>% names
						}else{
							if(method == "wilcox"){
								res1 <- suppressWarnings(wilcox.test(formu, data = table_compare, ...))
								max_group_select <- tapply(table_compare$Value, table_compare[, group], median) %>% {.[which.max(.)]} %>% names
							}else{
								if(method == "KW" & length(use_comp_group_num) == 1){
									res1 <- kruskal.test(formu, data = table_compare, ...)
									max_group_select <- tapply(table_compare$Value, table_compare[, group], median) %>% {.[which.max(.)]} %>% names
								}else{
									next
								}
							}
						}
					}else{
						if(method == "KW"){
							res1 <- kruskal.test(formu, data = table_compare, ...)
							max_group_select <- tapply(table_compare$Value, table_compare[, group], median) %>% {.[which.max(.)]} %>% names
						}else{
							next
						}
					}
					res2 <- res1$p.value
					res_list$comnames %<>% c(., paste0(as.character(all_name[,j]), collapse = " - "))
					res_list$p_value %<>% c(., res2)
					res_list$measure_use %<>% c(., measure)
					res_list$test_method %<>% c(., use_method)
					res_list$max_group %<>% c(., max_group_select)
					if(is.null(group_by)){
						res_list$group_by %<>% c(., NA)
					}else{
						res_list$group_by %<>% c(., group_by)
					}
				}
			}
			res_list
		},
		kdunn_test = function(input_table = NULL, group = NULL, measure = NULL, group_by = NULL, ...){
			use_method <- "Dunn's Kruskal-Wallis Multiple Comparisons"
			input_table[, group] %<>% factor(., levels = unique(as.character(.)))
			formu <- reformulate(group, "Value")
			dunnTest_raw <- FSA::dunnTest(formu, data = input_table, ...)
			max_group <- lapply(dunnTest_raw$res$Comparison, function(x){
				group_select <- unlist(strsplit(x, split = " - "))
				table_compare_select <- input_table[as.character(input_table[, group]) %in% group_select, ]
				tapply(table_compare_select$Value, table_compare_select[, group], median) %>% {.[which.max(.)]} %>% names
			}) %>% unlist
			if(is.null(group_by)){
				dunnTest_res <- data.frame(Measure = measure, Test_method = use_method, Group = max_group, dunnTest_raw$res)
			}else{
				dunnTest_res <- data.frame(Measure = measure, by_group = group_by, Test_method = use_method, Group = max_group, dunnTest_raw$res)
			}
			dunnTest_res
		},
		anova_test = function(input_table = NULL, group = NULL, measure = NULL, compare_result = NULL, ...){
			model <- aov(reformulate(group, "Value"), input_table)
			out <- agricolae::duncan.test(model, group, main = measure, ...)
			res2 <- out$groups[, "groups", drop = FALSE]
			res2$groups <- as.character(res2$groups)
			res2 <- data.frame(rownames(res2), res2, stringsAsFactors = FALSE, check.names = FALSE)
			colnames(res2) <- c("name", measure)
			if(is.null(compare_result)){
				compare_result <- res2
			}else{
				compare_result <- dplyr::full_join(compare_result, res2, by = c("name" = "name"))
			}
			compare_result
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
