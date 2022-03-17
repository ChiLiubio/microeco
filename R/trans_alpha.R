#' @title Create trans_alpha object for alpha diveristy statistics and plotting.
#'
#' @description
#' This class is a wrapper for a series of alpha diveristy related analysis, including the statistics and plotting based on 
#' An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and Paul et al. (2013) <doi:10.1371/journal.pone.0061217>.
#'
#' @export
trans_alpha <- R6Class(classname = "trans_alpha",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param group default NULL; the sample column used for the statistics; If provided, can return alpha_stat.
		#' @param order_x default NULL; sample_table column name or a vector containg sample names; if provided, order samples by using factor.
		#' @return alpha_data and alpha_stat stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Group")
		#' }
		initialize = function(dataset = NULL, group = NULL, order_x = NULL) {
			self$group <- group
			if(is.null(dataset)){
				message("Parameter dataset not provided. Please run the functions with your other customized data!")
				self$alpha_data <- NULL
			}else{
				if(is.null(dataset$alpha_diversity)){
					message("Alpha diversity not found. Calculate it automatically ...")
					dataset$cal_alphadiv()
				}
				alpha_data <- dataset$alpha_diversity %>% 
					cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
					.[, !grepl("^se", colnames(.))] %>%
					# to long format
					reshape2::melt(id.vars = "Sample") %>%
					`colnames<-`(c("Sample", "Measure", "Value")) %>%
					dplyr::left_join(., rownames_to_column(dataset$sample_table), by = c("Sample" = "rowname"))
				if(!is.null(order_x)){
					if(length(order_x == 1)){
						alpha_data$Sample %<>% factor(., levels = unique(dataset$sample_table[, order_x]))
					} else {
						alpha_data$Sample %<>% factor(., levels = order_x)
					}
				}
				self$alpha_data <- alpha_data
				message('The transformed diversity data is stored in object$alpha_data ...')
			}

			if(!is.null(group)){
				if(is.null(dataset)){
					stop("Parameter dataset not provided, but group is provided!")
				}
				self$alpha_stat <- microeco:::summarySE_inter(alpha_data, measurevar = "Value", groupvars = c(self$group, "Measure"))
				message('The group statistics are stored in object$alpha_stat ...')
			}else{
				self$alpha_stat <- NULL
			}
		},
		#' @description
		#' Test the difference of alpha diversity across groups.
		#'
		#' @param method default "KW"; "KW_dunn" or "anova"; KW: Kruskal-Wallis Rank Sum Test (groups > 2) or Wilcoxon Rank Sum and Signed 
		#'   Rank Tests (groups = 2); KW_dunn: Dunn's Kruskal-Wallis Multiple Comparisons, see dunnTest function in FSA package; 
		#'   anova:  Duncan's multiple range test for anova;
		#' @param measures default NULL; a vector; if null, all indexes will be calculated; see names of microtable$alpha_diversity, 
		#' 	 e.g. Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage, PD.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of p.adjust function for available options.
		#' @param anova_set default NULL; specified group set for anova, such as 'block + N*P*K', see \code{\link{aov}}.
		#' @param ... parameters passed to kruskal.test or wilcox.test function (method = "KW") or dunnTest function of FSA package (method = "KW_dunn") or
		#'   agricolae::duncan.test (method = "anova").
		#' @return res_alpha_diff in object. A data.frame generally. A list for anova when anova_set is assigned.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(method = "KW")
		#' t1$cal_diff(method = "KW_dunn")
		#' t1$cal_diff(method = "anova")
		#' }
		cal_diff = function(method = c("KW", "KW_dunn", "anova")[1], measures = NULL, p_adjust_method = "fdr", anova_set = NULL, ...){
			group <- self$group
			alpha_data <- self$alpha_data
			if(is.null(measures)){
				measures <- unique(as.character(alpha_data$Measure))
			}else{
				if(! all(measures %in% as.character(alpha_data$Measure))){
					stop("One or more measures input not in the alpha_data! Please check the input!")
				}
			}
			if(method == "KW" & length(unique(as.character(alpha_data[, group]))) > 5){
				stop("There are too many groups to do paired comparisons using KW method, please use anova!")
			}
			if(!is.null(anova_set)){
				method <- "anova"
			}
			if(method == "KW"){
				comnames <- c()
				p_value <- c()
				measure_use <- c()
				test_method <- c()
				for(k in measures){
					div_table <- alpha_data[alpha_data$Measure == k, c(group, "Value")]
					groupvec <- as.character(div_table[, group])
					use_comp_group_num <- unique(c(2, length(unique(groupvec))))
					for(i in use_comp_group_num){
						all_name <- combn(unique(groupvec), i)
						use_method <- ifelse(i == 2, "Wilcoxon Rank Sum Test", "Kruskal-Wallis Rank Sum Test")
						for(j in 1:ncol(all_name)){
							table_compare <- div_table[groupvec %in% as.character(all_name[, j]), ]
							table_compare[, group] %<>% factor(., levels = unique(as.character(.)))
							formu <- reformulate(group, "Value")
							if(i == 2){
								res1 <- wilcox.test(formu, data = table_compare, ...)
							}else{
								res1 <- kruskal.test(formu, data = table_compare, ...)
							}
							res2 <- res1$p.value
							comnames %<>% c(., paste0(as.character(all_name[,j]), collapse = " vs "))
							p_value %<>% c(., res2)
							measure_use %<>% c(., k)
							test_method %<>% c(., use_method)
						}
					}
				}
				p_value_adjust <- p.adjust(p_value, method = p_adjust_method)
				significance_label <- cut(p_value_adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
				compare_result <- data.frame(comnames, measure_use, test_method, p_value, p_value_adjust, significance_label)
				colnames(compare_result) <- c("Groups", "Measure", "Test_method", "p.value", "p.adjust", "Significance")
			}
			if(method == "KW_dunn"){
				if(length(unique(alpha_data[, group])) == 2){
					stop("There are only 2 groups. Please select other method instead of KW_dunn !")
				}
				use_method <- "Dunn's Kruskal-Wallis Multiple Comparisons"
				compare_result <- data.frame()
				for(k in measures){
					div_table <- alpha_data[alpha_data$Measure == k, c(group, "Value")]
					table_compare <- div_table
					table_compare[, group] %<>% factor(., levels = unique(as.character(.)))
					formu <- reformulate(group, "Value")
					dunnTest_raw <- FSA::dunnTest(formu, data = table_compare, ...)
					dunnTest_res <- data.frame(Measure = k, Test_method = use_method, dunnTest_raw$res)
					compare_result %<>% rbind(., dunnTest_res)
				}
			}
			if(method == "anova"){
				# library(agricolae)
				if(is.null(anova_set)){
					compare_result <- NULL
				}else{
					compare_result <- list()
				}
				
				for(i in measures){
					use_data <- alpha_data[alpha_data$Measure == i, ]
					if(is.null(anova_set)){
						model <- aov(reformulate(group, "Value"), use_data)
						out <- agricolae::duncan.test(model, group, main = i, ...)
						res2 <- out$groups[, "groups", drop = FALSE]
						res2$groups <- as.character(res2$groups)
						res2 <- data.frame(rownames(res2), res2, stringsAsFactors = FALSE, check.names = FALSE)
						colnames(res2) <- c("name", i)
						if(is.null(compare_result)){
							compare_result <- res2
						} else {
							compare_result <- dplyr::full_join(compare_result, res2, by = c("name" = "name"))
						}
					}else{
						model <- aov(reformulate(anova_set, "Value"), use_data)
						compare_result[[i]] <- summary(model)
					}
				}
				if(is.null(anova_set)){
					compare_result %<>% `row.names<-`(.[,1]) %>% .[,-1]
				}
			}
			self$res_alpha_diff <- compare_result
			message('The result is stored in object$res_alpha_diff ...')
		},
		#' @description
		#' Plotting the alpha diveristy.
		#'
		#' @param color_values colors used for presentation.
		#' @param measure default Shannon; alpha diveristy measurement; see names of alpha_diversity of dataset, 
		#'   e.g. Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage, PD.
		#' @param group default NULL; group name used for the plot.
		#' @param add_letter default FALSE; If TRUE, the letters of duncan test will be added in the plot.
		#' @param use_boxplot default TRUE; TRUE: boxplot, FALSE: mean_se plot.
		#' @param boxplot_color default TRUE; TRUE: use color_values, FALSE: use "black".
		#' @param boxplot_add default "jitter"; points type, see the add parameter in ggpubr::ggboxplot.
		#' @param order_x_mean default FALSE; whether order x axis by the means of groups from large to small.
		#' @param pair_compare default FALSE; whether perform paired comparisons.
		#' @param pair_compare_filter default ""; groups that need to be removed in the comparisons.
		#' @param pair_compare_method default wilcox.test; wilcox.test, kruskal.test, t.test or anova.
		#' @param xtext_angle default NULL; number (e.g. 30) used to make x axis text generate angle.
		#' @param xtext_size default 10, x axis text size.
		#' @param ytitle_size default 17, y axis title size.
		#' @param base_font default "sans", font in the plot.
		#' @param ... parameters pass to ggpubr::ggboxplot function.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_alpha(measure = "Shannon", group = "Group", pair_compare = TRUE)
		#' }
		plot_alpha = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			measure = "Shannon",
			group = NULL,
			add_letter = FALSE,
			use_boxplot = TRUE,
			boxplot_color = TRUE,
			boxplot_add = "jitter",
			order_x_mean = TRUE,
			pair_compare = FALSE,
			pair_compare_filter = "",
			pair_compare_method = "wilcox.test",
			xtext_angle = NULL,
			xtext_size = 10,
			ytitle_size = 17,
			base_font = "sans",
			...
			){
			if(is.null(group)){
				group <- self$group
			}
			use_data <- self$alpha_data[self$alpha_data$Measure == measure, ]
			
			if(order_x_mean){
				use_data[, group] <- factor(use_data[, group], levels = names(sort(tapply(use_data$Value, use_data[, group], mean), decreasing = TRUE)))
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
				
				if(add_letter){
					order_groups <- names(tapply(use_data$Value, use_data[, group], max))
					group_position <- tapply(use_data$Value, use_data[, group], function(x) {res <- max(x); ifelse(is.na(res), x, res)}) %>% 
						{. + max(.)/30}
					textdata <- data.frame(
						x = order_groups, 
						y = group_position[order_groups], 
						add = self$res_alpha_diff[order_groups, measure], 
						stringsAsFactors = FALSE
						)
					p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = 7)
				}
				if(pair_compare){
					# construct and filter the paired comparisons list
					comparisons_list <- unique(as.character(self$alpha_data[, group])) %>% 
						combn(., 2) %>% 
						{.[, unlist(lapply(as.data.frame(.), function(x) any(grepl(pair_compare_filter, x)))), drop = FALSE]} %>% 
						{lapply(seq_len(ncol(.)), function(x) .[, x])}
					
					p <- p + ggpubr::stat_compare_means(
							comparisons = comparisons_list,
							method = pair_compare_method, 
							tip.length=0.01,
							label = "p.signif",
							symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
				}
			}else{
				p <- ggplot(use_data, aes_string(x = group, y = "Value")) + 
					theme_minimal() +
					stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2) +
					stat_summary(fun.y=mean, geom="point", size = rel(3))

				if(add_letter){
					order_groups <- names(tapply(use_data$Value, use_data[, group], max))
					group_position <- tapply(use_data$Value, use_data[, group], function(x) {res <- mean_se(x)$ymax; ifelse(is.na(res), x, res)}) %>% {. + max(.)/50}
					textdata <- data.frame(x = order_groups, y = group_position[order_groups], add = self$res_alpha_diff[order_groups, measure], stringsAsFactors = FALSE)
					p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = 7)
				}
				p <- p + theme(
					axis.title = element_text(face = "bold",size = rel(1.8)),
					axis.line.x = element_line(colour="black"),
					axis.line.y = element_line(colour="black"),
					axis.ticks = element_line(),
					panel.grid.major = element_line(colour="#f0f0f0"),
					panel.grid.minor = element_blank(),
					plot.margin=unit(c(10,5,5,5),"mm")
				)
			}
			
			p <- p + ylab(measure) + xlab("") + theme(legend.position="none")
			p <- p + theme(
					axis.text.x = element_text(colour = "black", size = xtext_size),
					axis.title.y= element_text(size=ytitle_size),
					axis.text.y = element_text(size = rel(1.1)),
					axis.title.x = element_blank()
					)
			if(!is.null(xtext_angle)){
				p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", vjust = 1, hjust = 1, size = xtext_size))
			}
			if(!is.null(base_font)){
				p <- p + theme(text=element_text(family=base_font))
			}
			p
		},
		#' @description
		#' Print the trans_alpha object.
		print = function() {
			cat("trans_alpha class:\n")
			cat(paste("alpha_data have", ncol(self$alpha_data), "columns: ", paste0(colnames(self$alpha_data), collapse = ", "), "\n"))
			if(!is.null(self$alpha_stat)) cat(paste("alpha_stat have", ncol(self$alpha_stat), "columns: ", paste0(colnames(self$alpha_stat), collapse = ", "), "\n"))
			invisible(self)
		}
		),
	lock_objects = FALSE,
	lock_class = FALSE
)
