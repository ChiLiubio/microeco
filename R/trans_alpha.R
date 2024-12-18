#' @title Create \code{trans_alpha} object for alpha diversity statistics and visualization.
#'
#' @description
#' This class is a wrapper for a series of alpha diversity analysis, including the statistics and visualization.
#'
#' @export
trans_alpha <- R6Class(classname = "trans_alpha",
	public = list(
		#' @param dataset \code{\link{microtable}} object.
		#' @param group default NULL; a column name of \code{sample_table} in the input microtable object used for the statistics across groups.
		#' @param by_group default NULL; a column name of \code{sample_table} used to perform the differential test 
		#'   among groups (from \code{group} parameter) for each group (from \code{by_group} parameter) separately.
		#' @param by_ID default NULL; a column name of \code{sample_table} used to perform paired T test or paired Wilcoxon test for the paired data,
		#'   such as continuous sampling of individual animals or plant compartments for different plant species (ID). 
		#'   So \code{by_ID} in sample_table should be the smallest unit of sample collection without any repetition in it.
		#'   When the \code{by_ID} parameter is provided, the function can automatically perform paired test, and no more parameters is required.
		#' @param order_x default NULL; a column name of \code{sample_table} or a vector with sample names. If provided, sort samples using \code{factor}.
		#' @return \code{data_alpha} and \code{data_stat} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Group")
		#' }
		initialize = function(dataset = NULL, group = NULL, by_group = NULL, by_ID = NULL, order_x = NULL) {
			if(is.null(dataset)){
				message("Parameter dataset not provided. Please run the functions with your other customized data ...")
				self$data_alpha <- NULL
				self$data_stat <- NULL
				if(! is.null(group)){
					stop("Parameter dataset is not provided, but group is provided!")
				}
			}else{
				dataset$tidy_dataset()
				if(is.null(dataset$alpha_diversity)){
					message("The alpha_diversity in dataset not found! Calculate it automatically ...")
					dataset$cal_alphadiv()
				}
				data_alpha <- dataset$alpha_diversity %>% .[, !grepl("^se", colnames(.)), drop = FALSE]
				tmp_filter <- c()
				for(i in colnames(data_alpha)){
					if(any(is.nan(data_alpha[, i]))){
						tmp_filter %<>% c(., i)
					}
				}
				if(length(tmp_filter) > 0){
					message("NaN is found for index ", paste(tmp_filter, collapse = " "), "! Filtering ...")
					data_alpha %<>% .[, ! colnames(.) %in% tmp_filter, drop = FALSE]
				}
				data_alpha %<>% cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
					reshape2::melt(id.vars = "Sample") %>%
					`colnames<-`(c("Sample", "Measure", "Value")) %>%
					dplyr::left_join(., rownames_to_column(dataset$sample_table), by = c("Sample" = "rowname"))
				if(!is.null(by_group)){
					check_table_variable(data_alpha, by_group, "by_group", "dataset$sample_table")
					if(is.factor(dataset$sample_table[, by_group])){
						data_alpha[, by_group] %<>% factor(., levels = levels(dataset$sample_table[, by_group]))
					}
				}
				if(! is.null(group)){
					check_table_variable(data_alpha, group, "group", "dataset$sample_table")
					self$data_stat <- microeco:::summarySE_inter(data_alpha, measurevar = "Value", groupvars = c(group, by_group, "Measure"))
					message('The group statistics are stored in object$data_stat ...')
					if(is.factor(dataset$sample_table[, group])){
						data_alpha[, group] %<>% factor(., levels = levels(dataset$sample_table[, group]))
					}
				}else{
					self$data_stat <- NULL
				}
				check_table_variable(data_alpha, by_ID, "by_ID", "dataset$sample_table")
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
			
			self$group <- group
			self$by_group <- by_group
			self$by_ID <- by_ID
		},
		#' @description
		#' Differential test on alpha diversity.
		#'
		#' @param measure default NULL; character vector; If NULL, all indexes will be used; see names of \code{microtable$alpha_diversity}, 
		#' 	 e.g. \code{c("Observed", "Chao1", "Shannon")}.
		#' @param method default "KW"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'KW'}}{Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons <10.1080/00401706.1964.10490181> based on \code{dunnTest} function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum Test for all paired groups
		#'     	  When \code{by_ID} parameter is provided in creating the object of the class, paired Wilcoxon test will be performed.}
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups. 
		#'     	  When \code{by_ID} parameter is provided in creating the object of the class, paired t-test will be performed.}
		#'     \item{\strong{'anova'}}{Variance analysis. For one-way anova, the default post hoc test is Duncan's new multiple range test.
		#'     	  Please use \code{anova_post_test} parameter to change the post hoc method.
		#'     	  For multi-way anova, Please use \code{formula} parameter to specify the model and see \code{\link{aov}} for more details}
		#'     \item{\strong{'scheirerRayHare'}}{Scheirer-Ray-Hare test (nonparametric test) for a two-way factorial experiment; 
		#'     	  see \code{scheirerRayHare} function of \code{rcompanion} package}
		#'     \item{\strong{'lm'}}{Linear Model based on the \code{lm} function}
		#'     \item{\strong{'lme'}}{Linear Mixed Effect Model based on the \code{lmerTest} package}
		#'     \item{\strong{'betareg'}}{Beta Regression for Rates and Proportions based on the \code{betareg} package}
		#'     \item{\strong{'glmm'}}{Generalized linear mixed model (GLMM) based on the \code{glmmTMB} package.
		#'     	  A family function can be provided using parameter passing, such as: \code{family = glmmTMB::lognormal(link = "log")}}
		#'     \item{\strong{'glmm_beta'}}{Generalized linear mixed model (GLMM) with a family function of beta distribution. 
		#'     	  This is an extension of the GLMM model in \code{'glmm'} option. 
		#'     	  The only difference is in \code{glmm_beta} the family function is fixed with the beta distribution function, 
		#'     	  facilitating the fitting for proportional data (ranging from 0 to 1). The link function is fixed with \code{"logit"}.}
		#'   }
		#' @param p_adjust_method default "fdr" (for "KW", "wilcox", "t.test" methods) or "holm" (for "KW_dunn"); P value adjustment method; 
		#' 	  For \code{method = 'KW', 'wilcox' or 't.test'}, please see \code{method} parameter of \code{p.adjust} function for available options;
		#' 	  For \code{method = 'KW_dunn'}, please see \code{dunn.test::p.adjustment.methods} for available options.
		#' @param formula default NULL; applied to two-way or multi-factor analysis when 
		#'   method is \code{"anova"}, \code{"scheirerRayHare"}, \code{"lm"}, \code{"lme"}, \code{"betareg"} or \code{"glmm"}; 
		#'   specified set for independent variables, i.e. the latter part of a general formula, 
		#'   such as \code{'block + N*P*K'}.
		#' @param KW_dunn_letter default TRUE; For \code{method = 'KW_dunn'}, \code{TRUE} denotes significances are presented by letters;
		#'   \code{FALSE} means significances are shown by asterisk for paired comparison.
		#' @param alpha default 0.05; Significant level; used for generating significance letters when method is 'anova' or 'KW_dunn'.
		#' @param anova_post_test default "duncan.test". The post hoc test method for one-way anova. 
		#'   The default option represents the Duncan's new multiple range test.
		#'   Other available options include "LSD.test" (LSD post hoc test) and "HSD.test" (HSD post hoc test). 
		#'   All those are the function names from \code{agricolae} package.
		#' @param anova_varequal_test default FALSE; whether conduct Levene's Test for equality of variances.
		#'   Only available for one-way anova. Significant P value means the variance among groups is not equal.
		#' @param return_model default FALSE; whether return the original "lm", "lmer" or "glmm" model list in the object.
		#' @param ... parameters passed to \code{kruskal.test} (when \code{method = "KW"}) or \code{wilcox.test} function (when \code{method = "wilcox"}) or 
		#'   \code{dunnTest} function of \code{FSA} package (when \code{method = "KW_dunn"}) or 
		#'   \code{agricolae::duncan.test}/\code{agricolae::LSD.test}/\code{agricolae::HSD.test} (when \code{method = "anova"}, one-way anova) or 
		#'   \code{rcompanion::scheirerRayHare} (when \code{method = "scheirerRayHare"}) or 
		#'   \code{stats::lm} (when \code{method = "lm"}) or 
		#'   \code{lmerTest::lmer} (when \code{method = "lme"}) or 
		#'   \code{betareg::betareg} (when \code{method = "betareg"}) or 
		#'   \code{glmmTMB::glmmTMB} (when \code{method = "glmm"}).
		#' @return \code{res_diff}, stored in object with the format \code{data.frame}.\cr
		#'     When method is "betareg", "lm", "lme" or "glmm", 
		#'     "Estimate" and "Std.Error" columns represent the fitted coefficient and its standard error, respectively.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(method = "KW")
		#' t1$cal_diff(method = "anova")
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Type", by_group = "Group")
		#' t1$cal_diff(method = "anova")
		#' }
		cal_diff = function(
			measure = NULL,
			method = c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lm", "lme", "betareg", "glmm", "glmm_beta")[1], 
			formula = NULL,
			p_adjust_method = "fdr", 
			KW_dunn_letter = TRUE,
			alpha = 0.05,
			anova_post_test = "duncan.test",
			anova_varequal_test = FALSE,
			return_model = FALSE,
			...
			){
			method <- match.arg(method, c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lm", "lme", "betareg", "glmm", "glmm_beta"))
			group <- self$group
			
			if(method %in% c("scheirerRayHare", "lm", "lme", "betareg", "glmm", "glmm_beta") & is.null(formula)){
				if(is.null(formula)){
					stop("The formula parameter is NULL! It is necessary for the method: ", method, " !")
				}
			}
			if(!method %in% c("anova", "scheirerRayHare", "lm", "lme", "betareg", "glmm", "glmm_beta")){
				if(is.null(group)){
					stop("For the method: ", method, " , group is necessary! Please recreate the object and set the group parameter!")
				}
			}
			if(method == "anova"){
				if(is.null(group) & is.null(formula)){
					stop("Please provide either formula parameter or group parameter (in the object creation) for ANOVA method!")
				}
			}
			# 'by_group' for test inside each by_group instead of all groups in 'group'
			by_group <- self$by_group
			data_alpha <- self$data_alpha
			by_ID <- self$by_ID
			
			if(! method %in% c("lm", "lme", "betareg", "glmm", "glmm_beta")){
				data_alpha %<>% dropallfactors
			}
			if(is.null(measure)){
				measure <- unique(as.character(data_alpha$Measure))
			}else{
				if(! all(measure %in% as.character(data_alpha$Measure))){
					stop("One or more measures input not in the object$data_alpha! Please check the input!")
				}
			}
			if(is.null(by_group)){
				if(method %in% c("wilcox", "t.test") & length(unique(as.character(data_alpha[, group]))) > 10){
					stop("There are too many groups to do paired comparisons! please use method = 'KW' or 'KW_dunn' or 'anova' !")
				}
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
						if(is.null(by_ID)){
							div_table <- data_alpha[data_alpha$Measure == k, c(group, "Value")]
						}else{
							div_table <- data_alpha[data_alpha$Measure == k, c(group, by_ID, "Value")]
						}
						res_list <- private$kwwitt_test(method = method, input_table = div_table, group = group, res_list = res_list, measure = k, by_ID = by_ID, ...)
					}else{
						for(each_group in unique_bygroups){
							if(is.null(by_ID)){
								div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, c(by_group, group, "Value")]
							}else{
								div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, c(by_group, group, by_ID, "Value")]
							}
							if(length(unique(as.character(div_table[, group]))) < 2){
								next
							}
							res_list <- private$kwwitt_test(method = method, input_table = div_table, group = group, res_list = res_list, 
								group_by = each_group, measure = k, by_ID = by_ID, ...)
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
					colnames(compare_result) <- c("Comparison", "Measure", "Method", "Group", "P.unadj", "P.adj")
				}else{
					compare_result <- data.frame(res_list$comnames, res_list$group_by, res_list$measure_use, res_list$test_method, 
						res_list$max_group, res_list$p_value, p_value_adjust)
					colnames(compare_result) <- c("Comparison", "by_group", "Measure", "Method", "Group", "P.unadj", "P.adj")
				}
			}
			if(method == "KW_dunn"){
				if(is.null(by_group)){
					if(length(unique(data_alpha[, group])) == 2){
						stop("There are only 2 groups. Please select other method instead of KW_dunn !")
					}
				}
				if(is.null(p_adjust_method)){
					p_adjust_method <- "none"
				}else{
					if(p_adjust_method == "fdr"){
						p_adjust_method <- "holm"
					}else{
						if(!p_adjust_method %in% c('none', 'bonferroni', 'sidak', 'holm', 'hs', 'hochberg', 'bh', 'by')){
							stop("For KW_dunn method, p_adjust_method must be one of 'none', 'bonferroni', 'sidak', 'holm', 'hs', 'hochberg', 'bh' and 'by'!")
						}
					}
				}
				if(p_adjust_method != "none"){
					message("P value adjustment method: ", p_adjust_method, " ...")
				}
				compare_result <- data.frame()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, c(group, "Value")]
						if(private$check_skip_comp(div_table, k, by_group)){
							next
						}
						tmp_res <- private$kdunn_test(input_table = div_table, group = group, measure = k, KW_dunn_letter = KW_dunn_letter, 
							p_adjust_method = p_adjust_method, alpha = alpha, ...)
						compare_result %<>% rbind(., tmp_res)
					}else{
						for(each_group in unique_bygroups){
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, c(by_group, group, "Value")]
							if(length(unique(as.character(div_table[, group]))) < 3){
								message("Skip the by_group: ", each_group, " as groups number < 3!")
								next
							}
							if(private$check_skip_comp(div_table, k, by_group, each_group)){
								next
							}
							tmp_res <- private$kdunn_test(input_table = div_table, group = group, measure = k, KW_dunn_letter = KW_dunn_letter, 
								p_adjust_method = p_adjust_method, alpha = alpha, ...)
							tmp_res <- cbind.data.frame(by_group = each_group, tmp_res)
							compare_result %<>% rbind(., tmp_res)
						}
					}
				}
			}
			if(method == "anova" & is.null(formula)){
				if(!require("agricolae")){
					stop("Please first install the agricolae package!")
				}
				anova_post_test <- match.arg(anova_post_test, c("duncan.test", "LSD.test", "HSD.test"))
				message("Perform post hoc test with the method: ", anova_post_test, " ...")
				compare_result <- data.frame()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, ]
						tmp_res <- private$anova_test(input_table = div_table, group = group, measure = k, post_test = anova_post_test, 
							alpha = alpha, varequal_test = anova_varequal_test, ...)
						compare_result %<>% rbind(., tmp_res)
					}else{
						for(each_group in unique_bygroups){
							test <- data_alpha[all_bygroups == each_group, ]
							if(length(unique(as.character(test[, group]))) < 2){
								message("Skip the by_group: ", each_group, " as groups number < 2!")
								next
							}
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, ]
							check_res <- tryCatch(tmp_res <- private$anova_test(input_table = div_table, group = group, measure = k, post_test = anova_post_test, 
								alpha = alpha, varequal_test = anova_varequal_test, ...),
								error = function(e) {skip_to_next <- TRUE})
							if(rlang::is_true(check_res)){
								next
							}
							tmp_res <- cbind.data.frame(by_group = each_group, tmp_res)
							compare_result %<>% rbind.data.frame(., tmp_res)
						}
					}
				}
			}
			if(method %in% c("anova", "scheirerRayHare", "betareg") & !is.null(formula)){
				# to make multi-factor distinct from one-way anova
				compare_result <- data.frame()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, ]
						tmp_res <- private$formula_anova_sche_betareg_test(method = method, input_table = div_table, formula = formula, measure = k, ...)
						if(is.null(tmp_res)){
							next
						}else{
							compare_result %<>% rbind(., tmp_res)
						}
					}else{
						for(each_group in unique_bygroups){
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, ]
							tmp_res <- private$formula_anova_sche_betareg_test(method = method, input_table = div_table, formula = formula, measure = k, ...)
							if(is.null(tmp_res)){
								next
							}else{
								tmp_res <- cbind.data.frame(by_group = each_group, tmp_res)
								compare_result %<>% rbind(., tmp_res)
							}
						}
					}
				}
				method <- paste0(method, "-formula")
			}
			if(method %in% c("lm", "lme", "glmm", "glmm_beta")){
				compare_result <- data.frame()
				if(return_model){
					res_model <- list()
				}
				for(k in measure){
					div_table <- data_alpha[data_alpha$Measure == k, ]
					if(method == "lm"){
						tmp <- stats::lm(reformulate(formula, "Value"), data = div_table, ...)
						if(return_model){
							res_model[[k]] <- tmp
						}
						tmp_summary <- summary(tmp)
						tmp_coefficients <- as.data.frame(tmp_summary$coefficients, check.names = FALSE)
						tmp_model_p <- anova(tmp)
						R2data <- private$R2_extract(tmp, nrow(tmp_model_p) + nrow(tmp_coefficients))

						tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), 
							Measure = k, 
							Factors = c("Model", rownames(tmp_model_p), rownames(tmp_coefficients)), 
							R2data, 
							Estimate  = c(NA, rep(NA, nrow(tmp_model_p)), tmp_coefficients$Estimate), 
							Std.Error = c(NA, rep(NA, nrow(tmp_model_p)), tmp_coefficients$`Std. Error`), 
							P.unadj = c(NA, tmp_model_p$`Pr(>F)`, tmp_coefficients$`Pr(>|t|)`)
						)
					}else{
						if(method == "lme"){
							tmp <- lmerTest::lmer(reformulate(formula, "Value"), data = div_table, ...)
							if(return_model){
								res_model[[k]] <- tmp
							}
							tmp_summary <- summary(tmp)
							tmp_coefficients <- as.data.frame(tmp_summary$coefficients, check.names = FALSE)
							tmp_model_p <- anova(tmp)
							tmp_random_p <- lmerTest::ranova(tmp)
							R2data <- private$R2_extract(tmp, nrow(tmp_model_p) + nrow(tmp_random_p) + nrow(tmp_coefficients))
							
							tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), 
								Measure = k, 
								Factors = c("Model", rownames(tmp_model_p), rownames(tmp_random_p), rownames(tmp_coefficients)), 
								R2data, 
								Estimate  = c(NA, rep(NA, nrow(tmp_model_p) + nrow(tmp_random_p)), tmp_coefficients$Estimate), 
								Std.Error = c(NA, rep(NA, nrow(tmp_model_p) + nrow(tmp_random_p)), tmp_coefficients$`Std. Error`), 
								P.unadj   = c(NA, tmp_model_p$`Pr(>F)`, tmp_random_p$`Pr(>Chisq)`, tmp_coefficients$`Pr(>|t|)`)
							)
						}else{
							if(method == "glmm"){
								tmp <- glmmTMB::glmmTMB(reformulate(formula, "Value"), data = div_table, ...)
							}else{
								tmp <- glmmTMB::glmmTMB(reformulate(formula, "Value"), data = div_table, family = glmmTMB::beta_family(link = "logit"), ...)
							}
							if(return_model){
								res_model[[k]] <- tmp
							}
							tmp_summary <- summary(tmp)
							tmp_coefficients <- as.data.frame(tmp_summary$coefficients$cond, check.names = FALSE)
							tmp_model_p <- car::Anova(tmp)
							R2data <- private$R2_extract(tmp, nrow(tmp_model_p) + nrow(tmp_coefficients))

							tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), 
								Measure = k, 
								Factors = c("Model", rownames(tmp_model_p), rownames(tmp_coefficients)), 
								R2data, 
								Estimate = c(NA, rep(NA, nrow(tmp_model_p)), tmp_coefficients$Estimate), 
								Std.Error = c(NA, rep(NA, nrow(tmp_model_p)), tmp_coefficients$`Std. Error`), 
								P.unadj = c(NA, tmp_model_p$`Pr(>Chisq)`, tmp_coefficients$`Pr(>|z|)`)
							)
						}
					}
					compare_result %<>% rbind(., tmp_res)
				}
				if(return_model){
					self$res_model <- res_model
					message("The original ", method, " models list is stored in object$res_model ...")
				}
				method <- paste0(method, "-formula")
			}
			if("P.adj" %in% colnames(compare_result)){
				compare_result$Significance <- generate_p_siglabel(compare_result$P.adj, nonsig = "ns") %>% as.character
			}else{
				if("P.unadj" %in% colnames(compare_result)){
					compare_result$Significance <- generate_p_siglabel(compare_result$P.unadj, nonsig = "ns") %>% as.character
				}
			}
			self$res_diff <- compare_result
			self$cal_diff_method <- method
			self$measure <- measure
			message('The result is stored in object$res_diff ...')
			invisible(self)
		},
		#' @description
		#' Plot the alpha diversity. 
		#'   Box plot (and others for visualizing data in groups of single factor) is used for the visualization of alpha diversity when the \code{group} is found in the object.
		#'   When the formula is found in the \code{res_diff} table in the object, 
		#'   heatmap is employed automatically to show the significances of differential test for multiple indexes, 
		#' 	 and errorbar (coefficient and standard errors) can be used for single index.
		#'
		#' @param plot_type default "ggboxplot"; plot type; available options include "ggboxplot", "ggdotplot", "ggviolin", 
		#'   "ggstripchart", "ggerrorplot", "errorbar" and "barerrorbar".
		#'   The options starting with "gg" are function names coming from \code{ggpubr} package.
		#'   All those methods with \code{ggpubr} package use the \code{data_alpha} table in the object. 
		#'   "errorbar" represents Mean±SD or Mean±SE plot based on \code{ggplot2} package by invoking the \code{data_stat} table in the object.
		#'   "barerrorbar" denotes "bar plot + error bar". It is similar with "errorbar" and has a bar plot.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color pallete for groups.
		#' @param measure default "Shannon"; one alpha diversity index in the object.
		#' @param group default NULL; group name used for the plot.
		#' @param add default NULL; add another plot element; passed to the \code{add} parameter of the function (e.g., \code{ggboxplot}) from \code{ggpubr} package 
		#'   when \code{plot_type} starts with "gg" (functions coming from ggpubr package).
		#' @param add_sig default TRUE; whether add significance label using the result of \code{cal_diff} function, i.e. \code{object$res_diff};
		#'   This is manily designed to add post hoc test of anova or other significances to make the label mapping easy.
		#' @param add_sig_label default "Significance"; select a colname of \code{object$res_diff} for the label text when 'Letter' is not in the table, 
		#'   such as 'P.adj' or 'Significance'.
		#' @param add_sig_text_size default 3.88; the size of text in added label.
		#' @param add_sig_label_num_dec default 4; reserved decimal places when the parameter \code{add_sig_label} use numeric column, like 'P.adj'.
		#' @param order_x_mean default FALSE; whether order x axis by the means of groups from large to small.
		#' @param y_start default 0.1; the y axis value from which to add the significance asterisk label; 
		#' 	  the default 0.1 means \code{max(values) + 0.1 * (max(values) - min(values))}.
		#' @param y_increase default 0.05; the increasing y axia space to add the label (asterisk or letter); the default 0.05 means \code{0.05 * (max(values) - min(values))}; 
		#' 	  this parameter is also used to label the letters of anova result with the fixed space.
		#' @param xtext_angle default 30; number (e.g. 30). Angle of text in x axis.
		#' @param xtext_size default 13; x axis text size. NULL means the default size in ggplot2.
		#' @param ytitle_size default 17; y axis title size.
		#' @param bar_width default 0.9; the bar width when \code{plot_type = "barerrorbar"}.
		#' @param bar_alpha default 0.8; the alpha of bar color when \code{plot_type = "barerrorbar"}.
		#' @param dodge_width default 0.9; the dodge width used in \code{position_dodge} function of ggplot2 package when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param plot_SE default TRUE; TRUE: the errorbar is \eqn{mean±se}; FALSE: the errorbar is \eqn{mean±sd}. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param errorbar_size default 1; errorbar size. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param errorbar_width default 0.2; errorbar width. Available when \code{plot_type} is "errorbar" or "barerrorbar" and \code{by_group} is NULL.
		#' @param errorbar_addpoint default TRUE; whether add point for mean. Available when \code{plot_type} is "errorbar" or "barerrorbar" and \code{by_group} is NULL.
		#' @param errorbar_color_black default FALSE; whether use black for the color of errorbar when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param point_size default 3; point size for taxa. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param point_alpha default 0.8; point transparency. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param add_line default FALSE; whether add line. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param line_size default 0.8; line size when \code{add_line = TRUE}. Available when \code{plot_type} is "errorbar" or "barerrorbar".
		#' @param line_type default 2; an integer; line type when \code{add_line = TRUE}. The available case is same with \code{line_size}.
		#' @param line_color default "grey50"; line color when \code{add_line = TRUE}. Available when \code{by_group} is NULL. Other available case is same with \code{line_size}.
		#' @param line_alpha default 0.5; line transparency when \code{add_line = TRUE}. The available case is same with \code{line_size}.
		#' @param heatmap_cell default "P.unadj"; the column of \code{res_diff} table for the cell of heatmap when formula with multiple factors is found in the method.
		#' @param heatmap_sig default "Significance"; the column of \code{res_diff} for the significance label of heatmap.
		#' @param heatmap_x default "Factors"; the column of \code{res_diff} for the x axis of heatmap.
		#' @param heatmap_y default "Taxa"; the column of \code{res_diff} for the y axis of heatmap.
		#' @param heatmap_lab_fill default "P value"; legend title of heatmap.
		#' @param coefplot_sig_pos default 2; Significance label position in the coefficient point and errorbar plot. 
		#' 	  The formula is \code{Estimate + coefplot_sig_pos * Std.Error}.
		#' 	  This plot is used when there is only one measure found in the table, 
		#' 	  and 'Estimate' and 'Std.Error' are both in the column names (such as for \code{lm} and \code{lme methods}).
		#' 	  The x axis is 'Estimate', and y axis denotes 'Factors'.
		#' 	  When coefplot_sig_pos is a negative value, the label is in the left of the errorbar.
		#' 	  Errorbar size and width in the coefficient point plot can be adjusted with the parameters \code{errorbar_size} and \code{errorbar_width}. 
		#' 	  Point size and alpha can be adjusted with parameters \code{point_size} and \code{point_alpha}. 
		#' 	  The significance label size can be adjusted with parameter \code{add_sig_text_size}.
		#' 	  Furthermore, the vertical line around 0 can be adjusted with parameters \code{line_size}, \code{line_type}, \code{line_color} and \code{line_alpha}. 
		#' @param ... parameters passing to \code{ggpubr::ggboxplot} function (or other functions shown by \code{plot_type} parameter when it starts with "gg") or 
		#' 	  \code{plot_cor} function in \code{\link{trans_env}} class for the heatmap of multiple factors when formula is found in the \code{res_diff} of the object.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Group")
		#' t1$cal_diff(method = "wilcox")
		#' t1$plot_alpha(measure = "Shannon", add_sig = TRUE)
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Type", by_group = "Group")
		#' t1$cal_diff(method = "wilcox")
		#' t1$plot_alpha(measure = "Shannon", add_sig = TRUE)
		#' }
		plot_alpha = function(
			plot_type = "ggboxplot",
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			measure = "Shannon",
			group = NULL,
			add = NULL,
			add_sig = TRUE,
			add_sig_label = "Significance",
			add_sig_text_size = 3.88,
			add_sig_label_num_dec = 4,
			order_x_mean = FALSE,
			y_start = 0.1,
			y_increase = 0.05,
			xtext_angle = 30,
			xtext_size = 13,
			ytitle_size = 17,
			bar_width = 0.9,
			bar_alpha = 0.8,
			dodge_width = 0.9,
			plot_SE = TRUE,
			errorbar_size = 1,
			errorbar_width = 0.2,
			errorbar_addpoint = TRUE,
			errorbar_color_black = FALSE,
			point_size = 3,
			point_alpha = 0.8,
			add_line = FALSE,
			line_size = 0.8, 
			line_type = 2,
			line_color = "grey50",
			line_alpha = 0.5, 
			heatmap_cell = "P.unadj",
			heatmap_sig = "Significance",
			heatmap_x = "Factors",
			heatmap_y = "Measure",
			heatmap_lab_fill = "P value",
			coefplot_sig_pos = 2,
			...
			){
			# first determine visualization way
			if(is.null(self$res_diff)){
				use_heatmap <- FALSE
			}else{
				if(any(grepl("formula", self$res_diff[, "Method"]))){
					use_heatmap <- TRUE
				}else{
					use_heatmap <- FALSE
				}
			}
			if(use_heatmap){
				# For only one measure case, use errorplot
				if(all(c("Estimate", "Std.Error") %in% colnames(self$res_diff))){
					if(length(unique(self$res_diff$Measure)) == 1){
						message("For one measure, employ coefficient point and errorbar instead of heatmap ...")
						use_single_errorplot <- TRUE
					}else{
						use_single_errorplot <- FALSE
					}
				}else{
					use_single_errorplot <- FALSE
				}
				if(use_single_errorplot){
					tmp_data <- self$res_diff
					tmp_data %<>% .[!is.na(.$Estimate), ]
					
					tmp_data$sig_pos <- tmp_data$Estimate + coefplot_sig_pos * tmp_data$Std.Error
					use_color_values <- expand_colors(color_values, length(unique(tmp_data$Factors)))
					
					p <- ggplot(tmp_data, aes(x = Estimate, y = Factors, color = Factors)) +
						theme_bw() +
						geom_point(size = point_size, alpha = point_alpha) + 
						geom_errorbar(aes(xmin = Estimate - Std.Error, xmax = Estimate + Std.Error), width = errorbar_width, size = errorbar_size) +
						scale_color_manual(values = use_color_values) + 
						geom_text(aes(x = sig_pos, y = Factors, label = Significance), data = tmp_data, inherit.aes = FALSE, size = add_sig_text_size) +
						geom_vline(xintercept = 0, linetype = line_type, size = line_size, color = line_color, alpha = line_alpha) + 
						theme(panel.grid = element_blank(), legend.position = "none") +
						ylab("") +
						xlab("Coefficient")
					
				}else{
					tmp <- self$res_diff
					if("by_group" %in% colnames(tmp)){
						if(!is.factor(tmp[, "by_group"])){
							if(is.factor(self$data_alpha[, self$by_group])){
								tmp[, "by_group"] %<>% factor(., levels = levels(self$data_alpha[, self$by_group]))
							}else{
								tmp[, "by_group"] %<>% as.factor
							}
						}
					}
					tmp_trans_env <- convert_diff2transenv(tmp, heatmap_x, heatmap_y, heatmap_cell, heatmap_sig, heatmap_lab_fill)
					p <- tmp_trans_env$plot_cor(keep_full_name = TRUE, keep_prefix = TRUE, xtext_angle = xtext_angle, xtext_size = xtext_size, ...)
				}
			}else{
				
				cal_diff_method <- self$cal_diff_method
				if(is.null(group)){
					if(is.null(self$group)){
						stop("No group provided !")
					}else{
						group <- self$group
					}
				}
				by_group <- self$by_group
				if(add_sig){
					if(is.null(self$res_diff)){
						message("No res_diff found! Only plot the data ...")
						add_sig <- FALSE
					}
					if(is.null(self$group)){
						add_sig <- FALSE
					}else{
						if(group != self$group){
							add_sig <- FALSE
						}
					}
				}
				if(is.null(measure)){
					if(!is.null(self$measure)){
						message("Use the measure in the object ...")
						measure <- self$measure[1]
					}else{
						message("Input measure is NULL! Select the first one ...")
						measure <- unique(self$data_alpha$Measure)[1]
					}
				}else{
					if(length(measure) > 1){
						stop("Input measure must be one!")
					}else{
						if(! measure %in% self$data_alpha$Measure){
							stop("Please provide a correct measure parameter !")
						}	
					}
				}
				
				use_data <- self$data_alpha[self$data_alpha$Measure == measure, ]
				if(plot_type %in% c("errorbar", "barerrorbar")){
					use_data_plot <- self$data_stat[self$data_stat$Measure == measure, ]
				}
				
				if(order_x_mean){
					if(is.null(by_group)){
						mean_orders <- names(sort(tapply(use_data$Value, use_data[, group], mean), decreasing = TRUE))
						use_data[, group] %<>% factor(., levels = mean_orders)
						if(plot_type %in% c("errorbar", "barerrorbar")){
							use_data_plot[, group] %<>% factor(., levels = mean_orders)
						}
					}else{
						# first sort by_group
						mean_orders_bygroup <- names(sort(tapply(use_data$Value, use_data[, by_group], mean), decreasing = TRUE))
						use_data[, by_group] %<>% factor(., levels = mean_orders_bygroup)
					}
				}
				if(!is.factor(use_data[, group])){
					use_data[, group] %<>% as.factor
					if(plot_type %in% c("errorbar", "barerrorbar")){
						use_data_plot[, group] %<>% as.factor
					}
				}
				if(!is.null(by_group)){
					if(!is.factor(use_data[, by_group])){
						use_data[, by_group] %<>% as.factor
					}
				}

				color_values <- expand_colors(color_values, length(unique(use_data[, group])))
				
				if(! plot_type %in% c("ggboxplot", "ggdotplot", "ggviolin", "ggstripchart", "ggerrorplot", "errorbar", "barerrorbar")){
					stop("Unknown plot_type: ", plot_type, "!")
				}
				if(is.null(add)){
					add <- "none"
				}else{
					require(ggpubr)
				}
				if(plot_type %in% c("ggboxplot", "ggdotplot", "ggviolin", "ggstripchart", "ggerrorplot")){
					use_ggpubr_function <- eval(parse(text = paste0("ggpubr::", plot_type)))
					if(is.null(by_group)){
						p <- use_ggpubr_function(use_data, x = group, y = "Value", color = group, palette = color_values, add = add, ...)
					}else{
						p <- use_ggpubr_function(use_data, x = by_group, y = "Value", color = group, palette = color_values, add = add, ...)
					}
				}
				if(plot_type %in% c("errorbar", "barerrorbar")){
					colnames(use_data_plot)[colnames(use_data_plot) == "Mean"] <- "Value"
					if(plot_SE){
						colnames(use_data_plot)[colnames(use_data_plot) == "SE"] <- "errvalue"
					}else{
						colnames(use_data_plot)[colnames(use_data_plot) == "SD"] <- "errvalue"
					}
					if(is.null(by_group)){
						if(plot_type == "barerrorbar"){
							p <- ggplot(use_data_plot, aes(x = .data[[group]], y = .data[["Value"]], color = .data[[group]], fill = .data[[group]], group = 1))
							p <- p + geom_bar(stat = "identity", position = position_dodge(width = dodge_width), width = bar_width, alpha = bar_alpha)
						}else{
							p <- ggplot(use_data_plot, aes(x = .data[[group]], y = .data[["Value"]], color = .data[[group]], group = 1))
						}
						if(errorbar_color_black){
							p <- p + geom_errorbar(aes(ymin = Value - errvalue, ymax = Value + errvalue), linewidth = errorbar_size, width = errorbar_width, color = "black")
						}else{
							p <- p + geom_errorbar(aes(ymin = Value - errvalue, ymax = Value + errvalue), linewidth = errorbar_size, width = errorbar_width)
						}
						
						if(errorbar_addpoint){
							p <- p + geom_point(size = point_size, alpha = point_alpha)
						}
					}else{
						if(plot_type == "barerrorbar"){
							p <- ggplot(use_data_plot, aes(x = .data[[by_group]], y = .data[["Value"]], color = .data[[group]], fill = .data[[group]], group = .data[[group]]))
							p <- p + geom_bar(stat = "identity", position = position_dodge(width = dodge_width), width = bar_width, alpha = bar_alpha)
						}else{
							p <- ggplot(use_data_plot, aes(x = .data[[by_group]], y = .data[["Value"]], color = .data[[group]], group = .data[[group]]))
						}
						if(errorbar_color_black){
							p <- p + geom_errorbar(aes(ymin = Value - errvalue, ymax = Value + errvalue), position = position_dodge(width = dodge_width), 
								linewidth = errorbar_size, width = errorbar_width, color = "black")
						}else{
							p <- p + geom_errorbar(aes(ymin = Value - errvalue, ymax = Value + errvalue), position = position_dodge(width = dodge_width), 
								linewidth = errorbar_size, width = errorbar_width)						
						}
						if(errorbar_addpoint){
							p <- p + geom_point(size = point_size, alpha = point_alpha, position = position_dodge(width = dodge_width))
						}
					}
					if(add_line){
						if(is.null(by_group)){
							p <- p + geom_line(colour = line_color, linewidth = line_size, alpha = line_alpha, 
								linetype = line_type)
						}else{
							p <- p + geom_line(linewidth = line_size, alpha = line_alpha, 
								linetype = line_type, position = position_dodge(width = dodge_width))
						}
					}
					p <- p + scale_color_manual(values = color_values) + theme_bw()
					if(plot_type == "barerrorbar"){
						p <- p + scale_fill_manual(values = color_values)
					}
				}
				
				if(add_sig){
					if(!is.null(by_group)){
						# identify the position for the label when by_group is not NULL
						if(! plot_type %in% c("errorbar", "barerrorbar")){
							bygroup_width <- 0.8
						}else{
							bygroup_width <- dodge_width
						}
					}
					diff_res <- self$res_diff
					if(cal_diff_method %in% c("KW", "KW_dunn", "wilcox", "t.test", "anova")){
						if("Letter" %in% colnames(diff_res)){
							if(is.null(by_group)){
								tmp <- diff_res[diff_res$Measure == measure, ]
								rownames(tmp) <- tmp$Group
								x_axis_order <- levels(use_data[, group])
								add_letter_text <- tmp[x_axis_order, "Letter"]
								group_position <- tapply(use_data$Value, use_data[, group], function(x) {res <- max(x); ifelse(is.na(res), x, res)}) %>% 
									{. + y_increase * (max(use_data$Value) - min(use_data$Value))}
								textdata <- data.frame(
									x = x_axis_order, 
									y = group_position[x_axis_order], 
									add = add_letter_text, 
									stringsAsFactors = FALSE
									)
							}else{
								# for each group
								x_mid <- c()
								annotations <- c()
								y_position <- c()
								all_by_groups <- levels(use_data[, by_group])
								all_groups <- levels(use_data[, group])
								tmp_use_data <- dropallfactors(use_data)
								y_range_use <- max(tmp_use_data$Value) - min(tmp_use_data$Value)
								
								for(j in all_by_groups){
									select_tmp_data <- tmp_use_data %>% .[.[, by_group] == j, ]
									x_axis_order <- all_groups %>% .[. %in% select_tmp_data[, group]]
									# identify the start x position for each by_group
									start_bar_mid <- 1 - (bygroup_width/2 - bygroup_width/(length(x_axis_order) * 2))
									increase_bar_mid <- bygroup_width/length(x_axis_order)
									for(i in x_axis_order){
										# first determine the bar range
										mid_num <- match(j, all_by_groups) - 1
										check_anno <- diff_res[diff_res$by_group == j & diff_res$Group == i & diff_res$Measure == measure, "Letter"]
										if(length(check_anno) == 0){
											annotations %<>% c(., "")
										}else{
											annotations %<>% c(., check_anno)
										}
										x_mid %<>% c(., mid_num + (start_bar_mid + (match(i, x_axis_order) - 1) * increase_bar_mid))
									}
									tmp_position <- tapply(select_tmp_data$Value, select_tmp_data[, group], function(x) {res <- max(x); ifelse(is.na(res), x, res)}) %>% 
										{. + y_increase * y_range_use} %>% 
										.[x_axis_order]
									y_position %<>% c(., tmp_position)
								}
								textdata <- data.frame(
									x = x_mid, 
									y = y_position, 
									add = annotations, 
									stringsAsFactors = FALSE
								)
							}
							p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size, inherit.aes = FALSE)
						}else{
							use_diff_data <- diff_res %>% .[.$Measure == measure, ]
							# check group numbers
							if(cal_diff_method == "KW"){
								comp_group_num <- sapply(use_diff_data$Comparison, function(x){length(unlist(gregexpr(" - ", x)))})
								paired_add_sig <- ifelse(any(comp_group_num > 1), FALSE, TRUE)
							}else{
								paired_add_sig <- TRUE
							}
							if(paired_add_sig){
								if(! add_sig_label %in% colnames(use_diff_data)){
									stop("Please provide a correct add_sig_label parameter! Must be a colname of object$res_diff !")
								}
								if(is.numeric(use_diff_data[, add_sig_label])){
									use_diff_data[, add_sig_label] %<>% round(., add_sig_label_num_dec) %>% as.character
								}
								annotations <- c()
								y_position <- c()
								x_min <- c()
								x_max <- c()
								if(is.null(by_group)){
									x_axis_order <- levels(use_data[, group])
									y_range <- max(use_data$Value) - min(use_data$Value)
									y_start <- max(use_data$Value) + y_start * y_range
									for(i in seq_len(nrow(use_diff_data))){
										annotations %<>% c(., as.character(use_diff_data[i, add_sig_label]))
										x_min %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\1", use_diff_data[i, "Comparison"]), x_axis_order))
										x_max %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\2", use_diff_data[i, "Comparison"]), x_axis_order))
										y_position %<>% c(., y_start + i * y_increase * y_range)
									}
								}else{
									all_by_groups <- levels(use_data[, by_group])
									all_groups <- levels(use_data[, group])
									tmp_use_data <- dropallfactors(use_data)
									
									diff_by_groups <- use_diff_data$by_group %>% unique
									# same y_range_use for all elements in by_group instead of each one
									y_range_use <- max(tmp_use_data$Value) - min(tmp_use_data$Value)

									for(j in diff_by_groups){
										select_tmp_data <- tmp_use_data %>% .[.[, by_group] == j, ]
										# y axix starting position
										y_start_use <- max(select_tmp_data$Value) + y_start * y_range_use
										x_axis_order <- all_groups %>% .[. %in% select_tmp_data[, group]]
										# identify the start x position for each by_group
										start_bar_mid <- 1 - (bygroup_width/2 - bygroup_width/(length(x_axis_order) * 2))
										increase_bar_mid <- bygroup_width/length(x_axis_order)
										# loop for each group of use_diff_data
										tmp_diff_res <- use_diff_data[use_diff_data$by_group == j, ]
										for(i in seq_len(nrow(tmp_diff_res))){
											annotations %<>% c(., as.character(tmp_diff_res[i, add_sig_label]))
											# first determine the bar range
											mid_num <- match(j, all_by_groups) - 1
											x_min %<>% c(., mid_num + 
												(start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\1", tmp_diff_res[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
											x_max %<>% c(., mid_num + 
												(start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\2", tmp_diff_res[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
											y_position %<>% c(., y_start_use + i * y_increase * y_range_use)
										}
									}
								}
								p <- p + ggsignif::geom_signif(
									annotations = annotations,
									y_position = y_position, 
									xmin = x_min, 
									xmax = x_max,
									textsize = add_sig_text_size,
									color = "black"
								)
							}else{
								if(is.null(by_group)){
									comp_group_num <- sapply(use_diff_data$Comparison, function(x){length(unlist(gregexpr(" - ", x)))})
									y_position <- max(use_data$Value) + y_start * (max(use_data$Value) - min(use_data$Value))
									annotations <- as.character(use_diff_data[, add_sig_label])
									textdata <- data.frame(
										x = comp_group_num/2, 
										y = y_position, 
										add = annotations, 
										stringsAsFactors = FALSE
									)
									p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size, inherit.aes = FALSE)
								}else{
									x_mid <- c()
									annotations <- c()
									all_by_groups <- levels(use_data[, by_group])
									diff_by_groups <- use_diff_data$by_group %>% unique
									for(j in diff_by_groups){
										x_mid %<>% c(., match(j, all_by_groups))
										annotations %<>% c(., as.character(use_diff_data[use_diff_data$by_group == j, add_sig_label]))
									}
									textdata <- data.frame(
										x = x_mid, 
										y = max(use_data$Value) + y_start * (max(use_data$Value) - min(use_data$Value)), 
										add = annotations, 
										stringsAsFactors = FALSE
									)
									p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size, inherit.aes = FALSE)
								}
							}
						}
					}
				}
				p <- p + ylab(measure) + xlab("")
				p <- p + theme(
						axis.title.y = element_text(size = ytitle_size),
						axis.text.y = element_text(size = rel(1.1))
						)
				p <- p + ggplot_xtext_anglesize(xtext_angle, xtext_size)
				if(is.null(by_group)){
					p <- p + theme(legend.position = "none")
				}
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
		kwwitt_test = function(method = NULL, input_table = NULL, group = NULL, res_list = NULL, group_by = NULL, measure = NULL, by_ID = NULL, ...){
			if(!is.null(by_ID)){
				if(method == "KW"){
					stop("Please use wilcox instead of KW when by_ID is provided!")
				}
			}
			groupvec <- as.character(input_table[, group])
			use_comp_group_num <- unique(c(2, length(unique(groupvec))))
			for(i in use_comp_group_num){
				all_name <- combn(unique(groupvec), i)
				use_method <- switch(method, KW = "Kruskal-Wallis Rank Sum Test", wilcox = "Wilcoxon Rank Sum Test", t.test = "t.test")
				for(j in 1:ncol(all_name)){
					table_compare <- input_table[groupvec %in% as.character(all_name[, j]), ]
					if(is.null(by_ID)){
						table_compare[, group] %<>% factor(., levels = unique(as.character(.)))
						formu <- reformulate(group, "Value")
					}else{
						table_compare %<>% dropallfactors
						tmp <- table(table_compare[, by_ID])
						# check the missing item for paired comparison
						if(length(unique(tmp)) > 1){
							# delete the missing item
							tmp %<>% .[. != min(.)]
							table_compare %<>% .[.[, by_ID] %in% names(tmp), ]
						}
						order_ID <- unique(table_compare[, by_ID])
						order_group <- unique(table_compare[, group])
						x_value_table <- table_compare[table_compare[, group] == order_group[1], ]
						y_value_table <- table_compare[table_compare[, group] == order_group[2], ]
						x_value <- x_value_table[match(x_value_table[, by_ID], order_ID), "Value"]
						y_value <- y_value_table[match(y_value_table[, by_ID], order_ID), "Value"]
					}
					if(i == 2){
						if(method == "t.test"){
							if(is.null(by_ID)){
								res1 <- t.test(formu, data = table_compare, ...)
							}else{
								res1 <- t.test(x = x_value, y = y_value, paired = TRUE, ...)
							}
							max_group_select <- private$group_value_compare(table_compare$Value, table_compare[, group], mean)
						}else{
							if(method == "wilcox"){
								if(is.null(by_ID)){
									res1 <- suppressWarnings(wilcox.test(formu, data = table_compare, ...))
								}else{
									res1 <- suppressWarnings(wilcox.test(x = x_value, y = y_value, paired = TRUE, ...))
								}
								max_group_select <- private$group_value_compare(table_compare$Value, table_compare[, group], median)
							}else{
								if(method == "KW" & length(use_comp_group_num) == 1){
									res1 <- kruskal.test(formu, data = table_compare, ...)
									max_group_select <- private$group_value_compare(table_compare$Value, table_compare[, group], median)
								}else{
									next
								}
							}
						}
					}else{
						if(method == "KW"){
							res1 <- kruskal.test(formu, data = table_compare, ...)
							max_group_select <- private$group_value_compare(table_compare$Value, table_compare[, group], median)
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
		kdunn_test = function(input_table = NULL, group = NULL, measure = NULL, KW_dunn_letter = TRUE, p_adjust_method = NULL, alpha = 0.05, ...){
			use_method <- "Dunn's Kruskal-Wallis Multiple Comparisons"
			raw_groups <- input_table[, group]
			if(any(grepl("-", raw_groups))){
				input_table[, group] %<>% gsub("-", "sub&hyphen&sub", ., fixed = TRUE)
			}
			if(KW_dunn_letter){
				if(any(grepl("\\s", raw_groups))){
					input_table[, group] %<>% gsub(" ", "&space&", ., fixed = TRUE)
				}
			}
			orderd_groups <- tapply(input_table[, "Value"], input_table[, group], median) %>% sort(decreasing = TRUE) %>% names
			input_table[, group] %<>% factor(., levels = orderd_groups)
			formu <- reformulate(group, "Value")
			dunnTest_raw <- FSA::dunnTest(formu, data = input_table, method = p_adjust_method, ...)
			dunnTest_table <- dunnTest_raw$res
			# adjust the orders in each comparison
			tmp <- strsplit(dunnTest_table$Comparison, split = " - ")
			tmp <- lapply(tmp, function(x){
				matchord <- match(x, orderd_groups)
				if(matchord[1] > matchord[2]){
					rev(x)
				}else{
					x
				}
			})
			orderd_groups %<>% .[. %in% unique(unlist(tmp))]
			dunnTest_table$Comparison <- lapply(tmp, function(x){paste0(x, collapse = " - ")}) %>% unlist
			# generate ordered combined paired groups
			tmp <- combn(orderd_groups, 2) %>% t %>% as.data.frame %>% apply(., 1, function(x){paste0(x, collapse = " - ")})
			dunnTest_table <- dunnTest_table[match(tmp, dunnTest_table$Comparison), ]
			if(KW_dunn_letter){
				dunntest_final <- rcompanion::cldList(P.adj ~ Comparison, data = dunnTest_table, threshold = alpha,
					remove.space = TRUE, remove.equal = FALSE, remove.zero = FALSE, swap.colon = FALSE, swap.vs = FALSE)
				if(any(grepl("-", raw_groups))){
					dunntest_final$Group %<>% gsub("sub&hyphen&sub", "-", ., fixed = TRUE)
				}
				if(any(grepl("\\s", raw_groups))){
					dunntest_final$Group %<>% gsub("&space&", " ", .)
				}
				dunnTest_res <- data.frame(Measure = measure, Method = use_method, dunntest_final)
			}else{
				max_group <- lapply(dunnTest_table$Comparison, function(x){
					group_select <- unlist(strsplit(x, split = " - "))
					table_compare_select <- input_table[as.character(input_table[, group]) %in% group_select, ]
					private$group_value_compare(table_compare_select$Value, table_compare_select[, group], median)
				}) %>% unlist
				if(any(grepl("-", raw_groups))){
					dunnTest_table$Comparison %<>% gsub("sub&hyphen&sub", "-", ., fixed = TRUE)
					max_group %<>% gsub("sub&hyphen&sub", "-", ., fixed = TRUE)
				}
				dunnTest_res <- data.frame(Measure = measure, Method = use_method, Group = max_group, dunnTest_table)
			}
			dunnTest_res
		},
		anova_test = function(input_table, group, measure, post_test, alpha, varequal_test, ...){
			model <- aov(reformulate(group, "Value"), input_table)
			post_test_function <- get(post_test)
			out <- post_test_function(model, group, main = measure, alpha = alpha, ...)
			res1 <- out$groups[, "groups", drop = FALSE]
			res1$groups <- as.character(res1$groups)
			res1 <- data.frame(rownames(res1), res1, stringsAsFactors = FALSE, check.names = FALSE)
			colnames(res1) <- c("Group", "Letter")
			rownames(res1) <- NULL
			if(varequal_test){
				if(!is.factor(input_table[, group])){
					input_table[, group] %<>% as.factor
				}
				ve_res <- car::leveneTest(reformulate(group, "Value"), input_table)
				ve_p <- c(ve_res$`Pr(>F)`[1], rep(NA, nrow(res1) - 1))
				res <- data.frame(Measure = measure, Method = "anova", res1, leveneTest_Pvalue = ve_p)
			}else{
				res <- data.frame(Measure = measure, Method = "anova", res1)
			}
			res
		},
		group_value_compare = function(value, group, ...){
			group %<>% as.character
			group_values <- tapply(value, group, ...)
			if(any(is.na(group_values))){
				"NA"
			}else{
				group_values %>% {.[which.max(.)]} %>% names
			}
		},
		check_skip_comp = function(longdata, variable, by_group, by_group_var = NULL){
			if(length(unique(longdata$Value)) == 1){
				if(is.null(by_group)){
					message("Only one value is found! Skip ", variable, " ...")
				}else{
					message("Only one value is found! Skip ", variable, " for ", by_group, ": ", by_group_var," ...")
				}
				TRUE
			}else{
				FALSE
			}
		},
		formula_anova_sche_betareg_test = function(method, input_table, formula, measure, ...){
			if(method == "anova"){
				model <- aov(reformulate(formula, "Value"), input_table, ...)
				tmp <- summary(model)[[1]]
				tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), Measure = measure, Factors = gsub("\\s+$", "", rownames(tmp)), 
					Df = tmp$Df, Fvalue = tmp$`F value`, P.unadj = tmp$`Pr(>F)`)
			}
			if(method == "scheirerRayHare"){
				invisible(capture.output(tmp <- rcompanion::scheirerRayHare(reformulate(formula, "Value"), input_table, ...)))
				tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), Measure = measure, Factors = rownames(tmp), 
					Df = tmp$Df, Fvalue = tmp$H, P.unadj = tmp$p.value)
			}
			if(method == "betareg"){
				check_res <- tryCatch(tmp <- betareg::betareg(reformulate(formula, "Value"), data = input_table, ...), error = function(e) { skip_to_next <- TRUE})
				if(rlang::is_true(check_res)) {
					message("Model fitting failed for ", measure, " ! Skip ...")
					tmp_res <- NULL
				}else{
					tmp <- betareg::betareg(reformulate(formula, "Value"), data = input_table, ...)
					# extract the first element: coefficients
					tmp_coefficients <- summary(tmp)[[1]]
					tmp_mean <- tmp_coefficients$mean %>% as.data.frame
					tmp_precision <- tmp_coefficients$precision %>% as.data.frame
					tmp_res <- data.frame(Method = paste0(method, " formula for ", formula), Measure = measure, 
						Factors = c(rownames(tmp_mean), rownames(tmp_precision)), 
						Estimate = c(tmp_mean$Estimate, tmp_precision$Estimate), 
						Std.Error = c(tmp_mean$`Std. Error`, tmp_precision$`Std. Error`), 
						Zvalue = c(tmp_mean$`z value`, tmp_precision$`z value`), 
						P.unadj = c(tmp_mean$`Pr(>|z|)`, tmp_precision$`Pr(>|z|)`))
				}
			}
			tmp_res
		},
		# extract R2 from models
		R2_extract = function(model, number_supp){
			tmp_model_R2 <- try(performance::r2(model), silent = TRUE)
			if(inherits(tmp_model_R2, "try-error")) {
				R2data <- data.frame(R2 = NA)
			}else{
				if(all(is.na(tmp_model_R2))){
					R2data <- data.frame(R2 = NA)
				}else{
					R2data <- data.frame(R2 = c(
						paste0(unlist(lapply(tmp_model_R2, function(x){paste0(names(x), ": ", x)})), collapse = "; "), 
						rep(NA, number_supp)))
				}
			}
			R2data
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
