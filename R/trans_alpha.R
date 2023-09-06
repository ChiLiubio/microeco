#' @title Create \code{trans_alpha} object for alpha diversity statistics and plot.
#'
#' @description
#' This class is a wrapper for a series of alpha diversity related analysis, including the statistics and visualization based on 
#' An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035> and Paul et al. (2013) <doi:10.1371/journal.pone.0061217>.
#'
#' @export
trans_alpha <- R6Class(classname = "trans_alpha",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param group default NULL; a column of \code{sample_table} used for the statistics; If provided, can return \code{data_stat}.
		#' @param by_group default NULL; a column of \code{sample_table} used to perform the differential test 
		#'   among groups (\code{group} parameter) for each group (\code{by_group} parameter). So \code{by_group} has a higher level than \code{group} parameter.
		#' @param by_ID default NULL; a column of \code{sample_table} used to perform paired t test or paired wilcox test for the paired data,
		#'   such as the data of plant compartments for different plant species (ID). 
		#'   So \code{by_ID} in sample_table should be the smallest unit of sample collection without any repetition in it.
		#' @param order_x default NULL; a \code{sample_table} column name or a vector containg sample names; if provided, order samples by using \code{factor}.
		#' @return \code{data_alpha} and \code{data_stat} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Group")
		#' }
		initialize = function(dataset = NULL, group = NULL, by_group = NULL, by_ID = NULL, order_x = NULL) {
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
			if(! is.null(by_ID)){
				if(! by_ID %in% colnames(data_alpha)){
					stop("Provided by_ID: ", by_ID, " is not found in dataset$sample_table!")
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
			self$by_ID <- by_ID
		},
		#' @description
		#' Differential test of alpha diversity.
		#'
		#' @param method default "KW"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons, see \code{dunnTest} function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum Test for all paired groups}
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'anova'}}{Duncan's new multiple range test for one-way anova; see \code{duncan.test} function of \code{agricolae} package.
		#'     	  For multi-factor anova, see \code{\link{aov}}}
		#'     \item{\strong{'scheirerRayHare'}}{Scheirer Ray Hare test (nonparametric test) for a two-way factorial experiment; 
		#'     	  see \code{scheirerRayHare} function of \code{rcompanion} package}
		#'     \item{\strong{'lme'}}{Linear Mixed Effect Model based on the \code{lmerTest} package}
		#'     \item{\strong{'betareg'}}{Beta Regression for Rates and Proportions based on the \code{betareg} package}
		#'     \item{\strong{'glmm'}}{Generalized linear mixed model (GLMM) based on the \code{glmmTMB} package}
		#'   }
		#' @param measure default NULL; a vector; If NULL, all indexes will be calculated; see names of \code{microtable$alpha_diversity}, 
		#' 	 e.g. Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage and PD.
		#' @param p_adjust_method default "fdr" (for "KW", "wilcox", "t.test") or "holm" (for "KW_dunn"); P value adjustment method; 
		#' 	  For \code{method = 'KW', 'wilcox' or 't.test'}, please see method parameter of \code{p.adjust} function for available options;
		#' 	  For \code{method = 'KW_dunn'}, please see \code{dunn.test::p.adjustment.methods} for available options.
		#' @param formula default NULL; applied to two-way or multi-factor anova when 
		#'   method = \code{"anova"} or \code{"scheirerRayHare"} or \code{"lme"} or \code{"betareg"} or \code{"glmm"}; 
		#'   specified set for independent variables, i.e. the latter part of a general formula, 
		#'   such as \code{'block + N*P*K'}.
		#' @param KW_dunn_letter default TRUE; For \code{method = 'KW_dunn'}, \code{TRUE} denotes paired significances are presented by letters;
		#'   \code{FALSE} means significances are shown by asterisk for paired comparison.
		#' @param alpha default 0.05; Significant level; used for generating significance letters when method is 'anova' or 'KW_dunn'.
		#' @param return_model default FALSE; whether return the original lmer or glmm model list in the object.
		#' @param ... parameters passed to \code{kruskal.test} (\code{method = "KW"}) or \code{wilcox.test} function (\code{method = "wilcox"}) or 
		#'   \code{dunnTest} function of \code{FSA} package (\code{method = "KW_dunn"}) or 
		#'   \code{agricolae::duncan.test} (\code{method = "anova"}, one-way) or 
		#'   \code{rcompanion::scheirerRayHare} (\code{method = "scheirerRayHare"}) or 
		#'   \code{lmerTest::lmer} (\code{method = "lme"}) or 
		#'   \code{betareg::betareg} (\code{method = "betareg"}) or 
		#'   \code{glmmTMB::glmmTMB} (\code{method = "glmm"}).
		#' @return \code{res_diff}, stored in object with the format \code{data.frame}.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(method = "KW")
		#' t1$cal_diff(method = "anova")
		#' t1 <- trans_alpha$new(dataset = dataset, group = "Type", by_group = "Group")
		#' t1$cal_diff(method = "anova")
		#' }
		cal_diff = function(
			method = c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lme", "betareg", "glmm")[1], 
			measure = NULL, 
			p_adjust_method = "fdr", 
			formula = NULL,
			KW_dunn_letter = TRUE,
			alpha = 0.05,
			return_model = FALSE,
			...
			){
			method <- match.arg(method, c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lme", "betareg", "glmm"))
			group <- self$group
			
			if(method %in% c("scheirerRayHare", "lme", "betareg", "glmm") & is.null(formula)){
				if(is.null(formula)){
					stop("formula is necessary! Please provide formula parameter!")
				}
			}
			if(!method %in% c("anova", "scheirerRayHare", "lme", "betareg", "glmm")){
				if(is.null(group)){
					stop("For the method: ", method, " , group is necessary! Please recreate the object!")
				}
			}
			if(method == "anova"){
				if(is.null(group) & is.null(formula)){
					stop("Both formula and group is NULL! Please provide either formula or group in the object creation!")
				}
			}
			# 'by_group' for test inside each by_group instead of all groups in 'group'
			by_group <- self$by_group
			data_alpha <- self$data_alpha
			by_ID <- self$by_ID

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
					colnames(compare_result) <- c("Comparison", "Measure", "Test_method", "Group", "P.unadj", "P.adj")
				}else{
					compare_result <- data.frame(res_list$comnames, res_list$group_by, res_list$measure_use, res_list$test_method, 
						res_list$max_group, res_list$p_value, p_value_adjust)
					colnames(compare_result) <- c("Comparison", "by_group", "Measure", "Test_method", "Group", "P.unadj", "P.adj")
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
							tmp_res <- private$kdunn_test(input_table = div_table, group = group, measure = k, KW_dunn_letter = KW_dunn_letter, 
								p_adjust_method = p_adjust_method, alpha = alpha, ...)
							tmp_res <- cbind.data.frame(by_group = each_group, tmp_res)
							compare_result %<>% rbind(., tmp_res)
						}
					}
				}
			}
			if(method == "anova" & is.null(formula)){
				compare_result <- data.frame()
				for(k in measure){
					if(is.null(by_group)){
						div_table <- data_alpha[data_alpha$Measure == k, ]
						tmp_res <- private$anova_test(input_table = div_table, group = group, measure = k, alpha = alpha, ...)
						compare_result %<>% rbind(., tmp_res)
					}else{
						for(each_group in unique_bygroups){
							test <- data_alpha[all_bygroups == each_group, ]
							if(length(unique(as.character(test[, group]))) < 2){
								message("Skip the by_group: ", each_group, " as groups number < 2!")
								next
							}
							div_table <- data_alpha[data_alpha$Measure == k & all_bygroups == each_group, ]
							tmp_res <- private$anova_test(input_table = div_table, group = group, measure = k, alpha = alpha, ...)
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
					div_table <- data_alpha[data_alpha$Measure == k, ]
					if(method == "anova"){
						model <- aov(reformulate(formula, "Value"), div_table, ...)
						tmp <- summary(model)[[1]]
						tmp1 <- data.frame(method = paste0(method, " formula for ", formula), Measure = k, Factors = rownames(tmp), 
							Df = tmp$Df, Fvalue = tmp$`F value`, P.unadj = tmp$`Pr(>F)`)
					}
					if(method == "scheirerRayHare"){
						invisible(capture.output(tmp <- rcompanion::scheirerRayHare(reformulate(formula, "Value"), div_table, ...)))
						tmp1 <- data.frame(method = paste0(method, " formula for ", formula), Measure = k, Factors = rownames(tmp), 
							Df = tmp$Df, Fvalue = tmp$H, P.unadj = tmp$p.value)
					}
					if(method == "betareg"){
						check_res <- tryCatch(tmp <- betareg::betareg(reformulate(formula, "Value"), data = div_table, ...), error = function(e) { skip_to_next <- TRUE})
						if(rlang::is_true(check_res)) {
							message("Model fitting failed for ", k, " ! Skip ...")
							next
						}else{
							tmp <- betareg::betareg(reformulate(formula, "Value"), data = div_table, ...)
							# extract the first element: coefficients
							tmp_coefficients <- summary(tmp)[[1]]
							tmp_mean <- tmp_coefficients$mean %>% as.data.frame
							tmp_precision <- tmp_coefficients$precision %>% as.data.frame
							tmp1 <- data.frame(method = paste0(method, " formula for ", formula), Measure = k, 
								Factors = c(rownames(tmp_mean), rownames(tmp_precision)), 
								Estimate = c(tmp_mean$Estimate, tmp_precision$Estimate), 
								Zvalue = c(tmp_mean$`z value`, tmp_precision$`z value`), 
								P.unadj = c(tmp_mean$`Pr(>|z|)`, tmp_precision$`Pr(>|z|)`))
						}
					}
					compare_result %<>% rbind(., tmp1)
				}
				compare_result$Significance <- cut(compare_result$P.unadj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", "ns"))
				method <- paste0(method, "-formula")
			}
			if(method %in% c("lme", "glmm")){
				compare_result <- data.frame()
				if(return_model){
					res_model <- list()
				}
				for(k in measure){
					div_table <- data_alpha[data_alpha$Measure == k, ]
					tmp_res <- data.frame()
					if(method == "lme"){
						tmp <- lmerTest::lmer(reformulate(formula, "Value"), data = div_table, ...)
						if(return_model){
							res_model[[k]] <- tmp
						}
						tmp_summary <- summary(tmp)
						tmp_coefficients <- as.data.frame(tmp_summary$coefficients, check.names = FALSE)
						tmp_model_R2 <- performance::r2(tmp)
						tmp_model_p <- anova(tmp)
						tmp_random_p <- lmerTest::ranova(tmp)
						tmp_res <- data.frame(method = paste0(method, " formula for ", formula), 
							Measure = k, 
							Factors = c("Model", rownames(tmp_model_p), rownames(tmp_random_p), rownames(tmp_coefficients)), 
							Conditional_R2 = c(tmp_model_R2$R2_conditional, rep(NA, nrow(tmp_model_p) + nrow(tmp_random_p) + length(rownames(tmp_coefficients)))),
							Marginal_R2 = c(tmp_model_R2$R2_marginal, rep(NA, nrow(tmp_model_p) + nrow(tmp_random_p) + length(rownames(tmp_coefficients)))),
							Estimate = c(NA, rep(NA, nrow(tmp_model_p) + nrow(tmp_random_p)), tmp_coefficients$Estimate), 
							P.unadj = c(NA, tmp_model_p$`Pr(>F)`, tmp_random_p$`Pr(>Chisq)`, tmp_coefficients$`Pr(>|t|)`)
						)
					}else{
						tmp <- glmmTMB::glmmTMB(reformulate(formula, "Value"), data = div_table, ...)
						if(return_model){
							res_model[[k]] <- tmp
						}
						tmp_summary <- summary(tmp)
						tmp_coefficients <- as.data.frame(tmp_summary$coefficients$cond, check.names = FALSE)
						tmp_model_p <- car::Anova(tmp)
						tmp_model_R2 <- performance::r2(tmp)
						test <- try(tmp_model_R2$R2_conditional, silent = TRUE)
						if(inherits(test, "try-error")) {
							message("R2 unavailable for ", k, " !")
							tmp_model_R2 <- list(R2_conditional = NA, R2_marginal = NA)
						}
						tmp_res <- data.frame(method = paste0(method, " formula for ", formula), 
							Measure = k, 
							Factors = c("Model", rownames(tmp_model_p), rownames(tmp_coefficients)), 
							Conditional_R2 = c(tmp_model_R2$R2_conditional, rep(NA, nrow(tmp_model_p) + length(rownames(tmp_coefficients)))),
							Marginal_R2 = c(tmp_model_R2$R2_marginal, rep(NA, nrow(tmp_model_p) + length(rownames(tmp_coefficients)))),
							Estimate = c(NA, rep(NA, nrow(tmp_model_p)), tmp_coefficients$Estimate), 
							P.unadj = c(NA, tmp_model_p$`Pr(>Chisq)`, tmp_coefficients$`Pr(>|z|)`)
						)
					}
					compare_result %<>% rbind(., tmp_res)
				}
				if(return_model){
					self$res_model <- res_model
					message("The original ", method, " models list is stored in object$res_model ...")
				}
			}
			if(! method %in% c("anova", paste0(c("anova", "scheirerRayHare", "betareg"), "-formula"), "lme", "glmm")){
				if("P.adj" %in% colnames(compare_result)){
					compare_result$Significance <- cut(compare_result$P.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", "ns"))
					compare_result$Significance %<>% as.character
				}
			}
			self$res_diff <- compare_result
			self$cal_diff_method <- method
			self$measure <- measure
			message('The result is stored in object$res_diff ...')
		},
		#' @description
		#' Plot the alpha diversity.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color pallete for groups.
		#' @param measure default Shannon; one alpha diversity measurement; see names of alpha_diversity of dataset, 
		#'   e.g., Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher, Coverage, PD.
		#' @param group default NULL; group name used for the plot.
		#' @param add_sig default TRUE; wheter add significance label using the result of \code{cal_diff} function, i.e. \code{object$res_diff};
		#'   This is manily designed to add post hoc test of anova or other significances to make the label mapping easy.
		#' @param add_sig_label default "Significance"; select a colname of \code{object$res_diff} for the label text when 'Letter' is not in the table, 
		#'   such as 'P.adj' or 'Significance'.
		#' @param add_sig_text_size default 3.88; the size of text in added label.
		#' @param use_boxplot default TRUE; TRUE: boxplot; FALSE: mean-se plot.
		#' @param boxplot_add default "jitter"; points type, see the add parameter in \code{ggpubr::ggboxplot}.
		#' @param order_x_mean default FALSE; whether order x axis by the means of groups from large to small.
		#' @param y_start default 0.1; the y axis value from which to add the significance asterisk label; 
		#' 	  the default 0.1 means \code{max(values) + 0.1 * (max(values) - min(values))}.
		#' @param y_increase default 0.05; the increasing y axia space to add the label (asterisk or letter); the default 0.05 means \code{0.05 * (max(values) - min(values))}; 
		#' 	  this parameter is also used to label the letters of anova result with the fixed space.
		#' @param xtext_angle default NULL; number (e.g. 30) used to make x axis text generate angle.
		#' @param xtext_size default 15; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @param barwidth default 0.9; the bar width in plot; applied when by_group is not NULL.
		#' @param ... parameters pass to \code{ggpubr::ggboxplot} function.
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
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			measure = "Shannon",
			group = NULL,
			add_sig = TRUE,
			add_sig_label = "Significance",
			add_sig_text_size = 3.88,
			use_boxplot = TRUE,
			boxplot_add = "jitter",
			order_x_mean = FALSE,
			y_start = 0.1,
			y_increase = 0.05,
			xtext_angle = NULL,
			xtext_size = 15,
			ytitle_size = 17,
			barwidth = 0.9,
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
			if(order_x_mean){
				if(is.null(by_group)){
					mean_orders <- names(sort(tapply(use_data$Value, use_data[, group], mean), decreasing = TRUE))
				}else{
					# first sort by_group
					mean_orders_bygroup <- names(sort(tapply(use_data$Value, use_data[, by_group], mean), decreasing = TRUE))
					use_data[, by_group] %<>% factor(., levels = mean_orders_bygroup)
					# sort each group by by_group levels
					mean_orders <- lapply(mean_orders_bygroup, function(x){
						tmp_data <- use_data[use_data[, by_group] == x, ]
						names(sort(tapply(tmp_data$Value, tmp_data[, group], mean), decreasing = TRUE))
					}) %>% unlist
				}
				use_data[, group] %<>% factor(., levels = mean_orders)
			}else{
				if(!is.factor(use_data[, group])){
					use_data[, group] %<>% as.factor
				}
				if(!is.null(by_group)){
					if(!is.factor(use_data[, by_group])){
						use_data[, by_group] %<>% as.factor
					}
				}
			}
			color_values <- expand_colors(color_values, length(unique(use_data[, group])))
			if(use_boxplot){
				if(is.null(by_group)){
					p <- ggpubr::ggboxplot(
						use_data, x = group, y = "Value", color = group, 
						palette = color_values, add = boxplot_add, outlier.colour = "white", 
						...
						)
				}else{
					p <- ggpubr::ggboxplot(
						use_data, x = by_group, y = "Value", color = group, 
						palette = color_values, add = boxplot_add, outlier.colour = "white", 
						...
						)
				}
			}else{
				p <- ggplot(use_data, aes_meco(x = group, y = "Value")) + 
					theme_minimal() +
					stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
					stat_summary(fun = mean, geom = "point", size = rel(3)) + 
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
			if(add_sig){
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
							for(j in all_by_groups){
								select_tmp_data <- tmp_use_data %>% .[.[, by_group] == j, ]
								x_axis_order <- all_groups %>% .[. %in% select_tmp_data[, group]]
								# identify the start x position for each by_group
								start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
								increase_bar_mid <- barwidth/length(x_axis_order)
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
									{. + y_increase * (max(select_tmp_data$Value) - min(select_tmp_data$Value))} %>% 
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
						p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size)
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
								use_diff_data[, add_sig_label] %<>% round(., 4)
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
									annotations %<>% c(., use_diff_data[i, add_sig_label])
									x_min %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\1", use_diff_data[i, "Comparison"]), x_axis_order))
									x_max %<>% c(., match(gsub("(.*)\\s-\\s(.*)", "\\2", use_diff_data[i, "Comparison"]), x_axis_order))
									y_position %<>% c(., y_start + i * y_increase * y_range)
								}
							}else{
								all_by_groups <- levels(use_data[, by_group])
								all_groups <- levels(use_data[, group])
								tmp_use_data <- dropallfactors(use_data)
								
								diff_by_groups <- use_diff_data$by_group %>% unique
								for(j in diff_by_groups){
									select_tmp_data <- tmp_use_data %>% .[.[, by_group] == j, ]
									# y axix starting position
									y_range_use <- max(select_tmp_data$Value) - min(select_tmp_data$Value)
									y_start_use <- max(select_tmp_data$Value) + y_start * y_range_use
									x_axis_order <- all_groups %>% .[. %in% select_tmp_data[, group]]
									# identify the start x position for each by_group
									start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
									increase_bar_mid <- barwidth/length(x_axis_order)
									# loop for each group of use_diff_data
									tmp_diff_res <- use_diff_data[use_diff_data$by_group == j, ]
									for(i in seq_len(nrow(tmp_diff_res))){
										annotations %<>% c(., tmp_diff_res[i, add_sig_label])
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
								textsize = add_sig_text_size
							)
						}else{
							if(is.null(by_group)){
								comp_group_num <- sapply(use_diff_data$Comparison, function(x){length(unlist(gregexpr(" - ", x)))})
								y_position <- max(use_data$Value) + y_start * (max(use_data$Value) - min(use_data$Value))
								annotations <- use_diff_data[, add_sig_label]
								textdata <- data.frame(
									x = comp_group_num/2, 
									y = y_position, 
									add = annotations, 
									stringsAsFactors = FALSE
								)
								p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size)
							}else{
								x_mid <- c()
								annotations <- c()
								all_by_groups <- levels(use_data[, by_group])
								diff_by_groups <- use_diff_data$by_group %>% unique
								for(j in diff_by_groups){
									x_mid %<>% c(., match(j, all_by_groups))
									annotations %<>% c(., use_diff_data[use_diff_data$by_group == j, add_sig_label])
								}
								textdata <- data.frame(
									x = x_mid, 
									y = max(use_data$Value) + y_start * (max(use_data$Value) - min(use_data$Value)), 
									add = annotations, 
									stringsAsFactors = FALSE
								)
								p <- p + geom_text(aes(x = x, y = y, label = add), data = textdata, size = add_sig_text_size)
							}
						}
					}
				}
			}
			if(is.null(by_group)){
				p <- p + theme(legend.position="none")
			}
			p <- p + ylab(measure) + xlab("")
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
				input_table[, group] %<>% gsub("-", "sub&&&sub", ., fixed = TRUE)
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
				dunnTest_final <- rcompanion::cldList(P.adj ~ Comparison, data = dunnTest_table, threshold = alpha)
				if(any(grepl("-", raw_groups))){
					dunnTest_final$Group %<>% gsub("sub&&&sub", "-", ., fixed = TRUE)
				}
				dunnTest_res <- data.frame(Measure = measure, Test_method = use_method, dunnTest_final)
			}else{
				max_group <- lapply(dunnTest_table$Comparison, function(x){
					group_select <- unlist(strsplit(x, split = " - "))
					table_compare_select <- input_table[as.character(input_table[, group]) %in% group_select, ]
					private$group_value_compare(table_compare_select$Value, table_compare_select[, group], median)
				}) %>% unlist
				if(any(grepl("-", raw_groups))){
					dunnTest_table$Comparison %<>% gsub("sub&&&sub", "-", ., fixed = TRUE)
					max_group %<>% gsub("sub&&&sub", "-", ., fixed = TRUE)
				}
				dunnTest_res <- data.frame(Measure = measure, Test_method = use_method, Group = max_group, dunnTest_table)
			}
			dunnTest_res
		},
		anova_test = function(input_table = NULL, group = NULL, measure = NULL, alpha = 0.05, ...){
			model <- aov(reformulate(group, "Value"), input_table)
			out <- agricolae::duncan.test(model, group, main = measure, alpha = alpha, ...)
			res1 <- out$groups[, "groups", drop = FALSE]
			res1$groups <- as.character(res1$groups)
			res1 <- data.frame(rownames(res1), res1, stringsAsFactors = FALSE, check.names = FALSE)
			colnames(res1) <- c("Group", "Letter")
			rownames(res1) <- NULL
			res2 <- data.frame(Measure = measure, Test_method = "anova", res1)
			res2
		},
		group_value_compare = function(value, group, ...){
			group %<>% as.character
			group_values <- tapply(value, group, ...)
			if(any(is.na(group_values))){
				"NA"
			}else{
				group_values %>% {.[which.max(.)]} %>% names
			}
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
