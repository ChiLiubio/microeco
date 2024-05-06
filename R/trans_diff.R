#' @title 
#' Create \code{trans_diff} object for the differential analysis on the taxonomic abundance
#'
#' @description
#' This class is a wrapper for a series of differential abundance test and indicator analysis methods, including 
#'  LEfSe based on the Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>,
#'  random forest <doi:10.1016/j.geoderma.2018.09.035>, metastat based on White et al. (2009) <doi:10.1371/journal.pcbi.1000352>,
#'  non-parametric Kruskal-Wallis Rank Sum Test,
#'  Dunn's Kruskal-Wallis Multiple Comparisons based on the \code{FSA} package, Wilcoxon Rank Sum and Signed Rank Tests, t-test, anova, 
#'  Scheirer Ray Hare test, 
#'  R package \code{metagenomeSeq} Paulson et al. (2013) <doi:10.1038/nmeth.2658>, 
#'  R package \code{ANCOMBC} <doi:10.1038/s41467-020-17041-7>, R package \code{ALDEx2} <doi:10.1371/journal.pone.0067019; 10.1186/2049-2618-2-15>, 
#'  R package \code{MicrobiomeStat} <doi:10.1186/s13059-022-02655-5>, beta regression <doi:10.18637/jss.v034.i02>, R package \code{maaslin2},
#'  linear mixed-effects model and generalized linear mixed model.
#'  
#' @export
trans_diff <- R6Class(classname = "trans_diff",
	public = list(
		#' @param dataset default NULL; \code{\link{microtable}} object.
		#' @param method default "lefse"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'lefse'}}{LEfSe method based on Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>}
		#'     \item{\strong{'rf'}}{random forest and non-parametric test method based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>}
		#'     \item{\strong{'metastat'}}{Metastat method for all paired groups based on White et al. (2009) <doi:10.1371/journal.pcbi.1000352>}
		#'     \item{\strong{'metagenomeSeq'}}{zero-inflated log-normal model-based differential test method from \code{metagenomeSeq} package.}
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons when group number > 2; see dunnTest function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'anova'}}{ANOVA for one-way or multi-factor analysis; see \code{cal_diff} function of \code{trans_alpha} class}
		#'     \item{\strong{'scheirerRayHare'}}{Scheirer Ray Hare test for nonparametric test used for a two-way factorial experiment; 
		#'     	  see \code{scheirerRayHare} function of \code{rcompanion} package}
		#'     \item{\strong{'lm'}}{Linear Model based on the \code{lm} function}
		#'     \item{\strong{'ALDEx2_t'}}{runs Welch's t and Wilcoxon tests with \code{ALDEx2} package; see also the test parameter in \code{ALDEx2::aldex} function;
		#'     	  ALDEx2 uses the centred log-ratio (clr) transformation and estimates per-feature technical variation within each sample using Monte-Carlo instances 
		#'     	  drawn from the Dirichlet distribution; Reference: <doi:10.1371/journal.pone.0067019> and <doi:10.1186/2049-2618-2-15>; 
		#'     	  require \code{ALDEx2} package to be installed 
		#'     	  (\href{https://bioconductor.org/packages/release/bioc/html/ALDEx2.html}{https://bioconductor.org/packages/release/bioc/html/ALDEx2.html})}
		#'     \item{\strong{'ALDEx2_kw'}}{runs Kruskal-Wallace and generalized linear model (glm) test with \code{ALDEx2} package; 
		#'     	  see also the \code{test} parameter in \code{ALDEx2::aldex} function.}
		#'     \item{\strong{'DESeq2'}}{Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution based on the \code{DESeq2} package.}
		#'     \item{\strong{'edgeR'}}{The \code{exactTest} method of \code{edgeR} package is implemented.}
		#'     \item{\strong{'ancombc2'}}{Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC) 
		#'        based on the \code{ancombc2} function from \code{ANCOMBC} package.
		#'        If the \code{fix_formula} parameter is not provided, the function can automatically assign it by using group parameter.
		#'        For this method, the \code{group} parameter is directly passed to the group parameter of \code{ancombc2} function.
		#'        Reference: <doi:10.1038/s41467-020-17041-7><10.1038/s41592-023-02092-7>; Require \code{ANCOMBC} package to be installed 
		#'        (\href{https://bioconductor.org/packages/release/bioc/html/ANCOMBC.html}{https://bioconductor.org/packages/release/bioc/html/ANCOMBC.html})}
		#'     \item{\strong{'linda'}}{Linear Model for Differential Abundance Analysis of High-dimensional Compositional Data 
		#'     	  based on the \code{linda} function of \code{MicrobiomeStat} package. 
		#'     	  For linda method, please provide either the group parameter or the formula parameter.
		#'     	  When the formula parameter is provided, it should start with '~' as it is directly used by the linda function.
		#'     	  If the group parameter is used, the prefix '~' is not necessary as the function can automatically add it.
		#'     	  The parameter \code{feature.dat.type = 'count'} has been fixed. Other parameters can be passed to the \code{linda} function.
		#'     	  Reference: <doi:10.1186/s13059-022-02655-5>}
		#'     \item{\strong{'maaslin2'}}{finding associations between metadata and potentially high-dimensional microbial multi-omics data based on the Maaslin2 package.
		#'     	  Using this option can invoke the \code{trans_env$cal_cor} function with \code{cor_method = "maaslin2"}.}
		#'     \item{\strong{'betareg'}}{Beta Regression based on the \code{betareg} package. 
		#'     	  Please see the \code{beta_pseudo} parameter for the use of pseudo value when there is 0 or 1 in the data}
		#'     \item{\strong{'lme'}}{Linear Mixed Effect Model based on the \code{lmerTest} package.
		#'     	  In the return table, the significance of fixed factors are tested by function \code{anova}.
		#'     	  The significance of 'Estimate' in each term of fixed factors comes from the model.}
		#'     \item{\strong{'glmm'}}{Generalized linear mixed model (GLMM) based on the \code{glmmTMB} package.
		#'     	  The \code{formula} and \code{family} parameters are needed. 
		#'     	  Please refer to glmmTMB package to select the family function, e.g. \code{family = glmmTMB::lognormal(link = "log")}.
		#'     	  The usage of formula is similar with that in 'lme' method.
		#'     	  For more available parameters, please see \code{glmmTMB::glmmTMB} function and use parameter passing.
		#'     	  In the return table, Conditional_R2 and Marginal_R2 represent total variance (explained by both fixed and random effects) and the variance explained by 
		#'     	  fixed effects, respectively. The significance of fixed factors are tested by Chi-square test from function \code{car::Anova}.
		#'     	  The significance of 'Estimate' in each term of fixed factors comes from the model.}
		#'     \item{\strong{'glmm_beta'}}{Generalized linear mixed model with a family function of beta distribution, 
		#'     	  developed for the relative abundance (ranging from 0 to 1) of taxa specifically. 
		#'     	  This is an extension of the GLMM model in \code{'glmm'} option.
		#'     	  The only difference is in \code{glmm_beta} the family function is fixed with the beta distribution function, 
		#'     	  i.e. \code{family = glmmTMB::beta_family(link = "logit")}.
		#'     	  Please see the \code{beta_pseudo} parameter for the use of pseudo value when there is 0 or 1 in the data}
		#'   }
		#' @param group default NULL; sample group used for the comparision; a colname of input \code{microtable$sample_table};
		#' 	  It is necessary when method is not "anova" or method is "anova" but formula is not provided.
		#' 	  Once group is provided, the return res_abund will have mean and sd values for group.
		#' @param taxa_level default "all"; 'all' represents using abundance data at all taxonomic ranks; 
		#' 	  For testing at a specific rank, provide taxonomic rank name, such as "Genus".
		#' 	  If the provided taxonomic name is neither 'all' nor a colname in tax_table of input dataset, 
		#' 	  the function will use the features in input \code{microtable$otu_table} automatically.
		#' @param filter_thres default 0; the abundance threshold, such as 0.0005 when the input is relative abundance; only available when method != "metastat".
		#' 	  The features with abundances lower than filter_thres will be filtered.
		#' @param alpha default 0.05; significance threshold to select taxa when method is "lefse" or "rf"; 
		#'    or used to generate significance letters when method is 'anova' or 'KW_dunn' like the alpha parameter in \code{cal_diff} of \code{trans_alpha} class.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of \code{p.adjust} function for other available options; 
		#'    "none" means disable p value adjustment; So when \code{p_adjust_method = "none"}, P.adj is same with P.unadj.
		#' @param transformation default NULL; feature abundance transformation method in the class \code{\link{trans_norm}},
		#'    such as 'AST' for the arc sine square root transformation.
		#'    Only available when \code{method} is one of "KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "betareg" and "lme".
		#' @param remove_unknown default TRUE; whether remove unknown features that donot have clear classification information.
		#' @param lefse_subgroup default NULL; sample sub group used for sub-comparision in lefse; Segata et al. (2011) <doi:10.1186/gb-2011-12-6-r60>.
		#' @param lefse_min_subsam default 10; sample numbers required in the subgroup test.
		#' @param lefse_norm default 1000000; scale value in lefse.
		#' @param nresam default 0.6667; sample number ratio used in each bootstrap for method = "lefse" or "rf".
		#' @param boots default 30; bootstrap test number for method = "lefse" or "rf".
		#' @param rf_ntree default 1000; see ntree in randomForest function of randomForest package when method = "rf".
		#' @param group_choose_paired default NULL; a vector used for selecting the required groups for paired testing, only used for method = "metastat" or "metagenomeSeq".
		#' @param metagenomeSeq_count default 1; Filter features to have at least 'counts' counts.; see the count parameter in MRcoefs function of \code{metagenomeSeq} package.
		#' @param ALDEx2_sig default c("wi.eBH", "kw.eBH"); which column of the final result is used as the significance asterisk assignment;
		#'   applied to method = "ALDEx2_t" or "ALDEx2_kw"; the first element is provided to "ALDEx2_t"; the second is provided to "ALDEx2_kw";
		#'   for "ALDEx2_t", the available choice is "wi.eBH" (Expected Benjamini-Hochberg corrected P value of Wilcoxon test)
		#'   and "we.eBH" (Expected BH corrected P value of Welch's t test); for "ALDEx2_kw"; for "ALDEx2_t",
		#'   the available choice is "kw.eBH" (Expected BH corrected P value of Kruskal-Wallace test) and "glm.eBH" (Expected BH corrected P value of glm test).
		#' @param by_group default NULL; a column of sample_table used to perform the differential test 
		#'   among groups (\code{group} parameter) for each group (\code{by_group} parameter). So \code{by_group} has a higher level than \code{group} parameter.
		#'   Same with the \code{by_group} parameter in \code{trans_alpha} class. 
		#'   Only available when method is one of \code{c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare")}.
		#' @param by_ID default NULL; a column of sample_table used to perform paired t test or paired wilcox test for the paired data,
		#'   such as the data of plant compartments for different plant species (ID). 
		#'   So \code{by_ID} in sample_table should be the smallest unit of sample collection without any repetition in it.
		#'   Same with the \code{by_ID} parameter in trans_alpha class.
		#' @param beta_pseudo default .Machine$double.eps; the pseudo value used when the parameter \code{method} is \code{'betareg'} or \code{'glmm_beta'}.
		#'   As the beta distribution function limits 0 < response value < 1, a pseudo value will be added for the data that equal to 0.
		#'   The data that equal to 1 will be replaced by \code{1/(1 + beta_pseudo)}.
		#' @param ... parameters passed to \code{cal_diff} function of \code{trans_alpha} class when method is one of 
		#' 	 "KW", "KW_dunn", "wilcox", "t.test", "anova", "betareg", "lme", "glmm" or "glmm_beta";
		#' 	 passed to \code{ANCOMBC::ancombc2} function when method is "ancombc2" (except tax_level, global and fix_formula parameters);
		#' 	 passed to \code{ALDEx2::aldex} function when method = "ALDEx2_t" or "ALDEx2_kw";
		#' 	 passed to \code{DESeq2::DESeq} function when method = "DESeq2";
		#' 	 passed to \code{MicrobiomeStat::linda} function when method = "linda";
		#' 	 passed to \code{trans_env$cal_cor} function when method = "maaslin2".
		#' @return res_diff and res_abund.\cr
		#'   \strong{res_abund} includes mean abundance of each taxa (Mean), standard deviation (SD), standard error (SE) and sample number (N) in the group (Group).\cr
		#'   \strong{res_diff} is the detailed differential test result, may containing:\cr
		#'     \strong{"Comparison"}: The groups for the comparision, maybe all groups or paired groups. If this column is not found, all groups are used;\cr
		#'     \strong{"Group"}: Which group has the maximum median or mean value across the test groups; 
		#'        For non-parametric methods, median value; For t.test, mean value;\cr
		#'     \strong{"Taxa"}: which taxa is used in this comparision;\cr
		#'     \strong{"Method"}: Test method used in the analysis depending on the method input;\cr
		#'     \strong{"LDA" or "MeanDecreaseGini"}: LDA: linear discriminant score in LEfSe; MeanDecreaseGini: mean decreasing gini index in random forest;\cr
		#'     \strong{"P.unadj"}: original p value;\cr
		#'     \strong{"P.adj"}: adjusted p value;\cr
		#'     \strong{Others}: qvalue: qvalue in metastat analysis.
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
			method = c("lefse", "rf", "metastat", "metagenomeSeq", "KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lm",
				"ancombc2", "ALDEx2_t", "ALDEx2_kw", "DESeq2", "edgeR", "linda", "maaslin2", "betareg", "lme", "glmm", "glmm_beta")[1],
			group = NULL,
			taxa_level = "all",
			filter_thres = 0,
			alpha = 0.05,
			p_adjust_method = "fdr",
			transformation = NULL,
			remove_unknown = TRUE,
			lefse_subgroup = NULL,
			lefse_min_subsam = 10,
			lefse_norm = 1000000,
			nresam = 0.6667,
			boots = 30,
			rf_ntree = 1000,
			group_choose_paired = NULL,
			metagenomeSeq_count = 1,
			ALDEx2_sig = c("wi.eBH", "kw.eBH"),
			by_group = NULL,
			by_ID = NULL,
			beta_pseudo = .Machine$double.eps,
			...
			){
			if(is.null(p_adjust_method)){
				message("Redefine p_adjust_method = ", sQuote("fdr"), " instead of NULL. To disable p value adjustment, please set p_adjust_method = ", sQuote("none"), " ...")
				p_adjust_method <- "fdr"
			}
			if(is.null(dataset)){
				self$method <- NULL
				message("Input dataset is NULL. Please run the functions with customized data ...")
			}else{
				method <- match.arg(method, c("lefse", "rf", "metastat", "metagenomeSeq", "KW", "KW_dunn", "wilcox", "t.test", 
					"anova", "scheirerRayHare", "lm", "ancombc2", "ALDEx2_t", "ALDEx2_kw", "DESeq2", "edgeR", "linda", "maaslin2", "betareg", "lme", "glmm", "glmm_beta"))

				tmp_dataset <- clone(dataset)
				sampleinfo <- tmp_dataset$sample_table
				if(is.null(group) & ! method %in% c("anova", "scheirerRayHare", "lm", "betareg", "lme", "glmm", "glmm_beta", "maaslin2", "ancombc2", "linda", "DESeq2")){
					stop("The group parameter is necessary for differential test method: ", method, " !")
				}
				if(!is.null(group)){
					if(! method %in% c("linda", "DESeq2")){
						if(length(group) > 1){
							stop("Please provide only one colname of sample_table for group parameter!")
						}
						if(! group %in% colnames(sampleinfo)){
							stop("Please provide a correct colname of sample_table for group parameter!")
						}
						if(is.factor(sampleinfo[, group])){
							self$group_order <- levels(sampleinfo[, group])
							sampleinfo[, group] %<>% as.character
						}else{
							self$group_order <- unique(as.character(sampleinfo[, group]))
						}
					}
				}
				check_taxa_abund(tmp_dataset)
				
				if(grepl("all", taxa_level, ignore.case = TRUE)){
					abund_table <- do.call(rbind, unname(tmp_dataset$taxa_abund))
				}else{
					if(! taxa_level %in% names(tmp_dataset$taxa_abund)){
						message("Provided taxa_level: ", taxa_level, " not in tax_table of dataset; use features in otu_table ...")
						tmp_dataset$add_rownames2taxonomy(use_name = taxa_level)
						suppressMessages(tmp_dataset$cal_abund(rel = TRUE))
					}
					abund_table <- tmp_dataset$taxa_abund[[taxa_level]]
				}
				
				message(nrow(abund_table), " input features ...")
				filter_output <- filter_lowabund_feature(abund_table = abund_table, filter_thres = filter_thres)
				abund_table <- filter_output$abund_table
				filter_features <- filter_output$filter_features
				
				if(method %in% c("lefse", "rf", "KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lm", "betareg", "lme", "glmm", "glmm_beta")){
					if(remove_unknown){
						abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]}
						message(nrow(abund_table), " features are remained after removing unknown features ...")
						if(nrow(abund_table) == 0){
							stop("No available feature! Please set the parameter: remove_unknown = FALSE")
						}
					}
					if(method == "lefse"){
						abund_table %<>% {. * lefse_norm}
						self$lefse_norm <- lefse_norm
					}
					
					if(method %in% c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lm", "betareg", "lme", "glmm", "glmm_beta")){
						abund_foralpha <- abund_table
						if(!is.null(transformation)){
							message("Perform data transformation with method: ", transformation, " ...")
							tmp <- suppressMessages(trans_norm$new(as.data.frame(t(abund_foralpha))))
							abund_foralpha <- tmp$norm(method = transformation) %>% t %>% as.data.frame
						}
						if(method == "KW_dunn"){
							# filter taxa with 1 across all samples
							abund_foralpha %<>% .[apply(., 1, function(x){length(unique(x)) != 1}), ]
						}
						if(method %in% c("betareg", "glmm_beta")){
							if(any(abund_foralpha) > 1){
								abund_foralpha %<>% {./(max(.) + beta_pseudo)}
							}
							if(any(abund_foralpha == 0)){
								abund_foralpha[abund_foralpha == 0] <- beta_pseudo
							}
							if(any(abund_foralpha == 1)){
								abund_foralpha[abund_foralpha == 1] <- 1/(1 + beta_pseudo)
							}
						}
						tem_data <- clone(tmp_dataset)
						tem_data$alpha_diversity <- as.data.frame(t(abund_foralpha))
						tem_data1 <- suppressMessages(trans_alpha$new(dataset = tem_data, group = group, by_group = by_group, by_ID = by_ID))
						tem_data1$cal_diff(method = method, p_adjust_method = p_adjust_method, alpha = alpha, ...)
						output <- tem_data1$res_diff
						if(!is.null(tem_data1$res_model)){
							self$res_model <- tem_data1$res_model
							message("The original ", method, " models list is stored in object$res_model ...")
						}
						colnames(output)[colnames(output) == "Measure"] <- "Taxa"
						method <- tem_data1$cal_diff_method
						if("Letter" %in% colnames(output)){
							output <- cbind.data.frame(output, Significance = output$Letter)
						}
					}
				}
				
				if(method %in% c("lefse", "rf")){
					group_vec <- sampleinfo[, group] %>% as.factor
					comparisions <- paste0(levels(group_vec), collapse = " - ")
					message("Start Kruskal-Wallis rank sum test for ", group, " ...")
					res_class <- suppressWarnings(lapply(seq_len(nrow(abund_table)), function(x) private$test_mark(abund_table[x, ], group_vec, method = "kruskal.test")))
					
					pvalue_raw <- unlist(lapply(res_class, function(x) x$p_value))
					names(pvalue_raw) <- rownames(abund_table)
					pvalue_raw[is.nan(pvalue_raw)] <- 1
					message(sum(pvalue_raw < alpha), " taxa found significant ...")

					pvalue <- p.adjust(pvalue_raw, method = p_adjust_method)
					sel_taxa <- pvalue < alpha
					message("After P value adjustment, ", sum(sel_taxa), " taxa found significant ...")
					private$check_taxa_number(sel_taxa, p_adjust_method)
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
					output$Significance <- generate_p_siglabel(output$P.adj, nonsig = "ns")
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
					for(num in seq_len(boots)){
						res_lda_pair <- list()
						sample_names_resample <- colnames(abund_table_sub)[base::sample(1:ncol(abund_table_sub), size = ceiling(ncol(abund_table_sub) * nresam))]
						abund_table_sub_resample <- abund_table_sub[, sample_names_resample]
						sampleinfo_resample <- sampleinfo[sample_names_resample, , drop = FALSE]
						# make sure the groups and samples number available
						if(sum(table(as.character(sampleinfo_resample[, group])) > 1) < 2){
							res_lda[[num]] <- NA
							next
						}
						for(i in seq_len(ncol(all_class_pairs))){
							sel_samples <- sampleinfo_resample[, group] %in% all_class_pairs[, i]
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
								# consider subgroup as an independent variable
								abund1 <- cbind.data.frame(t(abund_table_sub_lda), Group = group_vec_lda, lefse_subgroup = subgroup_vec)
							}
							check_res <- tryCatch(mod1 <- MASS::lda(Group ~ ., abund1, tol = 1.0e-10), error = function(e) { skip_to_next <- TRUE})
							if(rlang::is_true(check_res)) {
								res_lda_pair[[i]] <- NA
								next
							}else{
								w <- mod1$scaling[,1]
								if(is.null(names(w))){
									names(w) <- rownames(mod1$scaling)
								}
								w_unit <- w/sqrt(sum(w^2))
								w_unit %<>% {.[!grepl("lefse_subgroup", names(.))]}
								ss <- abund1[, !colnames(abund1) %in% c("Group", "lefse_subgroup")]
								xy_matrix <- as.matrix(ss)
								LD <- xy_matrix %*% w_unit
								effect_size <- tapply(LD, group_vec_lda, mean) %>% as.vector %>% {.[1] - .[2]} %>% abs
								coeff <- abs(w_unit * effect_size)
								coeff[is.nan(coeff)] <- 0
								names(coeff) %<>% gsub("^`|`$", "", .)
								rres <- mod1$means %>% as.data.frame
								colnames(rres) %<>% gsub("^`|`$", "", .)
								rres <- rres[, rownames(abund_table_sub_lda), drop = FALSE]
								rres1 <- apply(rres, 2, function(x) abs(x[1] - x[2]))
								res_lda_pair[[i]] <- (rres1 + coeff[names(rres1)]) *0.5
							}
						}
						res_lda[[num]] <- res_lda_pair
					}
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
					output$Significance <- generate_p_siglabel(output$P.adj, nonsig = "ns")
				}
				
				if(method %in% c("metastat", "metagenomeSeq", "ALDEx2_t", "edgeR")){
					tmp_dataset$sample_table %<>% dropallfactors
					private$check_taxa_level_all(taxa_level)
					if(is.null(group_choose_paired)){
						all_name <- combn(unique(as.character(sampleinfo[, group])), 2)
					}else{
						all_name <- combn(unique(group_choose_paired), 2)
					}
					message("Total ", ncol(all_name), " paired group for test ...")
					output <- data.frame()
				}
				if(method == "metastat"){
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

					for(i in 1:ncol(all_name)) {
						message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ..."))
						use_data <- new_abund[ , unlist(lapply(as.character(all_name[,i]), function(x) which(as.character(sampleinfo[, group]) %in% x)))]
						use_data %<>% .[!grepl("__$", rownames(.)), ]
						use_data <- use_data[apply(use_data, 1, sum) != 0, ]
						g <- sum(as.character(sampleinfo[, group]) == as.character(all_name[1, i])) + 1
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
					output$Significance <- generate_p_siglabel(output$qvalue, nonsig = "ns")
					output <- data.frame(output, Group = max_group)
				}
				if(method == "metagenomeSeq"){
					if(!require("metagenomeSeq")){
						stop("metagenomeSeq package not installed !")
					}
					for(i in 1:ncol(all_name)) {
						message(paste0("Run ", i, " : ", paste0(as.character(all_name[, i]), collapse = " - "), " ...\n"))
						use_dataset <- clone(tmp_dataset)
						use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[, i]), , drop = FALSE]
						newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
						obj <- newMRexperiment(
							newdata$otu_table, 
							phenoData= AnnotatedDataFrame(newdata$sample_table)
	#						featureData = AnnotatedDataFrame(use_dataset$tax_table)
							)
						## Normalization and Statistical testing
						obj_1 <- cumNorm(obj)
						pd <- pData(obj)
						if(group != "Group"){
							if("Group" %in% colnames(pd)){
								pd <- pd[, colnames(pd) != "Group", drop = FALSE]
							}
							colnames(pd)[which(colnames(pd) == group)] <- "Group"
						}
						mod <- model.matrix(~1 + Group, data = pd)
						objres1 <- fitFeatureModel(obj_1, mod)
						tb <- data.frame(logFC = objres1@fitZeroLogNormal$logFC, se = objres1@fitZeroLogNormal$se)
						p <- objres1@pvalues
						if(p_adjust_method == "ihw-ubiquity" | p_adjust_method == "ihw-abundance"){
							padj <- MRihw(objres1, p, p_adjust_method, 0.1)
						}else{
							padj <- p.adjust(p, method = p_adjust_method)
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
				if(method == "DESeq2"){
					if(!require("DESeq2")){
						stop("DESeq2 package is not installed !")
					}
					use_dataset <- clone(tmp_dataset)
					newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
					if(is.null(group)){
						all_parameters <- c(as.list(environment()), list(...))
						if(! "formula" %in% names(all_parameters)){
							stop("Either group or formula parameter should be provided!")
						}
						group <- all_parameters[["formula"]]
					}
					if(!grepl("^~", group)){
						group <- paste0("~", group)
					}
					use_formula <- stats::as.formula(group)
					deseq_obj <- DESeqDataSetFromMatrix(
						countData = newdata$otu_table,
						colData = newdata$sample_table,
						design = use_formula
						)
					res_deseq <- DESeq(deseq_obj, ...)
					self$res_diff_raw <- res_deseq
					message('Original result is stored in object$res_diff_raw ...')
					message('Merge the output tables ...')
					group %<>% gsub("^~", "", .)
					if(grepl("+", group, fixed = TRUE)){
						method <- paste0("DESeq formula: ", group)
					}
					all_groups <- strsplit(group, split = "+", fixed = TRUE) %>% unlist %>% gsub("^\\s+|\\s+$", "", .)
					output <- data.frame()
					res_names <- resultsNames(res_deseq) %>% .[. != "Intercept"]
					for(i in res_names){
						res <- results(res_deseq, name = i)
						res <- data.frame(res)
						if(grepl("_vs_", i, fixed = TRUE)){
							num <- lapply(all_groups, function(j){pmatch(j, i)}) %>% unlist %>% {!is.na(.)} %>% which
							i <- gsub(paste0(all_groups[num], "_"), "", i, fixed = TRUE)
							res$Comparison <- gsub("_vs_", " - ", i, fixed = TRUE)
						}else{
							res$Comparison <- i
						}
						res$Taxa <- rownames(res)
						rownames(res) <- NULL
						output <- rbind.data.frame(output, res)
					}
					colnames(output)[colnames(output) == "padj"] <- "P.adj"
					colnames(output)[colnames(output) == "pvalue"] <- "P.unadj"					
					output <- output[, c("Comparison", "Taxa", colnames(output)[1:(ncol(output) - 2)])]
				}
				if(method == "edgeR"){
					if(!require("edgeR")){
						stop("edgeR package is not installed!")
					}
					for(i in 1:ncol(all_name)) {
						message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ...\n"))
						use_dataset <- clone(tmp_dataset)
						use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[,i]), , drop = FALSE]
						newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
						y <- DGEList(counts = newdata$otu_table, group = newdata$sample_table[, group])
						z <- calcNormFactors(y)
						z <- estimateTagwiseDisp(estimateCommonDisp(z))
						res_raw <- exactTest(z)
						st <- topTags(res_raw, n = nrow(newdata$otu_table), adjust.method = p_adjust_method, sort.by = "PValue")
						res <- st@.Data[[1]]
						colnames(res)[colnames(res) == "PValue"] <- "P.unadj"
						colnames(res)[colnames(res) == "FDR"] <- "P.adj"
						res <- cbind.data.frame(feature = rownames(res), res)
						rownames(res) <- NULL
						add_name <- paste0(as.character(all_name[, i]), collapse = " - ") %>% rep(., nrow(res))
						res <- cbind.data.frame(compare = add_name, res)
						output <- rbind.data.frame(output, res)
					}
				}
				if(method %in% c("ALDEx2_t", "ALDEx2_kw")){
					if(!require("ALDEx2")){
						stop("ALDEx2 package is not installed !")
					}
				}
				if(method == "ALDEx2_t"){
					for(i in 1:ncol(all_name)) {
						message(paste0("Run ", i, " : ", paste0(as.character(all_name[,i]), collapse = " - "), " ...\n"))
						use_dataset <- clone(tmp_dataset)
						use_dataset$sample_table %<>% .[.[, group] %in% as.character(all_name[,i]), , drop = FALSE]
						newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
						res_raw <- ALDEx2::aldex(newdata$otu_table, newdata$sample_table[, group], test = "t", ...)
						res <- cbind.data.frame(feature = rownames(res_raw), res_raw)
						rownames(res) <- NULL
						add_name <- paste0(as.character(all_name[, i]), collapse = " - ") %>% rep(., nrow(res))
						res <- cbind.data.frame(compare = add_name, res) %>%
							.[, !grepl("^rab\\.", colnames(.))]
						output <- rbind.data.frame(output, res)
					}
				}
				if(method == "ALDEx2_kw"){
					if(taxa_level == "all"){
						stop("Please provide the taxa_level instead of 'all', such as 'Genus' !")
					}
					use_dataset <- clone(tmp_dataset)
					use_dataset$sample_table %<>% dropallfactors
					newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
					res_raw <- ALDEx2::aldex(newdata$otu_table, newdata$sample_table[, group], test = "kw", ...)
					res <- cbind.data.frame(feature = rownames(res_raw), res_raw)
					comparisions <- paste0(unique(as.character(use_dataset$sample_table[, group])), collapse = " - ")
					output <- cbind.data.frame(compare = comparisions, res)
				}
				if(method %in% c("ALDEx2_t", "ALDEx2_kw")){
					if(! any(colnames(output) %in% ALDEx2_sig)){
						stop("ALDEx2_sig is not found in the result! Please check the input!")
					}else{
						output$P.adj <- output[, ALDEx2_sig %>% .[. %in% colnames(output)] %>% .[1]]
					}
				}
				if(method == "ancombc2"){
					if(!require("ANCOMBC")){
						stop("ANCOMBC package is not installed !")
					}
					if(!require("file2meco")){
						stop("Please install file2meco package! The function meco2phyloseq is required!")
					}
					use_dataset <- clone(tmp_dataset)
					newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
					newdata <- file2meco::meco2phyloseq(newdata)
					# check fix_formula parameter
					all_parameters <- c(as.list(environment()), list(...))
					if("fix_formula" %in% names(all_parameters)){
						res_raw <- ANCOMBC::ancombc2(newdata, tax_level = NULL, assay_name = "counts", group = group, ...)
					}else{
						res_raw <- ANCOMBC::ancombc2(newdata, tax_level = NULL, assay_name = "counts", group = group, fix_formula = group, ...)
					}
					self$res_diff_raw <- res_raw
					tmp <- res_raw$res
					if(is.null(tmp)){
						stop("The res in the ancombc2 results is NULL!")
					}
					message('Converting res to long format ...')
					res_convert <- data.frame()
					for(i in seq_len(nrow(tmp))){
						taxon_data <- tmp[i, ]
						raw_colnames <- colnames(taxon_data)
						all_factors <- raw_colnames %>% .[grepl("lfc_", .)] %>% gsub("lfc_", "", .)
						res_convert <- rbind(res_convert, data.frame(Taxa = tmp[i, "taxon"], 
							Factors = all_factors, 
							lfc = taxon_data %>% .[, grepl("^lfc_", raw_colnames)] %>% unlist,
							se = taxon_data %>% .[, grepl("^se_", raw_colnames)] %>% unlist,
							W = taxon_data %>% .[, grepl("^W_", raw_colnames)] %>% unlist,
							p = taxon_data %>% .[, grepl("^p_", raw_colnames)] %>% unlist,
							P.adj = taxon_data %>% .[, grepl("^q_", raw_colnames)] %>% unlist,
							diff = taxon_data %>% .[, grepl("^diff_", raw_colnames)] %>% unlist,
							passed_ss = taxon_data %>% .[, grepl("^passed_ss_", raw_colnames)] %>% unlist
							)
						)
					}
					output <- res_convert
				}
				if(method == "linda"){
					if(!require("MicrobiomeStat")){
						stop("MicrobiomeStat package is not installed!")
					}
					private$check_taxa_level_all(taxa_level)
					use_dataset <- clone(tmp_dataset)
					newdata <- private$generate_microtable_unrel(use_dataset, taxa_level, filter_thres, filter_features)
					if(is.null(group)){
						all_parameters <- c(as.list(environment()), list(...))
						if(! "formula" %in% names(all_parameters)){
							stop("Either group or formula parameter should be provided!")
						}
						group <- all_parameters[["formula"]]
						if(!grepl("^~", group)){
							stop("The input formula parameter should start with ~! Please read the document of formula parameter of MicrobiomeStat::linda function!")
						}
						res <- MicrobiomeStat::linda(as.matrix(newdata$otu_table), newdata$sample_table, feature.dat.type = 'count', ...)
					}else{
						if(!grepl("^~", group)){
							group <- paste0("~", group)
						}
						res <- MicrobiomeStat::linda(as.matrix(newdata$otu_table), newdata$sample_table, formula = group, feature.dat.type = 'count', ...)
					}
					self$res_diff_raw <- res
					message('Original result is stored in object$res_diff_raw ...')
					message('Merge the output tables ...')
					group %<>% gsub("^~", "", .)
					# different cases
					output <- data.frame()
					if(grepl("+", group, fixed = TRUE)){
						# multi-factor
						for(i in names(res$output)){
							tmp <- res$output[[i]]
							tmp <- data.frame(Taxa = rownames(tmp), Factors = i, tmp)
							output %<>% rbind(., tmp)
						}
						method <- paste0("linda formula: ", group)
					}else{
						# normal case
						if(group %in% colnames(newdata$sample_table)){
							all_groups <- newdata$sample_table[, group] %>% as.character %>% unique
							for(i in names(res$output)){
								tmp <- res$output[[i]]
								tmp_group_element1 <- gsub(group, "", i, fixed = TRUE)
								tmp_group_element2 <- paste0(group, all_groups) %>% .[!. %in% names(res$output)] %>% gsub(group, "", .)
								tmp <- data.frame(compare = paste0(tmp_group_element1, " - ", tmp_group_element2), Taxa = rownames(tmp), tmp)
								output %<>% rbind(., tmp)
							}
						}else{
							stop('The group is not found in the column names of sample_table! Skip the merging step!')
						}
					}
					colnames(output)[colnames(output) %in% c("pvalue", "padj")] <- c("P.unadj", "P.adj")
				}
				if(method == "maaslin2"){
					tmp_trans_env <- trans_env$new(dataset = tmp_dataset, env_cols = 1:ncol(tmp_dataset$sample_table))
					tmp_trans_env$cal_cor(use_data = taxa_level, cor_method = method, filter_thres = filter_thres,
						plot_heatmap = FALSE, plot_scatter = FALSE, ...)
					output <- tmp_trans_env$res_cor
					self$res_trans_env <- tmp_trans_env
					message('Raw trans_env object for maaslin2 is stored in object$res_trans_env ...')
				}
				
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
				# further calculate mean and sd of res_abund with group parameter
				if(!is.null(group)){
					if(group %in% colnames(sampleinfo)){
						colnames(res_abund) <- c("Taxa", "Sample", "Abund")
						res_abund <- suppressWarnings(dplyr::left_join(res_abund, rownames_to_column(sampleinfo), by = c("Sample" = "rowname")))
						res_abund <- microeco:::summarySE_inter(res_abund, measurevar = "Abund", groupvars = c("Taxa", group))
						colnames(res_abund)[colnames(res_abund) == group] <- "Group"
					}
				}
				if(method %in% c("metagenomeSeq", "ALDEx2_t", "ALDEx2_kw", "DESeq2", "edgeR", "linda")){
					output %<>% dropallfactors(unfac2num = TRUE)
					colnames(output)[1:2] <- c("Comparison", "Taxa")
					if(group %in% colnames(sampleinfo)){
						# filter the unknown taxa in output
						output %<>% .[.$Taxa %in% res_abund$Taxa, ]
						output$Group <- lapply(seq_along(output$Taxa), function(x){
							select_group_split <- strsplit(output[x, "Comparison"], split = " - ") %>% unlist
							res_abund[res_abund$Taxa == output[x, "Taxa"] & res_abund$Group %in% select_group_split, ] %>%
							{.[which.max(.$Mean), "Group"]}
						}) %>% unlist
					}
				}
				self$res_abund <- res_abund
				message('Taxa abundance table is stored in object$res_abund ...')
				if("Factors" %in% colnames(output)){
					output[, "Factors"] %<>% gsub("\\s+$", "", .)
				}
				if(!is.null(output)){
					if("P.adj" %in% colnames(output)){
						if(!"Significance" %in% colnames(output)){
							output$Significance <- generate_p_siglabel(output$P.adj, nonsig = "")
						}
					}
				}
				self$res_diff <- output
				if(method == "ancombc2"){
					message("Original ancombc2 results are stored in object$res_diff_raw ...")
					if(!is.null(output)){
						message("Converted test result is stored in object$res_diff ...")
					}
					method <- "ancombc2_formula"
				}else{
					message(method , " analysis result is stored in object$res_diff ...")
				}
				self$method <- method
				self$taxa_level <- taxa_level
				# save abund_table for the cladogram
				self$abund_table <- abund_table
			}
		},
		#' @description
		#' Plot the abundance of differential taxa
		#'
		#' @param use_number default 1:20; numeric vector; the taxa numbers (1:n) selected in the plot; 
		#'   If the n is larger than the number of total significant taxa, automatically use all the taxa.		
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette.
		#' @param select_group default NULL; this is used to select the paired groups. 
		#'   This parameter is especially useful when the comparision methods is applied to paired groups;
		#'   The input select_group must be one of \code{object$res_diff$Comparison}.
		#' @param select_taxa default NULL; character vector to provide taxa names. 
		#' 	 The taxa names should be same with the names shown in the plot, not the 'Taxa' column names in \code{object$res_diff$Taxa}.
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
		#' @param text_y_size default 10; the size for the y axis text, i.e. feature text.
		#' @param coord_flip default TRUE; whether flip cartesian coordinates so that horizontal becomes vertical, and vertical becomes horizontal.
		#' @param xtext_angle default 45; number ranging from 0 to 90; used to make x axis text generate angle to reduce text overlap; 
		#' 	  only available when coord_flip = FALSE.
		#' @param ... parameters passed to \code{ggsignif::stat_signif} when add_sig = TRUE.
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
			xtext_angle = 45,
			...
			){
			abund_data <- self$res_abund
			method <- self$method
			diff_data <- self$res_diff
			if(grepl("ancombc2", method)){
				stop("The function can not be applied to ancombc2!")
			}
			if(grepl("formula", method)){
				stop("The function can not be applied to multi-factor analysis!")
			}
			# first determine how to select compared groups
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
			# sort according to different columns
			if(method == "metastat"){
				message('Reorder taxa according to qvalue in res_diff from low to high ...')
				diff_data %<>% .[order(.$qvalue, decreasing = FALSE), ]
				# diff_data %<>% .[.$qvalue < 0.05, ]
			}else{
				# lefse and rf are ordered
				if(! method %in% c("lefse", "rf", "anova")){
					if("P.adj" %in% colnames(diff_data)){
						message('Reorder taxa according to P.adj in res_diff from low to high ...')
						diff_data %<>% .[order(.$P.adj, decreasing = FALSE), ]
					}
				}
			}
			if(nrow(diff_data) == 0){
				stop("No significant taxa can be used to plot the abundance!")
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
				message('Use provided select_taxa to filter and reorder taxa ...')
				diff_data %<>% .[.$Taxa %in% select_taxa, ]
				if(nrow(diff_data) == 0){
					stop("No significant taxa can be used to plot the abundance!")
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
			
			p <- ggplot(abund_data, aes(x = Taxa, y = Mean, color = Group, fill = Group)) +
				theme_bw() +
				geom_bar(stat="identity", position = position_dodge(), width = barwidth)
			if(use_se == T){
				p <- p + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.45, position=position_dodge(barwidth), color = "black")
			}else{
				p <- p + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.45, position=position_dodge(barwidth), color = "black")
			}
			
			if(add_sig){
				# assign labels by factor orders
				x_axis_order <- levels(abund_data$Group)
				if(! add_sig_label %in% colnames(diff_data)){
					stop("The add_sig_label parameter must be one of colnames of object$res_diff!")
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
				
				# for groups > 3, show the global comparision
				if((length(levels(abund_data$Group)) > 2 & method %in% c("lefse", "rf", "KW", "ALDEx2_kw")) | 
					(length(unlist(gregexpr(" - ", diff_data$Comparison[1]))) > 1 & method == "ancombc2")){
					add_letter_text <- diff_data[match(all_taxa, diff_data$Taxa), add_sig_label]
					textdf <- data.frame(
						x = all_taxa, 
						y = y_start_use, 
						add = add_letter_text, 
						stringsAsFactors = FALSE
						)
					p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
				}else{
					if(! "Letter" %in% colnames(diff_data)){
						if(any(grepl("\\s-\\s", x_axis_order))){
							stop("The group names have ' - ' characters, which can hinder the group recognition and mapping in the plot! Please rename groups and rerun!")
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
						x_mid <- c()
						annotations <- c()
						y_position <- c()

						start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
						increase_bar_mid <- barwidth/length(x_axis_order)

						for(j in all_taxa){
							select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Taxa == j, ]
							for(i in seq_len(nrow(select_use_diff_data))){
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
						p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
					}
				}
			}

			p <- p + scale_color_manual(values = color_values) +
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
				p <- p + guides(fill = guide_legend(reverse = FALSE, ncol = 1), color = "none")
				if(xtext_angle != 0){
					p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", vjust = 1, hjust = 1, size = text_y_size))
				}else{
					p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", size = text_y_size))
				}
				p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = text_y_size, color = "black")) +
					theme(plot.margin = unit(c(.1, .1, .1, 1), "cm"))
			}
			p
		},
		#' @description
		#' Bar plot for indicator index, such as LDA score and P value.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for different groups.
		#' @param color_group_map default FALSE; whether match the colors to groups in order to fix the color in each group when part of groups are not shown in the plot.
		#'    When \code{color_group_map = TRUE}, the group_order inside the object will be used as full groups set to guide the color extraction.
		#' @param use_number default 1:10; numeric vector; the taxa numbers used in the plot, i.e. 1:n.
		#' @param threshold default NULL; threshold value of indicators for selecting taxa, such as 3 for LDA score of LEfSe.
		#' @param select_group default NULL; this is used to select the paired group when multiple comparisions are generated;
		#'   The input select_group must be one of \code{object$res_diff$Comparison}.
		#' @param keep_full_name default FALSE; whether keep the taxonomic full lineage names.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix, such as "g__".
		#' @param group_order default NULL; a vector to order the legend and colors in plot; 
		#' 	  If NULL, the function can first determine whether the group column of \code{microtable$sample_table} is factor. If yes, use the levels in it.
		#' 	  If provided, this parameter can overwrite the levels in the group of \code{microtable$sample_table}.
		#' @param axis_text_y default 12; the size for the y axis text.
		#' @param coord_flip default TRUE; whether flip cartesian coordinates so that horizontal becomes vertical, and vertical becomes horizontal.
		#' @param xtext_angle default 45; number ranging from 0 to 90; used to make x axis text generate angle to reduce text overlap; 
		#' 	  only available when coord_flip = FALSE.
		#' @param xtext_size default 10; the text size of x axis.
		#' @param heatmap_cell default "P.unadj"; the column of data for the cell of heatmap when formula with multiple factors is found in the method.
		#' @param heatmap_sig default "Significance"; the column of data for the significance label of heatmap.
		#' @param heatmap_x default "Factors"; the column of data for the x axis of heatmap.
		#' @param heatmap_y default "Taxa"; the column of data for the y axis of heatmap.
		#' @param heatmap_lab_fill default "P value"; legend title of heatmap.
		#' @param ... parameters passing to \code{\link{geom_bar}} for the bar plot or 
		#' 	  \code{plot_cor} function in \code{\link{trans_env}} class for the heatmap of multiple factors when formula is found in the method.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_diff_bar(use_number = 1:20)
		#' }
		plot_diff_bar = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			color_group_map = FALSE,
			use_number = 1:10,
			threshold = NULL,
			select_group = NULL,
			keep_full_name = FALSE,
			keep_prefix = TRUE,
			group_order = NULL,
			axis_text_y = 12,
			coord_flip = TRUE,
			xtext_angle = 45,
			xtext_size = 10,
			heatmap_cell = "P.unadj",
			heatmap_sig = "Significance",
			heatmap_x = "Factors",
			heatmap_y = "Taxa",
			heatmap_lab_fill = "P value",
			...
			){
			use_data <- self$res_diff
			method <- self$method
			if(is.null(method)){
				if(is.null(use_data)){
					stop("The res_diff should be provided when method is NULL!")
				}else{
					if(!"Value" %in% colnames(use_data)){
						stop("No Value column found in the object$res_diff!")
					}
				}
			}
			if(keep_full_name == F){
				if(any(grepl("\\..__", use_data$Taxa))){
					use_data$Taxa %<>% gsub(".*(.__.*?$)", "\\1", .)
				}else{
					use_data$Taxa %<>% gsub(".*\\|", "", .)
				}
			}
			if(keep_prefix == F){
				use_data$Taxa %<>% gsub(".__", "", .)
			}

			if((! grepl("formula", method)) & (! method %in% c("maaslin2"))){
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
				if("Group" %in% colnames(use_data)){
					if(is.null(group_order)){
						if((!is.null(self$group_order)) & (length(unique(use_data$Group)) == length(self$group_order))){
							use_data$Group %<>% factor(., levels = self$group_order)
						}else{
							use_data$Group %<>% as.character %>% as.factor
						}
					}else{
						use_data$Group %<>% factor(., levels = group_order)
					}
					if(color_group_map){
						# fix colors for each group
						all_groups <- self$group_order
						# make sure colors length enough for selection
						color_values <- expand_colors(color_values, length(all_groups))
						use_groups <- levels(use_data$Group)
						color_values %<>% .[match(use_groups, all_groups)]
					}else{
						color_values <- expand_colors(color_values, length(levels(use_data$Group)))
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
				}else{
					use_data %<>% .[order(.$Value, decreasing = TRUE), ]
					if(coord_flip){
						use_data$Taxa %<>% factor(., levels = rev(.))
					}else{
						use_data$Taxa %<>% factor(., levels = .)
					}
					ylab_title <- "Value"
				}
				self$plot_diff_bar_taxa <- levels(use_data$Taxa) %>% rev
				
				if("Group" %in% colnames(use_data)){
					p <- ggplot(use_data, aes(x = Taxa, y = Value, color = Group, fill = Group, group = Group)) +
						scale_color_manual(values = color_values) +
						scale_fill_manual(values = color_values)
				}else{
					p <- ggplot(use_data, aes(x = Taxa, y = Value))
				}
				p <- p +
					geom_bar(stat = "identity", position = position_dodge(), ...) +
					theme_bw() +
					ylab(ylab_title) +
					xlab("") +
					theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = axis_text_y, color = "black")) +
					theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()) +
					theme(panel.border = element_blank()) +
					theme(axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square"))
				
				if(coord_flip){
					p <- p + coord_flip()
				}else{
					p <- p + ggplot_xtext_anglesize(xtext_angle = xtext_angle, xtext_size = xtext_size)
				}
				p
			}else{
				# heatmap for multi-factor
				message("Perform heatmap instead of bar plot as formula is found ...")
				tmp <- use_data
				tmp_trans_env <- convert_diff2transenv(tmp, heatmap_x, heatmap_y, heatmap_cell, heatmap_sig, heatmap_lab_fill)
				tmp_trans_env$plot_cor(keep_full_name = keep_full_name, keep_prefix = keep_prefix, ...)
			}
		},
		#' @description
		#' Plot the cladogram using taxa with significant difference.
		#'
		#' @param color default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color palette used in the plot.
		#' @param group_order default NULL; a vector to order the legend in plot; 
		#' 	  If NULL, the function can first check whether the group column of sample_table is factor. If yes, use the levels in it.
		#' 	  If provided, this parameter can overwrite the levels in the group of sample_table. 
		#' 	  If the number of provided group_order is less than the number of groups in \code{res_diff$Group}, the function will select the groups of group_order automatically.
		#' @param use_taxa_num default 200; integer; The taxa number used in the background tree plot; select the taxa according to the mean abundance .
		#' @param filter_taxa default NULL; The mean relative abundance used to filter the taxa with low abundance.
		#' @param use_feature_num default NULL; integer; The feature number used in the plot; 
		#'	  select the features according to the LDA score (method = "lefse") or MeanDecreaseGini (method = "rf") from high to low.
		#' @param clade_label_level default 4; the taxonomic level for marking the label with letters, root is the largest.
		#' @param select_show_labels default NULL; character vector; The features to show in the plot with full label names, not the letters.
		#' @param only_select_show default FALSE; whether only use the the select features in the parameter \code{select_show_labels}.
		#' @param sep default "|"; the seperate character in the taxonomic information.
		#' @param branch_size default 0.2; numberic, size of branch.
		#' @param alpha default 0.2; shading of the color.
		#' @param clade_label_size default 2; basic size for the clade label; please also see \code{clade_label_size_add} and \code{clade_label_size_log}.
		#' @param clade_label_size_add default 5; added basic size for the clade label; see the formula in \code{clade_label_size_log} parameter.
		#' @param clade_label_size_log default \code{exp(1)}; the base of \code{log} function for added size of the clade label; the size formula: 
		#'   \code{clade_label_size + log(clade_label_level + clade_label_size_add, base = clade_label_size_log)}; 
		#'   so use \code{clade_label_size_log}, \code{clade_label_size_add} and \code{clade_label_size}
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
				stop("This function is available only when taxa_level = 'all' !")
			}
			if(!is.null(use_feature_num)){
				if(use_feature_num > nrow(marker_table)){
					message("Input use_feature_num ", use_feature_num, " larger than available features number ", nrow(marker_table), " ! ", 
						"Use ", nrow(marker_table), " instead of it ...")
				}else{
					message("Select ", use_feature_num, " significant features ...")
					marker_table %<>% .[1:use_feature_num, ]
				}
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
			if(! all(color_groups %in% unique(marker_table$Group))){
				tmp_message <- color_groups[! color_groups %in% unique(marker_table$Group)]
				message("Part of groups in group_order, ", paste(tmp_message, collapse = " "), ", not found in the feature table ...")
				color_groups %<>% .[. %in% unique(marker_table$Group)]
			}
			if(! all(unique(marker_table$Group) %in% color_groups)){
				tmp_message <- unique(marker_table$Group)[! unique(marker_table$Group) %in% color_groups]
				message("Part of groups in the feature table, ", paste(tmp_message, collapse = " "), ", not found in the group_order ...")
				marker_table %<>% .[.$Group %in% color_groups, ]
			}
			# get the color palette
			if(length(color) < length(unique(marker_table$Group))){
				stop("Please provide enough color palette! There are ", length(unique(marker_table$Group)), 
					" groups, but only ", length(color), " colors provideed in color parameter!")
			}else{
				color <- color[1:length(unique(marker_table$Group))]
			}

			tree <- private$plot_backgroud_tree(abund_table = abund_table, use_taxa_num = use_taxa_num, filter_taxa = filter_taxa, sep = sep)
			
			# generate annotation
			annotation <- private$generate_cladogram_annotation(marker_table, tree = tree, color = color, color_groups = color_groups, sep = sep)
			# check again filtered groups
			if(!is.null(use_taxa_num)){
				if(! all(color_groups %in% unique(annotation$enrich_group))){
					tmp_message <- color_groups[! color_groups %in% unique(annotation$enrich_group)]
					message("Biomarkers in group(s), ", paste(tmp_message, collapse = " "), ", not found in the background tree!", 
						" Please try to enlarge parameter use_taxa_num from ", use_taxa_num, " to a larger one ...")
					color_groups %<>% .[. %in% unique(annotation$enrich_group)]
				}
			}

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
				geom_rect(aes(xmin = x, xmax = x, ymax = y, ymin = y, fill = enrich_group), data = hilights_df, inherit.aes = FALSE) +
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
		check_taxa_level_all = function(taxa_level){
			if(taxa_level == "all"){
				stop("The taxa_level parameter cannot be 'all' for this method! Please provide a taxonomic level, such as 'Genus' !")
			}
		},
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
		check_taxa_number = function(sel_taxa, p_adjust_method){
			if(sum(sel_taxa) == 0){
				if(p_adjust_method == "none"){
					stop('No significant feature found!')
				}else{
					stop('No significant feature found! To disable p value adjustment, please use p_adjust_method = "none"!')
				}
			}
			if(sum(sel_taxa) == 1){
				if(p_adjust_method == "none"){
					stop('Only one significant feature found! Stop running subsequent process!')
				}else{
					stop('Only one significant feature found! Stop running subsequent process! To disable p value adjustment, please use p_adjust_method = "none"!')
				}
			}
		},
		generate_microtable_unrel = function(use_dataset, taxa_level, filter_thres, filter_features){
			use_dataset$tidy_dataset()
			suppressMessages(use_dataset$cal_abund(rel = FALSE))
			use_feature_table <- use_dataset$taxa_abund[[taxa_level]]
			if(filter_thres > 0){
				use_feature_table %<>% .[! rownames(.) %in% names(filter_features), ]
			}
			message("Available feature number: ", nrow(use_feature_table))
			newdata <- microtable$new(otu_table = use_feature_table, sample_table = use_dataset$sample_table)
			newdata$tidy_dataset()
			newdata
		},
		# plot the background tree according to raw abundance table
		plot_backgroud_tree = function(abund_table, use_taxa_num = NULL, filter_taxa = NULL, sep = "|"){
			# filter the taxa with unidentified classification or with space, in case of the unexpected error in the following operations
			abund_table %<>% {.[!grepl("\\|.__\\|", rownames(.)), ]} %>%
				{.[!grepl("\\s", rownames(.)), ]} %>%
				# also filter uncleared classification to make it in line with the lefse above
				{.[!grepl("Incertae_sedis|unculture", rownames(.), ignore.case = TRUE), ]}
			if(nrow(abund_table) <= 2){
				stop("After filtering out non-standard taxonomy information, the abundance table only has ", nrow(abund_table), " feature(s)! ", 
					"Is there an issue with the taxonomy table? ", 
					"Please first use the tidy_taxonomy function to process the taxonomy information table before constructing the microtable object.")
			}
			if(!is.null(use_taxa_num)){
				if(use_taxa_num < nrow(abund_table)){
					message("Select ", use_taxa_num, " most abundant taxa as the background cladogram ...")
					abund_table %<>% .[names(sort(apply(., 1, mean), decreasing = TRUE)[1:use_taxa_num]), ]
				}else{
					message("Provided use_taxa_num: ", use_taxa_num, " >= ", " total effective taxa number. Skip the selection ...")
				}
			}
			if(!is.null(filter_taxa)){
				abund_table %<>% .[apply(., 1, mean) > (self$lefse_norm * filter_taxa), ]
			}
			abund_table %<>% .[sort(rownames(.)), ]
			tree_table <- data.frame(taxa = row.names(abund_table), abd = rowMeans(abund_table), stringsAsFactors = FALSE) %>%
				dplyr::mutate(taxa = paste("r__Root", .data$taxa, sep = sep), abd = .data$abd/max(.data$abd)*100)
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
			# levels used for extending clade label
			label_levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
				factor(levels = rev(unlist(lapply(taxa_split, function(x) gsub("(.)__.*", "\\1", x))) %>% .[!duplicated(.)]))

			# root must be a parent node
			nodes_parent <- purrr::map_chr(taxa_split, ~ .x[length(.x) - 1]) %>% c("root", .)

			# add index for nodes
			is_tip <- !nodes %in% nodes_parent
			index <- vector("integer", length(is_tip))
			index[is_tip] <- 1:sum(is_tip)
			index[!is_tip] <- (sum(is_tip) + 1):length(is_tip)

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
			tree
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
		get_offset = function(x){(x*0.2+0.2)^2},
		# dependent functions in utility
		calculate_metastat = function(inputdata, g, pflag = FALSE, threshold = NULL, B = NULL){
			trans_data <- load_frequency_matrix(input = inputdata)
			res <- detect_differentially_abundant_features(jobj = trans_data, g = g, pflag = pflag, threshold = threshold, B = B)
			res
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
