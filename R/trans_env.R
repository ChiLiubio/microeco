#' @title 
#' Create \code{trans_env} object for analyzing the effects of environmental factors on communities.
#'
#' @description
#' This class is a wrapper for a series of operations associated with environmental measurements, including redundancy analysis, 
#' mantel test, correlation analysis and linear fitting.
#'
#' @export
trans_env <- R6Class(classname = "trans_env",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param env_cols default NULL; either numeric vector or character vector to select columns in \code{microtable$sample_table}, i.e. dataset$sample_table. 
		#'    This parameter should be used in the case that all the required environmental data is in \code{sample_table} of your \code{microtable} object.
		#'    Otherwise, please use \code{add_data} parameter.
		#' @param add_data default NULL; \code{data.frame} format; provide the environmental data in the format \code{data.frame}; rownames should be sample names.
		#'   This parameter should be used when the \code{microtable$sample_table} object does not have environmental data. 
		#'   Under this circumstance, the \code{env_cols} parameter can not be used because no data can be selected.
		#' @param character2numeric default TRUE; whether convert the characters or factors to numeric attributes.
		#' @param complete_na default FALSE; Whether fill the NA (missing value) in the environmental data;
		#'   If TRUE, the function can run the interpolation with the \code{mice} package.
		#' @return \code{data_env} stored in the object.
		#' @examples
		#' data(dataset)
		#' data(env_data_16S)
		#' t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
		initialize = function(
			dataset = NULL, 
			env_cols = NULL, 
			add_data = NULL, 
			character2numeric = TRUE, 
			complete_na = FALSE
			){
			# support all the dataset, env_cols and add_data = NULL from v0.7.0
			if(is.null(dataset)){
				message("The dataset not provided. Remember to provide additional data in the correponding function ...")
			}
			if(is.null(add_data)){
				if(!is.null(env_cols)){
					env_data <- dataset$sample_table[, env_cols, drop = FALSE]
				}else{
					env_data <- NULL
					message("Both env_cols and add_data are NULL. Remember to provide additional data in the correponding function ...")
				}
			}else{
				if(is.null(dataset)){
					env_data <- add_data
				}else{
					# first check rownames of add_data
					if(!any(rownames(add_data) %in% rownames(dataset$sample_table))){
						stop("No valid rowname in add_data! Rownames must be sample names in the input add_data!")
					}
					env_data <- add_data[rownames(add_data) %in% rownames(dataset$sample_table), , drop = FALSE]
				}
			}
			if(!is.null(dataset)){
				use_dataset <- clone(dataset)
				if(!is.null(env_data)){
					inter_sum <- sum(rownames(use_dataset$sample_table) %in% rownames(env_data))
					if(inter_sum < nrow(use_dataset$sample_table)){
						message(nrow(use_dataset$sample_table) - inter_sum, " sample(s) not found in environmental data ...")
						filter_names <- use_dataset$sample_table %>% rownames %>% .[!. %in% rownames(env_data)]
						message("Filter sample(s): ", paste(filter_names, collapse = " "), " ...")
						use_dataset$sample_table %<>% base::subset(rownames(.) %in% rownames(env_data))
						use_dataset$tidy_dataset(main_data = FALSE)
						message("The pruned microtable is stored in object$dataset ...")
					}
					env_data %<>% .[rownames(use_dataset$sample_table), , drop = FALSE]
				}
				self$dataset <- use_dataset
			}else{
				self$dataset <- NULL
			}
			if(!is.null(env_data)){
				if(complete_na == T){
					env_data[env_data == ""] <- NA
					env_data %<>% dropallfactors(., unfac2num = TRUE)
					env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
					env_data %<>% mice::mice(print = FALSE) %>% mice::complete(., 1)
				}
				if(character2numeric == T){
					env_data %<>% dropallfactors(., unfac2num = TRUE, char2num = TRUE)
				}
			}
			self$data_env <- env_data
			message("Env data is stored in object$data_env ...")
		},
		#' @description
		#' Test the difference of environmental variable across groups.
		#'
		#' @param group default NULL; a colname of sample_table used to compare values across groups.
		#' @param by_group default NULL; perform the differential test among groups (\code{group} parameter) within each group (\code{by_group} parameter).
		#' @param method default "KW"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'KW'}}{KW: Kruskal-Wallis Rank Sum Test for all groups (>= 2)}
		#'     \item{\strong{'KW_dunn'}}{Dunn's Kruskal-Wallis Multiple Comparisons, see \code{dunnTest} function in \code{FSA} package}
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'anova'}}{Duncan's new multiple range test for one-way anova; see \code{duncan.test} function of \code{agricolae} package.
		#'     	  For multi-factor anova, see \code{aov}}
		#'     \item{\strong{'scheirerRayHare'}}{Scheirer Ray Hare test for nonparametric test used for a two-way factorial experiment; 
		#'     	  see \code{scheirerRayHare} function of \code{rcompanion} package}
		#'     \item{\strong{'lme'}}{lme: Linear Mixed Effect Model based on the \code{lmerTest} package}
		#'   }
		#' @param ... parameters passed to \code{cal_diff} function of \code{\link{trans_alpha}} class.
		#' @return \code{res_diff} in object.
		#'   In the data frame, 'Group' column means that the group has the maximum median or mean value across the test groups;
		#'   For non-parametric methods, maximum median value; For t.test, maximum mean value.
		#' @examples
		#' \donttest{
		#' t1$cal_diff(group = "Group", method = "KW")
		#' t1$cal_diff(group = "Group", method = "KW_dunn")
		#' t1$cal_diff(group = "Group", method = "anova")
		#' }
		cal_diff = function(group = NULL, by_group = NULL, method = c("KW", "KW_dunn", "wilcox", "t.test", "anova", "scheirerRayHare", "lme")[1], ...){
			if(is.null(group) & ! method %in% c("anova", "scheirerRayHare")){
				stop("The group parameter is necessary for the method: ", method, "!")
			}
			if(!is.null(group)){
				if(!group %in% colnames(self$dataset$sample_table)){
					stop("Provided parameter group must be one colname of sample_table! Please check it!")
				}
			}
			if(is.null(self$data_env)){
				stop("The data_env is NULL! Please check the data input when creating the object !")
			}else{
				env_data <- self$data_env
			}
			tem_data <- clone(self$dataset)
			# use test method in trans_alpha
			tem_data$alpha_diversity <- env_data
			tem_data1 <- suppressMessages(trans_alpha$new(dataset = tem_data, group = group, by_group = by_group))
			suppressMessages(tem_data1$cal_diff(method = method, ...))
			self$res_diff <- tem_data1$res_diff
			self$res_diff_tmp <- tem_data1
			self$group <- group
			message('The result is stored in object$res_diff ...')
		},
		#' @description
		#' Plot environmental variables across groups and add the significance label.
		#'
		#' @param ... parameters passed to \code{plot_alpha} in \code{\link{trans_alpha}} class. 
		#' 	 Please see \code{plot_alpha} function of \code{\link{trans_alpha}} for all the available parameters.
		plot_diff = function(...){
			if(is.null(self$res_diff_tmp)){
				stop("Please first run cal_diff function!")
			}
			res_diff_tmp <- self$res_diff_tmp
			res_diff_tmp$plot_alpha(...)
		},
		#' @description
		#' Calculate the autocorrelations among environmental variables and plot the result.
		#'
		#' @param group default NULL; a colname of sample_table; used to perform calculations for different groups.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette.
		#' @param alpha default 0.8; the alpha value to add transparency in colors; useful when group is not NULL.
		#' @param ... default parameters passed to \code{GGally::ggpairs}.
		#' @return \code{ggmatrix}.
		#' @examples
		#' \dontrun{
		#' # Spearman correlation
		#' t1$cal_autocor(upper = list(continuous = GGally::wrap("cor", method= "spearman")))
		#' }
		cal_autocor = function(group = NULL, color_values = RColorBrewer::brewer.pal(8, "Dark2"), alpha = 0.8, ...){
			if(!requireNamespace("GGally", quietly = TRUE)){
				stop("Please first install GGally with the command: install.packages('GGally') !")
			}
			if(is.null(self$data_env)){
				stop("The data_env is NULL! Please check the data input when creating the object !")
			}else{
				env_data <- self$data_env
			}
			if(is.null(group)){
				g <- GGally::ggpairs(env_data, ...)
			}else{
				sample_table <- self$dataset$sample_table
				if(! group %in% colnames(sample_table)){
					stop("Please provide a correct group name!")
				}
				merge_data <- cbind.data.frame(sample_table[, group, drop = FALSE], env_data[rownames(sample_table), ])
				g <- GGally::ggpairs(merge_data, aes_meco(colour = group, alpha = alpha),  ...)
				# Loop through each plot changing relevant scales 
				for(i in 1:g$nrow){
					for(j in 1:g$ncol){
						g[i, j] <- g[i, j] + 
							scale_fill_manual(values = color_values) +
							scale_color_manual(values = color_values)
					}
				}
			}
			g
		},
		#' @description
		#' Redundancy analysis (RDA) and Correspondence Analysis (CCA) based on the \code{vegan} package.
		#'
		#' @param method default c("RDA", "dbRDA", "CCA")[1]; the ordination method.
		#' @param feature_sel default FALSE; whether perform the feature selection based on forward selection method.
		#' @param taxa_level default NULL; If use RDA or CCA, provide the taxonomic rank, such as "Phylum" or "Genus";
		#'   If use otu_table; please set \code{taxa_level = "OTU"}.
		#' @param taxa_filter_thres default NULL; The relative abundance threshold used to filter taxa when method is "RDA" or "CCA".
		#' @param use_measure default NULL; a name of beta diversity matrix; only useful when parameter \code{method = "dbRDA"};
		#' 	 If not provided, use the first beta diversity matrix in the \code{microtable$beta_diversity} automatically.
		#' @param add_matrix default NULL; additional distance matrix provided, when the user does not want to use the beta diversity matrix within the dataset;
		#'   only available when method = "dbRDA".
		#' @param ... paremeters passed to \code{dbrda}, \code{rda} or \code{cca} function according to the input of \code{method} parameter.
		#' @return \code{res_ordination} and \code{res_ordination_R2} in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination(method = "dbRDA", use_measure = "bray")
		#' t1$cal_ordination(method = "RDA", taxa_level = "Genus")
		#' t1$cal_ordination(method = "CCA", taxa_level = "Genus")
		#' }
		cal_ordination = function(
			method = c("RDA", "dbRDA", "CCA")[1],
			feature_sel = FALSE,
			taxa_level = NULL,
			taxa_filter_thres = NULL,
			use_measure = NULL,
			add_matrix = NULL,
			...
			){
			if(length(method) > 1){
				stop("The method parameter should have only one option!")
			}
			if(! method %in% c("RDA", "dbRDA", "CCA")){
				stop("The method should be one of 'RDA', 'dbRDA' and 'CCA' !")
			}
			if(is.null(self$data_env)){
				stop("The data_env is NULL! Please check the data input when creating the object !")
			}else{
				env_data <- self$data_env
			}
			
			if(method == "dbRDA"){
				# add_matrix has the priority
				if(!is.null(add_matrix)){
					use_matrix <- add_matrix
					message("Use the additional matrix provided for dbRDA ...")
				}else{
					if(is.null(self$dataset)){
						stop("No dataset is found in the object! please check the input!")
					}else{
						if(is.null(self$dataset$beta_diversity)){
							message("The beta_diversity in dataset is NULL; try to calculate it ...")
							self$dataset$cal_betadiv(unifrac = FALSE)
							message("Calculating done ...")
						}
						if(!is.null(use_measure)){
							if(! use_measure %in% names(self$dataset$beta_diversity)){
								stop("Please make sure use_measure: ", use_measure, " is in the dataset$beta_diversity!")
							}else{
								use_matrix <- self$dataset$beta_diversity[[use_measure]]
								message("Use ", use_measure, " in dataset$beta_diversity for dbRDA ...")
							}
						}else{
							use_matrix <- self$dataset$beta_diversity[[1]]
							message("Parameter use_measure not provided; use the first matrix in dataset$beta_diversity ...")
						}
					}
				}
				use_data <- use_matrix[rownames(env_data), rownames(env_data)] %>% as.dist
			}
			if(method %in% c("RDA", "CCA")){
				if(is.null(self$dataset)){
					stop("The dataset in the object is NULL! Please provide dataset when creating the object!")
				}
				if(is.null(taxa_level)){
					if("Genus" %in% colnames(self$dataset$tax_table)){
						taxa_level <- "Genus"
						message("No taxa_level provided; use Genus level automatically !")
					}else{
						taxa_level <- "OTU"
						message("No taxa_level provided; use otu_table as the feature abundance input!")
					}
				}
				if(taxa_level == "OTU"){
					use_abund <- self$dataset$otu_table
					cat("The taxa_level is OTU; Use otu_table as abundance table input ...\n")
					# add OTU as one column to make the operations concordant
					if(! taxa_level %in% colnames(self$dataset$tax_table)){
						self$dataset$add_rownames2taxonomy(use_name = "OTU")
					}
				}else{
					if(! taxa_level %in% colnames(self$dataset$tax_table)){
						stop("The taxa_level provided is neither OTU nor one of the column names of dataset$tax_table!")
					}else{
						newdat <- self$dataset$merge_taxa(taxa_level)
						use_abund <- newdat$otu_table
					}
				}
				if(!is.null(taxa_filter_thres)){
					use_abund <- use_abund[apply(use_abund, 1, sum)/sum(use_abund) > taxa_filter_thres, ]
				}
				use_data <- as.data.frame(t(use_abund))
			}
			if(feature_sel == T){
				message('Start forward selection ...')
				if(method == "dbRDA"){
					mod0 <- dbrda(use_data ~ 1, env_data, ...)
					mod1 <- dbrda(use_data ~ ., env_data, ...)
				}
				if(method == "RDA"){
					mod0 <- rda(use_data ~ 1, env_data, ...)
					mod1 <- rda(use_data ~ ., env_data, ...)					
				}
				if(method == "CCA"){
					mod0 <- cca(use_data ~ 1, env_data, ...)
					mod1 <- cca(use_data ~ ., env_data, ...)					
				}
				forward_res <- ordiR2step(mod0, scope = formula(mod1), direction="forward", perm.max = 999)
				res_sign <- gsub("+ ", "", rownames(data.frame(forward_res$anova)), fixed = TRUE)
				if(length(res_sign) == 0){
					stop("Non variables obtained after selection according to model. Check method and data!")
				}
				res_sign <- res_sign[1:(length(res_sign) - 1)]
				env_data <- env_data[, c(res_sign), drop=FALSE]
			}
			self$ordination_method <- method
			self$taxa_level <- taxa_level
			if(method == "dbRDA"){
				res_ordination <- dbrda(use_data ~ ., env_data, ...)
			}
			if(method == "RDA"){
				res_ordination <- rda(use_data ~ ., env_data, ...)
			}
			if(method == "CCA"){
				res_ordination <- cca(use_data ~ ., env_data, ...)
			}
			self$res_ordination <- res_ordination
			message('The original ordination result is stored in object$res_ordination ...')
			self$res_ordination_R2 <- unlist(RsquareAdj(res_ordination))
			message('The R2 is stored in object$res_ordination_R2 ...')
		},
		#' @description
		#' Use anova to test the significance of the terms and axis in ordination.
		#'
		#' @param ... the parameters passed to \code{anova} function.
		#' @return \code{res_ordination_terms and res_ordination_axis} in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination_anova()
		#' }
		cal_ordination_anova = function(...){
			if(is.null(self$res_ordination)){
				stop("Please first run cal_ordination function to obtain the ordination result!")
			}else{
				# test for sig.environ.variables
				self$res_ordination_terms <- anova(self$res_ordination, by = "terms", permu = 1000, ...)
				message('The terms anova result is stored in object$res_ordination_terms ...')
				self$res_ordination_axis <- anova(self$res_ordination, by = "axis", perm.max = 1000, ...)
				message('The axis anova result is stored in object$res_ordination_axis ...')
			}
		},
		#' @description
		#' Fits each environmental vector onto the ordination to obtain the contribution of each variable.
		#'
		#' @param ... the parameters passed to \code{vegan::envfit} function.
		#' @return \code{res_ordination_envfit} in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination_envfit()
		#' }
		cal_ordination_envfit = function(...){
			if(is.null(self$res_ordination)){
				stop("Please first run cal_ordination function to obtain the ordination result!")
			}else{
				self$res_ordination_envfit <- vegan::envfit(self$res_ordination, self$data_env, ...)
				message('Result is stored in object$res_ordination_envfit ...')
			}
		},
		#' @description
		#' Transform ordination results for the following plot.
		#'
		#' @param show_taxa default 10; taxa number shown in the plot.
		#' @param adjust_arrow_length default FALSE; whether adjust the arrow length to be clearer.
		#' @param min_perc_env default 0.1; used for scaling up the minimum of env arrow; multiply by the maximum distance between samples and origin.
		#' @param max_perc_env default 0.8; used for scaling up the maximum of env arrow; multiply by the maximum distance between samples and origin.
		#' @param min_perc_tax default 0.1; used for scaling up the minimum of tax arrow; multiply by the maximum distance between samples and origin.
		#' @param max_perc_tax default 0.8; used for scaling up the maximum of tax arrow; multiply by the maximum distance between samples and origin.
		#' @return \code{res_ordination_trans} in the object.
		#' @examples
		#' \donttest{
		#' t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.1, max_perc_env = 1)
		#' }
		trans_ordination = function(
			show_taxa = 10, 
			adjust_arrow_length = FALSE, 
			min_perc_env = 0.1, 
			max_perc_env = 0.8, 
			min_perc_tax = 0.1, 
			max_perc_tax = 0.8
			){
			if(is.null(self$res_ordination)){
				stop("Please first run cal_ordination function !")
			}
			res_ordination <- self$res_ordination
			scrs <- scores(res_ordination)
			scrs$biplot <- scores(res_ordination, choices = c(1, 2), "bp", scaling = "sites")
			df_sites <- cbind.data.frame(scrs$sites, self$dataset$sample_table[rownames(scrs$sites), ])
			colnames(df_sites)[1:2] <- c("x","y")
			
			multiplier <- vegan:::ordiArrowMul(scrs$biplot)
			df_arrows<- scrs$biplot * multiplier
			colnames(df_arrows)<-c("x","y")
			df_arrows <- as.data.frame(df_arrows)
			eigval <- res_ordination$CCA$eig/sum(res_ordination$CCA$eig)
			eigval <- round(100 * eigval, 1)
			eigval[1] <- paste0("RDA1", " [", eigval[1], "%]")
			eigval[2] <- paste0("RDA2", " [", eigval[2], "%]")

			if(self$ordination_method != "dbRDA"){
				scrs$biplot_spe <- scores(res_ordination, choices=c(1, 2), "sp", scaling="species")
				df_species <- scrs$species
				colnames(df_species)[1:2] <- c("x","y")
				multiplier_spe <- vegan:::ordiArrowMul(scrs$biplot_spe)
				df_arrows_spe <- scrs$biplot_spe * multiplier_spe
				colnames(df_arrows_spe) <- c("x","y")
				df_arrows_spe <- dropallfactors(cbind.data.frame(
					df_arrows_spe, 
					self$dataset$tax_table[rownames(df_arrows_spe), self$taxa_level, drop = FALSE]
					))
				df_arrows_spe %<>% .[!grepl("__$|__uncultured|sp$", .[, 3]), ]
				df_arrows_spe <- df_arrows_spe %>% 
					{.[,1]^2 + .[,2]^2} %>% 
					`names<-`(rownames(df_arrows_spe)) %>% 
					sort(., decreasing = TRUE) %>% 
					.[1:show_taxa] %>% 
					names %>% 
					df_arrows_spe[., ]
			}else{
				df_species <- NULL
				df_arrows_spe <- NULL
			}
			if(adjust_arrow_length == T){
				df_arrows[, 1:2] <- private$stand_fun(arr = df_arrows[, 1:2], ref = df_sites, min_perc = min_perc_env, max_perc = max_perc_env)
				if(self$ordination_method != "dbRDA"){
					df_arrows_spe[, 1:2] <- private$stand_fun(arr = df_arrows_spe[, 1:2], ref = df_sites, min_perc = min_perc_tax, max_perc = max_perc_tax)
				}
			}
			self$res_ordination_trans <- list(
				df_sites = df_sites, 
				df_arrows = df_arrows, 
				eigval = eigval, 
				df_species = df_species, 
				df_arrows_spe = df_arrows_spe
				)
			message('The result list is stored in object$res_ordination_trans ...')
		},
		#' @description
		#' plot ordination result.
		#'
		#' @param plot_color default NULL; a colname of \code{sample_table} to assign colors to different groups.
		#' @param plot_shape default NULL; a colname of \code{sample_table} to assign shapes to different groups.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color pallete for different groups.
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); a vector for point shape types of groups, see ggplot2 tutorial.
		#' @param env_text_color default "black"; environmental variable text color.
		#' @param env_arrow_color default "grey30"; environmental variable arrow color.
		#' @param taxa_text_color default "firebrick1"; taxa text color.
		#' @param taxa_arrow_color default "firebrick1"; taxa arrow color.
		#' @param env_text_size default 3.7; environmental variable text size.
		#' @param taxa_text_size default 3; taxa text size.
		#' @param taxa_text_italic default TRUE; "italic"; whether use "italic" style for the taxa text in the plot.
		#' @param plot_type default "point"; one or more elements of "point", "ellipse", "chull" and "centroid".
		#'   \describe{
		#'     \item{\strong{'point'}}{add point}
		#'     \item{\strong{'ellipse'}}{add confidence ellipse for points of each group}
		#'     \item{\strong{'chull'}}{add convex hull for points of each group}
		#'     \item{\strong{'centroid'}}{add centroid line of each group}
		#'   }
		#' @param point_size default 3; point size in plot when "point" is in plot_type.
		#' @param point_alpha default .8; point transparency in plot when "point" is in plot_type.
		#' @param centroid_segment_alpha default 0.6; segment transparency in plot when "centroid" is in plot_type.
		#' @param centroid_segment_size default 1; segment size in plot when "centroid" is in plot_type.
		#' @param centroid_segment_linetype default 3; an integer; the line type related with centroid in plot when "centroid" is in plot_type.
		#' @param ellipse_chull_fill default TRUE; whether fill colors to the area of ellipse or chull.
		#' @param ellipse_chull_alpha default 0.1; color transparency in the ellipse or convex hull depending on whether "ellipse" or "centroid" is in plot_type.
		#' @param ellipse_level default .9; confidence level of ellipse when "ellipse" is in plot_type.
		#' @param ellipse_type default "t"; ellipse type when "ellipse" is in plot_type; see type in \code{\link{stat_ellipse}}.
		#' @param add_sample_label default NULL; the column name in sample table, if provided, show the point name in plot.
		#' @param env_nudge_x default NULL; numeric vector to adjust the env text x axis position; passed to nudge_x parameter of \code{ggrepel::geom_text_repel} function;
		#'   default NULL represents automatic adjustment; the length must be same with the row number of \code{object$res_ordination_trans$df_arrows}. For example, 
		#'   if there are 5 env variables, env_nudge_x should be something like \code{c(0.1, 0, -0.2, 0, 0)}. 
		#'   Note that this parameter and env_nudge_y is generally used when the automatic text adjustment is not very well.
		#' @param env_nudge_y default NULL; numeric vector to adjust the env text y axis position; passed to nudge_y parameter of ggrepel::geom_text_repel function;
		#'   default NULL represents automatic adjustment; the length must be same with the row number of \code{object$res_ordination_trans$df_arrows}. For example, 
		#'   if there are 5 env variables, env_nudge_y should be something like \code{c(0.1, 0, -0.2, 0, 0)}.
		#' @param taxa_nudge_x default NULL; numeric vector to adjust the taxa text x axis position; passed to nudge_x parameter of ggrepel::geom_text_repel function;
		#'   default NULL represents automatic adjustment; the length must be same with the row number of \code{object$res_ordination_trans$df_arrows_spe}. For example, 
		#'   if 3 taxa are shown, taxa_nudge_x should be something like \code{c(0.3, -0.2, 0)}.
		#' @param taxa_nudge_y default NULL; numeric vector to adjust the taxa text y axis position; passed to nudge_y parameter of ggrepel::geom_text_repel function;
		#'   default NULL represents automatic adjustment; the length must be same with the row number of \code{object$res_ordination_trans$df_arrows_spe}. For example, 
		#'   if 3 taxa are shown, taxa_nudge_y should be something like \code{c(-0.2, 0, 0.4)}.
		#' @param ... paremeters passed to \code{geom_point} for controlling sample points.
		#' @return ggplot object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination(method = "RDA")
		#' t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
		#' t1$plot_ordination(plot_color = "Group")
		#' t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "chull"))
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "centroid"), 
		#' 	  centroid_segment_linetype = 1)
		#' t1$plot_ordination(plot_color = "Group", env_nudge_x = c(0.4, 0, 0, 0, 0, -0.2, 0, 0), 
		#' 	  env_nudge_y = c(0.6, 0, 0.2, 0.5, 0, 0.1, 0, 0.2))
		#' }
		plot_ordination = function(
			plot_color = NULL,
			plot_shape = NULL,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			env_text_color = "black",
			env_arrow_color = "grey30",
			taxa_text_color = "firebrick1",
			taxa_arrow_color = "firebrick1",
			env_text_size = 3.7,
			taxa_text_size = 3,
			taxa_text_italic = TRUE,
			plot_type = "point",
			point_size = 3,
			point_alpha = 0.8,
			centroid_segment_alpha = 0.6,
			centroid_segment_size = 1,
			centroid_segment_linetype = 3,
			ellipse_chull_fill = TRUE,
			ellipse_chull_alpha = 0.1,
			ellipse_level = 0.9,
			ellipse_type = "t",
			add_sample_label = NULL,
			env_nudge_x = NULL,
			env_nudge_y = NULL,
			taxa_nudge_x = NULL,
			taxa_nudge_y = NULL,
			...
			){
			if(is.null(self$res_ordination_trans)){
				stop("Please first run trans_ordination function !")
			}
			if(is.null(plot_color)){
				if(any(c("ellipse", "chull", "centroid") %in% plot_type)){
					stop("Plot ellipse or chull or centroid need groups! Please provide plot_color parameter!")
				}
			}
			if(! all(plot_type %in% c("point", "ellipse", "chull", "centroid"))){
				message("There maybe a typo in your plot_type input! plot_type should be one or more from 'point', 'ellipse', 'chull' and 'centroid'!")
			}
			df_sites <- self$res_ordination_trans$df_sites
			
			p <- ggplot()
			p <- p + theme_bw()
			p <- p + theme(panel.grid=element_blank())
			p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey80")
			p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey80")
			if("point" %in% plot_type){
				p <- p + geom_point(
					data = df_sites, 
					aes_meco("x", "y", colour = plot_color, shape = plot_shape), 
					size = point_size,
					alpha = point_alpha,
					...
					)
			}
			# plot arrows
			env_text_data <- self$res_ordination_trans$df_arrows %>% dplyr::mutate(label = gsub("`", "", rownames(.)))
			p <- p + geom_segment(
				data = env_text_data, 
				aes(x = 0, y = 0, xend = x, yend = y), 
				arrow = arrow(length = unit(0.2, "cm")), 
				color = env_arrow_color
				)
			if(is.null(env_nudge_x) & is.null(env_nudge_y)){
				p <- p + ggrepel::geom_text_repel(
					data = env_text_data, 
					aes(x, y, label = label), 
					size = env_text_size, 
					color = env_text_color, 
					segment.color = "white"
				)
			}else{
				if(is.null(env_nudge_x)){
					env_nudge_x <- rep(0, nrow(env_text_data))
				}else{
					if(length(env_nudge_x) != nrow(env_text_data)){
						stop("The length of env_nudge_x not equal to the number of environmental variables !")
					}
					if(! inherits(env_nudge_x, "numeric")){
						stop("The env_nudge_x must be numeric !")
					}
				}
				if(is.null(env_nudge_y)){
					env_nudge_y <- rep(0, nrow(env_text_data))
				}else{
					if(length(env_nudge_y) != nrow(env_text_data)){
						stop("The length of env_nudge_y not equal to the number of environmental variables !")
					}
					if(! inherits(env_nudge_y, "numeric")){
						stop("The env_nudge_y must be numeric !")
					}
				}
				for(i in seq_len(nrow(env_text_data))){
					p <- p + ggrepel::geom_text_repel(data = env_text_data[i, ], aes(x, y, label = label), size = env_text_size, 
						color = env_text_color, segment.color = "white", nudge_x = env_nudge_x[i], nudge_y = env_nudge_y[i])
				}
			}
			if("centroid" %in% plot_type){
				centroid_xy <- data.frame(group = df_sites[, plot_color], x = df_sites[, "x"], y = df_sites[, "y"]) %>%
					dplyr::group_by(group) %>%
					dplyr::summarise(cx = mean(x), cy = mean(y)) %>%
					as.data.frame()
				combined_centroid_xy <- merge(df_sites, centroid_xy, by.x = plot_color, by.y = "group")
				p <- p + geom_segment(
					data = combined_centroid_xy, 
					aes_meco(x = "x", xend = "cx", y = "y", yend = "cy", colour = plot_color),
					alpha = centroid_segment_alpha, 
					size = centroid_segment_size, 
					linetype = centroid_segment_linetype
				)
			}
			if(any(c("ellipse", "chull") %in% plot_type)){
				if(ellipse_chull_fill){
					ellipse_chull_fill_color <- plot_color
				}else{
					ellipse_chull_fill_color <- NULL
					ellipse_chull_alpha <- 0
				}
				mapping <- aes_meco(x = "x", y = "y", group = plot_color, colour = plot_color, fill = ellipse_chull_fill_color)
				if("ellipse" %in% plot_type){
					p <- p + ggplot2::stat_ellipse(
						mapping = mapping, 
						data = df_sites, 
						level = ellipse_level, 
						type = ellipse_type, 
						alpha = ellipse_chull_alpha, 
						geom = "polygon"
						)
				}
				if("chull" %in% plot_type){
					p <- p + ggpubr::stat_chull(
						mapping = mapping, 
						data = df_sites, 
						alpha = ellipse_chull_alpha,
						geom = "polygon"
						)
				}
				if(ellipse_chull_fill){
					p <- p + scale_fill_manual(values = color_values)
				}
			}
			if(!is.null(plot_color)){
				p <- p + scale_color_manual(values = color_values)
			}
			if(!is.null(plot_shape)){
				p <- p + scale_shape_manual(values = shape_values)
			}
			p <- p + xlab(self$res_ordination_trans$eigval[1]) + ylab(self$res_ordination_trans$eigval[2])
			
			if(self$ordination_method != "dbRDA"){
				df_arrows_spe1 <- self$res_ordination_trans$df_arrows_spe
				p <- p + geom_segment(
					data = df_arrows_spe1, 
					aes(x = 0, y = 0, xend = x, yend = y), 
					arrow = arrow(length = unit(0.2, "cm")), 
					color = taxa_arrow_color, 
					alpha = .6
					)
				df_arrows_spe1[, self$taxa_level] %<>% gsub(".*__", "", .) %>% gsub("Candidatus ", "", .) 
				if(taxa_text_italic == T){
					df_arrows_spe1[, self$taxa_level] %<>% paste0("italic('", .,"')")
				}
				if(is.null(taxa_nudge_x) & is.null(taxa_nudge_y)){
					p <- p + ggrepel::geom_text_repel(
						data = df_arrows_spe1, 
						aes_meco("x", "y", label = self$taxa_level), 
						size = taxa_text_size, 
						color = taxa_text_color, 
						segment.alpha = .01, 
						parse = TRUE
					)
				}else{
					if(is.null(taxa_nudge_x)){
						taxa_nudge_x <- rep(0, nrow(df_arrows_spe1))
					}else{
						if(length(taxa_nudge_x) != nrow(df_arrows_spe1)){
							stop("The length of taxa_nudge_x not equal to the number of environmental variables !")
						}
						if(! inherits(taxa_nudge_x, "numeric")){
							stop("The taxa_nudge_x must be numeric !")
						}
					}
					if(is.null(taxa_nudge_y)){
						taxa_nudge_y <- rep(0, nrow(df_arrows_spe1))
					}else{
						if(length(taxa_nudge_y) != nrow(df_arrows_spe1)){
							stop("The length of taxa_nudge_y not equal to the number of environmental variables !")
						}
						if(! inherits(taxa_nudge_y, "numeric")){
							stop("The taxa_nudge_y must be numeric !")
						}
					}
					for(j in seq_len(nrow(df_arrows_spe1))){
						p <- p + ggrepel::geom_text_repel(data = df_arrows_spe1[j, ], aes_meco("x", "y", label = self$taxa_level), size = taxa_text_size, 
							color = taxa_text_color, segment.alpha = .01, parse = TRUE, nudge_x = taxa_nudge_x[j], nudge_y = taxa_nudge_y[j])
					}
				}
			}
			if(!is.null(add_sample_label)){
				p <- p + ggrepel::geom_text_repel(
					data = df_sites,
					mapping = aes_meco(x = "x", y = "y", label = add_sample_label)
					)
			}
			p
		},
		#' @description
		#' Mantel test between beta diversity matrix and environmental data.
		#'
		#' @param select_env_data default NULL; numeric or character vector to select columns in \code{data_env} of the object; 
		#' 	 if not provided, automatically select the columns with numeric attributes.
		#' @param partial_mantel default FALSE; whether use partial mantel test; If TRUE, use other all measurements as the zdis in each calculation.
		#' @param add_matrix default NULL; additional distance matrix provided when the beta diversity matrix in the dataset is not used.
		#' @param use_measure default NULL; a name of beta diversity matrix. If necessary and not provided, use the first beta diversity matrix.
		#' @param method default "pearson"; one of "pearson", "spearman" and "kendall"; correlation method; see method parameter in \code{vegan::mantel} function.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of \code{p.adjust} function for available options.
		#' @param by_group default NULL; one column name or number in sample_table; used to perform mantel test for different groups separately.
		#' @param ... paremeters passed to \code{\link{mantel}} of vegan package.
		#' @return \code{res_mantel} in object.
		#' @examples
		#' \donttest{
		#' t1$cal_mantel(use_measure = "bray")
		#' t1$cal_mantel(partial_mantel = TRUE, use_measure = "bray")
		#' }
		cal_mantel = function(
			select_env_data = NULL, 
			partial_mantel = FALSE, 
			add_matrix = NULL, 
			use_measure = NULL, 
			method = "pearson", 
			p_adjust_method = "fdr",
			by_group = NULL,
			...
			){
			if(is.null(self$dataset) & is.null(add_matrix)){
				stop("No beta diversity data found! please see add_matrix parameter or provide dataset when creating the object!")
			}
			if(is.null(add_matrix)){
				if(is.null(self$dataset$beta_diversity)){
					message("The beta_diversity in dataset is NULL; try to calculate it ...")
					self$dataset$cal_betadiv(unifrac = FALSE)
					message("Calculating done ...")
				}else{
					if(!is.null(use_measure)){
						use_matrix <- self$dataset$beta_diversity[[use_measure]]
					}else{
						use_matrix <- self$dataset$beta_diversity[[1]]
					}
				}
			}else{
				use_matrix <- add_matrix
			}
			if(is.null(self$data_env)){
				stop("The data_env is NULL! Please check the data input when creating the object !")
			}else{
				env_data <- self$data_env
			}
			# if no selection, automatically select the numeric columns
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric))]
			}else{
				env_data <- env_data[, select_env_data]
			}
			if(is.null(by_group)){
				veg_dist <- as.dist(use_matrix[rownames(env_data), rownames(env_data)])
				res_mantel <- private$mantel_test(env = env_data, veg = veg_dist, partial_mantel = partial_mantel, method = method, 
					p_adjust_method = p_adjust_method, by_group = by_group, ...)
			}else{
				res_mantel <- data.frame()
				all_groups <- self$dataset$sample_table %>% dropallfactors %>% .[, by_group] %>% unique
				for(k in all_groups){
					use_sample_names <- self$dataset$sample_table %>% .[.[, by_group] == k, ] %>% rownames
					use_env_data <- env_data[use_sample_names, ]
					use_veg_dist <- as.dist(use_matrix[use_sample_names, use_sample_names])
					tmp_res <- private$mantel_test(env = use_env_data, veg = use_veg_dist, partial_mantel = partial_mantel, 
						method = method, p_adjust_method = p_adjust_method, by_group = k, ...)
					res_mantel %<>% rbind(., tmp_res)
				}
			}
			self$res_mantel <- res_mantel
			message('The result is stored in object$res_mantel ...')
		},
		#' @description
		#' Calculate the correlations between taxa abundance and environmental variables.
		#' Actually, it can also be applied for calculating other correlation between any two variables from two tables.
		#'
		#' @param use_data default "Genus"; "Genus", "all" or "other"; "Genus" or other taxonomic name: use genus or other taxonomic abundance table in \code{taxa_abund}; 
		#'    "all": use all merged taxa abundance table; "other": provide additional taxa name with \code{other_taxa} parameter which is necessary.
		#' @param select_env_data default NULL; numeric or character vector to select columns in data_env; if not provided, automatically select the columns with numeric attributes.
		#' @param cor_method default "pearson"; "pearson", "spearman" or "kendall"; correlation method.
		#' @param p_adjust_method default "fdr"; p.adjust method; see method parameter of \code{p.adjust} function for available options.
		#' 	  \code{p_adjust_method = "none"} can disable the p value adjustment.
		#' @param p_adjust_type default "All"; "All", "Type", "Taxa" or "Env"; P value adjustment type.
		#' 	  "Env": adjustment for each environmental variable separately; 
		#' 	  "Taxa": adjustment for each taxon separately; 
		#' 	  "Type": adjustment according to the groups provided. If \code{by_group} is NULL, adjustment is performed for all the data together.
		#' 	  If \code{by_group} is provided, for each group in it separately.
		#' 	  These three options are the first three colnames of return table \code{res_cor}.
		#' 	  "All": adjustment for all the data together no matter whether \code{by_group} is provided. If \code{by_group} is NULL, it is same with the "Type" option.
		#' @param add_abund_table default NULL; additional data table to be used. Samples must be rows.
		#' @param by_group default NULL; one column name or number in sample_table; calculate correlations for different groups separately.
		#' @param use_taxa_num default NULL; integer; a number used to select high abundant taxa; only useful when \code{use_data} parameter is a taxonomic level, e.g., "Genus".
		#' @param other_taxa default NULL; character vector containing a series of taxa names; used when use_data = "other"; 
		#' 	  the provided names should be standard full names used to select taxa from all the tables in taxa_abund list of the microtable object;
		#' 	  please see the example.
		#' @param group_use default NULL; numeric or character vector to select one column in sample_table for selecting samples; together with group_select.
		#' @param group_select default NULL; the group name used; remain samples within the group.
		#' @param taxa_name_full default TRUE; Whether use the complete taxonomic name of taxa.
		#' @return \code{res_cor} in object.
		#' @examples
		#' \donttest{
		#' t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", rf_taxa_level = "Genus")
		#' t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
		#' t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_diff$Taxa[1:40])
		#' t1$cal_cor(use_data = "other", p_adjust_type = "Env", other_taxa = t2$res_diff$Taxa[1:40])
		#' }
		cal_cor = function(
			use_data = c("Genus", "all", "other")[1],
			select_env_data = NULL,
			cor_method = c("pearson", "spearman", "kendall")[1],
			p_adjust_method = "fdr",
			p_adjust_type = c("All", "Type", "Taxa", "Env")[1],
			add_abund_table = NULL,
			by_group = NULL,
			use_taxa_num = NULL,
			other_taxa = NULL,
			group_use = NULL,
			group_select = NULL,
			taxa_name_full = TRUE
			){
			if(is.null(self$data_env)){
				stop("The data_env is NULL! Please check the data input when creating the object !")
			}else{
				env_data <- self$data_env
			}
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric)), drop = FALSE]
			}
			if(!is.null(add_abund_table)){
				abund_table <- add_abund_table
			}else{
				if(is.null(self$dataset$taxa_abund)){
					message('Invoke cal_abund function to calculate taxonomic abundance automatically ...')
					self$dataset$cal_abund()
				}
				if(use_data %in% names(self$dataset$taxa_abund)){
					abund_table <- self$dataset$taxa_abund[[use_data]]
				}else{
					if(grepl("all|other", use_data, ignore.case = TRUE)){
						abund_table <- do.call(rbind, unname(self$dataset$taxa_abund))
						if(use_data == "other"){
							if(is.null(other_taxa)){
								stop("You select other, but no other_taxa provided!")
							}
							abund_table <- abund_table[other_taxa, ]
						}
					}else{
						stop("Unknown use_data parameter provided!")
					}
				}
				abund_table %<>% .[!grepl("__$", rownames(.)), ]
				if(nrow(abund_table) == 0){
					stop("No available feature! Please check the input data!")
				}
				if(use_data %in% names(self$dataset$taxa_abund) & !is.null(use_taxa_num)){
					if(nrow(abund_table) > use_taxa_num){
						abund_table %<>% .[1:use_taxa_num, ] 
					}
				}
				abund_table <- as.data.frame(t(abund_table))
			}
			# filter samples by one group
			if(!is.null(group_use)){
				if(is.null(group_select)){
					stop("You select group_use parameter, but no group_select parameter provided!")
				}
				sel_sample_names <- self$dataset$sample_table %>% 
					.[.[, group_use] %in% group_select, ] %>% 
					rownames
				abund_table <- abund_table[sel_sample_names, ]
			}
			env_data %<>% .[rownames(.) %in% rownames(abund_table), , drop = FALSE]
			abund_table <- abund_table[rownames(env_data), , drop = FALSE]
			if(is.null(by_group)){
				groups <- rep("All", nrow(env_data))
			}else{
				groups <- self$dataset$sample_table[, by_group] %>% as.character
				message("Perform correlation by the groups in ", by_group, " of sample_table, respectively ...")
			}
			comb_names <- expand.grid(unique(groups), colnames(abund_table), colnames(env_data)) %>% 
				t %>% 
				as.data.frame(stringsAsFactors = FALSE)
			res <- sapply(comb_names, function(x){
				suppressWarnings(cor.test(abund_table[groups == x[1], x[2]], env_data[groups == x[1], x[3]], method = cor_method)) %>%
				{c(x, Correlation = unname(.$estimate), Pvalue = unname(.$p.value))}
				}
			)
			res %<>% t %>% as.data.frame(stringsAsFactors = FALSE)
			colnames(res) <- c("Type", "Taxa", "Env", "Correlation","Pvalue")
			res$Pvalue %<>% as.numeric
			res$Correlation %<>% as.numeric
			res$AdjPvalue <- rep(0, nrow(res))
			if(p_adjust_type == "All"){
				res$AdjPvalue <- p.adjust(res$Pvalue, method = p_adjust_method)
			}else{
				choose_col <- which(c("Type", "Taxa", "Env") %in% p_adjust_type)
				comb_names2 <- comb_names[choose_col, ] %>% t %>% as.data.frame %>% unique %>% t %>% as.data.frame(stringsAsFactors = FALSE)
				# p value adjustment separately
				for(i in seq_len(ncol(comb_names2))){
					x <- comb_names2[, i]
					row_sel <- unlist(lapply(as.data.frame(t(res[, choose_col, drop = FALSE])), function(y) all(y %in% x)))
					res$AdjPvalue[row_sel] <- p.adjust(res[row_sel, "Pvalue"], method = p_adjust_method)
				}
			}
			res$Significance <- cut(res$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
			res <- res[complete.cases(res), ]
			res$Env %<>% factor(., levels = unique(as.character(.)))
			if(taxa_name_full == F){
				res$Taxa %<>% gsub(".*__(.*?)$", "\\1", .)
			}
			self$res_cor <- res
			message('The correlation result is stored in object$res_cor ...')
			self$cor_method <- cor_method
		},
		#' @description
		#' Plot correlation heatmap.
		#'
		#' @param color_vector default \code{c("#053061", "white", "#A50026")}; colors with only three values representing low, middle and high values.
		#' @param color_palette default NULL; a customized palette with more color values; if provided, use it instead of color_vector.
		#' @param pheatmap default FALSE; whether use pheatmap package to plot the heatmap.
		#' @param filter_feature default NULL; character vector; used to filter features that only have significance labels in the filter_feature vector. 
		#'   For example, filter_feature = "" can be used to filter features that only have "", no any "*".
		#' @param ylab_type_italic default FALSE; whether use italic type for y lab text.
		#' @param keep_full_name default FALSE; whether use the complete taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param text_y_order default NULL; character vector; provide customized text order for y axis; shown in the plot from the top down.
		#' @param text_x_order default NULL; character vector; provide customized text order for x axis.
		#' @param font_family default NULL; font family used in \code{ggplot2}; only available when pheatmap = FALSE.
		#' @param cluster_ggplot default "none"; add clustering dendrogram for \code{ggplot2} based heatmap; 
		#'   available options: "none", "row", "col" or "both". "none": no any clustering used;
		#'   "row": add clustering for rows; "col": add clustering for columns; "both":  add clustering for both rows and columns.
		#'   Only available when pheatmap = FALSE.
		#' @param cluster_height_rows default 0.2, the dendrogram plot height for rows; available when cluster_ggplot != "none".
		#' @param cluster_height_cols default 0.2, the dendrogram plot height for columns; available cluster_ggplot != "none".
		#' @param text_y_position default "right"; "left" or "right"; the y axis text position; ggplot2 based heatmap.
		#' @param mylabels_x default NULL; provide x axis text labels additionally; only available when pheatmap = TRUE.
		#' @param ... paremeters passed to \code{ggplot2::geom_tile} or pheatmap, depending on the \code{pheatmap} parameter is FALSE or TRUE.
		#' @return plot.
		#' @examples
		#' \donttest{
		#' t1$plot_cor(pheatmap = FALSE)
		#' }
		plot_cor = function(
			color_vector = c("#053061", "white", "#A50026"),
			color_palette = NULL,
			pheatmap = FALSE,
			filter_feature = NULL,
			ylab_type_italic = FALSE,
			keep_full_name = FALSE,
			keep_prefix = TRUE,
			text_y_order = NULL,
			text_x_order = NULL,
			font_family = NULL,
			cluster_ggplot = "none",
			cluster_height_rows = 0.2,
			cluster_height_cols = 0.2,
			text_y_position = "right",
			mylabels_x = NULL,
			...
			){
			if(is.null(self$res_cor)){
				stop("Please first run cal_cor function to get plot data !")
			}
			if(length(color_vector) != 3){
				stop("color_vector parameter must have three values !")
			}
			cluster_ggplot <- match.arg(cluster_ggplot, c("none", "row", "col", "both"))
			use_data <- self$res_cor
			
			# filter features
			if(!is.null(filter_feature)){
				x1 <- unlist(lapply(unique(use_data$Taxa), function(x){
					t2 <- use_data %>% .[.$Taxa == x, "Significance"] %>% {all(. %in% filter_feature)}
					if(t2 == F){
						x
					}
				}))
				use_data %<>% .[.$Taxa %in% x1, ]
			}
			if(keep_full_name == F){
				use_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				use_data$Taxa %<>% gsub(".*__", "", .)
			}
			if(pheatmap & length(unique(use_data$Type)) > 1){
				use_data$Env <- paste0(use_data$Type, ": ", use_data$Env)
			}
			if(length(unique(use_data$Type)) == 1 | pheatmap){
				clu_data <- reshape2::dcast(use_data, Taxa~Env, value.var = "Correlation") %>% 
					`row.names<-`(.[,1]) %>% 
					.[, -1, drop = FALSE]
				clu_data[is.na(clu_data)] <- 0
				sig_data <- reshape2::dcast(use_data, Taxa~Env, value.var = "Significance") %>% 
					`row.names<-`(.[,1]) %>% 
					.[, -1, drop = FALSE]
				sig_data[is.na(sig_data)] <- ""
				if(length(unique(use_data$Type)) == 1 & pheatmap == F){
					row_cluster <- hclust(dist(clu_data)) 
					lim_y <- row_cluster %>% {.$labels[.$order]}
					col_cluster <- hclust(dist(t(clu_data)))
					lim_x <- col_cluster %>% {.$labels[.$order]}
				}
			}
			# the input text_y_order or text_x_order has priority
			if(!is.null(text_y_order) | !is.null(text_x_order)){
				if(cluster_ggplot != "none"){
					cluster_ggplot <- "none"
					message("Change cluster_ggplot to none, as text_y_order and/or text_x_order provided!")
				}
				if(!is.null(text_y_order) & !is.null(text_x_order)){
					lim_y <- rev(text_y_order)
					lim_x <- text_x_order
				}else{
					if(!is.null(text_y_order) & is.null(text_x_order)){
						lim_y <- rev(text_y_order)
					}else{
						lim_x <- text_x_order
					}
				}
			}
			
			# check whether the cor values all larger or smaller than 0
			if(all(use_data$Correlation >= 0) | all(use_data$Correlation <= 0)){
				color_palette <- color_vector
			}
			if(pheatmap == T){
				if(!require("pheatmap")){
					stop("pheatmap package not installed")
				}
				# check sd for each feature, if 0, delete
				if(any(apply(clu_data, MARGIN = 1, FUN = function(x) sd(x) == 0))){
					select_rows <- apply(clu_data, MARGIN = 1, FUN = function(x) sd(x) != 0)
					clu_data %<>% {.[select_rows, ]}
					sig_data %<>% {.[select_rows, ]}
				}
				if(ylab_type_italic == T){
					eval(parse(text = paste0(
						"mylabels_y <- c(", paste0("expression(italic(", paste0('"', rownames(clu_data),'"'), "))", collapse = ","),")", 
						collapse = "")))
				}else{
					mylabels_y <- rownames(clu_data)
				}
				# create the palette
				if(is.null(color_palette)){
					minpiece <- max(abs(min(clu_data)), max(clu_data))/100
					pos <- seq(from = 0, to = max(clu_data), by = minpiece)[-1]
					neg <- rev(seq(from = 0, to = abs(min(clu_data)), by = minpiece)[-1]) * (-1)
					posColor <- colorRampPalette(c(color_vector[2], color_vector[3]))(100)[1:length(pos)]
					negColor <- colorRampPalette(c(color_vector[2], color_vector[1]))(100)[1:length(neg)] %>% rev
					# make sure color_vector_use and myBreaks have same length
					color_vector_use <- c(negColor, color_vector[2], posColor)
					myBreaks <- c(neg, 0, pos)		
				}else{
					color_vector_use <- colorRampPalette(color_palette)(100)
					myBreaks <- NA
				}
				
				p <- pheatmap(
					clu_data, 
					clustering_distance_row = "correlation", 
					clustering_distance_cols= "correlation",
					border_color = NA, 
					scale = "none",
					labels_row = mylabels_y, 
					labels_col = mylabels_x, 
					display_numbers = sig_data, 
					number_color = "black", 
					color = color_vector_use,
					breaks = myBreaks,
					...
					)
				p$gtable
			}else{
				p <- ggplot(aes(x = Env, y = Taxa, fill = Correlation), data = use_data) +
					theme_bw() + 
					geom_tile(...)
				if(is.null(color_palette)){
					p <- p + scale_fill_gradient2(low = color_vector[1], high = color_vector[3], mid = color_vector[2])
				} else {
					p <- p + scale_fill_gradientn(colours = color_palette)
				}
				
				p <- p + geom_text(aes(label = Significance), color="black", size=4) + 
					labs(y = NULL, x = "Measure", fill = self$cor_method) +
					theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 10)) +
					theme(strip.background = element_rect(fill = "grey85", colour = "white"), axis.title = element_blank()) +
					theme(strip.text=element_text(size=11), panel.border = element_blank(), panel.grid = element_blank())
				if(length(unique(use_data$Type)) == 1){
					if(cluster_ggplot == "none"){
						p <- p + scale_y_discrete(limits = lim_y, position = text_y_position) + scale_x_discrete(limits = lim_x)
					}else{
						p <- p + scale_y_discrete(limits = lim_y, position = text_y_position) + scale_x_discrete(limits = lim_x)
						if(cluster_ggplot %in% c("row", "both")){
							row_plot <- ggtree::ggtree(row_cluster, hang = 0)
							p %<>% aplot::insert_left(row_plot, width = cluster_height_rows)
						}
						if(cluster_ggplot %in% c("col", "both")){
							col_plot <- ggtree::ggtree(col_cluster, hang = 0) + ggtree::layout_dendrogram()							
							p %<>% aplot::insert_top(col_plot, height = cluster_height_cols)
						}
					}
				}else{
					p <- p + facet_grid(. ~ Type, drop = TRUE, scale = "free", space = "free_x")
				}
				if(ylab_type_italic == T){
					p <- p + theme(axis.text.y = element_text(face = 'italic'))
				}
				if(!is.null(font_family)){
					p <- p + theme(text = element_text(family = font_family))
				}
				p
			}
		},
		#' @description
		#' Scatter plot with fitted line. The most important thing is to make sure that the input x and y
		#'  have correponding sample orders. If one of x and y is a matrix, the other will be also transformed to matrix with Euclidean distance.
		#'  Then, both of them are transformed to be vectors. If x or y is a vector with a single value, x or y will be
		#'  assigned according to the column selection of the \code{data_env} inside.
		#'
		#' @param x default NULL; a single numeric or character value or a vector or a distance matrix used for the x axis.
		#'     If x is a single value, it will be used to select the column of \code{data_env} inside.
		#'     If x is a distance matrix, it will be transformed to be a vector.
		#' @param y default NULL; a single numeric or character value or a vector or a distance matrix used for the y axis.
		#'     If y is a single value, it will be used to select the column of \code{data_env} inside.
		#'     If y is a distance matrix, it will be transformed to be a vector.
		#' @param group default NULL; a character vector; if length is 1, must be a colname of dataset$sample_table;
		#'    Otherwise, group should be a vector with same length of x/y (for vector) or ncol of x/y (for matrix).
		#' @param group_order default NULL; a vector to order groups, i.e. reorder the legend and colors in plot when group is not NULL; 
		#' 	  If group_order is NULL and group is provided, the function can first check whether the group column of dataset$sample_table is factor. 
		#' 	  If provided, overlook the levels in the group of dataset$sample_table.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); color pallete for different groups.
		#' @param shape_values default NULL; a numeric vector for point shape types of groups when group is not NULL, see ggplot2 tutorial.
		#' @param type default c("cor", "lm")[1]; "cor": correlation; "lm" for regression.
		#' @param cor_method default "pearson"; one of "pearson", "kendall" and "spearman"; correlation method.
		#' @param label_sep default ";"; the separator string between different label parts.
		#' @param label.x.npc default "left"; can be numeric or character vector of the same length as the number of groups and/or panels. If too short they will be recycled.
		#' \describe{
		#'   \item{numeric}{value should be between 0 and 1. Coordinates to be used for positioning the label, expressed in "normalized parent coordinates"}
		#'   \item{character}{allowed values include: i) one of c('right', 'left', 'center', 'centre', 'middle') for x-axis; ii) and one of 
		#'      c( 'bottom', 'top', 'center', 'centre', 'middle') for y-axis.}
		#' }
		#' @param label.y.npc default "top"; same usage with label.x.npc; see also label.y.npc parameter of \code{ggpubr::stat_cor} package.
		#' @param label.x default NULL; x axis absolute position for adding the statictic label.
		#' @param label.y default NULL; x axis absolute position for adding the statictic label.
		#' @param x_axis_title default ""; the title of x axis.
		#' @param y_axis_title default ""; the title of y axis.
		#' @param point_size default 5; point size value.
		#' @param point_alpha default 0.6; alpha value for the point color transparency.
		#' @param line_size default 0.8; line size value.
		#' @param line_se default TRUE; Whether show the confidence interval for the fitting.
		#' @param line_se_color default "grey70"; the color to fill the confidence interval when \code{line_se = TRUE}.
		#' @param line_alpha default 1; alpha value for the line color transparency.
		#' @param line_color default "black"; fitted line color only useful when group = NULL.
		#' @param pvalue_trim default 4; trim the decimal places of p value.
		#' @param cor_coef_trim default 3; trim the decimal places of correlation coefficient.
		#' @param lm_equation default TRUE; whether include the equation in the label when \code{type = "lm"}.
		#' @param lm_fir_trim default 2; trim the decimal places of regression first coefficient.
		#' @param lm_sec_trim default 2; trim the decimal places of regression second coefficient.
		#' @param lm_squ_trim default 2; trim the decimal places of regression R square.
		#' @param ... other arguments passed to \code{geom_text} or \code{geom_label}.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_scatterfit(x = 1, y = 2, type = "cor")
		#' t1$plot_scatterfit(x = 1, y = 2, type = "lm", point_alpha = .3)
		#' t1$plot_scatterfit(x = "pH", y = "TOC", type = "lm", group = "Group", line_se = FALSE)
		#' t1$plot_scatterfit(x = 
		#' 	 dataset$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)], y = "pH")
		#' }
		plot_scatterfit = function(
			x = NULL, 
			y = NULL, 
			group = NULL,
			group_order = NULL,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			shape_values = NULL,
			type = c("cor", "lm")[1],
			cor_method = "pearson",
			label_sep = ";",
			label.x.npc = "left", 
			label.y.npc = "top",
			label.x = NULL, 
			label.y = NULL, 
			x_axis_title = "", 
			y_axis_title = "",
			point_size = 5,
			point_alpha = 0.6,
			line_size = 0.8,
			line_alpha = 1,
			line_color = "black",
			line_se = TRUE,
			line_se_color = "grey70",
			pvalue_trim = 4, 
			cor_coef_trim = 3,
			lm_equation = TRUE,
			lm_fir_trim = 2,
			lm_sec_trim = 2,
			lm_squ_trim = 2,
			...
			){
			if(!(is.vector(x) | is.matrix(x))){
				stop("The input x is neither a vector nor a matrix !")
			}
			if(!(is.vector(y) | is.matrix(y))){
				stop("The input y is neither a vector nor a matrix !")
			}
			type <- match.arg(type, c("cor", "lm"))

			if(length(x) == 1){
				x <- self$data_env[, x]
			}
			if(length(y) == 1){
				y <- self$data_env[, y]
			}
			# first get the raw data length according to x
			x_raw_length <- ifelse(is.vector(x), length(x), ncol(x))
			if(is.null(group)){
				group_vector <- rep("all", x_raw_length)
			}else{
				if(length(group) == 1){
					if(is.null(self$dataset)){
						stop("No dataset$sample_table found! Please check the dataset input when creading the object!")
					}else{
						if(! group %in% colnames(self$dataset$sample_table)){
							stop("The group must be one of colnames of dataset$sample_table!")
						}
						# the sample_table names have been same with data_env when creading the object
						group_vector <- self$dataset$sample_table[, group] %>% as.character
					}
				}
			}
			if(length(group_vector) != x_raw_length){
				stop("The group length is not same with x data length! Please check the input!")
			}
			if(is.vector(x) & is.vector(y)){
				if(length(x) != length(y)){
					stop("x length is not equal to y length!")
				}
				use_data <- data.frame(x = x, y = y, Group = group_vector)
			}
			# Matrix is transformed to be a vector
			if(is.matrix(x) | is.matrix(y)){
				if(is.vector(x)){
					x1 <- as.data.frame(as.numeric(as.character(x)))
					x <- as.matrix(vegdist(x1, "euclidean"))
				}
				# has make sure x is a matrix
				x <- lapply(unique(group_vector), function(i){
						data.frame("x" = as.vector(as.dist(x[group_vector == i, group_vector == i])), "Group" = i)
					}) %>% 
					do.call(rbind, .)
				if(is.vector(y)){
					y1 <- as.data.frame(as.numeric(as.character(y)))
					y <- as.matrix(vegdist(y1, "euclidean"))
				}
				y <- lapply(unique(group_vector), function(i){
						data.frame("y" = as.vector(as.dist(y[group_vector == i, group_vector == i])), "Group" = i)
					}) %>% 
					do.call(rbind, .)
				if(nrow(x) != nrow(y)){
					stop("Converted x length is not equal to y length! Please check the input vector or matrix !")
				}else{
					use_data <- data.frame(x = x[, "x"], y = y[, "y"], Group = x[, "Group"])
				}
			}
			if(!is.null(group_order)){
				if(!all(unique(as.character(use_data$Group)) %in% group_order)){
					stop("Please check the group_order input! Part of group(s) not found in group_order!")
				}
			}
			if(!is.null(group)){
				if(length(group) == 1){
					if(is.null(group_order)){
						if(is.factor(self$dataset$sample_table[, group])){
							use_data$Group %<>% factor(., levels = levels(self$dataset$sample_table[, group]))
						}else{
							use_data$Group %<>% as.character %>% as.factor
						}
					}else{
						use_data$Group %<>% factor(., levels = group_order)
					}
				}else{
					if(is.null(group_order)){
						use_data$Group %<>% as.character %>% as.factor
					}else{
						use_data$Group %<>% factor(., levels = group_order)
					}
				}
			}
			if(!is.null(shape_values)){
				plot_shape <- "Group"
			}else{
				plot_shape <- NULL
			}
			if(is.null(group)){
				p <- ggplot(use_data, aes_meco(x = "x", y = "y"))
			}else{
				p <- ggplot(use_data, aes_meco(x = "x", y = "y", colour = "Group", shape = plot_shape))
			}
			p <- p + geom_point(size = point_size, alpha = point_alpha)
			if(is.null(group)){
				p <- p + geom_smooth(method = "lm", size = line_size, color = line_color, alpha = line_alpha, se = line_se, fill = line_se_color)
			}else{
				p <- p + geom_smooth(method = "lm", size = line_size, alpha = line_alpha, se = line_se, fill = line_se_color)
			}
			p <- p + stat_corlm(
				type = type, 
				cor_method = cor_method, 
				label_sep = label_sep, 
				pvalue_trim = pvalue_trim,
				cor_coef_trim = cor_coef_trim,
				lm_equation = lm_equation,
				lm_fir_trim = lm_fir_trim, 
				lm_sec_trim = lm_sec_trim,
				lm_squ_trim = lm_squ_trim,
				label.x.npc = label.x.npc, 
				label.y.npc = label.y.npc, 
				label.x = label.x, 
				label.y = label.y,
				...
				)
			p <- p + xlab(x_axis_title) + ylab(y_axis_title)
			if(!is.null(group)){
				p <- p + scale_color_manual(values = color_values)
				if(!is.null(shape_values)){
					p <- p + scale_shape_manual(values = shape_values)
				}
			}
			p
		},
		#' @description
		#' Print the trans_env object.
		print = function(){
			cat("trans_env object:\n")
			if(!is.null(self$data_env)){
				cat(paste0("Env table have ", ncol(self$data_env), " variables: ", paste0(colnames(self$data_env), collapse = ","), "\n"))
			}else{
				cat("No environmental variable table stored in the object.\n")
			}
		}
	),
	private = list(
		# transformation function
		stand_fun = function(arr, ref, min_perc = NULL, max_perc = NULL) {
			# arr and ref must be a two column data.frame or matrix
			if(is.null(min_perc)){
				min_perc <- 0.1
			}
			if(is.null(max_perc)){
				max_perc <- 0.7
			}
			ref_dis <- ref[, 1]^2 + ref[, 2]^2
			max_ref_dis <- max(ref_dis)
			a <- max_ref_dis * min_perc
			b <- max_ref_dis * max_perc
			arr_dis <- arr[, 1]^2 + arr[, 2]^2
			maxdis <- max(arr_dis)
			mindis <- min(arr_dis)
			k <- (b-a)/(maxdis - mindis) 
			norDis <- a + k * (arr_dis - mindis)
			per <- abs(arr[,1]/arr[,2])
			newx <- (((norDis*per^2) / (per^2 + 1)) ^ (1/2)) * sapply(arr[, 1], function(y) ifelse(y > 0, 1, -1))
			newy <- (newx/per) * sapply(arr[, 2], function(y) ifelse(y > 0, 1, -1))
			res <- data.frame(newx, newy)
			colnames(res) <- colnames(arr)
			res
		},
		mantel_test = function(env, veg, partial_mantel, method, p_adjust_method, by_group, ...){
			variable_name <- c()
			corr_res <- c()
			p_res <- c()
			# with Scaling and Centering normalization
			for(i in 1:ncol(env)){
				env_dist <- vegdist(scale(env[, i, drop = FALSE]), "euclid")
				if(partial_mantel == T){
					zdis <- vegdist(scale(env[, -i, drop = FALSE]), "euclid")
					man1 <- mantel.partial(veg, env_dist, zdis, method = method, ...)
				}else{
					man1 <- mantel(veg, env_dist, method = method, ...)
				}
				variable_name <- c(variable_name, colnames(env)[i])
				corr_res <- c(corr_res, man1$statistic)
				p_res <- c(p_res, man1$signif)
			}
			if(partial_mantel == T){
				mantel_type <- rep("partial mantel test", length(p_res))
			}else{
				mantel_type <- rep("mantel test", length(p_res))
			}
			cor_method <- rep(method, length(p_res))
			p_adjusted <- p.adjust(p_res, method = p_adjust_method)
			significance <- cut(p_adjusted, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))
			if(is.null(by_group)){
				by_group <- "All"
			}
			res_mantel <- data.frame(by_group, variable_name, mantel_type, cor_method, corr_res, p_res, p_adjusted, significance)
			colnames(res_mantel) <- c("by_group", "Variables", "mantel type", "Correlation method", "Correlation coefficient", "p.value", "p.adjusted", "Significance")
			res_mantel
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
