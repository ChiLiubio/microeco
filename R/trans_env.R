#' @title 
#' Create trans_env object for the analysis of the effects of environmental factors on communities.
#'
#' @description
#' This class is a wrapper for a series of operations associated with environmental measurements, including redundancy analysis, 
#' mantel test and correlation analysis based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>.
#'
#' @export
trans_env <- R6Class(classname = "trans_env",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param env_cols default NULL; a vector to select columns in sample_table, when the environmental data is in sample_table. 
		#'   Either numeric vector or character vector of colnames.
		#' @param add_data default NULL; data.frame format; provide the environmental data frame individually; rownames should be sample names.
		#' @param character2numeric default TRUE; whether transform the characters or factors to numeric attributes.
		#' @param complete_na default FALSE; Whether fill the NA (missing value) in the environmental data.
		#' @return env_data in trans_env object.
		#' @examples
		#' data(dataset)
		#' data(env_data_16S)
		#' t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S)
		initialize = function(
			dataset = NULL, 
			env_cols = NULL, 
			add_data = NULL, 
			character2numeric = TRUE, 
			complete_na = FALSE
			){
			# support all NULL from v0.7.0
			# first judge dataset input
			if(is.null(dataset)){
				cat("The dataset not provided. Remember to provide additional data in the function ...\n")
			}
			if(is.null(add_data)){
				if(!is.null(env_cols)){
					env_data <- dataset$sample_table[, env_cols, drop = FALSE]
				}else{
					env_data <- NULL
				}
			}else{
				if(is.null(dataset)){
					env_data <- add_data
				}else{
					env_data <- add_data[rownames(add_data) %in% rownames(dataset$sample_table), , drop = FALSE]
				}
			}
			
			if(!is.null(dataset)){
				dataset1 <- clone(dataset)
				if(!is.null(env_data)){
					inter_sum <- sum(rownames(dataset1$sample_table) %in% rownames(env_data))
					if(inter_sum == 0){
						stop("No sample names of sample_table found in env_data! Please check the names of env_data!")
					}
					if(inter_sum < nrow(dataset1$sample_table)){
						message(nrow(dataset1$sample_table)-inter_sum, " samples not found in env_data and removed!")
						dataset1$sample_table %<>% base::subset(rownames(.) %in% rownames(env_data))
						dataset1$tidy_dataset(main_data = FALSE)
					}
					env_data %<>% .[rownames(dataset1$sample_table), , drop = FALSE]
				}
				self$dataset <- dataset1
			}else{
				self$dataset <- NULL
			}
			if(complete_na == T){
				env_data[env_data == ""] <- NA
				env_data %<>% dropallfactors(., unfac2num = TRUE)
				env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
				env_data %<>% mice::mice(print = FALSE) %>% mice::complete(., 1)
			}
			if(character2numeric == T){
				if(!is.null(env_data)){
					env_data %<>% dropallfactors(., unfac2num = TRUE, char2num = TRUE)
				}
			}
			self$env_data <- env_data
		},
		#' @description
		#' Redundancy analysis (RDA) and Correspondence Analysis (CCA) based on the vegan package.
		#'
		#' @param method default c("RDA", "dbRDA", "CCA")[1]; the ordination method.
		#' @param feature_sel default FALSE; whether perform the feature selection.
		#' @param taxa_level default NULL; If use RDA or CCA, provide the taxonomic rank, such as "Phylum" or "Genus";
		#'   If use otu_table; please input "OTU".
		#' @param taxa_filter_thres default NULL; If want to filter taxa, provide the relative abundance threshold.
		#' @param use_measure default NULL; name of beta diversity matrix. If necessary and not provided, use the first beta diversity matrix.
		#' @param add_matrix default NULL; additional distance matrix provided, when you do not want to use the beta diversity matrix within the dataset;
		#'   only available when method = "dbRDA".
		#' @param ... paremeters pass to dbrda or rda or cca function according to the input of method.
		#' @return res_ordination, res_ordination_R2, res_ordination_terms and res_ordination_axis in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination(method = "RDA", use_measure = "bray")
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
			if(is.null(self$env_data)){
				stop("The env_data is NULL! Please check the env input when creating the object !")
			}else{
				env_data <- self$env_data
			}
			
			if(method == "dbRDA"){
				# add_matrix is a priority
				if(!is.null(add_matrix)){
					use_matrix <- add_matrix
					message("Use the additional matrix provided for dbRDA ...")
				}else{
					if(is.null(self$dataset)){
						stop("No distance matrix provided and the dataset in the object is NULL! please check the input!")
					}else{
						if(is.null(self$dataset$beta_diversity)){
							message("The beta_diversity in dataset is NULL; try to calculate it ...")
							self$dataset$cal_betadiv(unifrac = FALSE)
							message("Calculating done ...")
						}
						if(!is.null(use_measure)){
							use_matrix <- self$dataset$beta_diversity[[use_measure]]
							message("Use ", use_measure, " in dataset$beta_diversity for dbRDA ...")
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
					self$dataset$add_rownames2taxonomy(use_name = "OTU")
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
				env_data <- env_data[,c(res_sign),drop=FALSE]
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
			# test for sig.environ.variables
			self$res_ordination_terms <- anova(res_ordination, by = "terms", permu = 1000)
			message('The terms anova result is stored in object$res_ordination_terms ...')
			self$res_ordination_axis <- anova(res_ordination, by = "axis", perm.max = 1000)
			message('The axis anova result is stored in object$res_ordination_axis ...')
		},
		#' @description
		#' Fits each environmental vector onto the ordination to obtain the contribution of each variable.
		#'
		#' @param ... the parameters passing to vegan::envfit function.
		#' @return res_ordination_envsquare in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ordination_envsquare()
		#' }
		cal_ordination_envsquare = function(...){
			if(is.null(self$res_ordination)){
				stop("Please first run cal_ordination function to obtain the ordination result!")
			}else{
				self$res_ordination_envsquare <- vegan::envfit(self$res_ordination, self$env_data, ...)
				message('Result is stored in object$res_ordination_envsquare ...')
			}
		},
		#' @description
		#' transform ordination result for the following plotting.
		#'
		#' @param show_taxa default 10; taxa number shown in the plot.
		#' @param adjust_arrow_length default FALSE; whether adjust the arrow length to be clear
		#' @param min_perc_env default 0.1; used for scaling up the minimum of env arrow; multiply by the maximum distance between samples and origin.
		#' @param max_perc_env default 0.8; used for scaling up the maximum of env arrow; multiply by the maximum distance between samples and origin.
		#' @param min_perc_tax default 0.1; used for scaling up the minimum of tax arrow; multiply by the maximum distance between samples and origin.
		#' @param max_perc_tax default 0.8; used for scaling up the maximum of tax arrow; multiply by the maximum distance between samples and origin.
		#' @return res_ordination_trans in object.
		#' @examples
		#' \donttest{
		#' t1$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 10)
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
			scrs <- scores(res_ordination ,choices = c(1, 2), display = c("sp", "wa", "cn"))
			scrs$biplot <- scores(res_ordination, choices=c(1, 2), "bp", scaling="sites")
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
				colnames(df_arrows_spe)<-c("x","y")
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
				df_arrows[,1:2] <- private$stand_fun(arr = df_arrows[,1:2], ref = df_sites, min_perc = min_perc_env, max_perc = max_perc_env)
				if(self$ordination_method != "dbRDA"){
					df_arrows_spe[,1:2] <- private$stand_fun(arr = df_arrows_spe[,1:2], ref = df_sites, min_perc = min_perc_tax, max_perc = max_perc_tax)
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
		#' @param plot_color default NULL; group used for color.
		#' @param plot_shape default NULL; group used for shape.
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); color pallete.
		#' @param shape_values default see the function; vector used in the shape, see ggplot2 tutorial.
		#' @param taxa_text_color default "firebrick1"; taxa text colors.
		#' @param taxa_text_italic default TRUE; "italic"; whether use "italic" style for the taxa text in the plot.
		#' @return ggplot object.
		#' @examples
		#' \donttest{
		#' t1$plot_ordination(plot_color = "Group")
		#' }
		plot_ordination = function(
			plot_color = NULL,
			plot_shape = NULL,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			taxa_text_color = "firebrick1",
			taxa_text_italic = TRUE
			){
			if(is.null(self$res_ordination_trans)){
				stop("Please first run trans_ordination function !")
			}
			
			p <- ggplot()
			p <- p + theme_bw()
			p <- p + theme(panel.grid=element_blank())
			p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey80")
			p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey80")
			p <- p + geom_point(
				data=self$res_ordination_trans$df_sites, 
				aes_string("x", "y", colour = plot_color, shape = plot_shape), 
				size= 3.5)
			# plot arrows
			p <- p + geom_segment(
				data=self$res_ordination_trans$df_arrows, 
				aes(x = 0, y = 0, xend = x, yend = y), 
				arrow = arrow(length = unit(0.2, "cm")), 
				color = "grey30")
			p <- p + ggrepel::geom_text_repel(
				data = as.data.frame(self$res_ordination_trans$df_arrows * 1), 
				aes(x, y, label = gsub("`", "", rownames(self$res_ordination_trans$df_arrows))), 
				size=3.7, 
				color = "black", 
				segment.color = "white"
				)
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
					data=df_arrows_spe1, 
					aes(x = 0, y = 0, xend = x, yend = y), 
					arrow = arrow(length = unit(0.2, "cm")), 
					color = "firebrick1", 
					alpha = .6
					)
				df_arrows_spe1[, self$taxa_level] %<>% gsub(".*__", "", .) %>% gsub("Candidatus ", "", .) 
				if(taxa_text_italic == T){
					df_arrows_spe1[, self$taxa_level] %<>% paste0("italic('", .,"')")
				}
				p <- p + ggrepel::geom_text_repel(
					data = df_arrows_spe1, 
					aes_string("x", "y", label = self$taxa_level), 
					size=3, 
					color = taxa_text_color, 
					segment.alpha = .01, 
					parse = TRUE
				)
			}
			p
		},
		#' @description
		#' Mantel test between beta diversity matrix and environmental data.
		#'
		#' @param select_env_data default NULL; numeric or character vector to select columns in env_data; if not provided, automatically select the columns with numeric attributes.
		#' @param partial_mantel default FALSE; whether use partial mantel test; If TRUE, use other measurements as the zdis.
		#' @param add_matrix default NULL; additional distance matrix provided, if you donot want to use the beta diversity matrix in the dataset.
		#' @param use_measure default NULL; name of beta diversity matrix. If necessary and not provided, use the first beta diversity matrix.
		#' @param method default "pearson"; one of "pearson", "spearman" and "kendall"; correlation method.
		#' @param ... paremeters pass to \code{\link{mantel}}.
		#' @return res_mantel in object.
		#' @examples
		#' \donttest{
		#' t1$cal_mantel(use_measure = "bray")
		#' }
		cal_mantel = function(
			select_env_data = NULL, 
			partial_mantel = FALSE, 
			add_matrix = NULL, 
			use_measure = NULL, 
			method = "pearson", 
			...
			){
			if(is.null(self$dataset$beta_diversity) & is.null(add_matrix)){
				stop("No distance matrix provided; please use set add_matrix parameter or use provide dataset in creating Class !")
			}
			if(is.null(add_matrix)){
				if(!is.null(use_measure)){
					use_matrix <- self$dataset$beta_diversity[[use_measure]]
				}else{
					use_matrix <- self$dataset$beta_diversity[[1]]
				}
			}else{
				use_matrix <- add_matrix
			}
			env_data <- self$env_data
			# if no selection, automatically select the numeric columns
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric))]
			}else{
				env_data <- env_data[, select_env_data]
			}
			veg.dist <- as.dist(use_matrix[rownames(env_data), rownames(env_data)])
			variable_name <- c()
			corr_res <- c()
			p_res <- c()
			# with Scaling and Centering normalization
			for(i in 1:ncol(env_data)){
				env.dist <- vegdist(scale(env_data[, i, drop=FALSE]), "euclid")
				if(partial_mantel == T){
					zdis <- vegdist(scale(env_data[, -i, drop=FALSE]), "euclid")
					man1 <- mantel.partial(veg.dist, env.dist, zdis, ...)
				}else{
					man1 <- mantel(veg.dist, env.dist, ...)
				}
				variable_name <- c(variable_name, colnames(env_data)[i])
				corr_res <- c(corr_res, man1$statistic)
				p_res <- c(p_res, man1$signif)
			}
			cor_method <- rep(method, length(p_res))
			significance <- cut(p_res, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
			res_mantel <- data.frame(variable_name, cor_method, corr_res, p_res, significance)
			if(partial_mantel == T){
				self$res_mantel_partial <- res_mantel
			}else{
				self$res_mantel <- res_mantel
			}
			message('The result is stored in object$res_mantel or object$res_mantel_partial ...')
		},
		#' @description
		#' Calculating the correlations between taxa abundance and environmental variables.
		#' Indeed, it can also be used for calculating other correlation between any two variables from two tables.
		#'
		#' @param use_data default "Genus"; "Genus", "all" or "other"; "Genus" or other taxonomic name: use genus or other taxonomic abundance table in taxa_abund; 
		#'    "all": use all merged taxa abundance table; "other": provide additional taxa name with other_taxa parameter which is necessary.
		#' @param select_env_data default NULL; numeric or character vector to select columns in env_data; if not provided, automatically select the columns with numeric attributes.
		#' @param cor_method default "pearson"; "pearson", "spearman" or "kendall"; correlation method.
		#' @param p_adjust_method default "fdr"; p.adjust method.
		#' @param p_adjust_type default "Env"; "Type", "Taxa" or "Env"; p.adjust type; Env: environmental data; Taxa: taxa data; Type: group used.
		#' @param add_abund_table default NULL; additional data table to be used. Samples must be rows.
		#' @param by_group default NULL; one column name or number in sample_table; calculate correlations for different groups separately.
		#' @param use_taxa_num default NULL; integer; a number used to select high abundant taxa; only useful when use_data parameter is a taxonomic level, e.g. "Genus".
		#' @param other_taxa default NULL; provide additional taxa, see use_data parameter.
		#' @param group_use default NULL; numeric or character vector to select one column in sample_table for selecting samples; together with group_select.
		#' @param group_select default NULL; the group name used; will retain samples within the group.
		#' @param taxa_name_full default TRUE; Whether retain the complete taxonomic name of taxa.
		#' @return res_cor in object.
		#' @examples
		#' \donttest{
		#' t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", rf_taxa_level = "Genus")
		#' t1 <- trans_env$new(dataset = dataset, add_data = env_data_16S[, 4:11])
		#' t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_rf$Taxa[1:40])
		#' }
		cal_cor = function(
			use_data = c("Genus", "all", "other")[1],
			select_env_data = NULL,
			cor_method = c("pearson", "spearman", "kendall")[1],
			p_adjust_method = "fdr",
			p_adjust_type = c("Type", "Taxa", "Env")[3],
			add_abund_table = NULL,
			by_group = NULL,
			use_taxa_num = NULL,
			other_taxa = NULL,
			group_use = NULL,
			group_select = NULL,
			taxa_name_full = TRUE
			){
			env_data <- self$env_data
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric)), drop = FALSE]
			}
			if(!is.null(add_abund_table)){
				abund_table <- add_abund_table
			}else{
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
				abund_table %<>% .[!grepl("__$|__uncultured$", rownames(.)), ]
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
			abund_table <- abund_table[rownames(env_data), ]
			if(is.null(by_group)){
				groups <- rep("All", nrow(env_data))
			}else{
				groups <- self$dataset$sample_table[, by_group] %>% as.character
				message("Calculate the corr by the groups in ", by_group, " of sample_table, respectively ...")
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
			choose_col <- which(c("Type", "Taxa", "Env") %in% p_adjust_type)
			comb_names2 <- comb_names[choose_col, ] %>% 
				t %>% 
				as.data.frame %>% 
				unique %>% 
				t %>% 
				as.data.frame(stringsAsFactors = FALSE)
			# p value adjust by groups
			for(i in seq_len(ncol(comb_names2))){
				x <- comb_names2[, i]
				row_sel <- unlist(lapply(as.data.frame(t(res[, choose_col, drop = FALSE])), function(y) all(y %in% x)))
				res$AdjPvalue[row_sel] <- p.adjust(res[row_sel, "Pvalue"], method = p_adjust_method)
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
		#' @param color_vector default c("#053061", "white", "#A50026"); colors with only three values representing low, middle and high value.
		#' @param color_palette default NULL; a customized palette with more color values; if provided, use it instead of color_vector.
		#' @param pheatmap default FALSE; whether use pheatmap package to plot the heatmap.
		#' @param filter_feature default NULL; character vector; used to filter features that only have significance labels in the filter_feature vector. 
		#'   For example, filter_feature = "" can be used to filter features that only have "", no any "*".
		#' @param ylab_type_italic default FALSE; whether use italic type for y lab text.
		#' @param keep_full_name default FALSE; whether use the complete taxonomic name.
		#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
		#' @param text_y_order default NULL; character vector; provide customized text order for y axis; shown in the plot the top down.
		#' @param text_x_order default NULL; character vector; provide customized text order for x axis.
		#' @param font_family default NULL; font family used in ggplot2; only available when pheatmap = FALSE.
		#' @param cluster_ggplot default "none"; add clustering dendrogram for ggplot2 based heatmap; 
		#'   available options: "none", "row", "col" or "both". "none": no any clustering used;
		#'   "row": add clustering for rows; "col": add clustering for columns; "both":  add clustering for both rows and columns.
		#'   Only available when pheatmap = FALSE.
		#' @param cluster_height_rows default 0.2, the dendrogram plot height for rows; available when cluster_ggplot != "none".
		#' @param cluster_height_cols default 0.2, the dendrogram plot height for columns; available cluster_ggplot != "none".
		#' @param text_y_position default "right"; "left" or "right"; the y axis text position; ggplot2 based heatmap.
		#' @param fontsize default 9; base fontsize for the plot; see fontsize in pheatmap.
		#' @param mylabels_x default NULL; provide x axis text labels additionally; only available when pheatmap = TRUE.
		#' @param ... paremeters pass to ggplot2::geom_tile or pheatmap, depending on the pheatmap = FALSE or TRUE.
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
			fontsize = 9,
			mylabels_x = NULL,
			...
			){
			if(is.null(self$res_cor)){
				stop("Please first use cal_cor to get plot data !")
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
					fontsize = fontsize, 
					labels_row = mylabels_y, 
					labels_col = mylabels_x, 
					display_numbers = sig_data, 
					number_color = "black", 
					fontsize_number = fontsize * 1.2,
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
							row_plot <- factoextra::fviz_dend(
								row_cluster, 
								show_labels = FALSE, 
								labels_track_height = 0, 
								horiz = TRUE, 
								main = "",
								xlab = "",
								ylab = "",
								ggtheme = theme_classic()) + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
							p %<>% aplot::insert_left(row_plot, width = cluster_height_rows)
						}
						if(cluster_ggplot %in% c("col", "both")){
							col_plot <- factoextra::fviz_dend(
								col_cluster, 
								show_labels = FALSE, 
								labels_track_height = 0, 
								horiz = FALSE, 
								main = "",
								xlab = "",
								ylab = "",
								ggtheme = theme_classic()) + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
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
		#' Scatter plot and add fitted line. The most important thing is to make sure that the input x and y
		#'  have correponding sample orders. If one of x and y is a matrix, the other will be also transformed to matrix with Euclidean distance.
		#'  Then, both of them are transformed to be vectors. If x or y is a vector with a single value, x or y will be
		#'  assigned according to the column selection of the env_data inside.
		#'
		#' @param x default NULL; a single numeric or character value or a vector or a distance matrix used for the x axis.
		#'     If x is a single value, it will be used to select the column of env_data inside.
		#'     If x is a distance matrix, it will be transformed to be a vector.
		#' @param y default NULL; a single numeric or character value or a vector or a distance matrix used for the y axis.
		#'     If y is a single value, it will be used to select the column of env_data inside.
		#'     If y is a distance matrix, it will be transformed to be a vector.
		#' @param use_cor default TRUE; TRUE for correlation; FALSE for regression.
		#' @param cor_method default "pearson"; one of "pearson", "kendall" and "spearman".
		#' @param add_line default TRUE; whether add the fitted line in the plot.
		#' @param use_se default TRUE; Whether show the confidence interval for the fitting.
		#' @param text_x_pos default NULL; the central x axis position of the fitting text.
		#' @param text_y_pos default NULL; the central y axis position of the fitting text.
		#' @param x_axis_title default ""; the title of x axis.
		#' @param y_axis_title default ""; the title of y axis.
		#' @param pvalue_trim default 4; trim the decimal places of p value.
		#' @param cor_coef_trim default 3; trim the decimal places of correlation coefficient.
		#' @param lm_fir_trim default 2; trim the decimal places of regression first coefficient.
		#' @param lm_sec_trim default 2; trim the decimal places of regression second coefficient.
		#' @param lm_squ_trim default 2; trim the decimal places of regression R square.
		#' @param ... the parameters passing to ggplot2::geom_point function.
		#' @return plot.
		#' @examples
		#' \donttest{
		#' t1$plot_scatterfit(x = 1, y = 2, alpha = .5)
		#' }
		plot_scatterfit = function(
			x = NULL, 
			y = NULL, 
			use_cor = TRUE,
			cor_method = "pearson",
			add_line = TRUE,
			use_se = TRUE,
			text_x_pos = NULL, 
			text_y_pos = NULL, 
			x_axis_title = "", 
			y_axis_title = "",
			pvalue_trim = 4, 
			cor_coef_trim = 3,
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
			if(length(x) == 1){
				x <- self$env_data[, x]
			}
			if(length(y) == 1){
				y <- self$env_data[, y]
			}
			# Matrix is transformed to be a vector
			if(is.matrix(x) | is.matrix(y)){
				if(is.matrix(x) & is.matrix(y)){
					x <- as.vector(as.dist(x))
					y <- as.vector(as.dist(y))
				}else{
					if(is.matrix(x)){
						x <- as.vector(as.dist(x))
						y1 <- as.data.frame(as.numeric(as.character(y)))
						y <- as.vector(vegdist(y1, "euclidean"))
					}else{
						y <- as.vector(as.dist(y))
						x1 <- as.data.frame(as.numeric(as.character(x)))
						x <- as.vector(vegdist(x1, "euclidean"))
					}
				}
			}
			if(length(x) != length(y)){
				stop("The length of x axis vector is not equal to the length of y axis vector!")
			}
			use_data <- data.frame(x, y)
			if(use_cor == T){
				fit <- cor.test(x, y, method = cor_method)
			}else{
				fit <- lm(y ~ x)
			}
			# default position max * .8
			if(is.null(text_x_pos)){
				text_x_pos <- max(use_data$x) * 0.8
			}
			if(is.null(text_y_pos)){
				text_y_pos <- max(use_data$y) * 0.8
			}

			p <- ggplot(use_data, aes(x = x, y = y)) + 
				theme_bw() + 
				geom_point(shape = 20, ...) +
				theme(panel.grid = element_blank())
			if(add_line == T){
				p <- p + geom_smooth(method = "lm", size = .8, colour = "black", se = use_se)
			}
			p <- p + annotate("text", 
					x = text_x_pos, 
					y = text_y_pos, 
					label = private$fit_equat(
						fit, 
						use_cor = use_cor, 
						pvalue_trim = pvalue_trim, 
						cor_coef_trim = cor_coef_trim, 
						lm_fir_trim = lm_fir_trim, 
						lm_sec_trim = lm_sec_trim, 
						lm_squ_trim = lm_squ_trim
						), 
					parse = TRUE) +
#				scale_x_continuous(limits = c(min(x2) - 0.2 * (range(x2)[2] - range(x2)[1]), NA)) +
				xlab(x_axis_title) + 
				ylab(y_axis_title) +
				theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
			
			p
		},
		#' @description
		#' Print the trans_env object.
		print = function(){
			cat("trans_env class:\n")
			if(!is.null(self$env_data)){
				cat(paste0("Env table have ", ncol(self$env_data), " variables: ", paste0(colnames(self$env_data), collapse = ","), "\n"))
			}else{
				cat("No env table stored in the object.\n")
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
		# parse the equation and add the statistics
		fit_equat = function(
			equat, 
			use_cor = TRUE, 
			pvalue_trim = 4, 
			cor_coef_trim = 3, 
			lm_fir_trim = 2, 
			lm_sec_trim = 2, 
			lm_squ_trim = 2
			){
			if(class(equat) == "lm" & use_cor == T){
				stop("Input is lm class, but use_cor is TRUE! Please check the use_cor parameter!")
			}
			pvalue <- ifelse(use_cor == T, equat$p.value, anova(equat)$`Pr(>F)`[1])
			if(use_cor == T){
				estimate = equat$estimate
				all_coef <- list(
					R1 = unname(round(estimate, digits = cor_coef_trim)), 
					P1 = ifelse(pvalue < 0.0001, " < 0.0001", paste0(" = ", round(pvalue, digits = pvalue_trim)))
					)
				res <- substitute(italic(R)~"="~R1*";"~~italic(P)*P1, all_coef)
			}else{
				inte <- round(unname(coef(equat))[1], digits = lm_sec_trim)
				lm_coef <- list(a = ifelse(inte < 0, paste0(" - ", abs(inte)), paste0(" + ", as.character(inte))),
							  b = round(unname(coef(equat))[2], digits = lm_fir_trim),
							  r2 = round(summary(equat)$r.squared, digits = lm_squ_trim),
							  p1 = ifelse(pvalue < 0.0001, " < 0.0001", paste0(" = ", round(pvalue, digits = pvalue_trim)))
							  )
				res <- substitute(italic(y) == b %.% italic(x)*a*","~~italic(R)^2~"="~r2*","~~italic(P)*p1, lm_coef)
			}
			as.character(as.expression(res))
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)

