#' @title Create \code{trans_beta} object for beta-diversity analysis
#'
#' @description
#' This class is a wrapper for a series of beta-diversity related analysis, 
#' including ordination analysis based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>, group distance comparision, 
#' clustering, perMANOVA based on Anderson al. (2008) <doi:10.1111/j.1442-9993.2001.01070.pp.x>, ANOSIM and PERMDISP.
#' Note that the beta diversity analysis methods related with environmental variables are encapsulated within the \code{trans_env} class.
#'
#' @export
trans_beta <- R6Class(classname = "trans_beta",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param measure default NULL; a matrix name stored in \code{microtable$beta_diversity} list, such as "bray" or "jaccard", or a customized matrix; 
		#' 	 used for ordination, manova, group distance comparision, etc.;
		#' 	 Please see \code{cal_betadiv} function of \code{\link{microtable}} class for more details.
		#' @param group default NULL; sample group used for manova, betadisper or group distance comparision.
		#' @return measure, group and dataset stored in the object.
		#' @examples
		#' data(dataset)
		#' t1 <- trans_beta$new(dataset = dataset, measure = "bray", group = "Group")
		initialize = function(
			dataset = NULL, 
			measure = NULL, 
			group = NULL
			){
			check_microtable(dataset)
			if(!is.null(measure)){
				if(is.vector(measure)){
					if(length(measure) > 1){
						stop("The input measure should only have one element! Please check it!")
					}
					if(is.null(dataset$beta_diversity)){
						stop("No beta_diversity list found in the input dataset! Please first use cal_betadiv function in microtable class to calculate it!")
					}
					if(is.character(measure)){
						if(!measure %in% names(dataset$beta_diversity)){
							stop("Input measure: ", measure, " should be one of beta_diversity distance matrix names in dataset! ",
								"Please use names(dataset$beta_diversity) to check it!")
						}
					}else{
						if(is.numeric(measure)){
							if(measure > length(dataset$beta_diversity)){
								stop("Input measure: ", measure, " is larger than total beta_diversity distance matrixes number! Please check it")
							}
							measure %<>% round
						}else{
							stop("Unknown format of input measure parameter!")
						}
					}
					use_matrix <- dataset$beta_diversity[[measure]]
				}else{
					if(is.matrix(measure)){
						if(! any(rownames(measure) %in% rownames(dataset$sample_table))){
							stop("Provided measure is a matrix. The row names should be sample names!")
						}
						if(! any(colnames(measure) %in% rownames(dataset$sample_table))){
							stop("Provided measure is a matrix. The column names should be sample names!")
						}
						if(! all(dataset$sample_names() %in% rownames(measure))){
							stop("Some sample names are not found in the matrix of provided measure!")
						}
						use_matrix <- measure[dataset$sample_names(), dataset$sample_names()]
					}else{
						stop("Input measure parameter should be either a vector or a matrix!")
					}
				}
				self$use_matrix <- use_matrix
			}
			if(!is.null(group)){
				if(! group %in% colnames(dataset$sample_table)){
					stop("Provided group must be one of colnames in sample_table of dataset!")
				}
			}
			self$sample_table <- dataset$sample_table
			self$measure <- measure
			self$group <- group
			use_dataset <- clone(dataset)
			use_dataset$phylo_tree <- NULL
			use_dataset$rep_fasta <- NULL
			use_dataset$taxa_abund <- NULL
			use_dataset$alpha_diversity <- NULL
			self$dataset <- use_dataset
		},
		#' @description
		#' Unconstrained ordination.
		#'
		#' @param method default "PCoA"; "PCoA", "NMDS", "PCA", "DCA", "PLS-DA" or "OPLS-DA".
		#' 	  PCoA: principal coordinates analysis; NMDS: non-metric multidimensional scaling, PCA: principal component analysis; DCA: detrended correspondence analysis; 
		#' 	  PLS-DA: partial least squares discriminant analysis; OPLS-DA: orthogonal partial least squares discriminant analysis.
		#' 	  For the methods details, please refer to the papers <doi:10.1111/j.1574-6941.2007.00375.x> (for PCoA, NMDS, PCA and DCA) and 
		#' 	  <doi:10.1186/s12859-019-3310-7> (for PLS-DA or OPLS-DA).
		#' @param ncomp default 2; dimensions in the result. 
		#' 	  For the method "NMDS", this argument will be passed to the \code{k} parameter in the \code{vegan::metaMDS} function.
		#' @param taxa_level default NULL; available for PCA, DCA or NMDS (\code{NMDS_matrix = TRUE}).
		#' 	  Default NULL means using the \code{otu_table} in the microtable object.
		#' 	  For other options, please	provide the taxonomic rank names in \code{tax_table}, such as "Phylum" or "Genus".
		#' 	  In such cases, the data will be merged according to the provided taxonomic levels to generated a new abundance table.		
		#' @param NMDS_matrix default TRUE; For the NMDS method, whether use a distance matrix as input like PCoA. If it is FALSE, the input will be the abundance table like PCA.
		#' @param trans default FALSE; whether species abundance will be square root transformed; only available when \code{method} is "PCA" or "DCA".
		#' 	  For method "NMDS" and \code{NMDS_matrix = FALSE}, please set the \code{autotransform} parameter, which will be passed to \code{vegan::metaMDS} function directly.
		#' @param scale_species default FALSE; whether species loading in PCA, DCA or NMDS (\code{NMDS_matrix = FALSE}) is scaled.
		#' @param scale_species_ratio default 0.8; the ratio to scale up the loading; multiply by the maximum distance between samples and origin. 
		#' 	  Only available when \code{scale_species = TURE}.
		#' @param orthoI default NA; number of orthogonal components (for OPLS-DA only). Default NA means the number of orthogonal components is automatically computed.
		#' 	  Please also see \code{orthoI} parameter in \code{opls} function of ropls package.
		#' @param ordination deprecated. Please use \code{method} argument instead.
		#' @param ... parameters passed to \code{vegan::rda} function when \code{method = "PCA"}, 
		#' 	  or \code{vegan::decorana} function when \code{method = "DCA"}, 
		#' 	  or \code{ape::pcoa} function when \code{method = "PCoA"}, 
		#' 	  or \code{vegan::metaMDS} function when \code{method = "NMDS"},
		#' 	  or \code{ropls::opls} function when \code{method = "PLS-DA"} or \code{method = "OPLS-DA"} .
		#' @return \code{res_ordination} list stored in the object.
		#' 	  In the list, \code{model} is the original analysis results; \code{scores} is the sample scores table; \code{loading} is the feature loading table.
		#' @examples
		#' t1$cal_ordination(method = "PCoA")
		cal_ordination = function(
			method = "PCoA",
			ncomp = 2,
			taxa_level = NULL,
			NMDS_matrix = TRUE,
			trans = FALSE, 
			scale_species = FALSE,
			scale_species_ratio = 0.8,
			orthoI = NA,
			ordination = deprecated(),
			...
			){
			
			if(lifecycle::is_present(ordination)) {
				lifecycle::deprecate_warn("1.8.0", "cal_ordination(ordination)", "cal_ordination(method)")
				method <- ordination
			}
			
			if(is.null(method)){
				stop("Input method should not be NULL !")
			}
			if(!method %in% c("PCA", "PCoA", "NMDS", "DCA", "PLS-DA", "OPLS-DA")){
				stop("Input method should be one of 'PCoA', 'NMDS', 'PCA', 'DCA', 'PLS-DA' and 'OPLS-DA' !")
			}
			use_data <- self$dataset
			if(! is.null(taxa_level)){
				check_tax_level(taxa_level, use_data)
				use_data <- use_data$merge_taxa(taxa_level)
			}
			if(method %in% c("PCA", "DCA", "PLS-DA", "OPLS-DA")){
				plot.x <- switch(method, PCA = "PC1", DCA = "DCA1", 'PLS-DA' = "p1", 'OPLS-DA' = "p1")
				plot.y <- switch(method, PCA = "PC2", DCA = "DCA2", 'PLS-DA' = "p2", 'OPLS-DA' = "o1")
				if(trans == T){
					abund <- sqrt(use_data$otu_table)
				}else{
					abund <- use_data$otu_table
				}
				abund %<>% t
				if(method == "PCA"){
					model <- rda(abund, ...)
					expla <- round(model$CA$eig/model$CA$tot.chi*100, 1)
				}
				if(method == "DCA"){
					model <- decorana(abund, ...)
					expla <- round(eigenvals(model)/model$totchi * 100, 1)
				}
				if(method %in% c("PLS-DA", "OPLS-DA")){
					if(is.null(self$group)){
						stop("For the current method, group is necessary! Please recreate the trans_beta object and provide the group parameter!")
					}
					use_group <- self$sample_table[, self$group]
					if(method == "PLS-DA"){
						model <- ropls::opls(abund, use_group, orthoI = 0, ...)
					}else{
						model <- ropls::opls(abund, use_group, orthoI = orthoI, ...)
					}
					expla <- model@modelDF[, "R2X"] * 100
					names(expla) <- rownames(model@modelDF)
				}
				if(method %in% c("PCA", "DCA")){
					scores_sites <- scores(model, choices = 1:ncomp, display = "sites")
				}else{
					if(method == "PLS-DA"){
						model_score <- model@scoreMN
					}else{
						model_score <- cbind.data.frame(model@scoreMN, model@orthoScoreMN)
					}
					scores_sites <- model_score[, 1:ncomp, drop = FALSE]
				}
				combined <- cbind.data.frame(scores_sites, use_data$sample_table[rownames(scores_sites), , drop = FALSE])
				if(method %in% c("PCA", "DCA")){
					loading <- scores(model, choices = 1:ncomp, display = "species")
				}else{
					if(method == "PLS-DA"){
						model_loading <- model@loadingMN
					}else{
						model_loading <- cbind.data.frame(model@loadingMN, model@orthoLoadingMN)
					}
					loading <- model_loading[, 1:ncomp, drop = FALSE]
				}
				loading %<>% as.data.frame
				if(scale_species){
					for(i in 1:ncomp){
						maxratio <- max(abs(scores_sites[, i]))/max(abs(loading[, i]))
						loading[, i] %<>% {. * maxratio * scale_species_ratio}
					}
				}
				loading[, "dist"] <- apply(loading, 1, function(x){sum(x^2)})
				loading <- loading[order(loading[, "dist"], decreasing = TRUE), ]
				outlist <- list(model = model, scores = combined, loading = loading, eig = expla)
			}
			if(method == "PCoA"){
				if(is.null(self$use_matrix)){
					stop("Please recreate the object and set the parameter measure !")
				}
				model <- ape::pcoa(as.dist(self$use_matrix), ...)
				combined <- cbind.data.frame(model$vectors[,1:ncomp], use_data$sample_table)
				colnames(combined)[1:ncomp] <- paste0("PCo", 1:ncomp)
				expla <- round(model$values[,1]/sum(model$values[,1])*100, 1)
				names(expla) <- paste0("PCo", 1:length(expla))
				outlist <- list(model = model, scores = combined, eig = expla)
			}
			if(method == "NMDS"){
				if(NMDS_matrix){
					if(is.null(self$use_matrix)){
						stop("Please recreate the object and set the parameter measure !")
					}
					model <- vegan::metaMDS(as.dist(self$use_matrix), k = ncomp, ...)
				}else{
					abund <- use_data$otu_table %>% t
					model <- vegan::metaMDS(abund, k = ncomp, ...)
				}
				combined <- cbind.data.frame(model$points, use_data$sample_table)
				outlist <- list(model = model, scores = combined)
				if(! NMDS_matrix){
					loading <- scores(model, choices = 1:ncomp, display = "species") %>% as.data.frame
					if(scale_species){
						for(i in 1:ncomp){
							maxratio <- max(abs(combined[, i]))/max(abs(loading[, i]))
							loading[, i] %<>% {. * maxratio * scale_species_ratio}
						}
					}
					loading[, "dist"] <- apply(loading, 1, function(x){sum(x^2)})
					loading <- loading[order(loading[, "dist"], decreasing = TRUE), ]
					outlist$loading <- loading
				}
			}
			outlist$ncomp <- ncomp
			self$res_ordination <- outlist
			if(! is.null(taxa_level)){
				self$res_ordination$taxa_level <- taxa_level
			}
			message('The result is stored in object$res_ordination ...')
			self$ordination_method <- method
			invisible(self)
		},
		#' @description
		#' Plot the ordination result.
		#'
		#' @param plot_type default "point"; one or more elements of "point", "ellipse", "chull" and "centroid".
		#'   \describe{
		#'     \item{\strong{'point'}}{add sample points}
		#'     \item{\strong{'ellipse'}}{add confidence ellipse for points of each group}
		#'     \item{\strong{'chull'}}{add convex hull for points of each group}
		#'     \item{\strong{'centroid'}}{add centroid line of each group}
		#'   }
		#' @param choices default c(1, 2); selected axis for the visualization; must be numeric vector.
		#'   The maximum value must not exceed the parameter \code{ncomp} in the \code{cal_ordination} function.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for different groups.
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); a vector for point shape types of groups, see \code{ggplot2} tutorial.
		#' @param plot_color default NULL; a colname of \code{sample_table} to assign colors to different groups in plot.
		#' @param plot_shape default NULL; a colname of \code{sample_table} to assign shapes to different groups in plot.
		#' @param plot_group_order default NULL; a vector used to order the groups in the legend of plot.
		#' @param add_sample_label default NULL; a column name in \code{sample_table}; If provided, show the point name in plot.
		#' @param point_size default 3; point size when "point" is in \code{plot_type} parameter.
		#'   \code{point_size} can also be a variable name in \code{sample_table}, such as "pH".
		#' @param point_alpha default .8; point transparency in plot when "point" is in \code{plot_type} parameter.
		#' @param centroid_segment_alpha default 0.6; segment transparency in plot when "centroid" is in \code{plot_type} parameter.
		#' @param centroid_segment_size default 1; segment size in plot when "centroid" is in \code{plot_type} parameter.
		#' @param centroid_segment_linetype default 3; the line type related with centroid in plot when "centroid" is in \code{plot_type} parameter.
		#' @param ellipse_chull_fill default TRUE; whether fill colors to the area of ellipse or chull.
		#' @param ellipse_chull_alpha default 0.1; color transparency in the ellipse or convex hull depending on whether "ellipse" or "centroid" is in \code{plot_type} parameter.
		#' @param ellipse_level default .9; confidence level of ellipse when "ellipse" is in \code{plot_type} parameter.
		#' @param ellipse_type default "t"; ellipse type when "ellipse" is in \code{plot_type} parameter; see type in \code{stat_ellipse}.
		#' @param NMDS_stress_pos default c(1, 1); a numerical vector with two values used to represent the insertion position of the stress text. 
		#'   The first one denotes the x-axis, while the second one corresponds to the y-axis. 
		#'   The assigned position is determined by multiplying the respective value with the maximum point on the corresponding coordinate axis. 
		#'   Thus, the x-axis position is equal to \code{max(points of x axis) * NMDS_stress_pos[1]}, 
		#'   and the y-axis position is equal to \code{max(points of y axis) * NMDS_stress_pos[2]}. Negative values can also be utilized for the negative part of the axis.
		#'   \code{NMDS_stress_pos = NULL} denotes no stress text to show.
		#' @param NMDS_stress_text_prefix default ""; If NMDS_stress_pos is not NULL, this parameter can be used to add text in front of the stress value.
		#' @param loading_arrow default FALSE; whether show the loading using arrow.
		#' @param loading_taxa_num default 10; the number of taxa used for the loading. Only available when \code{loading_arrow = TRUE}.
		#' @param loading_text_taxlevel default NULL; which level of taxonomic table will be used.
		#'   Default NULL means using the \code{taxa_level} parameter in the previous \code{cal_ordination} function.
		#' @param loading_text_color default "black"; the color of taxa text. Only available when \code{loading_arrow = TRUE}.
		#' @param loading_arrow_color default "grey30"; the color of taxa arrow. Only available when \code{loading_arrow = TRUE}.
		#' @param loading_text_size default 3; the size of taxa text. Only available when \code{loading_arrow = TRUE}.
		#' @param loading_text_prefix default FALSE; whether show the prefix (e.g., g__) in the taxa text. Only available when \code{loading_arrow = TRUE}.
		#' @param loading_text_italic default FALSE; whether using italic for the taxa text. Only available when \code{loading_arrow = TRUE}.
		#' @return \code{ggplot}.
		#' @examples
		#' t1$plot_ordination(plot_type = "point")
		#' t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "point")
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "ellipse"))
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "centroid"), 
		#' 	  centroid_segment_linetype = 1)
		plot_ordination = function(
			plot_type = "point", 
			choices = c(1, 2), 
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			plot_color = NULL,
			plot_shape = NULL,
			plot_group_order = NULL,
			add_sample_label = NULL,
			point_size = 3,
			point_alpha = 0.8,
			centroid_segment_alpha = 0.6,
			centroid_segment_size = 1,
			centroid_segment_linetype = 3,
			ellipse_chull_fill = TRUE,
			ellipse_chull_alpha = 0.1,
			ellipse_level = 0.9,
			ellipse_type = "t",
			NMDS_stress_pos = c(1, 1),
			NMDS_stress_text_prefix = "",
			loading_arrow = FALSE,
			loading_taxa_num = 10, 
			loading_text_taxlevel = NULL,
			loading_text_color = "black",
			loading_arrow_color = "grey30",
			loading_text_size = 3,
			loading_text_prefix = FALSE,
			loading_text_italic = FALSE
			){
			ordination_method <- self$ordination_method
			if(is.null(ordination_method)){
				stop("Please first run cal_ordination function !")
			}
			if(is.null(plot_color)){
				if(any(c("ellipse", "chull", "centroid") %in% plot_type)){
					stop("Plot ellipse or chull or centroid need groups! Please provide plot_color parameter!")
				}
			}
			if(! all(plot_type %in% c("point", "ellipse", "chull", "centroid"))){
				message("There maybe a typo in the input plot_type! plot_type should be one or more of 'point', 'ellipse', 'chull' and 'centroid'!")
			}
			combined <- self$res_ordination$scores
			eig <- self$res_ordination$eig
			model <- self$res_ordination$model
			if(length(choices) > 2){
				stop("The maximum input length of choices parameter should be 2!")
			}
			if(max(choices) > self$res_ordination$ncomp){
				stop("Maximum of input choices is larger than the dimension, please try to enlarge the ncomp parameter in the cal_ordination function!")
			}
			plot_x <- colnames(self$res_ordination$scores)[choices[1]]
			plot_y <- colnames(self$res_ordination$scores)[choices[2]]
			if(!is.null(plot_group_order)){
				combined[, plot_color] %<>% factor(., levels = plot_group_order)
			}
			if(!is.null(plot_color)){
				color_values <- expand_colors(color_values, length(unique(combined[, plot_color])))
			}			
			
			if(is.numeric(point_size)){
				p <- ggplot(combined, aes_meco(x = plot_x, y = plot_y, colour = plot_color, shape = plot_shape))
				if("point" %in% plot_type){
					p <- p + geom_point(alpha = point_alpha, size = point_size)
				}
			}else{
				check_table_variable(combined, point_size, "point_size", "res_ordination$scores")
				p <- ggplot(combined, aes_meco(x = plot_x, y = plot_y, colour = plot_color, shape = plot_shape, size = point_size))
				if("point" %in% plot_type){
					p <- p + geom_point(alpha = point_alpha)
				}
			}
			if(ordination_method %in% c("PCA", "PCoA", "DCA", "PLS-DA", "OPLS-DA")){
				p <- p + xlab(paste(plot_x, " [", eig[plot_x],"%]", sep = "")) + 
					ylab(paste(plot_y, " [", eig[plot_y],"%]", sep = ""))
			}
			if(!is.null(NMDS_stress_pos)){
				if(ordination_method == "NMDS"){
					p <- p + annotate("text", x = max(combined[, 1]) * NMDS_stress_pos[1], y = max(combined[, 2]) * NMDS_stress_pos[2], 
						label = paste0(NMDS_stress_text_prefix, round(model$stress, 2)), parse = TRUE)
				}
			}
			if("centroid" %in% plot_type){
				centroid_xy <- data.frame(group = combined[, plot_color], x = combined[, plot_x], y = combined[, plot_y]) %>%
					dplyr::group_by(group) %>%
					dplyr::summarise(cx = mean(x), cy = mean(y)) %>%
					as.data.frame()
				combined_centroid_xy <- merge(combined, centroid_xy, by.x = plot_color, by.y = "group")
				p <- p + geom_segment(
					data = combined_centroid_xy, 
					aes_meco(x = plot_x, xend = "cx", y = plot_y, yend = "cy", colour = plot_color),
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
				mapping <- aes_meco(x = plot_x, y = plot_y, group = plot_color, colour = plot_color, fill = ellipse_chull_fill_color)
				if("ellipse" %in% plot_type){
					p <- p + ggplot2::stat_ellipse(
						mapping = mapping, 
						data = combined, 
						level = ellipse_level, 
						type = ellipse_type, 
						alpha = ellipse_chull_alpha, 
						geom = "polygon"
						)
				}
				if("chull" %in% plot_type){
					p <- p + ggpubr::stat_chull(
						mapping = mapping, 
						data = combined, 
						alpha = ellipse_chull_alpha,
						geom = "polygon"
						)
				}
				if(ellipse_chull_fill){
					p <- p + scale_fill_manual(values = color_values)
				}
			}
			if(!is.null(add_sample_label)){
				p <- p + ggrepel::geom_text_repel(aes_meco(label = add_sample_label))
			}
			if(!is.null(plot_color)){
				p <- p + scale_color_manual(values = color_values)
			}
			if(!is.null(plot_shape)){
				p <- p + scale_shape_manual(values = shape_values)
			}
			if(loading_arrow){
				if(! is.null(self$res_ordination$loading)){
					df_arrows <- self$res_ordination$loading[1:loading_taxa_num, ]
					colnames(df_arrows)[choices] <- c("x", "y")
					p <- p + geom_segment(
						data = df_arrows, 
						aes(x = 0, y = 0, xend = x, yend = y), 
						arrow = arrow(length = unit(0.2, "cm")), 
						color = loading_arrow_color, 
						alpha = .6,
						inherit.aes = FALSE
						)
					if(is.null(loading_text_taxlevel)){
						if(is.null(self$res_ordination$taxa_level)){
							df_arrows$label <- rownames(df_arrows)
						}else{
							df_arrows$label <- self$dataset$tax_table[rownames(df_arrows), self$res_ordination$taxa_level]
						}
					}else{
						df_arrows$label <- self$dataset$tax_table[rownames(df_arrows), loading_text_taxlevel]
					}
					if(! loading_text_prefix){
						df_arrows$label %<>% gsub("^.__", "", .)
					}
					if(loading_text_italic){
						df_arrows$label %<>% paste0("italic('", .,"')")
						loading_text_parse <- TRUE
					}else{
						loading_text_parse <- FALSE
					}
					p <- p + ggrepel::geom_text_repel(
						data = df_arrows, 
						aes_meco("x", "y", label = "label"), 
						size = loading_text_size, 
						color = loading_text_color, 
						segment.alpha = .01, 
						parse = loading_text_parse,
						inherit.aes = FALSE
					)
				}
			}
			p
		},
		#' @description
		#' Calculate perMANOVA (Permutational Multivariate Analysis of Variance) based on the \code{adonis2} function of vegan package <doi:10.1111/j.1442-9993.2001.01070.pp.x>.
		#'
		#' @param manova_all default TRUE; TRUE represents test for all the groups, i.e. the overall test;
		#'    FALSE represents test for all the paired groups.
		#' @param manova_set default NULL; other specified group set for manova, such as \code{"Group + Type"} and \code{"Group*Type"}.
		#'    Please also see the \code{formula} parameter (only right-hand side) in \code{adonis2} function of vegan package.
		#'    The parameter manova_set has higher priority than manova_all parameter. If manova_set is provided; manova_all is disabled.
		#' @param group default NULL; a column name of \code{sample_table} used for manova. If NULL, search \code{group} variable stored in the object.
		#'    Available when \code{manova_set} is not provided.
		#' @param by_group default NULL; one column name in \code{sample_table}; used to perform paired comparisions within each group. 
		#'    Only available when \code{manova_all = FALSE} and \code{manova_set} is not provided.
		#' @param p_adjust_method default "fdr"; p.adjust method; available when \code{manova_all = FALSE}; 
		#'    see \code{method} parameter of \code{p.adjust} function for available options.
		#' @param by default "terms"; same with the \code{by} parameter in \code{adonis2} function of vegan package. 
		#' @param permutations default 999; same with the \code{permutations} parameter in \code{adonis2} function of vegan package. 
		#' @param ... parameters passed to \code{adonis2} function of \code{vegan} package.
		#' @return \code{res_manova} stored in object with \code{data.frame} class.
		#' @examples
		#' t1$cal_manova(manova_all = TRUE)
		cal_manova = function(
			manova_all = TRUE,
			manova_set = NULL,
			group = NULL,
			by_group = NULL,
			p_adjust_method = "fdr",
			by = "terms",
			permutations = 999,
			...
			){
			if(is.null(self$use_matrix)){
				stop("Please recreate the object and set the parameter measure !")
			}
			use_matrix <- self$use_matrix
			metadata <- self$sample_table
			if(!is.null(manova_set)){
				use_formula <- reformulate(manova_set, substitute(as.dist(use_matrix)))
				res <- adonis2(use_formula, data = metadata, by = by, permutations = permutations, ...)
			}else{
				if(is.null(group)){
					if(is.null(self$group)){
						stop("Please provide the group parameter!")
					}else{
						group <- self$group
					}
				}else{
					check_table_variable(metadata, group, "group", "sample_table")
				}
				if(manova_all){
					use_formula <- reformulate(group, substitute(as.dist(use_matrix)))
					res <- adonis2(use_formula, data = metadata, by = by, permutations = permutations, ...)
				}else{
					res <- private$paired_manova_anosim_bygroup(
						by_group = by_group,
						test = "permanova",
						sample_info_use = metadata, 
						use_matrix = use_matrix, 
						group = group, 
						measure = self$measure, 
						p_adjust_method = p_adjust_method,
						permutations = permutations,
						...
					)
				}
			}
			if(inherits(res, "anova")){
				res %<>% as.data.frame
				res$Significance <- generate_p_siglabel(res$`Pr(>F)`)
			}
			self$res_manova <- res
			message('The result is stored in object$res_manova ...')
			invisible(self)
		},
		#' @description
		#' Analysis of similarities (ANOSIM) based on the \code{anosim} function of vegan package.
		#'
		#' @param paired default FALSE; whether perform paired test between any two combined groups from all the input groups.
		#' @param group default NULL; a column name of \code{sample_table}. If NULL, search \code{group} variable stored in the object.
		#' @param by_group default NULL; one column name in \code{sample_table}; used to perform paired comparisions within each group. 
		#'    Only available when \code{paired = TRUE}.
		#' @param p_adjust_method default "fdr"; p.adjust method; available when \code{paired = TRUE}; see method parameter of \code{p.adjust} function for available options.
		#' @param permutations default 999; same with the \code{permutations} parameter in \code{anosim} function of vegan package. 
		#' @param ... parameters passed to \code{anosim} function of \code{vegan} package.
		#' @return \code{res_anosim} stored in object with \code{data.frame} class.
		#' @examples
		#' t1$cal_anosim()
		cal_anosim = function(
			paired = FALSE,
			group = NULL,
			by_group = NULL,
			p_adjust_method = "fdr",
			permutations = 999,
			...
			){
			if(is.null(self$use_matrix)){
				stop("Please recreate the object and set the parameter measure!")
			}
			use_matrix <- self$use_matrix
			metadata <- self$sample_table

			if(is.null(group)){
				if(is.null(self$group)){
					stop("Please provide the group parameter!")
				}else{
					group <- self$group
				}
			}else{
				check_table_variable(metadata, group, "group", "sample_table")
			}
			if(paired){
				res <- private$paired_manova_anosim_bygroup(
					by_group = by_group,
					test = "anosim",
					sample_info_use = metadata, 
					use_matrix = use_matrix, 
					group = group, 
					measure = self$measure, 
					p_adjust_method = p_adjust_method,
					permutations = permutations,
					...
				)
			}else{
				tmp <- anosim(
					x = use_matrix, 
					grouping = metadata[, group], 
					permutations = permutations,
					...
				)
				res <- data.frame(Test = "ANOSIM for all groups", permutations = tmp$permutations, statistic.R = tmp$statistic, p.value = tmp$signif)
				res$Significance <- generate_p_siglabel(res$p.value)
			}
			self$res_anosim <- res
			message('The original result is stored in object$res_anosim ...')
			invisible(self)
		},
		#' @description
		#' Multivariate homogeneity test of groups dispersions (PERMDISP) based on \code{betadisper} function in vegan package.
		#'
		#' @param ... parameters passed to \code{betadisper} function.
		#' @return \code{res_betadisper} stored in object.
		#' @examples
		#' t1$cal_betadisper()
		cal_betadisper = function(...){
			if(is.null(self$use_matrix)){
				stop("Please recreate the object and set the parameter measure !")
			}
			use_matrix <- self$use_matrix
			res1 <- betadisper(as.dist(use_matrix), self$sample_table[, self$group], ...)
			res2 <- permutest(res1, pairwise = TRUE)
			self$res_betadisper <- res2
			message('The result is stored in object$res_betadisper ...')
			invisible(self)
		},
		#' @description
		#' Convert symmetric distance matrix to distance table of paired samples that are within groups or between groups.
		#' 
		#' @param within_group default TRUE; whether obtain distance table of paired samples within groups; if FALSE, obtain distances of paired samples between any two groups.
		#' @param by_group default NULL; one colname name of \code{sample_table} in \code{microtable} object.
		#'   If provided, transform distances by the provided \code{by_group} parameter. This is especially useful for ordering and filtering values further.
		#'   When \code{within_group = TRUE}, the result of by_group parameter is the format of paired groups.
		#'   When \code{within_group = FALSE}, the result of by_group parameter is the format same with the group information in \code{sample_table}.
		#' @param ordered_group default NULL; a vector representing the ordered elements of \code{group} parameter; only useful when within_group = FALSE.
		#' @param sep default TRUE; a character string to separate the group names after merging them into a new name.
		#' @return \code{res_group_distance} stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_group_distance(within_group = TRUE)
		#' }
		cal_group_distance = function(within_group = TRUE, by_group = NULL, ordered_group = NULL, sep = " vs "){
			if(!is.null(by_group)){
				if(!all(by_group %in% colnames(self$sample_table))){
					stop("Input by_group parameter must be colnames of sample_table in the microtable object!")
				}
			}
			if(is.null(self$group)){
				stop("The group inside the object is NULL! ",
					"Please provide the group parameter when creating the trans_beta object!")
			}
			if(within_group){
				res <- private$within_group_distance(distance = self$use_matrix, sampleinfo = self$sample_table, type = self$group, by_group = by_group, sep = sep)
			}else{
				res <- private$between_group_distance(distance = self$use_matrix, sampleinfo = self$sample_table, type = self$group, by_group = by_group, 
					ordered_group = ordered_group, sep = sep)
			}
			colnames(res)[colnames(res) == "value"] <- "Value"
			self$res_group_distance <- res
			message('The result is stored in object$res_group_distance ...')
			invisible(self)
		},
		#' @description
		#' Differential test of converted distances across groups.
		#'
		#' @param group default NULL; a column name of \code{object$res_group_distance} used for the statistics; If NULL, use the \code{group} inside the object.
		#' @param by_group default NULL; a column of \code{object$res_group_distance} used to perform the differential test 
		#'   among elements in \code{group} parameter for each element in \code{by_group} parameter. So \code{by_group} has a larger scale than \code{group} parameter.
		#'   This \code{by_group} is very different from the \code{by_group} parameter in the \code{cal_group_distance} function.
		#' @param by_ID default NULL; a column of \code{object$res_group_distance} used to perform paired t test or paired wilcox test for the paired data,
		#'   such as the data of plant compartments for different plant species (ID). 
		#'   So \code{by_ID} should be the smallest unit of sample collection without any repetition in it.
		#' @param ... parameters passed to \code{cal_diff} function of \code{\link{trans_alpha}} class.
		#' @return \code{res_group_distance_diff} stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_group_distance_diff()
		#' }
		cal_group_distance_diff = function(group = NULL, by_group = NULL, by_ID = NULL, ...){
			res_group_distance <- self$res_group_distance
			if(is.null(res_group_distance)){
				stop("Please first run cal_group_distance function!")
			}
			# use method in trans_alpha
			temp1 <- suppressMessages(trans_alpha$new(dataset = NULL))
			res_group_distance$Measure <- "group_distance"
			temp1$data_alpha <- res_group_distance
			if(is.null(group)){
				temp1$group <- self$group
			}else{
				temp1$group <- group
			}
			if(! temp1$group %in% colnames(res_group_distance)){
				stop("The group parameter: ", group, " is not in object$res_group_distance!")
			}
			if(!is.null(by_group)){
				if(! by_group %in% colnames(res_group_distance)){
					stop("Provided by_group parameter: ", by_group, " is not in object$res_group_distance!")
				}
				temp1$by_group <- by_group
			}
			if(!is.null(by_ID)){
				if(! by_ID %in% colnames(res_group_distance)){
					stop("Provided by_ID parameter: ", by_ID, " is not in object$res_group_distance!")
				}
				temp1$by_ID <- by_ID
			}
			suppressMessages(temp1$cal_diff(...))
			self$res_group_distance_diff <- temp1$res_diff
			self$res_group_distance_diff_tmp <- temp1
			message('The result is stored in object$res_group_distance_diff ...')
			invisible(self)
		},
		#' @description
		#' Plot the distances of paired groups within or between groups.
		#'
		#' @param plot_group_order default NULL; a vector used to order the groups in the plot.
		#' @param ... parameters (except measure) passed to \code{plot_alpha} function of \code{\link{trans_alpha}} class.
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' t1$plot_group_distance()
		#' }
		plot_group_distance = function(plot_group_order = NULL, ...){
			if(is.null(self$res_group_distance_diff)){
				if(is.null(self$res_group_distance)){
					stop("Please first run cal_group_distance function!")
				}
				group_distance <- self$res_group_distance
				group <- self$group
			}else{
				group_distance <- self$res_group_distance_diff_tmp$data_alpha
				group <- self$res_group_distance_diff_tmp$group
			}
			if(self$measure %in% c("wei_unifrac", "unwei_unifrac", "bray", "jaccard")){
				titlename <- switch(self$measure, 
					wei_unifrac = "Weighted Unifrac", 
					unwei_unifrac = "Unweighted Unifrac", 
					bray = "Bray-Curtis", 
					jaccard = "Jaccard")
				ylabname <- paste0(titlename, " distance")
			}else{
				ylabname <- self$measure
			}
			if(!is.null(plot_group_order)) {
				group_distance[, group] %<>% factor(., levels = plot_group_order)
			}else{
				if(!is.factor(group_distance[, group])){
					group_distance[, group] %<>% as.factor
				}
			}
			message("The ordered groups are ", paste0(levels(group_distance[, group]), collapse = " "), " ...")
			
			if(is.null(self$res_group_distance_diff)){
				temp1 <- suppressMessages(trans_alpha$new(dataset = NULL))
				group_distance$Measure <- "group_distance"
				temp1$data_alpha <- group_distance
				temp1$group <- group
				p <- temp1$plot_alpha(add_sig = FALSE, measure = "group_distance", ...) + ylab(ylabname)
			}else{
				# reassign res_diff for the case of customized manipulation on the object
				self$res_group_distance_diff_tmp$res_diff <- self$res_group_distance_diff
				# reassign group_distance for factors
				self$res_group_distance_diff_tmp$data_alpha <- group_distance
				p <- self$res_group_distance_diff_tmp$plot_alpha(measure = "group_distance", ...) + ylab(ylabname)
			}
			p
		},
		#' @description
		#' Plot clustering result based on the \code{ggdendro} package.
		#'
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); color palette for the text.
		#' @param measure default NULL; beta diversity index; If NULL, using the measure when creating object
		#' @param group default NULL; if provided, use this group to assign color.
		#' @param replace_name default NULL; if provided, use this as label.
		#' @return \code{ggplot}.
		#' @examples
		#' t1$plot_clustering(group = "Group", replace_name = c("Saline", "Type"))
		plot_clustering = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			measure = NULL, 
			group = NULL, 
			replace_name = NULL
			){
			dataset <- self$dataset
			if(is.null(measure)){
				if(is.null(self$use_matrix)){
					measure_matrix <- dataset$beta_diversity[[1]]
					measure <- names(dataset$beta_diversity)[1]
				}else{
					measure_matrix <- self$use_matrix
					measure <- self$measure
				}
			}else{
				measure_matrix <- dataset$beta_diversity[[measure]]
			}
			hc_measure <- hclust(as.dist(measure_matrix))
			hc_d_measure <- ggdendro::dendro_data(as.dendrogram(hc_measure))
			titlename <- switch(measure, wei_unifrac = "Weighted Unifrac", unwei_unifrac = "Unweighted Unifrac", bray = "Bray-Curtis", jaccard = "Jaccard")
			ylabname <- paste0("Distance (", titlename, ")")

			g1 <- ggplot(data = ggdendro::segment(hc_d_measure)) + 
				geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color = "grey30")
			if(!is.null(group) | !is.null(replace_name)){
				data2 <- suppressWarnings(dplyr::left_join(hc_d_measure$label, rownames_to_column(self$sample_table), by = c("label" = "rowname")))
				if(length(replace_name) > 1){
					data2$replace_name_use <- apply(data2[, replace_name], 1, function(x){paste0(x, collapse = "-")})
				}
			}
			if(is.null(group)){
				if(is.null(replace_name)){
					g1 <- g1 + geom_text(data=hc_d_measure$label, aes(x=x, y=y, label=label, hjust=-0.1), size=4)
				}else{
					if(length(replace_name) > 1){
						g1 <- g1 + geom_text(data=data2, aes_meco(x="x", y="y", label = "replace_name_use", hjust=-0.1), size=4)
					}else{
						g1 <- g1 + geom_text(data=data2, aes_meco(x="x", y="y", label = replace_name, hjust=-0.1), size=4)
					}
				}
			} else {
				if(is.null(replace_name)){
					g1 <- g1 + geom_text(data=data2, aes_meco(x="x", y="y", label="label", hjust=-0.1, colour = group), size=4)
				}else{
					if(length(replace_name) > 1){
						g1 <- g1 + geom_text(data=data2, aes_meco(x="x", y="y", label="replace_name_use", hjust=-0.1, colour = group), size=4)
					}else{
						g1 <- g1 + geom_text(data=data2, aes_meco(x="x", y="y", label=replace_name, hjust=-0.1, colour = group), size=4)
					}
				}
				g1 <- g1 + scale_color_manual(values = color_values)
			}
			g1 <- g1 + theme(legend.position="none") + coord_flip() +
				scale_x_discrete(labels=ggdendro::label(hc_d_measure)$label) +
				ylab(ylabname) +
				scale_y_reverse(expand=c(0.3, 0)) + 
				xlim(min(ggdendro::segment(hc_d_measure)[,1]) - 0.3, max(ggdendro::segment(hc_d_measure)[,1]) + 0.3) +
				theme(axis.line.y=element_blank(),
					  axis.ticks.y=element_blank(),
					  axis.text.y=element_blank(),
					  axis.title.y=element_blank(),
					  panel.background=element_rect(fill="white"),
					  panel.grid=element_blank(), 
					  panel.border = element_blank()) +
				theme(axis.line.x = element_line(color = "black", linetype = "solid", lineend = "square"))
			g1
		}
		),
	private = list(
		within_group_distance = function(distance, sampleinfo, type, by_group = NULL, sep = " vs "){
			all_group <- as.character(sampleinfo[, type]) %>% unique
			res <- data.frame()
			for(i in all_group){
				select_sample <- sampleinfo[, type] == i
				# filter group with sample number < 2
				if(sum(select_sample) < 2){
					next
				}
				distance_res <- distance[select_sample, select_sample] %>% as.dist %>% as.vector
				distance_res <- data.frame(value = distance_res, group = i)
				if(!is.null(by_group)){
					for(j in by_group){
						tmp <- sampleinfo[select_sample, j]
						paired_comb <- combn(length(tmp), 2)
						merged_bygroup <- lapply(seq_len(ncol(paired_comb)), function(x){paste0(sort(tmp[paired_comb[, x]]), collapse = sep)}) %>% unlist
						distance_res %<>% cbind(., merged_bygroup)
					}
				}
				res %<>% rbind(., distance_res)
			}
			colnames(res)[2] <- type
			if(ncol(res) > 2){
				colnames(res)[3:ncol(res)] <- by_group
			}
			res
		},
		between_group_distance = function(distance, sampleinfo, type, by_group = NULL, ordered_group = NULL, sep = " vs "){
			all_group <- as.character(sampleinfo[, type]) %>% unique
			# first check ordered_group
			if(!is.null(ordered_group)){
				ordered_group %<>% as.character
				if(!all(all_group %in% ordered_group)){
					stop("Please check the ordered_group! Part of groups are missing!")
				}
			}
			com1 <- combn(all_group, 2)
			res <- data.frame()
			for(i in seq_len(ncol(com1))){
				if(is.null(by_group)){
					f_name <- rownames(sampleinfo[sampleinfo[, type] == com1[1, i], ])
					s_name <- rownames(sampleinfo[sampleinfo[, type] == com1[2, i], ])
					if(!is.null(ordered_group)){
						vsname <- paste0(ordered_group[ordered_group %in% c(com1[1, i], com1[2, i])], collapse = sep)
					}else{
						vsname <- paste0(com1[1, i], sep, com1[2, i])
					}
					distance_res <- as.vector(distance[f_name, s_name])
					distance_res <- data.frame(Value = distance_res, vsname)
				}else{
					if(length(by_group) > 1){
						stop("The length of by_group must be one!")
					}
					tmp <- sampleinfo %>% dropallfactors %>% .[, by_group] %>% unique
					distance_res <- data.frame()
					for(j in tmp){
						f_name <- rownames(sampleinfo[sampleinfo[, type] == com1[1, i] & sampleinfo[, by_group] == j, ])
						s_name <- rownames(sampleinfo[sampleinfo[, type] == com1[2, i] & sampleinfo[, by_group] == j, ])
						if(!is.null(ordered_group)){
							vsname <- paste0(ordered_group[ordered_group %in% c(com1[1, i], com1[2, i])], collapse = sep)
						}else{
							vsname <- paste0(com1[1, i], sep, com1[2, i])
						}
						tmp_dis <- as.vector(distance[f_name, s_name])
						if(length(tmp_dis) == 0){
							next
						}else{
							tmp_dis <- data.frame(Value = tmp_dis, name = vsname)
							tmp_dis %<>% cbind(., by_group = j)
							distance_res %<>% rbind(., tmp_dis)
						}
					}
				}
				res %<>% rbind(., distance_res)
			}
			if(nrow(res) == 0){
				if(!is.null(by_group)){
					stop("No result is obtained! Please check the input by_group parameter! ", 
						"This probably comes from the wrong names in ", by_group, " of the sample information table! ",
						"Wrong names can lead to no overlap among ", type, " for each element of ", by_group, "!")
				}
			}
			colnames(res)[2] <- type
			if(ncol(res) > 2){
				colnames(res)[3:ncol(res)] <- by_group
			}
			res
		},
		paired_manova_anosim_bygroup = function(by_group, test, sample_info_use, use_matrix, group, measure, p_adjust_method, permutations, ...){
			if(is.null(by_group)){
				res <- private$paired_group_manova_anosim(
					test = test,
					sample_info_use = sample_info_use, 
					use_matrix = use_matrix, 
					group = group, 
					measure = measure, 
					p_adjust_method = p_adjust_method,
					permutations = permutations,
					...
				)
			}else{
				res <- data.frame()
				check_table_variable(sample_info_use, by_group, "by_group", "sample_table")
				all_bygroups <- unique(sample_info_use[, by_group])
				for(i in all_bygroups){
					message("For by_group: ", i, " ...")
					sub_meta <- sample_info_use[sample_info_use[, by_group] == i, ]
					if(length(unique(sub_meta[, group])) < 2){
						message("Skip by_group: ", i, ", because groups number < 2 ...")
						next
					}
					tmp <- private$paired_group_manova_anosim(
						test = test,
						sample_info_use = sub_meta, 
						use_matrix = use_matrix, 
						group = group, 
						measure = measure, 
						p_adjust_method = p_adjust_method,
						permutations = permutations,
						...
					)
					tmp <- data.frame(by_group = i, tmp)
					res <- rbind(res, tmp)
				}
			}
			res
		},
		paired_group_manova_anosim = function(test, sample_info_use, use_matrix, group, measure, p_adjust_method, permutations, ...){
			comnames <- c()
			test <- match.arg(test, choices = c("permanova", "anosim"))
			if(test == "permanova"){
				F <- c()
				R2 <- c()
			}else{
				R <- c()
			}
			p_value <- c()
			matrix_total <- use_matrix[rownames(sample_info_use), rownames(sample_info_use)]
			groupvec <- as.character(sample_info_use[, group])
			all_name <- combn(unique(sample_info_use[, group]), 2)
			for(i in 1:ncol(all_name)) {
				matrix_compare <- matrix_total[groupvec %in% as.character(all_name[,i]), groupvec %in% as.character(all_name[,i])]
				sample_info_compare <- sample_info_use[groupvec %in% as.character(all_name[,i]), , drop = FALSE]
				comnames <- c(comnames, paste0(as.character(all_name[,i]), collapse = " vs "))
				if(test == "permanova"){
					tmp_result <- adonis2(reformulate(group, substitute(as.dist(matrix_compare))), data = sample_info_compare, permutations = permutations, ...)
					F %<>% c(., tmp_result$F[1])
					R2 %<>% c(., tmp_result$R2[1])
					p_value %<>% c(., tmp_result$`Pr(>F)`[1])
				}else{
					tmp_result <- anosim(x = matrix_compare, grouping = sample_info_compare[, group], permutations = permutations, ...)
					R %<>% c(., tmp_result$statistic[1])
					p_value %<>% c(., tmp_result$signif[1])
				}
			}
			p_adjusted <- p.adjust(p_value, method = p_adjust_method)
			significance_label <- generate_p_siglabel(p_adjusted)
			measure_vec <- rep(measure, length(comnames))
			if(test == "permanova"){
				compare_result <- data.frame(comnames, measure_vec, F, R2, p_value, p_adjusted, significance_label)
				colnames(compare_result) <- c("Groups", "measure", "F", "R2","p.value", "p.adjusted", "Significance")
			}else{
				compare_result <- data.frame(comnames, measure_vec, R, p_value, p_adjusted, significance_label)
				colnames(compare_result) <- c("Groups", "measure", "R", "p.value", "p.adjusted", "Significance")
			}
			compare_result
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
