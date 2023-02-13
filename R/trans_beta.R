#' @title Create \code{trans_beta} object for beta-diversity analysis based on the distance matrix
#'
#' @description
#' This class is a wrapper for a series of beta-diversity related analysis, 
#' including several ordination calculations and plotting based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>, group distance comparision, 
#' clustering, perMANOVA based on Anderson al. (2008) <doi:10.1111/j.1442-9993.2001.01070.pp.x> and PERMDISP.
#'
#' @export
trans_beta <- R6Class(classname = "trans_beta",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param measure default NULL; bray, jaccard, wei_unifrac or unwei_unifrac, or other name of matrix you add in \code{microtable$beta_diversity}; 
		#' 	 used for ordination, manova or group distance. The measure must be one of names of microtable$beta_diversity list. 
		#' 	 Please see \code{microtable$cal_betadiv} function for more details.
		#' @param group default NULL; sample group used for manova, betadisper or group distance.
		#' @return parameters stored in the object.
		#' @examples
		#' data(dataset)
		#' t1 <- trans_beta$new(dataset = dataset, measure = "bray", group = "Group")
		initialize = function(
			dataset = NULL, 
			measure = NULL, 
			group = NULL
			) {
			if(is.null(dataset)){
				stop("dataset is necessary !")
			}
			if(!is.null(measure)){
				if(is.null(dataset$beta_diversity)){
					stop("No dataset$beta_diversity found! Please first use microtable$cal_betadiv to calculate ", measure, " !")
				}else{
					if(length(measure) > 1){
						stop("The input measure should only have one element! Please check it!")
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
						}else{
							stop("Unknown format of input measure parameter!")
						}
					}
					self$use_matrix <- dataset$beta_diversity[[measure]]
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
		#' Ordination based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>.
		#'
		#' @param ordination default "PCoA"; "PCA", "PCoA" or "NMDS". PCA: principal component analysis; 
		#' 	  PCoA: principal coordinates analysis; NMDS: non-metric multidimensional scaling.
		#' @param ncomp default 3; the returned dimensions.
		#' @param trans_otu default FALSE; whether species abundance will be square transformed, used for PCA.
		#' @param scale_species default FALSE; whether species loading in PCA will be scaled.
		#' @return \code{res_ordination} stored in the object.
		#' @examples
		#' t1$cal_ordination(ordination = "PCoA")		
		cal_ordination = function(
			ordination = "PCoA",
			ncomp = 3,
			trans_otu = FALSE, 
			scale_species = FALSE		
			){
			if(is.null(ordination)){
				stop("Input ordination should not be NULL !")
			}
			if(!ordination %in% c("PCA", "PCoA", "NMDS")){
				stop("Input ordination should be one of 'PCA', 'PCoA' and 'NMDS' !")
			}
			dataset <- self$dataset
			if(ordination == "PCA"){
				plot.x <- "PC1"
				plot.y <- "PC2"
				if(trans_otu == T){
					abund1 <- sqrt(dataset$otu_table)
				}else{
					abund1 <- dataset$otu_table
				}
				model <- rda(t(abund1))
				expla <- round(model$CA$eig/model$CA$tot.chi*100,1)
				scores <- scores(model, choices = 1:ncomp)$sites
				combined <- cbind.data.frame(scores, dataset$sample_table)

				if(is.null(dataset$tax_table)){
					loading <- scores(model, choices = 1:ncomp)$species
				}else{
					loading <- cbind.data.frame(scores(model, choices = 1:ncomp)$species, dataset$tax_table)
				}
				loading <- cbind.data.frame(loading, rownames(loading))

				if(scale_species == T){
					maxx <- max(abs(scores[,plot.x]))/max(abs(loading[,plot.x]))
					loading[, plot.x] <- loading[, plot.x] * maxx * 0.8
					maxy <- max(abs(scores[,plot.y]))/max(abs(loading[,plot.y]))
					loading[, plot.y] <- loading[, plot.y] * maxy * 0.8
				}

				species <- cbind(loading, loading[,plot.x]^2 + loading[,plot.y]^2)
				colnames(species)[ncol(species)] <- "dist"
				species <- species[with(species, order(-dist)), ]
				outlist <- list(model = model, scores = combined, loading = species, eig = expla)
			}
			if(ordination %in% c("PCoA", "NMDS")){
				if(is.null(self$use_matrix)){
					stop("Please recreate the object and set the parameter measure !")
				}
			}
			if(ordination == "PCoA"){
				model <- ape::pcoa(as.dist(self$use_matrix))
				combined <- cbind.data.frame(model$vectors[,1:ncomp], dataset$sample_table)
				pco_names <- paste0("PCo", 1:10)
				colnames(combined)[1:ncomp] <- pco_names[1:ncomp]
				expla <- round(model$values[,1]/sum(model$values[,1])*100, 1)
				expla <- expla[1:ncomp]
				names(expla) <- pco_names[1:ncomp]
				outlist <- list(model = model, scores = combined, eig = expla)
			}
			if(ordination == "NMDS"){
				model <- vegan::metaMDS(as.dist(self$use_matrix))
				combined <- cbind.data.frame(model$points, dataset$sample_table)
				outlist <- list(model = model, scores = combined)
			}
			self$res_ordination <- outlist
			message('The ordination result is stored in object$res_ordination ...')
			self$ordination <- ordination
		},
		#' @description
		#' Plotting the ordination result based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>.
		#'
		#' @param plot_type default "point"; one or more elements of "point", "ellipse", "chull" and "centroid".
		#'   \describe{
		#'     \item{\strong{'point'}}{add point}
		#'     \item{\strong{'ellipse'}}{add confidence ellipse for points of each group}
		#'     \item{\strong{'chull'}}{add convex hull for points of each group}
		#'     \item{\strong{'centroid'}}{add centroid line of each group}
		#'   }
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for different groups.
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); a vector for point shape types of groups, see \code{ggplot2} tutorial.
		#' @param plot_color default NULL; a colname of sample_table to assign colors to different groups in plot.
		#' @param plot_shape default NULL; a colname of sample_table to assign shapes to different groups in plot.
		#' @param plot_group_order default NULL; a vector used to order the groups in the legend of plot.
		#' @param add_sample_label default NULL; the column name in sample table, if provided, show the point name in plot.
		#' @param point_size default 3; point size in plot when "point" is in plot_type.
		#' @param point_alpha default .8; point transparency in plot when "point" is in plot_type.
		#' @param centroid_segment_alpha default 0.6; segment transparency in plot when "centroid" is in plot_type.
		#' @param centroid_segment_size default 1; segment size in plot when "centroid" is in plot_type.
		#' @param centroid_segment_linetype default 3; the line type related with centroid in plot when "centroid" is in plot_type.
		#' @param ellipse_chull_fill default TRUE; whether fill colors to the area of ellipse or chull.
		#' @param ellipse_chull_alpha default 0.1; color transparency in the ellipse or convex hull depending on whether "ellipse" or "centroid" is in plot_type.
		#' @param ellipse_level default .9; confidence level of ellipse when "ellipse" is in plot_type.
		#' @param ellipse_type default "t"; ellipse type when "ellipse" is in plot_type; see type in \code{\link{stat_ellipse}}.
		#' @return \code{ggplot}.
		#' @examples
		#' t1$plot_ordination(plot_type = "point")
		#' t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "point")
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "ellipse"))
		#' t1$plot_ordination(plot_color = "Group", plot_type = c("point", "centroid"), 
		#' 	  centroid_segment_linetype = 1)
		plot_ordination = function(
			plot_type = "point",
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
			ellipse_type = "t"
			){
			ordination <- self$ordination
			if(is.null(ordination)){
				stop("Please first run cal_ordination function !")
			}
			if(is.null(plot_color)){
				if(any(c("ellipse", "chull", "centroid") %in% plot_type)){
					stop("Plot ellipse or chull or centroid need groups! Please provide plot_color parameter!")
				}
			}
			if(! all(plot_type %in% c("point", "ellipse", "chull", "centroid"))){
				message("There maybe a typo in your plot_type input! plot_type should be one or more from 'point', 'ellipse', 'chull' and 'centroid'!")
			}
			combined <- self$res_ordination$scores
			eig <- self$res_ordination$eig
			model <- self$res_ordination$model
			plot_x <- colnames(self$res_ordination$scores)[1]
			plot_y <- colnames(self$res_ordination$scores)[2]
			
			if(!is.null(plot_group_order)){
				combined[, plot_color] %<>% factor(., levels = plot_group_order)
			}
			
			p <- ggplot(combined, aes_meco(x = plot_x, y = plot_y, colour = plot_color, shape = plot_shape))
			if("point" %in% plot_type){
				p <- p + geom_point(alpha = point_alpha, size = point_size)
			}
			if(ordination %in% c("PCA", "PCoA")){
				p <- p + xlab(paste(plot_x, " [", eig[plot_x],"%]", sep = "")) + 
					ylab(paste(plot_y, " [", eig[plot_y],"%]", sep = ""))
			}
			if(ordination == "NMDS"){
				p <- p + annotate("text", x = max(combined[,1]), y = max(combined[,2]) + 0.05, label = round(model$stress, 2), parse=TRUE)
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
			p
		},
		#' @description
		#' Calculate perMANOVA based on <doi:10.1111/j.1442-9993.2001.01070.pp.x> and R vegan \code{adonis2} function.
		#'
		#' @param manova_all default TRUE; TRUE represents test for all the groups, i.e. the overall test;
		#'    FALSE represents test for all the paired groups.
		#' @param manova_set default NULL; other specified group set for manova, such as \code{"Group + Type"} and \code{"Group*Type"}; see also \code{\link{adonis2}}.
		#'    manova_set has higher priority than manova_all parameter. If manova_set is provided; manova_all is disabled.
		#' @param group default NULL; a column name of \code{sample_table} used for manova. If NULL, search \code{group} variable stored in the object.
		#' @param p_adjust_method default "fdr"; p.adjust method when \code{manova_all = FALSE}; see method parameter of \code{p.adjust} function for available options.
		#' @param ... parameters passed to \code{\link{adonis2}} function of \code{vegan} package.
		#' @return \code{res_manova} stored in object.
		#' @examples
		#' t1$cal_manova(manova_all = TRUE)
		cal_manova = function(
			manova_all = TRUE,
			manova_set = NULL,
			group = NULL,
			p_adjust_method = "fdr",
			...
			){
			if(is.null(self$use_matrix)){
				stop("Please recreate the object and set the parameter measure !")
			}
			use_matrix <- self$use_matrix
			metadata <- self$sample_table
			if(!is.null(manova_set)){
				use_formula <- reformulate(manova_set, substitute(as.dist(use_matrix)))
				self$res_manova <- adonis2(use_formula, data = metadata, ...)
			}else{
				if(is.null(group)){
					if(is.null(self$group)){
						stop("Please provide the group parameter!")
					}else{
						group <- self$group
					}
				}
				if(manova_all){
					use_formula <- reformulate(group, substitute(as.dist(use_matrix)))
					self$res_manova <- adonis2(use_formula, data = metadata, ...)
				}else{
					self$res_manova <- private$paired_group_manova(
						sample_info_use = metadata, 
						use_matrix = use_matrix, 
						group = group, 
						measure = self$measure, 
						p_adjust_method = p_adjust_method,
						...
					)
				}
			}
			message('The result is stored in object$res_manova ...')
		},
		#' @description
		#' A wrapper for \code{betadisper} function in vegan package for multivariate homogeneity test of groups dispersions.
		#'
		#' @param ... parameters passed to \code{\link{betadisper}} function.
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
		},
		#' @description
		#' Transform sample distances within groups or between groups.
		#'
		#' @param within_group default TRUE; whether transform sample distance within groups, if FALSE, transform sample distance between any two groups.
		#' @param by_group default NULL; one colname name of sample_table in \code{microtable} object.
		#'   If provided, transform distances by the provided by_group parameter. This is especially useful for ordering and filtering values further.
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
		},
		#' @description
		#' Differential test of distances among groups.
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
		},
		#' @description
		#' Plotting the distance between samples within or between groups.
		#'
		#' @param plot_group_order default NULL; a vector used to order the groups in the plot.
		#' @param ... parameters (except measure) passed to \code{plot_alpha} function of \code{\link{trans_alpha}} class.
		#' @return \code{ggplot}.
		#' @examples
		#' \donttest{
		#' t1$plot_group_distance()
		#' }
		plot_group_distance = function(plot_group_order = NULL, ...){
			group_distance <- self$res_group_distance
			if(!is.null(self$res_group_distance_diff)){
				group <- self$res_group_distance_diff_tmp$group
			}else{
				group <- self$group
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
			if (!is.null(plot_group_order)) {
				group_distance[, group] %<>% factor(., levels = plot_group_order)
			}else{
				group_distance[, group] %<>% as.factor
			}
			message("The ordered groups are ", paste0(levels(group_distance[, group]), collapse = " "), " ...")
			
			if(is.null(self$res_group_distance_diff)){
				temp1 <- suppressMessages(trans_alpha$new(dataset = NULL))
				group_distance$Measure <- "group_distance"
				temp1$data_alpha <- group_distance
				temp1$group <- "Group"
				p <- temp1$plot_alpha(add_sig = FALSE, measure = "group_distance", ...) + ylab(ylabname)
			}else{
				self$res_group_distance_diff_tmp$res_diff <- self$res_group_distance_diff
				p <- self$res_group_distance_diff_tmp$plot_alpha(measure = "group_distance", ...) + ylab(ylabname)
			}
			p
		},
		#' @description
		#' Plotting clustering result based on the \code{ggdendro} package.
		#'
		#' @param use_colors colors for presentation.
		#' @param measure default NULL; beta diversity index; If NULL, using the measure when creating object
		#' @param group default NULL; if provided, use this group to assign color.
		#' @param replace_name default NULL; if provided, use this as label.
		#' @return \code{ggplot}.
		#' @examples
		#' t1$plot_clustering(group = "Group", replace_name = c("Saline", "Type"))
		plot_clustering = function(
			use_colors = RColorBrewer::brewer.pal(8, "Dark2"), 
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
				g1 <- g1 + scale_color_manual(values = use_colors)
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
		paired_group_manova = function(sample_info_use, use_matrix, group, measure, p_adjust_method, ...){
			comnames <- c()
			F <- c()
			R2 <- c()
			p_value <- c()
			matrix_total <- use_matrix[rownames(sample_info_use), rownames(sample_info_use)]
			groupvec <- as.character(sample_info_use[ , group])
			all_name <- combn(unique(sample_info_use[ , group]), 2)
			for(i in 1:ncol(all_name)) {
				matrix_compare <- matrix_total[groupvec %in% as.character(all_name[,i]), groupvec %in% as.character(all_name[,i])]
				sample_info_compare <- sample_info_use[groupvec %in% as.character(all_name[,i]), ]
				ad <- adonis2(reformulate(group, substitute(as.dist(matrix_compare))), data = sample_info_compare, ...)
				comnames <- c(comnames, paste0(as.character(all_name[,i]), collapse = " vs "))
				F %<>% c(., ad$F[1])
				R2 %<>% c(., ad$R2[1])
				p_value %<>% c(., ad$`Pr(>F)`[1])
			}
			p_adjusted <- p.adjust(p_value, method = p_adjust_method)
			significance_label <- cut(p_adjusted, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))
			measure_vec <- rep(measure, length(comnames))
			compare_result <- data.frame(comnames, measure_vec, F, R2, p_value, p_adjusted, significance_label)
			colnames(compare_result) <- c("Groups", "measure", "F", "R2","p.value", "p.adjusted", "Significance")
			compare_result
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
