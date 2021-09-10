#' @title Create trans_beta object for the analysis of distance matrix of beta-diversity.
#'
#' @description
#' This class is a wrapper for a series of beta-diversity related analysis, 
#' including several ordination calculations and plotting based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>, group distance comparision, 
#' clustering, perMANOVA based on Anderson al. (2008) <doi:10.1111/j.1442-9993.2001.01070.pp.x> and PERMDISP.
#'
#' @export
trans_beta <- R6Class(classname = "trans_beta",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param measure default NULL; bray, jaccard, wei_unifrac or unwei_unifrac, or other name of matrix you add; 
		#' 	 beta diversity index used for ordination, manova or group distance.
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
				if(!measure %in% c(names(dataset$beta_diversity), 1:length(dataset$beta_diversity))){
					stop("Input measure should be one of beta_diversity distance in dataset !")
				}else{
					self$use_matrix <- dataset$beta_diversity[[measure]]
				}
			}
			message("Please also cite the original paper: An et al. (2019). Soil bacterial community structure in Chinese wetlands. Geoderma, 337, 290-299.")
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
		#' @param ordination default "PCoA"; "PCA", "PCoA" or "NMDS".
		#' @param ncomp default 3; the returned dimensions.
		#' @param trans_otu default FALSE; whether species abundance will be square transformed, used for PCA.
		#' @param scale_species default FALSE; whether species loading in PCA will be scaled.
		#' @return res_ordination stored in the object.
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
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors for presentation.
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); a vector used in the shape type, see ggplot2 tutorial.
		#' @param plot_color default NULL; the sample group name used for color in plot.
		#' @param plot_shape default NULL; the sample group name used for shape in plot.
		#' @param plot_group_order default NULL; a vector used to order the groups in the legend of plot.
		#' @param plot_point_size default 3; point size in plot.
		#' @param plot_point_alpha default .9; point transparency in plot.
		#' @param plot_sample_label default NULL; the column name in sample table, if provided, show the point name in plot.
		#' @param plot_group_centroid default FALSE; whether show the centroid in each group of plot.
		#' @param plot_group default NULL; the column name in sample table, generally used with plot_group_centroid and plot_group_ellipse.
		#' @param segment_alpha default .6; segment transparency in plot.
		#' @param centroid_linetype default 3; the line type related with centroid in plot.
		#' @param plot_group_ellipse default FALSE; whether show the confidence ellipse in each group of plot.
		#' @param ellipse_level default .9; confidence level of ellipse.
		#' @param ellipse_alpha default .1; color transparency in the ellipse.
		#' @param ellipse_type default t; see type in \code{\link{stat_ellipse}}.
		#' @return ggplot.
		#' @examples
		#' t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_group_ellipse = TRUE)
		plot_ordination = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			plot_color = NULL,
			plot_shape = NULL,
			plot_group_order = NULL,
			plot_point_size = 3,
			plot_point_alpha = .9,
			plot_sample_label = NULL,
			plot_group_centroid = FALSE,
			plot_group = NULL,
			segment_alpha = .6,
			centroid_linetype = 3,
			plot_group_ellipse = FALSE,
			ellipse_level = 0.9,
			ellipse_alpha = 0.1,
			ellipse_type = "t"
			){
			ordination <- self$ordination
			if(is.null(ordination)){
				stop("Please first run cal_ordination function !")
			}
			
			combined <- self$res_ordination$scores
			eig <- self$res_ordination$eig
			model <- self$res_ordination$model
			plot_x <- colnames(self$res_ordination$scores)[1]
			plot_y <- colnames(self$res_ordination$scores)[2]
			
			if (!is.null(plot_group_order)) {
				combined[, plot_color] %<>% factor(., levels = plot_group_order)
			}
			
			p <- ggplot(combined, aes_string(x = plot_x, y = plot_y, color = plot_color, shape = plot_shape))
			p <- p + geom_point(alpha = plot_point_alpha, size = plot_point_size)
			if(ordination %in% c("PCA", "PCoA")){
				p <- p + xlab(paste(plot_x, " [", eig[plot_x],"%]", sep = "")) + 
					ylab(paste(plot_y, " [", eig[plot_y],"%]", sep = ""))
			}
			if(ordination == "NMDS"){
				p <- p + annotate("text",x= max(combined[,1]),y=max(combined[,2]) + 0.05,label = round(model$stress, 2), parse=TRUE)
			}
			if (plot_group_centroid == T){
				if(is.null(plot_group)){
					plot_group <- plot_color
				}
				os <- data.frame(group = combined[, plot_group], x = combined[, plot_x], y = combined[, plot_y]) %>%
					dplyr::group_by(group) %>%
					dplyr::summarise(cx = mean(x), cy = mean(y)) %>%
					as.data.frame()
				os2 <- merge(combined, os, by.x= plot_group, by.y = "group")
				p <- p + geom_segment(
					data=os2, 
					aes_string(x = plot_x, xend = "cx", y = plot_y, yend = "cy", color = plot_group),
					alpha = segment_alpha, 
					size = 0.5, 
					linetype = centroid_linetype
				)
			}
			if (plot_group_ellipse == T) {
				if(is.null(plot_group)){
					plot_group <- plot_color
				}
				mapping <- aes_string(x = plot_x, y = plot_y, group = plot_group, fill = plot_group)
				p <- p + ggplot2::stat_ellipse(
					mapping = mapping, 
					data = combined, 
					level = ellipse_level, 
					type = ellipse_type, 
					alpha = ellipse_alpha, 
					geom = "polygon"
					) + 
					scale_fill_manual(values = color_values)
			}
			if(!is.null(plot_color)){
				p <- p + scale_color_manual(values = color_values)
			}
			if(!is.null(plot_shape)){
				p <- p + scale_shape_manual(values = shape_values)
			}
			if(!is.null(plot_sample_label)){
				p <- p + ggrepel::geom_text_repel(aes_string(label = plot_sample_label))
			}
			p
		},
		#' @description
		#' Calculate perMANOVA based on Anderson al. (2008) <doi:10.1111/j.1442-9993.2001.01070.pp.x> and R vegan adonis function.
		#'
		#' @param cal_manova_all default FALSE; whether manova is used for all data.
		#' @param cal_manova_paired default FALSE; whether manova is used for all the paired groups.
		#' @param cal_manova_set default NULL; specified group set for manova, see \code{\link{adonis}}.
		#' @param permutations default 999; see permutations in \code{\link{adonis}}.
		#' @return res_manova stored in object.
		#' @examples
		#' t1$cal_manova(cal_manova_all = TRUE)
		cal_manova = function(
			cal_manova_all = FALSE, 
			cal_manova_paired = FALSE, 
			cal_manova_set = NULL, 
			permutations = 999
			){
			if(is.null(self$use_matrix)){
				stop("Please recreate the object and set the parameter measure !")
			}
			use_matrix <- self$use_matrix
			if(!is.null(cal_manova_set)){
				self$res_manova <- adonis(reformulate(cal_manova_set, substitute(as.dist(use_matrix))), data = self$sample_table, permutations = permutations)
			}
			if(cal_manova_all == T){
				self$res_manova <- adonis(reformulate(self$group, substitute(as.dist(use_matrix))), data = self$sample_table, permutations = permutations)
			}
			if(cal_manova_paired == T){
				self$res_manova <- private$paired_group_manova(
					sample_info_use = self$sample_table, 
					use_matrix = use_matrix, 
					group = self$group, 
					measure = self$measure, 
					permutations = permutations
				)
			}
			message('The result is stored in object$res_manova ...')
		},
		#' @description
		#' A wrapper for betadisper function in vegan package for multivariate homogeneity test of groups dispersions.
		#'
		#' @param ... parameters passed to \code{\link{betadisper}} function.
		#' @return res_betadisper stored in object.
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
		#' @return res_group_distance stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_group_distance(within_group = TRUE)
		#' }
		cal_group_distance = function(within_group = TRUE){
			if(within_group == T){
				self$res_group_distance <- private$within_group_distance(distance = self$use_matrix, sampleinfo=self$sample_table, type = self$group)
			}else{
				self$res_group_distance <- private$between_group_distance(distance = self$use_matrix, sampleinfo=self$sample_table, type = self$group)
			}
			message('The result is stored in object$res_group_distance ...')
		},
		#' @description
		#' Plotting the distance between samples within or between groups.
		#'
		#' @param plot_group_order default NULL; a vector used to order the groups in the plot.
		#' @param color_values colors for presentation.
		#' @param distance_pair_stat default FALSE; whether do the paired comparisions.
		#' @param pair_compare_filter_match default NULL; if provided, remove the matched groups; use the regular express to match the paired groups.
		#' @param pair_compare_filter_select default NULL; numeric vector; if provided, only select those input groups.
		#'   This parameter must be a numeric vector used to select the paired combination of groups. For example, pair_compare_filter_select = c(1, 3) 
		#'   can be used to select "CW"-"IW" and "IW"-"TW" from all the three pairs "CW"-"IW", "CW"-"TW" and "IW"-"TW" of ordered groups ("CW", "IW", "TW").
		#'   The parameter pair_compare_filter_select and pair_compare_filter_match can not be both used together.
		#' @param pair_compare_method default wilcox.test; wilcox.test, kruskal.test, t.test or anova.
		#' @param plot_distance_xtype default NULL; number used to make x axis text generate angle.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_group_distance(distance_pair_stat = TRUE)
		#' }
		plot_group_distance = function(
			plot_group_order = NULL,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			distance_pair_stat = FALSE,
			pair_compare_filter_match = NULL,
			pair_compare_filter_select = NULL,
			pair_compare_method = "wilcox.test",
			plot_distance_xtype = NULL
			){
			group_distance <- self$res_group_distance
			group <- self$group
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
			
			p <- ggplot(group_distance, aes_string(x = group, y = "value", color = group)) +
				theme_bw() +
				theme(panel.grid=element_blank()) +
				geom_boxplot(outlier.size =1,width=.6,linetype=1) +
				stat_summary(fun="mean", geom="point", shape=20, size=3, fill="white") +
				xlab("") +
				ylab(ylabname) +
				theme(axis.text=element_text(size=12)) +
				theme(axis.title=element_text(size=17), legend.position = "none") +
				scale_color_manual(values=color_values)
			if(!is.null(plot_distance_xtype)){
				p <- p + theme(axis.text.x = element_text(angle = plot_distance_xtype, colour = "black", vjust = 1, hjust = 1, size = 10))
			}
			if(distance_pair_stat == T){
				# remove some groups
				comparisons_list <- levels(group_distance[, group]) %>% 
					combn(., 2)
				if(!is.null(pair_compare_filter_match) & !is.null(pair_compare_filter_select)){
					stop("The parameter pair_compare_filter_select and pair_compare_filter_match can not be both used together!")
				}
				if(!is.null(pair_compare_filter_match)){
					comparisons_list %<>% {.[, unlist(lapply(as.data.frame(.), function(x) any(grepl(pair_compare_filter_match, x)))), drop = FALSE]}
				}
				if(!is.null(pair_compare_filter_select)){
					if(!is.numeric(pair_compare_filter_select)){
						stop("The parameter pair_compare_filter_select must be numeric !")
					}
					messages_use <- unlist(lapply(as.data.frame(comparisons_list[, pair_compare_filter_select, drop = FALSE]), 
						function(x){paste0(x, collapse = "-")}))
					
					message("Selected groups are ", paste0(messages_use, collapse = " "), " ...")
					comparisons_list %<>% .[, pair_compare_filter_select, drop = FALSE]
				}
				
				# generate the list
				comparisons_list %<>% {lapply(seq_len(ncol(.)), function(x) .[, x])}
				
				p <- p + ggpubr::stat_compare_means(comparisons = comparisons_list, method = pair_compare_method, 
						tip.length=0.01, label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
						symbols = c("****", "***", "**", "*", "ns")))
			}
			
			p
		},
		#' @description
		#' Plotting clustering result. Require ggdendro package.
		#'
		#' @param use_colors colors for presentation.
		#' @param measure default NULL; beta diversity index; If NULL, using the measure when creating object
		#' @param group default NULL; if provided, use this group to assign color.
		#' @param replace_name default NULL; if provided, use this as label.
		#' @return ggplot.
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
						g1 <- g1 + geom_text(data=data2, aes_string(x="x", y="y", label = "replace_name_use", hjust=-0.1), size=4)
					}else{
						g1 <- g1 + geom_text(data=data2, aes_string(x="x", y="y", label = replace_name, hjust=-0.1), size=4)
					}
				}
			} else {
				if(is.null(replace_name)){
					g1 <- g1 + geom_text(data=data2, aes_string(x="x", y="y", label="label", hjust=-0.1, color = group), size=4)
				}else{
					if(length(replace_name) > 1){
						g1 <- g1 + geom_text(data=data2, aes_string(x="x", y="y", label="replace_name_use", hjust=-0.1, color = group), size=4)
					}else{
						g1 <- g1 + geom_text(data=data2, aes_string(x="x", y="y", label=replace_name, hjust=-0.1, color = group), size=4)
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
		},
		#' @description
		#' Print the trans_beta object.
		print = function() {
			cat("trans_beta class:\n")
			if(!is.null(self$ordination)) cat(paste(self$ordination, "is used for ordination \n"))
			if(!is.null(self$res_manova)) cat("PerMANOVA is used for finding significance \n")
			invisible(self)
		}
		),
	private = list(
		within_group_distance = function(distance, sampleinfo, type){
			all_group <- as.character(sampleinfo[,type]) %>% unique
			res <- list()
			for (i in all_group) {
				res[[i]] <- as.vector(as.dist(distance[sampleinfo[,type] == i, sampleinfo[,type] == i]))
			}
			res <- reshape2::melt(res) 
			colnames(res)[2] <- type
			res
		},
		between_group_distance = function(distance, sampleinfo, type) {
			all_group <- as.character(sampleinfo[,type]) %>% unique
			com1 <- combn(all_group,2)
			res <- list()
			for (i in seq_len(ncol(com1))) {
				f_name <- rownames(sampleinfo[sampleinfo[, type] == com1[1,i], ])
				s_name <- rownames(sampleinfo[sampleinfo[, type] == com1[2,i], ])
				vsname <- paste0(com1[1,i], " vs ", com1[2,i])
				res[[vsname]] <- as.vector(distance[f_name, s_name])
			}
			res <- reshape2::melt(res) 
			colnames(res)[2] <- type
			res
		},
		paired_group_manova = function(sample_info_use, use_matrix, group, measure, permutations){
			comnames <- c()
			R2 <- c()
			p.value <- c()
			measure_vec <- c()
			matrix_total <- use_matrix[rownames(sample_info_use), rownames(sample_info_use)]
			groupvec <- as.character(sample_info_use[ , group])
			all_name <- combn(unique(sample_info_use[ , group]), 2)
			for(i in 1:ncol(all_name)) {
				matrix_compare <- matrix_total[groupvec %in% as.character(all_name[,i]), groupvec %in% as.character(all_name[,i])]
				sample_info_compare <- sample_info_use[groupvec %in% as.character(all_name[,i]), ]
				ad <- adonis(reformulate(group, substitute(as.dist(matrix_compare))), data = sample_info_compare, permutations = permutations)
				comnames <- c(comnames, paste0(as.character(all_name[,i]), collapse = " vs "))
				R2 <- c(R2, ad$aov.tab[1,5])
				p.value <- c(p.value, ad$aov.tab[1,6])
				measure_vec <- c(measure_vec, measure)
			}
			significance_label <- cut(p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
			#p.adjusted = p.adjust(p.value, method='fdr')
			permutations_res <- rep(permutations, length(comnames))
			compare_result <- data.frame(comnames, measure_vec, permutations_res, R2, p.value, significance_label)
			colnames(compare_result) <- c("Groups", "measure", "permutations", "R2","p.value", "Significance")
			compare_result
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
