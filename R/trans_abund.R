#' @title
#' Create \code{trans_abund} object for plotting taxonomic abundance.
#'
#' @description
#' This class is a wrapper for the taxonomic abundance transformations and visualization.
#' The converted data style is the long-format for \code{ggplot2} plot.
#' The plotting methods include bar plot, boxplot, heatmap, pie chart and line chart.
#'
#' @export
trans_abund <- R6Class(classname = "trans_abund",
	public = list(
		#' @param dataset default NULL; the object of \code{\link{microtable}} class.
		#' @param taxrank default "Phylum"; taxonomic rank.
		#' @param show default 0; the relative abundance threshold for filtering the taxa with low abundance.
		#' @param ntaxa default 10; how many taxa are selected to show. Taxa are ordered by abundance from high to low. 
		#'   This parameter does not conflict with the parameter \code{show}. Both can be used. \code{ntaxa = NULL} means it is unavailable.
		#' @param groupmean default NULL; calculate mean abundance for each group. Select a column name in \code{microtable$sample_table}.
		#' @param group_morestats default FALSE; only available when \code{groupmean} parameter is provided; 
		#'   Whether output more statistics for each group, including min, max, median and quantile;
		#'   Thereinto, quantile25 and quantile75 denote 25\% and 75\% quantiles, respectively.
		#' @param delete_full_prefix default TRUE; whether delete both the prefix of taxonomy and the character in front of them.
		#' @param delete_part_prefix default FALSE; whether only delete the prefix of taxonomy.
		#' @param prefix default NULL; character string; can be used when \code{delete_full_prefix = T} or \code{delete_part_prefix = T}; 
		#'   default NULL reprensents using the "letter+__", e.g. "k__" for Phylum level;
		#'   Please alter this parameter when the prefix is not standard.
		#' @param use_percentage default TRUE; show the abundance percentage.
		#' @param input_taxaname default NULL; character vector; input taxa names for selecting some taxa.
		#' @param high_level default NULL; a taxonomic rank, such as "Phylum", used to add the taxonomic information of higher level.
		#'   It is necessary for the legend with nested taxonomic levels in the bar plot.
		#' @param high_level_fix_nsub default NULL; an integer, used to fix the number of selected abundant taxa in each taxon from higher taxonomic level.
		#'   If the total number under one taxon of higher level is less than the high_level_fix_nsub, the total number will be used.
		#'   When \code{high_level_fix_nsub} is provided, the taxa number of higher level is calculated as: \code{ceiling(ntaxa/high_level_fix_nsub)}.
		#'   Note that \code{ntaxa} means either the parameter \code{ntaxa} or the taxonomic number obtained by filtering according to the \code{show} parameter.
		#' @return \code{data_abund} stored in the object. The column 'all_mean_abund' reprensents mean relative abundance across all the samples.
		#'   So the values in one taxon are all same across all the samples.
		#'   If the sum of column 'Abundance' in one sample is larger than 1, the 'Abundance', 'SD' and 'SE' has been multiplied by 100.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
		#' }
		initialize = function(
			dataset = NULL, 
			taxrank = "Phylum", 
			show = 0, 
			ntaxa = 10, 
			groupmean = NULL,
			group_morestats = FALSE,
			delete_full_prefix = TRUE,
			delete_part_prefix = FALSE,
			prefix = NULL,
			use_percentage = TRUE, 
			input_taxaname = NULL,
			high_level = NULL,
			high_level_fix_nsub = NULL
			){
			check_microtable(dataset)
			check_taxa_abund(dataset)
			sample_table <- dataset$sample_table
			if("Sample" %in% colnames(sample_table)){
				colnames(sample_table)[colnames(sample_table) == "Sample"] <- "Sample_replace"
			}
			if(! taxrank %in% names(dataset$taxa_abund)){
				stop("The input parameter taxrank: ", taxrank, " is not found! Please check whether it is correct!")
			}
			abund_data <- dataset$taxa_abund[[taxrank]] %>% 
				rownames_to_column(var = "Taxonomy") %>% 
				reshape2::melt(id.vars = "Taxonomy") %>% 
				`colnames<-`(c("Taxonomy", "Sample", "Abundance"))
			check_nd <- grepl("__$", abund_data$Taxonomy)
			if(any(check_nd)){
				abund_data$Taxonomy[check_nd] %<>% paste0(., "unidentified")
			}
			if(delete_full_prefix == T | delete_part_prefix == T){
				if(is.null(prefix)){
					prefix <- ".__"
				}
				if(delete_part_prefix == T){
					delete_full_prefix <- FALSE
				}
				if(delete_full_prefix == T){
					abund_data$Taxonomy %<>% gsub(paste0(".*", prefix, "(.*)"), "\\1", .)
				}else{
					abund_data$Taxonomy %<>% gsub(prefix, "", .)
				}
			}
			abund_data %<>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
				dplyr::summarise(Abundance = sum(Abundance)) %>%
				as.data.frame(stringsAsFactors = FALSE)
			# sort according to the abundance
			abund_data$Taxonomy %<>% as.character
			mean_abund <- tapply(abund_data$Abundance, abund_data$Taxonomy, FUN = mean)
			# add the mean abundance of all samples
			all_mean_abund <- data.frame(Taxonomy = names(mean_abund), all_mean_abund = mean_abund)
			rownames(all_mean_abund) <- NULL
			# add sample table
			abund_data %<>% {suppressWarnings(dplyr::left_join(., rownames_to_column(sample_table), by = c("Sample" = "rowname")))}
			# calculate mean vlaues for each group
			if(!is.null(groupmean)){
				message(paste0(groupmean, " column is used to calculate mean abundance ..."))
				# abund_data[, groupmean] %<>% as.character
				abund_data <- microeco:::summarySE_inter(abund_data, measurevar = "Abundance", groupvars = c("Taxonomy", groupmean), more = group_morestats)
				colnames(abund_data)[colnames(abund_data) == "Mean"] <- "Abundance"
				colnames(abund_data)[colnames(abund_data) == groupmean] <- "Sample"
				if(is.factor(sample_table[, groupmean])){
					abund_data$Sample %<>% factor(., levels = levels(sample_table[, groupmean]))
				}
			}
			abund_data <- dplyr::left_join(abund_data, all_mean_abund, by = c("Taxonomy" = "Taxonomy"))
			if(!is.null(high_level)){
				if(length(high_level) > 1){
					warning("Input high_level has multiple elements! Only select the first one!")
					high_level <- high_level[1]
				}
				message("Add higher taxonomic level into the table ...")
				if(! high_level %in% colnames(dataset$tax_table)){
					stop("Provided high_level must be a colname of input dataset$tax_table!")
				}else{
					# delete the prefix to make two tables consistent
					if(is.null(prefix)){
						prefix <- ".__"
					}
					extract_tax_table <- dataset$tax_table[, c(high_level, taxrank)] %>% unique
					extract_tax_table[, taxrank] %<>% gsub(prefix, "", .)
					abund_data <- dplyr::left_join(abund_data, extract_tax_table, by = c("Taxonomy" = taxrank))
				}
			}
			# get ordered taxa
			use_taxanames <- as.character(rev(names(sort(mean_abund))))
			if(!is.null(ntaxa)){
				ntaxa_theshold <- ntaxa_use <- ntaxa
			}else{
				ntaxa_theshold <- ntaxa_use <- sum(mean_abund > show)
			}
			if(ntaxa_use > sum(mean_abund > show)){
				ntaxa_use <- sum(mean_abund > show)
			}
			# filter useless taxa
			use_taxanames %<>% .[!grepl("unidentified|unculture|Incertae.sedis", .)]
			# identify used taxa
			if(is.null(input_taxaname)){
				if(is.null(high_level_fix_nsub)){
					if(length(use_taxanames) > ntaxa_use){
						use_taxanames %<>% .[1:ntaxa_use]
					}
				}else{
					high_level_n <- ceiling(ntaxa_use/high_level_fix_nsub)
					high_level_ordered_taxa <- abund_data[match(names(mean_abund), abund_data$Taxonomy), high_level] %>% unique %>% .[1:high_level_n]
					use_taxanames <- lapply(high_level_ordered_taxa, function(x){
						tmp <- abund_data[abund_data[, high_level] == x, ]
						tmp <- names(mean_abund) %>% .[. %in% tmp$Taxonomy]
						tmp[1:ifelse(length(tmp) < high_level_fix_nsub, length(tmp), high_level_fix_nsub)]
					}) %>% unlist
				}
			}else{
				# make sure input_taxaname are in use_taxanames
				if(!any(input_taxaname %in% use_taxanames)){
					stop("The input_taxaname does not match to taxa names! Please check the input!")
				}else{
					use_taxanames <- input_taxaname[input_taxaname %in% use_taxanames]
				}
			}
			if(!is.null(high_level)){
				# sort the taxa in high levels according to the sum of abundances
				tmp <- abund_data[abund_data$Taxonomy %in% use_taxanames, c(high_level, "Abundance")]
				tmp <- tapply(tmp$Abundance, tmp[, high_level], FUN = sum)
				data_taxanames_highlevel <- as.character(names(sort(tmp, decreasing = TRUE)))
				self$data_taxanames_highlevel <- data_taxanames_highlevel
			}
			# ylab title for different cases; more clear
			if(ntaxa_theshold < sum(mean_abund > show) | show == 0){
				if(use_percentage == T){
					abund_data$Abundance %<>% {. * 100}
					if("SE" %in% colnames(abund_data)) abund_data$SE %<>% {. * 100}
					if("SD" %in% colnames(abund_data)) abund_data$SD %<>% {. * 100}
					ylabname <- "Relative abundance (%)"
				}else{
					ylabname <- "Relative abundance"
				}
			}else{
				ylabname <- paste0("Relative abundance (", taxrank, " > ", show*100, "%)")
			}
			self$use_percentage <- use_percentage
			self$ylabname <- ylabname
			self$taxrank <- taxrank
			self$data_abund <- abund_data
			self$data_taxanames <- use_taxanames
			self$high_level <- high_level
			message('The transformed abundance data is stored in object$data_abund ...')
		},
		#' @description
		#' Bar plot.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the bars.
		#' @param bar_type default "full"; "full" or "notfull"; if \code{"full"}, total abundance are summed to 1 or 100 percentage.
		#' @param others_color default "grey90"; the color for "others" taxa.
		#' @param facet default NULL; a character vector for the facet; group column name of \code{sample_table}, such as, \code{"Group"};
		#'    If multiple facets are needed, please provide ordered names, such as \code{c("Group", "Type")}.
		#'    The latter should have a finer scale than the former one;
		#'    Please adjust the facet orders in the plot by assigning factors in \code{sample_table} before creating \code{trans_abund} object or 
		#'    assigning factors in the \code{data_abund} table of \code{trans_abund} object.
		#'    When multiple facets are used, please first install package \code{ggh4x} using the command \code{install.packages("ggh4x")}.
		#' @param order_x default NULL; vector; used to order the sample names in x axis; must be the samples vector, such as \code{c("S1", "S3", "S2")}.
		#' @param x_axis_name NULL; a character string; a column name of sample_table in dataset; used to show the sample names in x axis.
		#' @param barwidth default NULL; bar width, see \code{width} in \code{\link{geom_bar}}.
		#' @param use_alluvium default FALSE; whether add alluvium plot. If \code{TRUE}, please first install \code{ggalluvial} package.
		#' @param clustering default FALSE; whether order samples by the clustering.
		#' @param clustering_plot default FALSE; whether add clustering plot.
		#'     If \code{clustering_plot = TRUE}, \code{clustering} will be also TRUE in any case for the clustering.
		#' @param cluster_plot_width default 0.2, the dendrogram plot width; available when \code{clustering_plot = TRUE}.
		#' @param facet_color default "grey95"; facet background color.
		#' @param strip_text default 11; facet text size.
		#' @param legend_text_italic default FALSE; whether use italic in legend.
		#' @param xtext_angle default 0; number ranging from 0 to 90; used to adjust x axis text angle to reduce text overlap; 
		#' @param xtext_size default 10; x axis text size.
		#' @param xtext_keep default TRUE; whether retain x text.
		#' @param xtitle_keep default TRUE; whether retain x title.
		#' @param ytitle_size default 17; y axis title size.
		#' @param coord_flip default FALSE; whether flip cartesian coordinates so that horizontal becomes vertical, and vertical becomes horizontal.
		#' @param ggnested default FALSE; whether use nested legend. Need \code{ggnested} package to be installed (https://github.com/gmteunisse/ggnested).
		#'   To make it available, please assign \code{high_level} parameter when creating the object.
		#' @param high_level_add_other default FALSE; whether add 'Others' (all the unknown taxa) in each taxon of higher taxonomic level.
		#'   Only available when \code{ggnested = TRUE}.
		#' @return ggplot2 object. 
		#' @examples
		#' \donttest{
		#' t1$plot_bar(facet = "Group", xtext_keep = FALSE)
		#' }
		plot_bar = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			bar_type = "full",
			others_color = "grey90",
			facet = NULL,
			order_x = NULL,
			x_axis_name = NULL,
			barwidth = NULL,
			use_alluvium = FALSE,
			clustering = FALSE,
			clustering_plot = FALSE,
			cluster_plot_width = 0.2,
			facet_color = "grey95",
			strip_text = 11,
			legend_text_italic = FALSE,
			xtext_angle = 0,
			xtext_size = 10,
			xtext_keep = TRUE,
			xtitle_keep = TRUE,
			ytitle_size = 17,
			coord_flip = FALSE,
			ggnested = FALSE,
			high_level_add_other = FALSE
			){
			plot_data <- self$data_abund
			# try to filter useless columns
			plot_data %<>% .[, ! colnames(.) %in% c("N", "SD", "SE", "Median", "Min", "Max", "quantile25", "quantile75", "all_mean_abund")]
			use_taxanames <- self$data_taxanames
			if(ggnested){
				if(is.null(self$high_level)){
					stop("The high_level is necessary when ggnested = TRUE! Please assign high_level parameter when creating the object!")
				}
				if(high_level_add_other){
					plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
					use_taxanames %<>% c(., "Others")
					new_data <- plot_data %>% dplyr::group_by(!!! syms(c(self$high_level, "Taxonomy", "Sample"))) %>% 
						dplyr::summarise(Abundance = sum(Abundance)) %>%
						as.data.frame(stringsAsFactors = FALSE)
					plot_data_merge <- plot_data[, ! colnames(plot_data) %in% c(self$high_level, "Taxonomy", "Abundance"), drop = FALSE] %>% unique
					plot_data <- dplyr::left_join(new_data, plot_data_merge, by = c("Sample" = "Sample"))
				}
				bar_type <- "notfull"
			}
			if(bar_type == "full"){
				# make sure that taxonomy info are all in selected use_taxanames in case of special data
				if(!all(plot_data$Taxonomy %in% use_taxanames)){
					plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
					new_data <- plot_data %>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
						dplyr::summarise(Abundance = sum(Abundance)) %>%
						as.data.frame(stringsAsFactors = FALSE)
					plot_data_merge <- plot_data[, ! colnames(plot_data) %in% c("Taxonomy", "Abundance"), drop = FALSE] %>% unique
					plot_data <- dplyr::left_join(new_data, plot_data_merge, by = c("Sample" = "Sample"))
					plot_data$Taxonomy %<>% factor(., levels = rev(c(use_taxanames, "Others")))
				}else{
					plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
				}
			}else{
				if(ggnested){
					plot_data %<>% .[.[, self$high_level] %in% self$data_taxanames_highlevel, ]
					plot_data[, self$high_level] %<>% factor(., levels = self$data_taxanames_highlevel)
				}
				plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
				# two legend ordering types depending on ggnested
				if(ggnested){
					plot_data$Taxonomy %<>% factor(., levels = use_taxanames)
				}else{
					plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
				}
			}
			# order x axis samples
			plot_data <- private$adjust_axis_facet(
				plot_data = plot_data, 
				x_axis_name = x_axis_name, 
				order_x = order_x
				)
			# arrange plot_data--Abundance according to the Taxonomy-group column factor-levels
			plot_data <- plot_data[unlist(lapply(levels(plot_data$Taxonomy), function(x) which(plot_data$Taxonomy == x))),]
			if(!ggnested){
				if(any(grepl("Others", as.character(plot_data$Taxonomy)))){
					bar_colors_use <- expand_colors(color_values, length(unique(plot_data$Taxonomy)) - 1)
					bar_colors_use <- c(bar_colors_use, others_color)
				}else{
					bar_colors_use <- expand_colors(color_values, length(unique(plot_data$Taxonomy)))
				}
			}else{
				# high_level determine the colors
				bar_colors_use <- expand_colors(color_values, length(unique(plot_data[, self$high_level])))
			}
			if(clustering | clustering_plot){
				data_clustering <- reshape2::dcast(plot_data, Sample ~ Taxonomy, value.var = "Abundance", fun.aggregate = sum) %>% 
					`row.names<-`(.[,1]) %>% .[, -1]
				tmp_hclust <- hclust(dist(data_clustering)) 
				order_x_clustering <- tmp_hclust %>% {.$labels[.$order]} %>% as.character
				plot_data$Sample %<>% factor(., levels = order_x_clustering)
			}
			if(use_alluvium){
				p <- ggplot(plot_data, aes(
						x = Sample, y = Abundance, 
						fill = Taxonomy, color = Taxonomy, 
						weight = Abundance, 
						alluvium = Taxonomy, stratum = Taxonomy
					)) +
					ggalluvial::geom_flow(alpha = .4, width = 3/15) +
					ggalluvial::geom_stratum(width = .2) +
					scale_color_manual(values = rev(bar_colors_use))
			}else{
				if(ggnested){
					p <- ggnested::ggnested(plot_data, aes_meco(x = "Sample", y = "Abundance", main_group = self$high_level, sub_group = "Taxonomy"), main_palette = bar_colors_use)
				}else{
					p <- ggplot(plot_data, aes_meco(x = "Sample", y = "Abundance", fill = "Taxonomy"))
				}
				if(bar_type == "full"){
					if(self$use_percentage == T){
						p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth)
					}else{
						p <- p + geom_bar(stat = "identity", position = "fill", show.legend = T, width = barwidth)
					}
				}else{
					p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth)
				}
			}
			if(!ggnested){
				p <- p + scale_fill_manual(values = rev(bar_colors_use))
			}
			p <- p + xlab("") + ylab(self$ylabname)
			if(!is.null(facet)){
				if(coord_flip){
					facet_formula <- reformulate(".", paste0(facet, collapse = " + "))
				}else{
					facet_formula <- reformulate(facet, ".")
				}
				if(length(facet) == 1){
					p <- p + facet_grid(facet_formula, scales = "free", space = "free")
				}else{
					p <- p + ggh4x::facet_nested(facet_formula, nest_line = element_line(linetype = 2), scales = "free", space = "free")
				}
				p <- p + theme(strip.background = element_rect(fill = facet_color, color = facet_color), strip.text = element_text(size=strip_text))
				p <- p + scale_y_continuous(expand = c(0, 0.01))
			}else{
				if(bar_type == "full" & self$use_percentage == FALSE){
					p <- p + scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
				}else{
					p <- p + scale_y_continuous(expand = c(0, 0))
				}
			}
			p <- p + theme(panel.grid = element_blank(), panel.border = element_blank()) + 
				theme(axis.line.y = element_line(color = "grey60", linetype = "solid", lineend = "square"))
			if(legend_text_italic == T) {
				p <- p + theme(legend.text = element_text(face = 'italic'))
			}
			if(clustering_plot){
				if(! coord_flip){
					message("Rotate the axis automatically to add the clustering plot ...")
					coord_flip <- TRUE
				}
			}
			p <- p + private$ggplot_xtext_type(xtext_angle = xtext_angle, xtext_size = xtext_size, xtext_keep = xtext_keep, coord_flip = coord_flip)
			p <- p + theme(axis.title.y = element_text(size = ytitle_size))
			if(xtitle_keep == F){
				p <- p + theme(axis.title.x = element_blank())
			}
			p <- p + guides(fill = guide_legend(title = self$taxrank))
			if(use_alluvium | ggnested){
				p <- p + guides(color = guide_legend(title = self$taxrank))
			}
			if(coord_flip){
				p <- p + coord_flip()
			}
			if(clustering_plot){
				left_plot <- ggtree::ggtree(tmp_hclust, hang = 0)
				p %<>% aplot::insert_left(left_plot, width = cluster_plot_width)
			}
			p
		},
		#' @description
		#' Plot the heatmap.
		#'
		#' @param color_values default rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")); 
		#' 	  colors palette for the plotting.
		#' @param facet default NULL; a character vector for the facet; a group column name of \code{sample_table}, such as, \code{"Group"};
		#'    If multiple facets are needed, please provide ordered names, such as \code{c("Group", "Type")}.
		#'    The latter should have a finer scale than the former one;
		#'    Please adjust the facet orders in the plot by assigning factors in \code{sample_table} before creating \code{trans_abund} object or 
		#'    assigning factors in the \code{data_abund} table of \code{trans_abund} object.
		#'    When multiple facets are used, please first install package \code{ggh4x} using the command \code{install.packages("ggh4x")}.
		#' @param x_axis_name NULL; a character string; a column name of sample_table used to show the sample names in x axis.
		#' @param order_x default NULL; vector; used to order the sample names in x axis; must be the samples vector, such as, c("S1", "S3", "S2").
		#' @param withmargin default TRUE; whether retain the tile margin.
		#' @param plot_numbers default FALSE; whether plot the number in heatmap.
		#' @param plot_text_size default 4; If plot_numbers TRUE, text size in plot.
		#' @param plot_breaks default NULL; The legend breaks.
		#' @param margincolor default "white"; If withmargin TRUE, use this as the margin color.
		#' @param plot_colorscale default "log10"; color scale.
		#' @param min_abundance default .01; the minimum abundance percentage in plot.
		#' @param max_abundance default NULL; the maximum abundance percentage in plot, NULL reprensent the max percentage.
		#' @param strip_text default 11; facet text size.
		#' @param xtext_size default 10; x axis text size.
		#' @param ytext_size default 11; y axis text size.
		#' @param xtext_keep default TRUE; whether retain x text.
		#' @param xtitle_keep default TRUE; whether retain x title.
		#' @param grid_clean default TRUE; whether remove grid lines.
		#' @param xtext_angle default 0; number ranging from 0 to 90; used to adjust x axis text angle to reduce text overlap; 
		#' @param legend_title default "\% Relative\\nAbundance"; legend title text.
		#' @param pheatmap default FALSE; whether use pheatmap package to plot the heatmap.
		#' @param ... paremeters pass to pheatmap when pheatmap = TRUE.
		#' @return ggplot2 object or grid object based on pheatmap.
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
		#' t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
		#' }
		plot_heatmap = function(
			color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")), 
			facet = NULL,
			x_axis_name = NULL,
			order_x = NULL,
			withmargin = TRUE,
			plot_numbers = FALSE,
			plot_text_size = 4,
			plot_breaks = NULL,
			margincolor = "white",
			plot_colorscale = "log10",
			min_abundance = 0.01,
			max_abundance = NULL,
			strip_text = 11,
			xtext_size = 10,
			ytext_size = 11,
			xtext_keep = TRUE,
			xtitle_keep = TRUE,
			grid_clean = TRUE,
			xtext_angle = 0,
			legend_title = "% Relative\nAbundance",
			pheatmap = FALSE,
			...
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames
			plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}

			if(pheatmap == FALSE){
				# order x axis samples
				plot_data <- private$adjust_axis_facet(plot_data = plot_data, x_axis_name = x_axis_name, order_x = order_x)
				if (is.null(min_abundance)){
					min_abundance <- ifelse(min(plot_data$Abundance) > 0.001, min(plot_data$Abundance), 0.001)
				}
				if (is.null(max_abundance)){
					max_abundance <- max(plot_data$Abundance)
				}
				plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))

				p <- ggplot(plot_data, aes(x = .data[["Sample"]], y = .data[["Taxonomy"]], label = .data[[formatC("Abundance", format = "f", digits = 1)]]))
				
				if(withmargin == T){
					p <- p + geom_tile(aes(fill = Abundance), colour = margincolor, size = 0.5)
				}else{
					p <- p + geom_tile(aes(fill = Abundance))
				}
				p <- p + theme(axis.text.y = element_text(size = 12)) + theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

				if (plot_numbers == T){
					abund <- plot_data
					abund$Abundance <- round(abund$Abundance, 1)
					p <- p + geom_text(data = abund, size = plot_text_size, colour = "grey10")  
				}
				if (is.null(plot_breaks)){
					p <- p + scale_fill_gradientn(colours = color_values, trans = plot_colorscale, na.value = "#00008B", limits = c(min_abundance, max_abundance))
				}else{
					p <- p + scale_fill_gradientn(colours = color_values, trans = plot_colorscale, breaks=plot_breaks, na.value = "#00008B",
						limits = c(min_abundance, max_abundance))
				}
				if(!is.null(facet)){
					if(length(facet) == 1){
						p <- p + facet_grid(reformulate(facet, "."), scales = "free", space = "free")
					}else{
						p <- p + ggh4x::facet_nested(reformulate(facet), nest_line = element_line(linetype = 2), scales = "free", space = "free")
					}
					p <- p + theme(strip.background = element_rect(color = "white", fill = "grey92"), strip.text = element_text(size=strip_text))
				}
				p <- p + labs(x = "", y = "", fill = legend_title)
				if (!is.null(ytext_size)){
					p <- p + theme(axis.text.y = element_text(size = ytext_size))
				}
				p <- p + private$ggplot_xtext_type(xtext_angle = xtext_angle, xtext_size = xtext_size, xtext_keep = xtext_keep)
				if(grid_clean){
					p <- p + theme(panel.border = element_blank(), panel.grid = element_blank())
				}
				p
			} else {
				# first to wide format
				wide_table <- reshape2::dcast(plot_data, Taxonomy ~ Sample, value.var = "Abundance") %>% 
					`row.names<-`(.[,1]) %>% 
					.[, -1, drop = FALSE]
				# check sd for each feature, if 0, delete
				if(any(apply(wide_table, MARGIN = 1, FUN = function(x) sd(x) == 0))){
					select_rows <- apply(wide_table, MARGIN = 1, FUN = function(x) sd(x) != 0)
					wide_table %<>% {.[select_rows, ]}
				}
				p <- pheatmap::pheatmap(
					wide_table,
					...
					)
				p$gtable
			}
		},
		#' @description
		#' Box plot.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the box.
		#' @param group default NULL; a column name of sample table to show abundance across groups.
		#' @param show_point default FALSE; whether show points in plot.
		#' @param point_color default "black"; If show_point TRUE; use the color
		#' @param point_size default 3; If show_point TRUE; use the size
		#' @param point_alpha default .3; If show_point TRUE; use the transparency.
		#' @param plot_flip default FALSE; Whether rotate plot.
		#' @param boxfill default TRUE; Whether fill the box with colors.
		#' @param middlecolor default "grey95"; The middle line color.
		#' @param middlesize default 1; The middle line size.
		#' @param xtext_angle default 0; number ranging from 0 to 90; used to adjust x axis text angle to reduce text overlap; 
		#' @param xtext_size default 10; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @param ... parameters pass to \code{\link{geom_boxplot}} function.
		#' @return ggplot2 object. 
		#' @examples
		#' \donttest{
		#' t1$plot_box(group = "Group")
		#' }
		plot_box = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			group = NULL,
			show_point = FALSE,
			point_color = "black",
			point_size = 3,
			point_alpha = .3,
			plot_flip = FALSE,
			boxfill = TRUE,
			middlecolor = "grey95",
			middlesize = 1,
			xtext_angle = 0,
			xtext_size = 10,
			ytitle_size = 17,
			...
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames

			plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
			plot_data$Taxonomy %<>% factor(., levels = use_taxanames)

			p <- ggplot(plot_data, aes(x = .data[["Taxonomy"]], y = .data[["Abundance"]])) 
			p <- p + ylab(self$ylabname) + guides(col = guide_legend(reverse = TRUE)) + xlab("")
			if (plot_flip == T){ 
				p <- p + coord_flip()
			}
			if(is.null(group)) {
				p <- p + geom_boxplot(color = color_values[1], ...)
			} else {
				color_values <- expand_colors(color_values, length(unique(plot_data[, group])))
				if(boxfill == T){
					p <- p + geom_boxplot(aes(color = .data[[group]], fill = .data[[group]]), ...)
					p <- p + scale_fill_manual(values = color_values)
					p <- p + scale_color_manual(values = color_values) + guides(color = "none")
					## Change the default middle line
					dat <- ggplot_build(p)$data[[1]]
					p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour = middlecolor, size=middlesize)
				} else {	 
					p <- p + geom_boxplot(aes(color = .data[[group]]), ...) + scale_color_manual(values = color_values)
				}
			}
			if(show_point == T){
				p <- p + geom_point(size = point_size, color = point_color, alpha = point_alpha, position = "jitter")
			}
			p <- p + private$ggplot_xtext_type(xtext_angle = xtext_angle, xtext_size = xtext_size)
			p <- p + theme(axis.title.y = element_text(size = ytitle_size)) + scale_y_continuous(expand = c(0, 0.01))

			if(!is.null(group)) {
				p <- p + guides(fill=guide_legend(title=group))
			}
			p
		},
		#' @description
		#' Plot the line chart.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the points and lines.
		#' @param plot_SE default TRUE; TRUE: the errorbar is \eqn{mean±se}; FALSE: the errorbar is \eqn{mean±sd}.
		#' @param position default position_dodge(0.1); Position adjustment, either as a string (such as "identity"), or the result of a call to a position adjustment function.
		#' @param errorbar_size default 1; errorbar size.
		#' @param errorbar_width default 0.1; errorbar width.
		#' @param point_size default 3; point size for taxa.
		#' @param point_alpha default 0.8; point transparency.
		#' @param line_size default 0.8; line size.
		#' @param line_alpha default 0.8; line transparency.
		#' @param line_type default 1; an integer; line type.
		#' @param xtext_angle default 0; number ranging from 0 to 90; used to adjust x axis text angle to reduce text overlap; 
		#' @param xtext_size default 10; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @return ggplot2 object. 
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 5)
		#' t1$plot_line(point_size = 3)
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 5, groupmean = "Group")
		#' t1$plot_line(point_size = 5, errorbar_size = 1, xtext_angle = 30)
		#' }
		plot_line = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			plot_SE = TRUE,
			position = position_dodge(0.1),
			errorbar_size = 1,
			errorbar_width = 0.1,
			point_size = 3,
			point_alpha = 0.8,
			line_size = 0.8, 
			line_alpha = 0.8, 
			line_type = 1,
			xtext_angle = 0,
			xtext_size = 10,
			ytitle_size = 17
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames
			plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
			plot_data$Taxonomy %<>% factor(., levels = use_taxanames)
			color_values <- expand_colors(color_values, length(use_taxanames))
			
			p <- ggplot(plot_data, aes(x = .data[["Sample"]], y = .data[["Abundance"]], color = .data[["Taxonomy"]], group = .data[["Taxonomy"]]))
			if(("SE" %in% colnames(plot_data)) & plot_SE){
				p <- p + geom_errorbar(aes(ymin = Abundance - SE, ymax = Abundance + SE), width = errorbar_width, position = position, size = errorbar_size)
			}else{
				if(("SD" %in% colnames(plot_data)) & plot_SE){
					p <- p + geom_errorbar(aes(ymin = Abundance - SD, ymax = Abundance + SD), width = errorbar_width, position = position, size = errorbar_size)
				}
			}
			p <- p + geom_point(size = point_size, alpha = point_alpha, position = position)
			p <- p + geom_line(size = line_size, alpha = line_alpha, linetype = line_type, position = position)
			p <- p + ylab(self$ylabname) + guides(col = guide_legend(title=self$taxrank, reverse = TRUE)) + xlab("")
			p <- p + private$ggplot_xtext_type(xtext_angle = xtext_angle, xtext_size = xtext_size)
			p <- p + theme(axis.title.y = element_text(size = ytitle_size)) + scale_y_continuous(expand = c(0, 0.01))
			p <- p + scale_color_manual(values = color_values)
			p
		},
		#' @description
		#' Pie chart.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for each section.
		#' @param facet_nrow default 1; how many rows in the plot.
		#' @param strip_text default 11; sample title size.
		#' @param add_label default FALSE; Whether add the percentage label in each section of pie chart.
		#' @param legend_text_italic default FALSE; whether use italic in legend.
		#' @return ggplot2 object. 
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
		#' t1$plot_pie(facet_nrow = 1)
		#' }
		plot_pie = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			facet_nrow = 1, 
			strip_text = 11, 
			add_label = FALSE,
			legend_text_italic = FALSE
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames
			# sum others to one
			if(any(!plot_data$Taxonomy %in% use_taxanames)){
				plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
				plot_data %<>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
					dplyr::summarise(Abundance = sum(Abundance)) %>%
					as.data.frame(stringsAsFactors = FALSE)
				plot_data$Taxonomy %<>% factor(., levels = c(use_taxanames, "Others"))
				color_values <- expand_colors(color_values, length(use_taxanames) + 1)
			}else{
				color_values <- expand_colors(color_values, length(use_taxanames))
			}
			plot_data$label <- paste0(round(plot_data$Abundance, 1), "%")
			p <- ggplot(plot_data, aes(x = '', y = Abundance, fill = Taxonomy, label = label)) + 
				geom_bar(width = 1, stat = "identity") +
				coord_polar("y", start = 0)
			if(add_label){
				p <- p + ggrepel::geom_label_repel(position = position_stack(vjust = 0.5), show.legend = FALSE)
			}
			p <- p + private$blank_theme +
				scale_fill_manual(values = color_values) +
				theme(axis.text.x = element_blank()) +
				facet_wrap(~Sample, nrow = facet_nrow) +
				theme(strip.text = element_text(size = strip_text)) +
				guides(fill = guide_legend(title = self$taxrank))
			if(legend_text_italic == T) {
				p <- p + theme(legend.text = element_text(face = 'italic'))
			}
			p
		},
		#' @description
		#' Donut chart based on the \code{ggpubr::ggdonutchart} function.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the donut.
		#' @param label default TRUE; whether show the percentage label.
		#' @param facet_nrow default 1; how many rows in the plot.
		#' @param legend_text_italic default FALSE; whether use italic in legend.
		#' @param ... parameters passed to \code{ggpubr::ggdonutchart}.
		#' @return combined ggplot2 objects list, generated by \code{ggpubr::ggarrange} function. 
		#' @examples
		#' \dontrun{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
		#' t1$plot_donut(label = TRUE)
		#' }
		plot_donut = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			label = TRUE,
			facet_nrow = 1, 
			legend_text_italic = FALSE,
			...
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames
			# sum others to one
			if(any(!plot_data$Taxonomy %in% use_taxanames)){
				plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
				plot_data %<>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
					dplyr::summarise(Abundance = sum(Abundance)) %>%
					as.data.frame(stringsAsFactors = FALSE)
				plot_data$Taxonomy %<>% factor(., levels = c(use_taxanames, "Others"))
				color_values <- expand_colors(color_values, length(use_taxanames) + 1)
			}else{
				color_values <- expand_colors(color_values, length(use_taxanames))
			}
			plot_data$label <- paste0(round(plot_data$Abundance, 1), "%")

			# use ggarrange, because facet can not seperate the labels of all the samples
			plot_list <- list()
			for(i in unique(plot_data$Sample)){
				tmp <- plot_data[plot_data$Sample == i, ]
				p <- ggpubr::ggdonutchart(tmp, "Abundance", fill = "Taxonomy", label = "label", color = "white", palette = color_values, ...) + 
					guides(fill = guide_legend(title=self$taxrank))
					theme(axis.text.y = element_blank())
				if(label == F){
					p <- p + theme(axis.text.x = element_blank())
				}
				if(legend_text_italic == T) {
					p <- p + theme(legend.text = element_text(face = 'italic'))
				}
				plot_list[[i]] <- p
			}
			facet_ncol <- ceiling(length(plot_list)/facet_nrow)
			ggpubr::ggarrange(plotlist = plot_list, nrow = facet_nrow, ncol = facet_ncol, labels = names(plot_list), common.legend = TRUE, legend = "bottom")
		},
		#' @description
		#' Radar chart based on the \code{ggradar} package (https://github.com/ricardo-bion/ggradar).
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for samples.
		#' @param ... parameters passed to \code{ggradar::ggradar} function except group.colours parameter.
		#' @return ggplot2 object. 
		#' @examples
		#' \dontrun{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
		#' t1$plot_radar()
		#' }
		plot_radar = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			...
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames
			color_values <- expand_colors(color_values, length(use_taxanames))
			plot_data <- plot_data[plot_data$Taxonomy %in% use_taxanames, ]
			if(self$use_percentage){
				plot_data$Abundance %<>% {./100}
			}
			tmp_data <- reshape2::dcast(plot_data, Taxonomy ~ Sample, value.var = "Abundance")
			colnames(tmp_data)[1] <- "group"
			tmp_data$group %<>% factor(., levels = use_taxanames)
			# https://github.com/ricardo-bion/ggradar
			ggradar::ggradar(tmp_data, group.colours = color_values, ...)
		},
		#' @description
		#' Ternary diagrams based on the \code{ggtern} package.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the samples.
		#' @param color_legend_guide_size default 4; The size of legend guide for color.
		#' @return ggplot2 object. 
		#' @examples
		#' \dontrun{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
		#' t1$plot_tern()
		#' }
		plot_tern = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			color_legend_guide_size = 4
			){
			data_abund <- self$data_abund
			use_taxanames <- self$data_taxanames
			if(length(unique(data_abund$Sample)) != 3){
				stop("Ternary diagrams can only be used for three samples!")
			}
			data_abund %<>% {.[.$Taxonomy %in% use_taxanames, ]}
			data_abund$Taxonomy %<>% factor(., levels = use_taxanames)
			plot_data <- reshape2::dcast(data_abund, Taxonomy + all_mean_abund ~ Sample, value.var = "Abundance")
			colnames(plot_data)[2] <- "Abundance"
			if(is.factor(data_abund$Sample)){
				sample_names <- levels(data_abund$Sample)
			}else{
				sample_names <- unique(data_abund$Sample)
			}
			color_values <- expand_colors(color_values, length(use_taxanames))

			p <- ggtern::ggtern(data = plot_data, aes_meco(x = sample_names[1], y = sample_names[2], z = sample_names[3])) +
				theme_bw() +
				geom_point(aes(color = Taxonomy, size = Abundance)) +
				scale_color_manual(values = color_values) +
				scale_fill_manual(values = color_values) +
				guides(color = guide_legend(override.aes = list(size = color_legend_guide_size)))
			p
		},
		#' @description
		#' Print the trans_abund object.
		print = function(){
			cat("trans_abund object:\n")
			cat(paste("data_abund have", ncol(self$data_abund), "columns: ", paste0(colnames(self$data_abund), collapse = ", "), "\n"))
			cat(paste("data_abund have total", length(unique(as.character(self$data_abund$Taxonomy))), "taxa\n"))
			cat(paste("taxrank: ", self$taxrank, "\n"))
			if(!is.null(self$data_taxanames)){
				if(length(self$data_taxanames) > 50){
					cat(paste("data_taxanames: ", length(self$data_taxanames), "taxa\n"))
				}else{
					cat(paste("data_taxanames: ", paste0(self$data_taxanames, collapse = ", "), "\n"))
				}
			}
			invisible(self)
		}
		),
	private = list(
		adjust_axis_facet = function(plot_data, x_axis_name, order_x){
			# order x axis samples and facet
			if(!is.null(x_axis_name)){
				colnames(plot_data)[colnames(plot_data) == "Sample"] <- "Sample_rownames_before"
				if(! x_axis_name %in% colnames(plot_data)){
					stop(paste("No", x_axis_name, "found in the column names of sample_table!"))
				}else{
					colnames(plot_data)[colnames(plot_data) == x_axis_name] <- "Sample"
				}
			}
			if(!is.null(order_x)){
				if(length(order_x) == 1){
					stop("This may be wrong. Only one sample used to order the samples!")
				}else{
					plot_data$Sample %<>% factor(., levels = order_x)
				}
			}
			plot_data
		},
		ggplot_xtext_type = function(xtext_angle, xtext_size, xtext_keep = TRUE, coord_flip = FALSE){
			if(coord_flip){
				if(xtext_keep){
					theme(axis.text.y = element_text(colour = "black", size = xtext_size))
				}else{
					theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
				}
			}else{
				if(xtext_keep){
					ggplot_xtext_anglesize(xtext_angle, xtext_size)
				}else{
					theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
				}
			}
		},
		blank_theme = 
			theme_minimal() +
			theme(
			axis.title = element_blank(),
			panel.border = element_blank(),
			panel.grid=element_blank(),
			axis.ticks = element_blank(),
			legend.position="right",
			plot.title=element_text(size=14, face="bold")
		)
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
