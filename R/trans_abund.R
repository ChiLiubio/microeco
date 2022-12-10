#' @title
#' Create \code{trans_abund} object for plotting taxonomic abundance.
#'
#' @description
#' This class is a wrapper for the taxonomic abundance transformations and plotting. 
#' The transformed data style is the long-format for \code{ggplot2} plotting.
#' The plotting methods include bar plot, boxplot, heatmap and pie chart.
#'
#' @export
trans_abund <- R6Class(classname = "trans_abund",
	public = list(
		#' @param dataset default NULL; the object of \code{\link{microtable}} class.
		#' @param taxrank default "Phylum"; taxonomic rank.
		#' @param show default 0; the relative abundance threshold used for filtering.
		#' @param ntaxa default 10; how many taxa will be used, ordered by abundance from high to low; 
		#'   this parameter does not conflict with the parameter show; both can be used.
		#' @param groupmean default NULL; calculating mean abundance for each group; select a group column name in \code{microtable$sample_table}.
		#' @param delete_full_prefix default TRUE; whether delete both the prefix of taxonomy and the character in front of them.
		#' @param delete_part_prefix default FALSE; whether only delete the prefix of taxonomy.
		#' @param prefix default NULL; character string; can be used when \code{delete_full_prefix = T} or \code{delete_part_prefix = T}; 
		#'   default NULL reprensents using the "letter+__", e.g. "k__" for Phylum level;
		#'  Please alter this parameter when the prefix is not standard.
		#' @param use_percentage default TRUE; show the abundance percentage.
		#' @param input_taxaname default NULL; character vector; input taxa names for selecting some taxa.
		#' @return \code{data_abund} stored in the object.
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
			delete_full_prefix = TRUE,
			delete_part_prefix = FALSE,
			prefix = NULL,
			use_percentage = TRUE, 
			input_taxaname = NULL
			){
			if(is.null(dataset)){
				stop("No dataset provided!")
			}
			if(is.null(dataset$taxa_abund)){
				message('The taxa_abund in dataset is NULL. Calculate it now ...')
				dataset$cal_abund()
			}
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
			
			# add sample table
			abund_data %<>% {suppressWarnings(dplyr::left_join(., rownames_to_column(sample_table), by = c("Sample" = "rowname")))}
			# calculate mean vlaues for each group
			if(!is.null(groupmean)){
				message(paste0(groupmean, " column is used to calculate mean abundance ..."))
				# abund_data[, groupmean] %<>% as.character
				abund_data <- microeco:::summarySE_inter(abund_data, measurevar = "Abundance", groupvars = c("Taxonomy", groupmean))
				colnames(abund_data)[colnames(abund_data) == "Mean"] <- "Abundance"
				colnames(abund_data)[colnames(abund_data) == groupmean] <- "Sample"
				if(is.factor(sample_table[, groupmean])){
					abund_data$Sample %<>% factor(., levels = levels(sample_table[, groupmean]))
				}
			}
			# sort according to the abundance
			abund_data$Taxonomy %<>% as.character
			mean_abund <- tapply(abund_data$Abundance, abund_data$Taxonomy, FUN = mean)
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
				if(length(use_taxanames) > ntaxa_use){
					use_taxanames %<>% .[1:ntaxa_use]
				}
			}else{
				# make sure input_taxaname are in use_taxanames
				if(!any(input_taxaname %in% use_taxanames)){
					stop("input_taxaname do not match to taxa names! Please check the input!")
				}else{
					use_taxanames <- input_taxaname[input_taxaname %in% use_taxanames]
				}
			}
			# generate ylab title and store it
			if(ntaxa_theshold < sum(mean_abund > show) | show == 0){
				if(use_percentage == T){
					abund_data$Abundance %<>% {. * 100}
					if("SE" %in% colnames(abund_data)) abund_data$SE %<>% {. * 100}
					if("SD" %in% colnames(abund_data)) abund_data$SD %<>% {. * 100}
					ylabname <- "Relative abundance (%)"
				}else{
					ylabname <- "Relative abundance"
				}
			} else {
				ylabname <- paste0("Relative abundance (", taxrank, " > ", show*100, "%)")
			}
			self$use_percentage <- use_percentage
			self$ylabname <- ylabname
			self$taxrank <- taxrank
			self$data_abund <- abund_data
			self$data_taxanames <- use_taxanames
			message('The transformed abundance data is stored in object$data_abund ...')
		},
		#' @description
		#' Bar plot.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(12, "Paired"); colors palette for the plotting.
		#' @param bar_type default "full"; "full" or "notfull"; if full, the total abundance sum to 1 or 100 percentage.
		#' @param others_color default "grey90"; the color for "others" taxa.
		#' @param facet default NULL; a character vector for the facet; a group column name of \code{sample_table}, such as, \code{"Group"};
		#'    If multiple facets are needed, please provide ordered names, such as \code{c("Group", "Type")}.
		#'    the latter facet should have a finer scale than the former one;
		#'    Please use factors in \code{sample_table} to adjust the facet orders in the plot;
		#'    When multiple facets are used, please first install package \code{ggh4x} using the command \code{install.packages("ggh4x")}.
		#' @param order_x default NULL; vector; used to order the sample names in x axis; must be the samples vector, such as, c("S1", "S3", "S2").
		#' @param x_axis_name NULL; a character string; a column name of sample_table in dataset; used to show the sample names in x axis.
		#' @param barwidth default NULL; bar width, see width in \code{\link{geom_bar}}.
		#' @param use_alluvium default FALSE; whether add alluvium plot
		#' @param clustering default FALSE; whether order samples by the clustering
		#' @param facet_color default "grey95"; facet background color.
		#' @param strip_text default 11; facet text size.
		#' @param legend_text_italic default FALSE; whether use italic in legend.
		#' @param xtext_type_hor default TRUE; x axis text horizontal, if FALSE; text slant.
		#' @param xtext_size default 10; x axis text size.
		#' @param xtext_keep default TRUE; whether retain x text.
		#' @param xtitle_keep default TRUE; whether retain x title.
		#' @param ytitle_size default 17; y axis title size.
		#' @param ylab_title default NULL; y axis title.
		#' @return ggplot2 plot. 
		#' @examples
		#' \donttest{
		#' t1$plot_bar(facet = "Group", xtext_keep = FALSE)
		#' }
		plot_bar = function(
			color_values = RColorBrewer::brewer.pal(12, "Paired"),
			bar_type = "full",
			others_color = "grey90",
			facet = NULL,
			order_x = NULL,
			x_axis_name = NULL,
			barwidth = NULL,
			use_alluvium = FALSE,
			clustering = FALSE,
			facet_color= "grey95",
			strip_text = 11,
			legend_text_italic = FALSE,
			xtext_type_hor = TRUE,
			xtext_size = 10,
			xtext_keep = TRUE,
			xtitle_keep = TRUE,
			ytitle_size = 17,
			ylab_title = NULL
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames

			if(bar_type == "full"){
				# make sure whether taxonomy info are all in selected use_taxanames in case of special data
				if(!all(plot_data$Taxonomy %in% use_taxanames)){
					plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
					new_data <- plot_data %>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
						dplyr::summarise(Abundance = sum(Abundance)) %>%
						as.data.frame(stringsAsFactors = FALSE)
					plot_data_merge <- plot_data[, ! colnames(plot_data) %in% c("Taxonomy", "Abundance", "N", "SD", "SE"), drop = FALSE] %>% unique
					plot_data <- dplyr::left_join(new_data, plot_data_merge, by = c("Sample" = "Sample"))
					plot_data$Taxonomy %<>% factor(., levels = rev(c(use_taxanames, "Others")))
				}else{
					plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
				}
			}else{
				plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
				plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
			}
			# order x axis samples
			plot_data <- private$adjust_axis_facet(
				plot_data = plot_data, 
				x_axis_name = x_axis_name, 
				order_x = order_x
				)

			# arrange plot_data--Abundance according to the Taxonomy-group column factor-levels
			plot_data <- plot_data[unlist(lapply(levels(plot_data$Taxonomy), function(x) which(plot_data$Taxonomy == x))),]
			bar_colors_use <- color_values[1:length(unique(plot_data$Taxonomy))]
			if(any(grepl("Others", as.character(plot_data$Taxonomy)))){
				bar_colors_use[length(bar_colors_use)] <- others_color
			}
			if(clustering){
				data_clustering <- reshape2::dcast(plot_data, Sample ~ Taxonomy, value.var = "Abundance", fun.aggregate = sum) %>% 
					`row.names<-`(.[,1]) %>% .[, -1]
				order_x_clustering <- hclust(dist(data_clustering)) %>% {.$labels[.$order]} %>% as.character
				plot_data$Sample %<>% factor(., levels = order_x_clustering)
			}
			if(use_alluvium){
				p <- ggplot(plot_data, aes(
						x = .data[["Sample"]], y = .data[["Abundance"]], 
						fill = .data[["Taxonomy"]], color = .data[["Taxonomy"]], 
						weight = .data[["Abundance"]], 
						alluvium = .data[["Taxonomy"]], stratum = .data[["Taxonomy"]]
					)) +
					ggalluvial::geom_flow(alpha = .4, width = 3/15) +
					ggalluvial::geom_stratum(width = .2) +
					scale_color_manual(values = rev(bar_colors_use))
			}else{
				p <- ggplot(plot_data, aes(x = .data[["Sample"]], y = .data[["Abundance"]], fill = .data[["Taxonomy"]]))
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
			p <- p + scale_fill_manual(values = rev(bar_colors_use)) + xlab("")
			if(!is.null(ylab_title)){
				p <- p + ylab(ylab_title)
			}else{
				p <- p + ylab(self$ylabname)		
			}
			if(!is.null(facet)){
				if(length(facet) == 1){
					p <- p + facet_grid(reformulate(facet, "."), scales = "free", space = "free")
				}else{
					p <- p + ggh4x::facet_nested(reformulate(facet), nest_line = element_line(linetype = 2), scales = "free", space = "free")
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
			
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size, xtext_keep = xtext_keep)

			p <- p + theme(axis.title.y = element_text(size = ytitle_size))
			if(xtitle_keep == F) {
				p <- p + theme(axis.title.x = element_blank())
			}
			p <- p + guides(fill = guide_legend(title = self$taxrank))
			if(use_alluvium){
				p <- p + guides(color = guide_legend(title = self$taxrank))
			}
			p
		},
		#' @description
		#' Plot the heatmap.
		#'
		#' @param color_values default rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")); 
		#' 	  colors palette for the plotting.
		#' @param facet default NULL; a character vector for the facet; a group column name of sample_table, such as, \code{"Group"};
		#'    If multiple facets are needed, please provide ordered names, such as \code{c("Group", "Type")}.
		#'    the latter facet should have a finer scale than the former one;
		#'    Please use factors in \code{sample_table} to adjust the facet orders in the plot;
		#'    use multiple facets, please first install package \code{ggh4x} using the command \code{install.packages("ggh4x")}.
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
		#' @param xtext_type_hor default TRUE; x axis text horizontal, if FALSE; text slant.
		#' @param legend_title default "\% Relative\\nAbundance"; legend title text.
		#' @param pheatmap default FALSE; whether use pheatmap package to plot the heatmap.
		#' @param ... paremeters pass to pheatmap when pheatmap = TRUE.
		#' @return ggplot2 plot or grid plot based on pheatmap.
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
			xtext_type_hor = TRUE,
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
				p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size, xtext_keep = xtext_keep)
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
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the plotting.
		#' @param group default NULL; a column name of sample table to show abundance across groups.
		#' @param show_point default FALSE; whether show points in plot.
		#' @param point_color default "black"; If show_point TRUE; use the color
		#' @param point_size default 3; If show_point TRUE; use the size
		#' @param point_alpha default .3; If show_point TRUE; use the transparency.
		#' @param plot_flip default FALSE; Whether rotate plot.
		#' @param boxfill default TRUE; Whether fill the box with colors.
		#' @param middlecolor default "grey95"; The middle line color.
		#' @param middlesize default 1; The middle line size.
		#' @param xtext_type_hor default TRUE; x axis text horizontal, if FALSE; text slant.
		#' @param xtext_size default 10; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @param ... parameters pass to \code{\link{geom_boxplot}}.
		#' @return ggplot2 plot. 
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
			xtext_type_hor = FALSE,
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
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size)
			p <- p + theme(axis.title.y = element_text(size = ytitle_size)) + scale_y_continuous(expand = c(0, 0.01))

			if(!is.null(group)) {
				p <- p + guides(fill=guide_legend(title=group))
			}
			p
		},
		#' @description
		#' Plot the line chart.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the plotting.
		#' @param plot_SE default TRUE; TRUE: plot the errorbar with mean±se; FALSE: plot the errorbar with mean±sd.
		#' @param position default position_dodge(0.1); Position adjustment, either as a string (such as "identity"), or the result of a call to a position adjustment function.
		#' @param errorbar_size default 1; errorbar size.
		#' @param errorbar_width default 0.1; errorbar width.
		#' @param point_size default 3; point size for taxa.
		#' @param point_alpha default 0.8; point transparency.
		#' @param line_size default 0.8; line size.
		#' @param line_alpha default 0.8; line transparency.
		#' @param line_type default 1; an integer; line type.
		#' @param xtext_type_hor default TRUE; x axis text horizontal, if FALSE; text slant.
		#' @param xtext_size default 10; x axis text size.
		#' @param ytitle_size default 17; y axis title size.
		#' @return ggplot2 plot. 
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 5)
		#' t1$plot_line(point_size = 3)
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 5, groupmean = "Group")
		#' t1$plot_line(point_size = 5, errorbar_size = 1, xtext_type_hor = TRUE)
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
			xtext_type_hor = FALSE,
			xtext_size = 10,
			ytitle_size = 17
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames

			plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
			plot_data$Taxonomy %<>% factor(., levels = use_taxanames)
						
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
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size)
			p <- p + theme(axis.title.y = element_text(size = ytitle_size)) + scale_y_continuous(expand = c(0, 0.01))
			p <- p + scale_color_manual(values = color_values)
			
			p
		},
		#' @description
		#' Plot pie chart.
		#'
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the plotting.
		#' @param facet_nrow default 1; how many rows in the plot.
		#' @param strip_text default 11; sample title size.
		#' @param legend_text_italic default FALSE; whether use italic in legend.
		#' @return ggplot2 plot. 
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
		#' t1$plot_pie(facet_nrow = 1)
		#' }
		plot_pie = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 
			facet_nrow = 1, 
			strip_text = 11, 
			legend_text_italic = FALSE
			){
			plot_data <- self$data_abund
			use_taxanames <- self$data_taxanames

			plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
			plot_data$Taxonomy %<>% factor(., levels = rev(c(use_taxanames, "Others")))

			p <- ggplot(plot_data, aes(x='', y=Abundance, fill = Taxonomy, label = percent(Abundance))) + 
				geom_bar(width = 1, stat = "identity") +
				coord_polar("y", start=0) +
				private$blank_theme +
				scale_fill_manual(values = color_values) +
				theme(axis.text.x=element_blank()) +
				facet_wrap(~Sample, nrow = facet_nrow) +
				theme(strip.text = element_text(size=strip_text)) +
				guides(fill = guide_legend(title=self$taxrank, reverse = TRUE))
			if(legend_text_italic == T) {
				p <- p + theme(legend.text = element_text(face = 'italic'))
			}
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
		ggplot_xtext_type = function(xtext_type_hor, xtext_size, xtext_keep = TRUE){
			if(xtext_keep){
				if(xtext_type_hor == T){
					theme(axis.text.x = element_text(colour = "black", size = xtext_size))
				}else{
					theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = xtext_size))
				}
			}else{
				theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
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
