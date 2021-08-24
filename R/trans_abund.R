#' @title
#' Create trans_abund object to transform taxonomic abundance for plotting.
#'
#' @description
#' This class is a wrapper for the taxonomic abundance transformations and plotting. 
#' The transformed data style is the long-format for ggplot2 plotting.
#' The plotting approaches include the bar plot, boxplot, heatmap and pie chart based on An et al. (2019) <doi:10.1016/j.geoderma.2018.09.035>.
#'
#' @export
trans_abund <- R6Class(classname = "trans_abund",
	public = list(
		#' @param dataset default NULL; microtable object.
		#' @param taxrank default "Phylum"; taxonomic rank.
		#' @param show default 0; the relative abundance threshold used for filtering.
		#' @param ntaxa default 10; how many taxa will be used, ordered by abundance from high to low; 
		#'   this parameter does not conflict with the parameter show; both can be used.
		#' @param groupmean default NULL; calculating mean abundance for each group, select a group column name in sample_table.
		#' @param delete_full_prefix default TRUE; whether delete both the prefix and the character in front of them.
		#' @param delete_part_prefix default FALSE; whether only delete the prefix.
		#' @param prefix default NULL; character string; can be used when delete_full_prefix = T or delete_part_prefix = T; 
		#'   default NULL reprensents using the "letter+__", e.g. "k__" for Phylum level.
		#' @param use_percentage default TRUE; show the abundance percentage.
		#' @param input_taxaname default NULL; character vector; if some taxa are selected, input taxa names.
		#' @return abund_data for plotting. 
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
				stop("The taxonomic abundance has not been calculated! Please first calculate it using cal_abund() function in microtable class!")
			}
			self$sample_table <- dataset$sample_table
			self$taxrank <- taxrank
			self$use_percentage <- use_percentage
			
			abund_data <- dataset$taxa_abund[[self$taxrank]] %>% 
				cbind.data.frame(Taxonomy = rownames(.), ., stringsAsFactors = FALSE)
			abund_data <- reshape2::melt(abund_data, id.vars = "Taxonomy")
			colnames(abund_data) <- c("Taxonomy", "Sample", "Abundance")
			if(any(grepl("__$", abund_data$Taxonomy))){
				abund_data$Taxonomy[grepl("__$", abund_data$Taxonomy)] <- paste0(abund_data$Taxonomy[grepl("__$", abund_data$Taxonomy)], "unidentified")
			}
			if(delete_full_prefix == T | delete_part_prefix == T){
				if(is.null(prefix)){
					prefix <- paste0(tolower(substr(self$taxrank, 1, 1)), "__")
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
			# calculate mean vlaues for each group
			if(!is.null(groupmean)){
				abund_data$Sample %<>% as.character
				mdf <- suppressWarnings(dplyr::left_join(abund_data, rownames_to_column(self$sample_table), by=c("Sample" = "rowname")))
				message(paste0(groupmean, " column is used to calculate mean abundance ..."))
				abund_data <- mdf %>% dplyr::group_by(!!! syms(c("Taxonomy", groupmean))) %>% 
					dplyr::summarise(Abundance = mean(Abundance)) %>% 
					as.data.frame
				colnames(abund_data)[colnames(abund_data) == groupmean] <- "Sample"
			}else{
				abund_data <- suppressWarnings(dplyr::left_join(abund_data, rownames_to_column(self$sample_table), by=c("Sample" = "rowname")))
			}
			abund_data$Taxonomy %<>% as.character
			# sort according to the abundance
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
			use_taxanames <- use_taxanames[!grepl("unidentified|unculture|Incertae.sedis", use_taxanames)]
			# identify used taxa
			if(is.null(input_taxaname)){
				if(length(use_taxanames) > ntaxa_use) use_taxanames <- use_taxanames[1:ntaxa_use]
			} else {
				use_taxanames <- use_taxanames[use_taxanames %in% input_taxaname]
			}
			
			if(ntaxa_theshold < sum(mean_abund > show) | show == 0){
				if(self$use_percentage == T){
					abund_data$Abundance %<>% {. * 100}
					ylabname <- "Relative abundance (%)"
				}else{
					ylabname <- "Relative abundance"
				}
			} else {
				ylabname <- paste0("Relative abundance (", self$taxrank," > ", show*100, "%)")
			}
			self$ylabname <- ylabname
			self$abund_data <- abund_data
			self$use_taxanames <- use_taxanames
		},
		#' @description
		#' Bar plot in trans_abund object.
		#'
		#' @param use_colors default RColorBrewer::brewer.pal(12, "Paired"); providing the plotting colors.
		#' @param bar_type default "full"; "full" or "notfull"; if full, the total abundance sum to 1 or 100 percentage.
		#' @param others_color default "grey90"; the color for "others" taxa.
		#' @param facet default NULL; a character string; if using facet, providing a group column name of sample_table, such as, "Group".
		#' @param order_facet NULL; vector; used to order the facet, such as, c("Group1", "Group3", "Group2").
		#' @param x_axis_name NULL; a character string; a column name of sample_table used to show the sample names in x axis.
		#' @param order_x default NULL; vector; used to order the sample names in x axis; must be the samples vector, such as, c("S1", "S3", "S2").
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
		#' @param base_font default NULL; ggplot font family in the plot.
		#' @param ylab_title default NULL; y axis title.
		#' @return ggplot2 plot. 
		#' @examples
		#' \donttest{
		#' t1$plot_bar(facet = "Group", xtext_keep = FALSE)
		#' }
		plot_bar = function(
			use_colors = RColorBrewer::brewer.pal(12, "Paired"),
			bar_type = "full",
			others_color = "grey90",
			facet = NULL,
			order_facet = NULL,
			x_axis_name = NULL,
			order_x = NULL,
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
			base_font = NULL,
			ylab_title = NULL
			){
			plot_data <- self$abund_data
			
			# order x axis samples and facet
			plot_data <- private$adjust_axis_facet(plot_data = plot_data, x_axis_name = x_axis_name, order_x = order_x, facet = facet, order_facet = order_facet)

			if(use_alluvium){
				bar_type <- "notfull"
			}
			if(bar_type == "full"){
				# make sure whether taxonomy info are all in selected use_taxanames in case of special data
				if(!all(plot_data$Taxonomy %in% self$use_taxanames)){
					plot_data$Taxonomy[!plot_data$Taxonomy %in% self$use_taxanames] <- "Others"
					plot_data$Taxonomy %<>% factor(., levels = rev(c(self$use_taxanames, "Others")))
				}else{
					plot_data$Taxonomy %<>% factor(., levels = rev(self$use_taxanames))
				}
			}else{
				plot_data %<>% {.[.$Taxonomy %in% self$use_taxanames, ]}
				plot_data$Taxonomy %<>% factor(., levels = rev(self$use_taxanames))
			}
			# arrange plot_data--Abundance according to the Taxonomy-group column factor-levels
			plot_data <- plot_data[unlist(lapply(levels(plot_data$Taxonomy), function(x) which(plot_data$Taxonomy == x))),]
			bar_colors_use <- use_colors[1:length(unique(plot_data$Taxonomy))]
			if(any(grepl("Others", as.character(plot_data$Taxonomy)))){
				bar_colors_use[length(bar_colors_use)] <- others_color
			}
			if(clustering){
				data_clustering <- reshape2::dcast(plot_data, Sample ~ Taxonomy, value.var = "Abundance") %>% 
					`row.names<-`(.[,1]) %>% .[, -1]
				order_x_clustering <- hclust(dist(data_clustering)) %>% {.$labels[.$order]} %>% as.character
				plot_data$Sample %<>% factor(., levels = order_x_clustering)
			}
			if(use_alluvium){
				p <- ggplot(plot_data, aes_string(x = "Sample", y = "Abundance", fill = "Taxonomy", color = "Taxonomy",  weight = "Abundance", 
					alluvium = "Taxonomy", stratum = "Taxonomy")) +
					ggalluvial::geom_flow(alpha = .4, width = 3/15) +
					ggalluvial::geom_stratum(width = .2) +
					scale_color_manual(values = rev(bar_colors_use)) +
					scale_y_continuous(expand = c(0, 0))
			}else{
				p <- ggplot(plot_data, aes_string(x = "Sample", y = "Abundance", fill = "Taxonomy"))
				if(bar_type == "full"){
					if(self$use_percentage == T){
						p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth) + 
							scale_y_continuous(expand = c(0, 0))
					}else{
						p <- p + geom_bar(stat = "identity", position = "fill", show.legend = T, width = barwidth)
						p <- p + scale_y_continuous(limits = c(0,1), expand = c(0, 0))
					}
				}else{
					p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth) + 
						scale_y_continuous(expand = c(0, 0))
				}
			}
			p <- p + scale_fill_manual(values = rev(bar_colors_use)) + xlab("")
			if(!is.null(ylab_title)){
				p <- p + ylab(ylab_title)
			}else{
				p <- p + ylab(self$ylabname)		
			}
			if(!is.null(facet)) {
				p <- p + facet_grid(reformulate(facet, "."), scales = "free", space = "free")
				p <- p + theme(strip.background = element_rect(fill = facet_color, color = facet_color), strip.text = element_text(size=strip_text))
				p <- p + scale_y_continuous(expand = c(0, 0.01))
			}
			p <- p + theme(panel.grid = element_blank(), panel.border = element_blank()) + 
				theme(axis.line.y = element_line(color = "grey60", linetype = "solid", lineend = "square"))
			if(legend_text_italic == T) {
				p <- p + theme(legend.text = element_text(face = 'italic'))
			}
			
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size)
			if(!xtext_keep){
				p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
			}
			p <- p + theme(axis.title.y= element_text(size=ytitle_size))
			if(xtitle_keep == F) {
				p <- p + theme(axis.title.x = element_blank())
			}
			if(!is.null(base_font)){
				p <- p + theme(text=element_text(family=base_font))
			}
			p <- p + guides(fill=guide_legend(title=self$taxrank))
			if(use_alluvium){
			p <- p + guides(color = guide_legend(title=self$taxrank))
			}
			p
		},
		#' @description
		#' Plot the heatmap in trans_abund object.
		#'
		#' @param use_colors default RColorBrewer::brewer.pal(12, "Paired"); providing the plotting colors.
		#' @param facet default NULL; a character string; if using facet, providing a group column name of sample_table, such as, "Group".
		#' @param order_facet NULL; vector; used to order the facet, such as, c("Group1", "Group3", "Group2").
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
		#' @param base_font default NULL; font in the plot.
		#' @return ggplot2 plot.
		#' @examples
		#' \donttest{
		#' t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
		#' t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
		#' }
		plot_heatmap = function(
			use_colors = c("#00008B", "#102D9B", "#215AAC", "#3288BD", "#66C2A5",  "#E6F598", "#FFFFBF", "#FED690", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"), 
			facet = NULL,
			order_facet = NULL,
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
			base_font =NULL
			){
			plot_data <- self$abund_data
			# order x axis samples and facet
			plot_data <- private$adjust_axis_facet(plot_data = plot_data, x_axis_name = x_axis_name, order_x = order_x, facet = facet, order_facet = order_facet)

			if (is.null(min_abundance)){
				min_abundance <- ifelse(min(plot_data$Abundance) > 0.001, min(plot_data$Abundance), 0.001)
			}
			if (is.null(max_abundance)){
				max_abundance <- max(plot_data$Abundance)
			}
			plot_data %<>% {.[.$Taxonomy %in% self$use_taxanames, ]}
			plot_data$Taxonomy %<>% factor(., levels = rev(self$use_taxanames))

			p <- ggplot(plot_data, aes_string(x = "Sample", y = "Taxonomy", label = formatC("Abundance", format = "f", digits = 1)))
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
				p <- p + scale_fill_gradientn(colours = use_colors, trans = plot_colorscale, na.value = "#00008B", limits = c(min_abundance, max_abundance))
			}else{
				p <- p + scale_fill_gradientn(colours = use_colors, trans = plot_colorscale, breaks=plot_breaks, na.value = "#00008B",
					limits = c(min_abundance, max_abundance))
			}
			if(!is.null(order_facet)) {
				plot_data[, facet] <- factor(plot_data[, facet], levels = unique(order_facet))
			}
			if(!is.null(facet)){
				p <- p + facet_grid(reformulate(facet, "."), scales = "free", space = "free")
				p <- p + theme(strip.background = element_rect(color = "white", fill = "grey92"), strip.text = element_text(size=strip_text))
			}
			p <- p + labs(x = "", y = "", fill = "% Relative\nAbundance")
			if (!is.null(ytext_size)){
				p <- p + theme(axis.text.y = element_text(size = ytext_size))
			}
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size)
			if(!xtext_keep) {
				p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
			}
			if(grid_clean){
				p <- p + theme(panel.border = element_blank(), panel.grid = element_blank())
			}			
			if(!is.null(base_font)){
				p <- p + theme(text=element_text(family=base_font))
			}
			p
		},
		#' @description
		#' Box plot in trans_abund object.
		#'
		#' @param use_colors default RColorBrewer::brewer.pal(12, "Paired"); providing the plotting colors.
		#' @param group default NULL; column name of sample table to show abundance across groups.
		#' @param show_point default FALSE; whether show points in plot.
		#' @param point_color default "black"; If show_point TRUE; use the color
		#' @param point_size default 3; If show_point TRUE; use the size
		#' @param point_alpha default .3; If show_point TRUE; use the transparency.
		#' @param plot_flip default FALSE; Whether rotate plot.
		#' @param boxfill default TRUE; Whether fill the box.
		#' @param middlecolor default "grey95"; The middle line color.
		#' @param middlesize default 1; The middle line size.
		#' @param xtext_type_hor default TRUE; x axis text horizontal, if FALSE; text slant.
		#' @param xtext_size default 10; x axis text size.
		#' @param xtext_keep default TRUE; whether retain x text.
		#' @param xtitle_keep default TRUE; whether retain x title.
		#' @param ytitle_size default 17; y axis title size.
		#' @param base_font default NULL; font in the plot.
		#' @param ... parameters pass to \code{\link{geom_boxplot}}.
		#' @return ggplot2 plot. 
		#' @examples
		#' \donttest{
		#' t1$plot_box(group = "Group")
		#' }
		plot_box = function(
			use_colors = RColorBrewer::brewer.pal(8, "Dark2"),
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
			xtext_keep = TRUE,
			xtitle_keep = TRUE,
			ytitle_size = 17,
			base_font =NULL,
			...
			){
			plot_data <- self$abund_data
			plot_data %<>% {.[.$Taxonomy %in% self$use_taxanames, ]}

			plot_data$Taxonomy %<>% factor(., levels = self$use_taxanames)

			p <- ggplot(plot_data, aes_string(x = "Taxonomy", y = "Abundance")) 
			p <- p + ylab(self$ylabname) + guides(col = guide_legend(reverse = TRUE)) + xlab("")
			if (plot_flip == T){ 
				p <- p + coord_flip()
			}
			if(is.null(group)) {
				p <- p + geom_boxplot(color = use_colors[1], ...)
			} else {
				if(boxfill == T){
					p <- p + geom_boxplot(aes_string(color = group, fill = group), ...)
					p <- p + scale_fill_manual(values = use_colors)
					p <- p + scale_color_manual(values = use_colors) + guides(color = FALSE)
					## Change the default middle line
					dat <- ggplot_build(p)$data[[1]]
					p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour = middlecolor, size=middlesize)
				} else {	 
					p <- p + geom_boxplot(aes_string(color = group), ...) + scale_color_manual(values = use_colors)
				}
			}
			if(show_point == T){
				p <- p + geom_point(size = point_size, color = point_color, alpha = point_alpha, position = "jitter")
			}			
			p <- p + private$ggplot_xtext_type(xtext_type_hor = xtext_type_hor, xtext_size = xtext_size)

			if(!xtext_keep){
				p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
			}
			p <- p + theme(axis.title.y = element_text(size = ytitle_size)) + scale_y_continuous(expand = c(0, 0.01))

			if(xtitle_keep == F) {
				p <- p + theme(axis.title.x = element_blank())
			}
			if(!is.null(base_font)){
				p <- p + theme(text=element_text(family=base_font))
			}
			if(!is.null(group)) {
				p <- p + guides(fill=guide_legend(title=group))
			}
			p
		},
		#' @description
		#' Plot pie chart in trans_abund class.
		#'
		#' @param use_colors default RColorBrewer::brewer.pal(8, "Dark2"); providing the plotting colors.
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
			use_colors = RColorBrewer::brewer.pal(8, "Dark2"), 
			facet_nrow = 1, 
			strip_text = 11, 
			legend_text_italic = FALSE
			){
			plot_data <- self$abund_data
			plot_data$Taxonomy[!plot_data$Taxonomy %in% self$use_taxanames] <- "Others"
			plot_data$Taxonomy %<>% factor(., levels = rev(c(self$use_taxanames, "Others")))

			p <- ggplot(plot_data, aes(x='', y=Abundance, fill = Taxonomy, label = percent(Abundance))) + 
				geom_bar(width = 1, stat = "identity") +
				coord_polar("y", start=0) +
				private$blank_theme +
				scale_fill_manual(values = use_colors) +
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
		print = function() {
			cat("trans_abund class:\n")
			cat(paste("abund_data have", ncol(self$abund_data), "columns: ", paste0(colnames(self$abund_data), collapse = ", "), "\n"))
			cat(paste("abund_data have", nrow(self$abund_data), "rows\n"))
			if(!is.null(self$use_taxanames)){
				if(length(self$use_taxanames) > 20){
					cat(paste("Filtered taxa: ", length(self$use_taxanames), "\n"))
				}else{
					cat(paste("Filtered taxa names: ", paste0(self$use_taxanames, collapse = ", "), "\n"))
				}
			}
			invisible(self)
		}
		),
	private = list(
		adjust_axis_facet = function(plot_data, x_axis_name, order_x, facet, order_facet){
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
			if(!is.null(order_facet)){
				if(is.null(facet)){
					stop("You provide order_facet. It is necessary to provide facet!")
				}else{
					plot_data[, facet] %<>% factor(., levels = order_facet)
				}
			}
			plot_data
		},
		ggplot_xtext_type = function(xtext_type_hor, xtext_size){
			if(xtext_type_hor == T){
				theme(axis.text.x = element_text(colour = "black", size = xtext_size))
			} else {
				theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = xtext_size))
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
