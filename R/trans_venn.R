#' @title
#' Create \code{trans_venn} object for the Venn diagram, petal plot and UpSet plot.
#'
#' @description
#' This class is a wrapper for a series of intersection analysis related methods, including 2- to 5-way venn diagram, 
#' more than 5-way petal or UpSet plot and intersection transformations based on David et al. (2012) <doi:10.1128/AEM.01459-12>.
#'
#' @export
trans_venn <- R6Class(classname = "trans_venn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class or a matrix-like table (data.frame or matrix object).
		#' 	 If dataset is a matrix-like table, features must be rows.
		#' @param ratio default NULL; NULL, "numratio" or "seqratio"; "numratio": calculate the percentage of feature number; 
		#' 	 "seqratio": calculate the percentage of feature abundance; NULL: no additional percentage.
		#' @param name_joint default "&"; the joint mark for generating multi-sample names.
		#' @return \code{data_details} and \code{data_summary} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- dataset$merge_samples("Group")
		#' t1 <- trans_venn$new(dataset = t1, ratio = "numratio")
		#' }
		initialize = function(dataset, ratio = NULL, name_joint = "&"
			){
			if(is.null(dataset)){
				stop("The input dataset must be provided!")
			}
			if(inherits(dataset, "R6")){
				use_dataset <- clone(dataset)
				# filter the feature with abundance 0
				use_dataset$tidy_dataset()
				abund <- use_dataset$otu_table
				self$tax_table <- use_dataset$tax_table				
			}else{
				if(!any(is.data.frame(dataset), is.matrix(dataset))){
					stop("Input table must be data.frame or matrix class!")
				}
				abund <- dataset
			}
			res_names <- colnames(abund)
			colnumber <- ncol(abund)
			abund1 <- cbind.data.frame(OTU = rownames(abund), abund)
			abund2 <- reshape2::melt(abund1, id.var = "OTU", value.name = "Abundance", variable.name = "SeqID")
			# Create intersection matrix (removes duplicates)
			setunion <- rownames(abund)
			setmatrix <- abund
			setmatrix[setmatrix >= 1] <- 1
			samplesum <- apply(setmatrix, 2, sum)
			# Create all possible sample combinations within requested complexity levels
			# modified from the method in systemPipeR package
			allcombl <- lapply(seq_len(colnumber), function(x) combn(colnames(setmatrix), m = x, simplify = FALSE)) %>% 
				unlist(recursive = FALSE)
			venn_list <- sapply(seq_along(allcombl), function(x) private$vennSets(setmatrix = setmatrix, allcombl = allcombl, index = x, setunion = setunion))
			names(venn_list) <- sapply(allcombl, paste, collapse = name_joint)
			venn_abund <- sapply(venn_list, function(x){
				subset(abund2, OTU %in% x) %>% 
				dplyr::summarise(Sum = sum(Abundance)) %>% 
				unlist %>% 
				unname
			})
			venn_count_abund <- data.frame(Counts = sapply(venn_list, length), Abundance = venn_abund)
			if(!is.null(ratio)){
				if(!ratio %in% c("seqratio", "numratio")){
					stop("Provided parameter ratio must be one of NULL, 'seqratio' or 'numratio' !")
				}
				if(ratio == "seqratio"){
					venn_count_abund[, 2] <- paste0(round(venn_count_abund[,2]/sum(venn_count_abund[, 2]), 3) * 100, "%")
				}else{
					venn_count_abund[, 2] <- paste0(round(venn_count_abund[,1]/sum(venn_count_abund[, 1]), 3) * 100, "%")
				}
			}
			# make the length of elements same
			venn_maxlen <- max(sapply(venn_list, length))
			venn_table <- lapply(venn_list, function(x) {
				fill_length <- venn_maxlen - length(x)
				c(x, rep("", fill_length))
			})
			venn_table <- as.data.frame(t(do.call(rbind, venn_table)))
			self$data_summary <- venn_count_abund
			self$data_details <- venn_table
			self$data_samplesum <- samplesum
			self$colnumber <- colnumber
			self$res_names <- res_names
			self$ratio <- ratio
			self$otu_table <- abund
			self$name_joint <- name_joint
			message('The details of each venn part is stored in object$data_details ...')
			message('The venn summary table used for plot is stored in object$data_summary ...')
		},
		#' @description
		#' Plot venn diagram.
		#'
		#' @param color_circle default \code{RColorBrewer::brewer.pal(8, "Dark2")}; color pallete.
		#' @param fill_color default TRUE; whether fill the area color.
		#' @param text_size default 4.5; text size in plot.
		#' @param text_name_size default 6; name size in plot.
		#' @param text_name_position default NULL; name position in plot.
		#' @param alpha default .3; alpha for transparency.
		#' @param linesize default 1.1; cycle line size.
		#' @param petal_plot default FALSE; whether use petal plot.
		#' @param petal_color default "#BEAED4"; color of the petals; If petal_color only has one color value, all the petals will be assigned with this color value.
		#'    If petal_color has multiple colors, and the number of color values is smaller than the petal number, 
		#'    the function can append more colors automatically with the color interpolation.
		#' @param petal_color_center default "#BEBADA"; color of the center in the petal plot.
		#' @param petal_a default 4; the length of the ellipse.
		#' @param petal_r default 1; scaling up the size of the ellipse.
		#' @param petal_use_lim default c(-12, 12); the width of the plot.
		#' @param petal_center_size default 40; petal center circle size.
		#' @param petal_move_xy default 4; the distance of text to circle.
		#' @param petal_move_k default 2.3; the distance of title to circle.
		#' @param petal_move_k_count default 1.3; the distance of data text to circle.
		#' @param petal_text_move default 40; the distance between two data text.
		#' @param other_text_show default NULL; other characters used to show in the plot.
		#' @param other_text_position default c(1, 1); the text position for text in \code{other_text_show}.
		#' @param other_text_size default 5; the text size for text in \code{other_text_show}.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_venn()
		#' }
		plot_venn = function(
			color_circle = RColorBrewer::brewer.pal(8, "Dark2"),
			fill_color = TRUE,
			text_size = 4.5,
			text_name_size = 6,
			text_name_position = NULL,
			alpha = 0.3,
			linesize = 1.1,
			petal_plot = FALSE,
			petal_color = "#BEAED4",
			petal_color_center = "#BEBADA",
			petal_a = 4,
			petal_r = 1,
			petal_use_lim = c(-12, 12),
			petal_center_size = 40,
			petal_move_xy = 4,
			petal_move_k = 2.3,
			petal_move_k_count = 1.3,
			petal_text_move = 40,
			other_text_show = NULL,
			other_text_position = c(2, 2),
			other_text_size = 5
			){
			colnumber <- self$colnumber
			ratio <- self$ratio
			res_names <- self$res_names
			switch_num <- colnumber - 1
			summary_table <- self$data_summary

			if(colnumber > 5 & petal_plot == F){
				message("The number of elements is larger than 5! Automatically change petal_plot = TRUE! An alternative way of visualization is to use plot_bar function ...")
				petal_plot <- TRUE
			}
			# text position in venn
			if(is.null(text_name_position)){
				text_name_position <- switch(switch_num, 
					data.frame(x = c(1.5, 8.5), y = c(6, 6)),
					data.frame(x = c(2, 8, 5), y = c(7.9, 7.9, 1.6)),
					data.frame(x = c(1, 2.6, 6.8, 9), y = c(7.4, 8.2, 8.2,7.4)),
					data.frame(x = c(4.8, 9.2, 8.8, 1.65, 0.72), y = c(10.6, 7.7, 0.3, 0.2, 7.05))
				)
			}
			if(colnumber %in% 2:5 & petal_plot == F){
				plot_data <- data.frame(summary_table, private$pos_fun(switch_num))
			}
			if(colnumber == 2) {
				p <- ggplot(data.frame(), aes(x = c(5, 5), y = 0)) + 
					xlim(1, 9) + 
					ylim(2.5, 9) + 
					private$main_theme
				
				if(fill_color == T){
					p <- p + 
						geom_polygon(data = private$plotcircle(center = c(4, 6)), aes(x = x, y = y), fill=color_circle[1], alpha = alpha) +
						geom_polygon(data = private$plotcircle(center = c(6, 6)), aes(x = x, y = y), fill=color_circle[2], alpha = alpha)
				} else {
					p <- p +
						annotate("path", x = private$plotcircle(center = c(4, 6))$x, y = private$plotcircle(center = c(4, 6))$y, 
							color = color_circle[1], size = linesize) +
						annotate("path", x = private$plotcircle(center = c(6, 6))$x, y = private$plotcircle(center = c(6, 6))$y, 
							color = color_circle[2], size = linesize)
				}
			}
			if(colnumber == 3) {
				p <- ggplot(data.frame(), aes(x = c(5, 5), y = 0)) +
					xlim(1, 9) +	
					ylim(1, 9) + 
					private$main_theme
				
				if(fill_color == T){
					p <- p + 
					 geom_polygon(data = private$plotcircle(center = c(4, 6)), aes(x = x, y = y), fill = color_circle[1], alpha = alpha) +
					 geom_polygon(data = private$plotcircle(center = c(6, 6)), aes(x = x, y = y), fill = color_circle[2], alpha = alpha) +
					 geom_polygon(data = private$plotcircle(center = c(5, 4)), aes(x = x, y = y), fill = color_circle[3], alpha = alpha)
				} else {
					p <- p +
					annotate("path", x = private$plotcircle(center = c(4, 6))$x, y = private$plotcircle(center = c(4, 6))$y, 
						color = color_circle[1], size = linesize) +
					annotate("path", x = private$plotcircle(center = c(6, 6))$x, y = private$plotcircle(center = c(6, 6))$y, 
						color = color_circle[2], size = linesize) +
					annotate("path", x = private$plotcircle(center = c(5, 4))$x, y = private$plotcircle(center = c(5, 4))$y, 
						color = color_circle[3], size = linesize)
				}
			}
			
			if(colnumber == 4) {
				p <- ggplot(data.frame(), aes(x = c(5,5), y = 0)) + 
					xlim(0, 10) + 
					ylim(0, 10) + 
					private$main_theme
				
				map_data_1 <- private$plotellipse(center = c(3.5, 3.6), rotate = -35)
				map_data_2 <- private$plotellipse(center = c(4.7, 4.4), rotate = -35)
				map_data_3 <- private$plotellipse(center = c(5.3, 4.4), rotate = 35)
				map_data_4 <- private$plotellipse(center = c(6.5, 3.6), rotate = 35)
			
				if(fill_color == T){
					p <- p + 
						geom_polygon(data = map_data_1, aes(x = x, y = y), fill = color_circle[1], alpha = alpha) +
						geom_polygon(data = map_data_2, aes(x = x, y = y), fill = color_circle[2], alpha = alpha) +
						geom_polygon(data = map_data_3, aes(x = x, y = y), fill = color_circle[3], alpha = alpha) +
						geom_polygon(data = map_data_4, aes(x = x, y = y), fill = color_circle[4], alpha = alpha)
				} else {
					p <- p + 
						annotate("path", x = map_data_1$x, y = map_data_1$y, color = color_circle[1], size = linesize) +
						annotate("path", x = map_data_2$x, y = map_data_2$y, color = color_circle[2], size = linesize) +
						annotate("path", x = map_data_3$x, y = map_data_3$y, color = color_circle[3], size = linesize) +
						annotate("path", x = map_data_4$x, y = map_data_4$y, color = color_circle[4], size = linesize)
				}
			}
			if(colnumber == 5 & petal_plot == F) {
				p <- ggplot(data.frame(), aes(x = c(5, 5), y = 0)) + 
					xlim(0, 10.4) + 
					ylim(-0.5, 10.8) + 
					private$main_theme
				
				map_data_1 <- private$plotellipse(center = c(4.83, 6.2), radius = c(1.43, 4.11), rotate = 0)
				map_data_2 <- private$plotellipse(center = c(6.25, 5.4), radius = c(1.7, 3.6), rotate = 66)
				map_data_3 <- private$plotellipse(center = c(6.1, 3.5), radius = c(1.55, 3.9), rotate = 150)
				map_data_4 <- private$plotellipse(center = c(4.48, 3.15), radius = c(1.55, 3.92), rotate = 210)
				map_data_5 <- private$plotellipse(center = c(3.7, 4.8), radius = c(1.7, 3.6), rotate = 293.5)
			
				if(fill_color == T){
					p <- p + 
						geom_polygon(data = map_data_1, aes(x = x, y = y), fill=color_circle[1], alpha = alpha)+
						geom_polygon(data = map_data_2, aes(x = x, y = y), fill=color_circle[2], alpha = alpha)+
						geom_polygon(data = map_data_3, aes(x = x, y = y), fill=color_circle[3], alpha = alpha)+
						geom_polygon(data = map_data_4, aes(x = x, y = y), fill=color_circle[4], alpha = alpha)+
						geom_polygon(data = map_data_5, aes(x = x, y = y), fill=color_circle[5], alpha = alpha)
				} else {
					p <- p + 
						annotate("path", x = map_data_1$x, y = map_data_1$y, color = color_circle[1], size = linesize) +
						annotate("path", x = map_data_2$x, y = map_data_2$y, color = color_circle[2], size = linesize) +
						annotate("path", x = map_data_3$x, y = map_data_3$y, color = color_circle[3], size = linesize) +
						annotate("path", x = map_data_4$x, y = map_data_4$y, color = color_circle[4], size = linesize) +
						annotate("path", x = map_data_5$x, y = map_data_5$y, color = color_circle[5], size = linesize)
				}
			}
			if(colnumber %in% 2:5 & petal_plot == F){
				p <- p + annotate("text", x = text_name_position$x, y = text_name_position$y, label = res_names, size = text_name_size)
				if(!is.null(ratio)){
					p <- p + annotate("text", 
							x = plot_data[, 3], 
							y = plot_data[, 4], 
							label = c(paste(plot_data[, 1], "\n(", plot_data[, 2],")", sep = "")), 
							size = text_size
							)
				}else{
					p <- p + annotate("text", 
							x = plot_data[, 3], 
							y = plot_data[, 4], 
							label = plot_data[, 1], 
							size = text_size
							)
				}
			}
			if(colnumber > 4 & petal_plot == T) {
				nPetals <- colnumber
				plot_data <- summary_table[c(1:nPetals, nrow(summary_table)), ]
				if(length(petal_color) == 1){
					petal_color_use <- rep(petal_color, nPetals)
				}else{
					petal_color_use <- expand_colors(petal_color, nPetals)
				}
				
				p <- ggplot(data.frame(), aes(x=c(0, 0), y = 0)) +
					  xlim(petal_use_lim[1], petal_use_lim[2]) +
					  ylim(petal_use_lim[1], petal_use_lim[2]) +
					  private$main_theme
					  
				for(i in 1:nPetals){
					rotate <- 90 - (i - 1) * 360/nPetals
					rotate2 <- rotate * pi/180
					if(rotate < -90){
						rotate <- rotate + 180
					}
					mx <- petal_move_xy * cos(rotate2)
					my <- petal_move_xy * sin(rotate2)
					petal_data <- private$petal(mx = mx, my = my, rotate = rotate, a = petal_a, r = petal_r)
					p <- p + geom_polygon(data = petal_data, aes(x = x, y = y), fill = petal_color_use[i], alpha = alpha)
					p <- p + annotate("text", x = petal_move_k * mx, y = petal_move_k * my, label = rownames(plot_data)[i], size = text_name_size)
					p <- p + annotate("text", x = petal_move_k_count * mx, y = petal_move_k_count * my, label = plot_data[i, 1], size = text_size)
					if(!is.null(ratio)){
						p <- p + annotate("text", 
							x = petal_move_k_count * mx, 
							y = petal_move_k_count * my - sum(abs(petal_use_lim))/petal_text_move, 
							label = plot_data[i, 2], 
							size = text_size)
					}
				}
				p <- p + geom_point(aes(x = 0, y = 0), shape = 16, size = petal_center_size, colour = petal_color_center)
				p <- p + annotate("text", 
					x = 0, 
					y = 0 + sum(abs(petal_use_lim))/(petal_text_move * 2), 
					label = plot_data[nrow(plot_data), 1], 
					size = text_size)
				if(!is.null(ratio)){
					p <- p + annotate("text", 
						x = 0, 
						y = 0 - sum(abs(petal_use_lim))/(petal_text_move * 2), 
						label = plot_data[nrow(plot_data), 2], 
						size = text_size)
				}
			}
			if(!is.null(other_text_show)){
				p <- p + annotate("text", 
					x = other_text_position[1], 
					y = other_text_position[2], 
					label = other_text_show, 
					size = other_text_size)				
			}
			p
		},
		#' @description
		#' Plot the intersections using histogram, i.e. UpSet plot. Especially useful when samples > 5.
		#'
		#' @param left_plot default TRUE; whether add the left bar plot to show the feature number of each sample.
		#' @param sort_samples default FALSE; \code{TRUE} is used to sort samples according to the number of features in each sample.
		#'   \code{FALSE} means the sample order is same with that in sample_table of the raw dataset.
		#' @param up_y_title default "Intersection set"; y axis title of upper plot.
		#' @param up_y_title_size default 15; y axis title size of upper plot.
		#' @param up_y_text_size default 4; y axis text size of upper plot.
		#' @param up_bar_fill default "grey70"; bar fill color of upper plot.
		#' @param up_bar_width default 0.9; bar width of upper plot.
		#' @param bottom_y_text_size default 12; y axis text size, i.e. sample name size, of bottom sample plot.
		#' @param bottom_height default 1; bottom plot height relative to the upper bar plot. 1 represents the height of bottom plot is same with the upper bar plot.
		#' @param bottom_point_size default 3; point size of bottom plot.
		#' @param bottom_point_color default "black"; point color of bottom plot.
		#' @param bottom_background_fill default "grey95"; fill color for the striped background in the bottom sample plot.
		#'    If the parameter length is 1, use "white" to distinguish the color stripes.
		#'    If the parameter length is greater than 1, use all provided colors.
		#' @param bottom_background_alpha default 1; the color transparency for the parameter \code{bottom_background_fill}.
		#' @param bottom_line_width default 0.5; the line width in the bottom plot.
		#' @param bottom_line_colour default "black"; the line color in the bottom plot.
		#' @param left_width default 0.3; left bar plot width relative to the right bottom plot.
		#' @param left_bar_fill default "grey70"; fill color for the left bar plot presenting feature number.
		#' @param left_bar_alpha default 1; the color transparency for the parameter \code{left_bar_fill}.
		#' @param left_bar_width default 0.9; bar width of left plot.
		#' @param left_x_text_size default 10; x axis text size of the left bar plot.
		#' @param left_background_fill default "white"; fill color for the striped background in the left plot.
		#'    If the parameter length is 1, use "white" to distinguish the color stripes.
		#'    If the parameter length is greater than 1, use all provided colors.
		#' @param left_background_alpha default 1; the color transparency for the parameter \code{left_background_fill}.
		#' @return a ggplot2 object.
		#' @examples
		#' \donttest{
		#' t2 <- t1$plot_bar()
		#' }
		plot_bar = function(
			left_plot = TRUE,
			sort_samples = FALSE,
			up_y_title = "Intersection size",
			up_y_title_size = 15,
			up_y_text_size = 8,
			up_bar_fill = "grey70",
			up_bar_width = 0.9,
			bottom_y_text_size = 12,
			bottom_height = 1,
			bottom_point_size = 3,
			bottom_point_color = "black",
			bottom_background_fill = "grey95",
			bottom_background_alpha = 1,
			bottom_line_width = 0.5,
			bottom_line_colour = "black",
			left_width = 0.3,
			left_bar_fill = "grey70",
			left_bar_alpha = 1,
			left_bar_width = 0.9,
			left_x_text_size = 10,
			left_background_fill = "white",
			left_background_alpha = 1
			){
			colnumber <- self$colnumber
			ratio <- self$ratio
			res_names <- self$res_names
			switch_num <- colnumber - 1
			summary_table <- self$data_summary
			name_joint <- self$name_joint
			samplesum <- self$data_samplesum
			
			if(any(grepl(name_joint, res_names, fixed = TRUE))){
				stop("Please change name_joint parameter when creating trans_venn object!")
			}
			if(sort_samples){
				sample_levels <- sort(samplesum, decreasing = TRUE) %>% names %>% rev
			}else{
				sample_levels <- res_names %>% rev
			}
			
			plot_data <- summary_table %>% 
				rownames_to_column %>% 
				.[order(.$Counts, decreasing = TRUE), ]
			plot_data[, 1] %<>% factor(., levels = .)
			
			g1 <- ggplot(plot_data, aes(x = rowname, y = Counts)) +
				theme_classic() +
				geom_col(color = up_bar_fill, fill = up_bar_fill, width = up_bar_width) +
				ylab(up_y_title) +
				theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
				theme(axis.text = element_text(size = up_y_text_size), axis.title = element_text(size = up_y_title_size))
			
			matrix_data <- matrix(nrow = colnumber, ncol = nrow(plot_data)) %>% as.data.frame
			rownames(matrix_data) <- sample_levels
			colnames(matrix_data) <- plot_data[, 1]
			for(i in colnames(matrix_data)){
				tmp <- strsplit(i, name_joint, fixed = TRUE) %>% unlist
				for(j in rownames(matrix_data)){
					if(j %in% tmp){
						matrix_data[j, i] <- 1
					}
				}
			}
			sample_long <- matrix_data %>% rownames_to_column %>% reshape2::melt(., id.vars = "rowname")
			sample_long$variable %<>% factor(., levels = levels(plot_data[, 1]))
			sample_long$rowname %<>% factor(levels = sample_levels)
			sample_ture <- sample_long[!is.na(sample_long$value), ]
			sample_all <- sample_long
			sample_all$value <- 1

			g2 <- ggplot(sample_long, aes(x = variable, y = rowname))
				#theme(plot.background = element_rect(fill = "white", colour = "white", size = 0.1))
			number_sample <- length(unique(sample_ture$rowname))
			if(length(bottom_background_fill) == 1){
				# color strip type
				use_bottom_background_fill <- rep(c(bottom_background_fill, "white"), ceiling(number_sample / 2))
			}else{
				use_bottom_background_fill <- expand_colors(bottom_background_fill, number_sample)
			}
			for(i in seq_len(number_sample)){
				g2 <- g2 + annotate("rect", ymin = i - 0.5, ymax = i + 0.5, xmin = -Inf, xmax = Inf, fill = use_bottom_background_fill[i], 
					alpha = bottom_background_alpha)
			}
			g2 <- g2 + 
				geom_point(aes(x = variable, y = rowname), data = sample_all, size = bottom_point_size, color = "grey92", inherit.aes = FALSE) +
				geom_point(aes(x = variable, y = rowname), data = sample_ture, size = bottom_point_size, color = bottom_point_color, inherit.aes = FALSE) +
				theme_bw() +
				theme(legend.position = "none") +
				theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
				theme(axis.text = element_text(size = bottom_y_text_size)) +
				theme(panel.border = element_blank()) +
				theme(panel.grid = element_blank())
			
			line_data <- matrix_data
			line_data[] <- lapply(line_data, function(x){x <- 1:length(x); x})
			line_data[is.na(matrix_data)] <- NA

			line_data2 <- data.frame(y = apply(line_data, 2, min, na.rm = TRUE), yend = apply(line_data, 2, max, na.rm = TRUE), x = 1:ncol(line_data), xend = 1:ncol(line_data))
			g2 <- g2 + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = line_data2, linewidth = bottom_line_width, colour = bottom_line_colour)

			if(left_plot){
				g3_data <- data.frame(number = samplesum, rowname = names(samplesum))
				g3_data$rowname %<>% factor(levels = sample_levels)

				g3 <- ggplot(g3_data, aes(x = number, y = rowname))
				if(length(left_background_fill) == 1){
					# color strip type
					use_left_background_fill <- rep(c(left_background_fill, "white"), ceiling(number_sample / 2))
				}else{
					use_left_background_fill <- expand_colors(left_background_fill, number_sample)
				}
				for(i in seq_len(number_sample)){
					g3 <- g3 + annotate("rect", ymin = i - 0.5, ymax = i + 0.5, xmin = -Inf, xmax = Inf, fill = use_left_background_fill[i], 
						alpha = left_background_alpha)
				}
				if(length(left_bar_fill) == 1){
					g3 <- g3 + geom_col(fill = left_bar_fill, width = left_bar_width, alpha = left_bar_alpha)
				}else{
					use_left_bar_fill <- expand_colors(left_bar_fill, number_sample)
					g3 <- g3 + geom_col(aes(fill = rowname), width = left_bar_width, alpha = left_bar_alpha) +
						scale_fill_manual(values = use_left_bar_fill)
				}
				g3 <- g3 + 
					theme_bw() +
					theme(legend.position = "none") +
					scale_x_reverse() +
					theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
					theme(axis.text = element_text(size = left_x_text_size)) +
					theme(panel.border = element_blank()) +
					theme(panel.grid = element_blank())
				
				p1 <- g2 %>% 
					aplot::insert_top(g1, height = 1/bottom_height) %>% 
					aplot::insert_left(g3, width = left_width)
				p1
			}else{
				p1 <- aplot::insert_bottom(g1, g2, height = bottom_height)
				p1
			}
		},
		#' @description
		#' Transform intersection result to community-like microtable object for further composition analysis.
		#'
		#' @param use_frequency default TRUE; whether only use OTUs occurrence frequency, i.e. presence/absence data; if FALSE, use abundance data.
		#' @return a new \code{\link{microtable}} class.
		#' @examples
		#' \donttest{
		#' t2 <- t1$trans_comm(use_frequency = TRUE)
		#' }
		trans_comm = function(use_frequency = TRUE){
			otudata <- self$otu_table
			venn_table <- self$data_details
			sampledata <- data.frame(SampleID = colnames(venn_table), Group = colnames(venn_table)) %>% 'rownames<-'(colnames(venn_table))
			taxdata <- self$tax_table
			sum_table <- data.frame(apply(otudata, 1, sum))
			tt <- dplyr::full_join(rownames_to_column(sum_table[venn_table[,1] %>% as.character %>% .[. != ""], ,drop = FALSE]),
				rownames_to_column(sum_table[venn_table[,2] %>% as.character %>% .[. != ""], , drop = FALSE]), 
				by=c("rowname" = "rowname"))
			for(i in 3:ncol(venn_table)){
				tt <- dplyr::full_join(tt, 
						rownames_to_column(sum_table[venn_table[, i] %>% as.character %>% .[. != ""], , drop=FALSE]), 
						by=c("rowname" = "rowname"))
			}
			tt[is.na(tt)] <- 0
			tt %<>% 'rownames<-'(.[, 1]) %>% .[, -1, drop = FALSE]
			colnames(tt) <- colnames(venn_table)
			if(use_frequency == T){
				tt[tt != 0] <- 1
			}
			microtable$new(sample_table = sampledata, otu_table = tt, tax_table = taxdata, auto_tidy = TRUE)
		},
		#' @description
		#' Print the trans_venn object.
		print = function() {
			print(self$data_summary)
			invisible(self)
		}
	),
	private = list(
		# modified from vennSets function in systemPipeR package
		vennSets = function(setmatrix, allcombl, index, setunion){
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		},
		# fix the position for 2-5 way
		pos_fun = function(num){
			switch(num, 
				data.frame(x = c(3.1, 7, 5), y = c(6, 6, 6)),
				data.frame(x = c(3, 7, 5, 5, 3.8, 6.3, 5), y = c(6.5, 6.5, 3, 6.8, 4.6, 4.6, 5.3)),
				data.frame(
					x = c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5, 5, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5), 
					y = c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6, 2.2, 5.9, 4, 1.4, 1.4, 4, 2.8)
				),
				data.frame(
					x = c(4.85, 8, 7.1, 3.5, 2, 5.9, 4.4, 4.6, 3.6, 7.2, 6.5, 3.2, 5.4, 6.65, 3.4, 5, 6.02, 3.6, 5.2, 4.03, 4.2, 
						6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 4, 5.2, 6.4, 5.1), 
					y = c(8.3, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.4, 4.4, 6, 4.6, 2.1, 3.4, 3.25, 6.43, 6.38, 5.1, 2.49, 6.25, 
						3.08, 5.3, 4, 3.8, 3.2, 5.95, 5.75, 3.75, 3, 4.5, 4.6)
				)
			)
		},
		# Circle function for 2 or 3-way
		plotcircle = function(center = c(1, 1), diameter = 4, segments_split = 360) {
			r <- diameter / 2
			tt <- seq(0, 2*pi, length.out = segments_split)
			xx <- center[1] + r * cos(tt)
			yy <- center[2] + r * sin(tt)
			data.frame(x = xx, y = yy)
		},
		# Ellipse function for 4 or 5-way
		plotellipse = function(center = c(1, 1), radius = c(2, 4), rotate = 1, segments_split = 360) {
			angles <- (0:segments_split) * 2 * pi/segments_split
			rotate <- rotate * pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind(ellipse[, 1] * cos(rotate) + ellipse[, 2] * sin(rotate), ellipse[, 2] * cos(rotate) - ellipse[, 1] * sin(rotate))
			ellipse <- cbind(center[1] + ellipse[, 1], center[2] + ellipse[, 2])
			colnames(ellipse) <- c("x", "y")
			as.data.frame(ellipse)
		},
		# inspired by the code from Xu brother
		petal = function(r = 1, n = 1000, a = 4, b = 1.2, mx = 0, my = 0, rotate = 0){
			ang <- seq(0, 360, len = n+1)
			ang <- ang[1:n]
			ang <- ang*pi/180
			x <- r * cos(ang)
			y <- r * sin(ang)
			xy <- rbind(x, y)
			m <- diag(c(a, b))
			xy <- m %*% xy
			rotate <- rotate*pi/180
			m <- c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate))
			dim(m) <- c(2, 2)
			xy <- m %*% xy
			xy[1, ] <- xy[1, ] + mx
			xy[2, ] <- xy[2, ] + my
			xy <- as.data.frame(t(xy))
			colnames(xy) <- c("x", "y")
			xy
		},
		main_theme = theme(panel.grid.major=element_blank(), 
			panel.grid.minor=element_blank(), 
			axis.text=element_blank(),
			axis.title=element_blank(),
			axis.ticks=element_blank(),
			panel.border=element_blank(),
			panel.background = element_blank(),
			legend.key = element_blank(),
			plot.margin = unit(c(0,0,0,0), "mm")
		)
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
