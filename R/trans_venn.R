#' @title
#' Create trans_venn object.
#'
#' @description
#' This class is a wrapper for a series of venn analysis related methods, including venn result, 2- to 5-way venn diagram, 
#' more than 5-way petal plot and venn result transformations based on David et al. (2012) <doi:10.1128/AEM.01459-12>.
#'
#' @export
trans_venn <- R6Class(classname = "trans_venn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param sample_names default NULL; if provided, filter the samples.
		#' @param ratio default numratio; NULL, "numratio" or "seqratio"; numratio: calculate number percentage; seqratio: calculate sequence percentage; NULL: no additional percentage.
		#' @return venn_table venn_count_abund stored in trans_venn object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- dataset$merge_samples(use_group = "Group")
		#' t1 <- trans_venn$new(dataset = t1, ratio = "numratio")
		#' }
		initialize = function(dataset = NULL, sample_names = NULL, ratio = "numratio"
			) {
			use_dataset <- clone(dataset)
			if(!is.null(sample_names)){
				use_dataset$sample_table %<>% .[rownames(.) %in% sample_names, ]
			}
			use_dataset$tax_table %<>% base::subset(use_dataset$taxa_sums() != 0)
			use_dataset$tidy_dataset()
			abund <- use_dataset$otu_table
			col_names <- colnames(abund)
			colnumber <- ncol(abund)
			abund1 <- cbind.data.frame(OTU = rownames(abund), abund)
			abund2 <- reshape2::melt(abund1, id.var = "OTU", value.name= "Abundance", variable.name = "SeqID")
			## Create intersect matrix (removes duplicates)
			setunion <- rownames(abund)
			setmatrix <- abund
			setmatrix[setmatrix >= 1] <- 1
			## Create all possible sample combinations within requested complexity levels
			allcombl <- lapply(1:colnumber, function(x) combn(colnames(setmatrix), m=x, simplify=FALSE)) %>% unlist(recursive=FALSE)
			venn_list <- sapply(seq_along(allcombl), function(x) private$vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x, setunion = setunion))
			names(venn_list) <- sapply(allcombl, paste, collapse= "-")
			venn_abund <- sapply(venn_list, function(x){
				subset(abund2, OTU %in% x) %>% 
				dplyr::summarise(Sum = sum(Abundance)) %>% 
				unlist() %>% 
				unname()
			})
			venn_count_abund <- data.frame(Counts = sapply(venn_list, length), Abundance = venn_abund)
			if(!is.null(ratio)) {
				if(ratio == "seqratio") {
					venn_count_abund[,2] <- paste0(round(venn_count_abund[,2]/sum(venn_count_abund[,2]), 3) * 100, "%")
				}
				if(ratio == "numratio") {
					venn_count_abund[,2] <- paste0(round(venn_count_abund[,1]/sum(venn_count_abund[,1]), 3) * 100, "%")
				}
			}
			venn_maxlen <- max(sapply(venn_list, length))
			venn_table <- lapply(venn_list, function (x) {
				fill_length <- venn_maxlen - length(x)
				c(x, rep("", fill_length))
			})
			venn_table <- as.data.frame(t(do.call(rbind, venn_table)))
			self$venn_count_abund <- venn_count_abund
			self$venn_table <- venn_table
			self$colnumber <- colnumber
			self$col_names <- col_names
			self$ratio <- ratio
			self$otu_table <- abund
			self$tax_table <- use_dataset$tax_table
		},
		#' @description
		#' Plot venn diagram.
		#'
		#' @param color_circle default RColorBrewer::brewer.pal(8, "Dark2"); color pallete
		#' @param fill_color default TRUE; whether fill the area color
		#' @param text_size default 4.5; text size in plot
		#' @param text_name_size default 6; name size in plot
		#' @param text_name_position default NULL; name position in plot
		#' @param alpha default .3; alpha for transparency
		#' @param linesize default 1.1; cycle line size
		#' @param petal_plot default FALSE; whether use petal plot.
		#' @param petal_color default "skyblue"; color of the petal
		#' @param petal_a default 4; the length of the ellipse
		#' @param petal_r default 1; scaling up the size of the ellipse
		#' @param petal_use_lim default c(-12, 12); the width of the plot
		#' @param petal_center_size default 40; petal center circle size
		#' @param petal_move_xy default 4; the distance of text to circle
		#' @param petal_move_k default 2.3; the distance of title to circle
		#' @param petal_move_k_count default 1.3; the distance of data text to circle
		#' @param petal_text_move default 40; the distance between two data text
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
			petal_color = "skyblue",
			petal_a = 4,
			petal_r = 1,
			petal_use_lim = c(-12, 12),
			petal_center_size = 40,
			petal_move_xy = 4,
			petal_move_k = 2.3,
			petal_move_k_count = 1.3,
			petal_text_move = 40
			){
			colnumber <- self$colnumber
			ratio <- self$ratio
			col_names <- self$col_names
			switch_num <- colnumber-1
			
			# text position in venn
			if(is.null(text_name_position)){
				text_name_position <- switch(switch_num, 
					data.frame(x = c(1.5, 8.5), y = c(6, 6)),
					data.frame(x = c(2, 8, 5), y = c(7.9, 7.9, 1.6)),
					data.frame(x=c(1, 2.6, 6.8, 9), y = c(7.4, 8.2, 8.2,7.4)),
					data.frame(x=c(4.8, 9.2, 8.8, 1.65, 0.72), y = c(10.6, 7.7, 0.3, 0.2, 7.05))
				)
			}
			if(colnumber %in% 2:5 & petal_plot == F){
				plot_data <- data.frame(self$venn_count_abund, private$pos_fun(switch_num))
			}
			if(colnumber == 2) {
				p <- ggplot(data.frame(), aes(x=c(5,5), y=0)) + xlim(1,9) + ylim(2.5,9) + private$main_theme
				if(fill_color == T){
				p <- p + 
					geom_polygon(data = private$plotcircle(center = c(4, 6)),aes(x = x,y = y), fill=color_circle[1], alpha = alpha) +
					geom_polygon(data = private$plotcircle(center = c(6, 6)),aes(x = x,y = y), fill=color_circle[2], alpha = alpha)
				} else {
				p <- p +
					annotate("path", x = private$plotcircle(center = c(4, 6))$x, y = private$plotcircle(center = c(4, 6))$y, 
						color = color_circle[1], size = linesize) +
					annotate("path", x = private$plotcircle(center = c(6, 6))$x, y = private$plotcircle(center = c(6, 6))$y, 
						color = color_circle[2], size = linesize)
				}
			}
			if(colnumber == 3) {
				p <- ggplot(data.frame(), aes(x=c(5,5), y=0)) + xlim(1,9) +	ylim(1,9) + private$main_theme
				if(fill_color == T){
					p <- p + 
					 geom_polygon(data = private$plotcircle(center = c(4, 6)),aes(x = x,y = y), fill=color_circle[1], alpha = alpha) +
					 geom_polygon(data = private$plotcircle(center = c(6, 6)),aes(x = x,y = y), fill=color_circle[2], alpha = alpha) +
					 geom_polygon(data = private$plotcircle(center = c(5, 4)),aes(x = x,y = y), fill=color_circle[3], alpha = alpha)
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
				p <- ggplot(data.frame(), aes(x=c(5,5), y=0)) + xlim(0,10) + ylim(0,10) + private$main_theme

				if(fill_color == T){
					p <- p + 
					 geom_polygon(data = private$plotellipse(center = c(3.5, 3.6), rotate = -35),aes(x = x,y = y), fill=color_circle[1], alpha = alpha)+
					geom_polygon(data = private$plotellipse(center = c(4.7, 4.4), rotate = -35),aes(x = x,y = y), fill=color_circle[2], alpha = alpha) +
					geom_polygon(data = private$plotellipse(center = c(5.3, 4.4), rotate = 35),aes(x = x,y = y), fill=color_circle[3], alpha = alpha) +
					geom_polygon(data = private$plotellipse(center = c(6.5, 3.6), rotate = 35),aes(x = x,y = y), fill=color_circle[4], alpha = alpha)
				} else {
					p <- p +
					annotate("path", x = private$plotellipse(center = c(3.5, 3.6), rotate = -35)$x, 
						y = private$plotellipse(center = c(3.5, 3.6), rotate = -35)$y, color = color_circle[1], size = linesize) +
					annotate("path", x = private$plotellipse(center = c(4.7, 4.4), rotate = -35)$x, 
						y = private$plotellipse(center = c(4.7, 4.4), rotate = -35)$y, color = color_circle[2], size = linesize) +
					annotate("path", x = private$plotellipse(center = c(5.3, 4.4), rotate = 35)$x, 
						y = private$plotellipse(center = c(5.3, 4.4), rotate = 35)$y, color = color_circle[3], size = linesize) +
					annotate("path", x = private$plotellipse(center = c(6.5, 3.6), rotate = 35)$x, 
						y = private$plotellipse(center = c(6.5, 3.6), rotate = 35)$y, color = color_circle[4], size = linesize)
				}
			}
			if(colnumber == 5 & petal_plot == F) {
				p <- ggplot(data.frame(), aes(x=c(5,5), y=0)) + xlim(0,10.4) + ylim(-0.5,10.8) + private$main_theme
				if(fill_color == T){
					p <- p + 
					 geom_polygon(data = private$plotellipse(center = c(4.83, 6.2), radius = c(1.43, 4.11), rotate = 0),
						aes(x = x,y = y), fill=color_circle[1], alpha = alpha)+
					 geom_polygon(data = private$plotellipse(center = c(6.25, 5.4), radius = c(1.7, 3.6), rotate = 66),
						aes(x = x,y = y), fill=color_circle[2], alpha = alpha)+
					 geom_polygon(data = private$plotellipse(center = c(6.1, 3.5), radius = c(1.55, 3.9), rotate = 150),
						aes(x = x,y = y), fill=color_circle[3], alpha = alpha)+
					 geom_polygon(data = private$plotellipse(center = c(4.48, 3.15), radius = c(1.55, 3.92), rotate = 210),
						aes(x = x,y = y), fill=color_circle[4], alpha = alpha)+
					 geom_polygon(data = private$plotellipse(center = c(3.7, 4.8), radius = c(1.7, 3.6), rotate = 293.5), 
						aes(x = x,y = y), fill=color_circle[5], alpha = alpha)
				} else {
					p <- p +
					annotate("path", x = private$plotellipse(center = c(4.83, 6.2),  radius = c(1.43, 4.11), rotate = 0)$x, 
						y = private$plotellipse(center = c(4.83, 6.2), 
						radius = c(1.43, 4.11), rotate = 0)$y, color = color_circle[1], size = linesize)    +
					annotate("path", x = private$plotellipse(center = c(6.25, 5.4),  radius = c(1.7, 3.6), rotate = 66)$x, 
						y = private$plotellipse(center = c(6.25, 5.4), 
						radius = c(1.7, 3.6), rotate = 66)$y, color = color_circle[2], size = linesize)    +
					annotate("path", x = private$plotellipse(center = c(6.1, 3.5), radius = c(1.55, 3.9), rotate = 150)$x, 
						y = private$plotellipse(center = c(6.1, 3.5), 
						radius = c(1.55, 3.9), rotate = 150)$y, color = color_circle[3], size = linesize)    +
					annotate("path", x = private$plotellipse(center = c(4.48, 3.15), radius = c(1.55, 3.92), rotate = 210)$x, 
						y = private$plotellipse(center = c(4.48, 3.15), 
						radius = c(1.55, 3.92), rotate = 210)$y, color = color_circle[4], size = linesize)    +
					annotate("path", x = private$plotellipse(center = c(3.7, 4.8),   radius = c(1.7, 3.6), rotate = 293.5)$x, 
						y = private$plotellipse(center = c(3.7, 4.8), 
						radius = c(1.7, 3.6), rotate = 293.5)$y, color = color_circle[5], size = linesize)
				}
			}
			if(colnumber %in% 2:5 & petal_plot == F){
				p <- p + annotate("text", x = text_name_position$x, y = text_name_position$y, label = col_names, size = text_name_size)
				if(!is.null(ratio)){
					p <- p + annotate("text", x = plot_data[,3], y = plot_data[,4], label = c(paste(plot_data[,1], "\n(", plot_data[,2],")", sep = "")), 
						size = text_size)
				}else{
					p <- p + annotate("text", x = plot_data[,3], y = plot_data[,4], label = plot_data[,1], size = text_size)
				}
			}
			if(colnumber > 4 & petal_plot == T) {
				nPetals <- colnumber
				plot_data <- self$venn_count_abund[c(1:nPetals, nrow(self$venn_count_abund)), ]
				petal_color_use <- rep(petal_color, nPetals)
				
				p <- ggplot(data.frame(), aes(x=c(0, 0), y=0)) +
					  xlim(petal_use_lim[1], petal_use_lim[2]) +
					  ylim(petal_use_lim[1], petal_use_lim[2]) +
					  private$main_theme
				for(i in 1:nPetals){
					rotate <- 90 - (i-1)*360/nPetals
					rotate2 <- rotate*pi/180
					if(rotate < -90){
						rotate <- rotate + 180
					}
					mx <- petal_move_xy * cos(rotate2)
					my <- petal_move_xy * sin(rotate2)
					petal_data <- private$petal(mx = mx, my = my, rotate = rotate, a = petal_a, r = petal_r)

					p <- p + geom_polygon(data = petal_data, aes(x = x, y = y), fill = petal_color_use[i], alpha = alpha)
					p <- p + annotate("text", x = petal_move_k * mx, y = petal_move_k * my, label = rownames(plot_data)[i], size = text_name_size)
					p <- p + annotate("text", x = petal_move_k_count * mx, y = petal_move_k_count *my, label = plot_data[i, 1], size = text_size)
					if(!is.null(ratio)){
						p <- p + annotate("text", x = petal_move_k_count * mx, y = petal_move_k_count *my - sum(abs(petal_use_lim))/petal_text_move, 
							label = plot_data[i, 2], size = text_size)
					}
				}
				p <- p + geom_point(aes(x = 0, y = 0), shape = 16, size = petal_center_size, colour = petal_color)
				p <- p + annotate("text", x = 0, y = 0 + sum(abs(petal_use_lim))/(petal_text_move*2), label = plot_data[nrow(plot_data), 1], size = text_size)
				if(!is.null(ratio)){
					p <- p + annotate("text", x = 0, y = 0 - sum(abs(petal_use_lim))/(petal_text_move*2), label = plot_data[nrow(plot_data), 2], size = text_size)
				}
			}
			return(p)
		},
		#' @description
		#' Transform venn result for the composition analysis.
		#'
		#' @param use_OTUs_frequency default TRUE; whether only use OTUs occurrence frequency, i.e. presence/absence data; if FALSE, use abundance data.
		#' @return a new \code{\link{microtable}} class.
		#' @examples
		#' \donttest{
		#' t2 <- t1$trans_venn_com(use_OTUs_frequency = TRUE)
		#' }
		trans_venn_com = function(use_OTUs_frequency = TRUE){
			otudata <- self$otu_table
			venn_table <- self$venn_table
			sampledata <- data.frame(SampleID = colnames(venn_table), Group = colnames(venn_table)) %>% 'rownames<-'(colnames(venn_table))
			taxdata <- self$tax_table
			sum_table <- data.frame(apply(otudata, 1, sum))
			tt <- dplyr::full_join(rownames_to_column(sum_table[venn_table[,1] %>% as.character %>% .[. != ""], ,drop=FALSE]),
				rownames_to_column(sum_table[venn_table[,2] %>% as.character %>% .[. != ""], , drop=FALSE]), by=c("rowname" = "rowname"))
			for(i in 3:ncol(venn_table)){
				tt <- dplyr::full_join(tt, rownames_to_column(sum_table[venn_table[, i] %>% 
					as.character %>% .[. != ""], , drop=FALSE]), by=c("rowname" = "rowname"))
			}
			tt[is.na(tt)] <- 0
			tt %<>% 'rownames<-'(.[,1]) %>% .[,-1,drop = FALSE]
			colnames(tt) <- colnames(venn_table)
			if(use_OTUs_frequency == T){
				tt[tt != 0] <- 1
			}
			microtable$new(sample_table = sampledata, otu_table = tt, tax_table = taxdata)
		},
		#' @description
		#' Print the trans_venn object.
		print = function() {
			print(self$venn_count_abund)
		}
		),
	private = list(
		vennSets = function(setmatrix, allcombl, index, setunion){
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		},
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
		# Circle function for 2 or 3-way Venn
		plotcircle = function(center = c(1, 1), diameter = 4, segments_split = 360) {
			r <- diameter / 2
			tt <- seq(0, 2*pi, length.out = segments_split)
			xx <- center[1] + r * cos(tt)
			yy <- center[2] + r * sin(tt)
			return(data.frame(x = xx, y = yy))
		},
		# Ellipse function for 4 or 5-way Venn
		plotellipse = function(center = c(1, 1), radius = c(2, 4), rotate = 1, segments_split = 360) {
			angles <- (0:segments_split) * 2 * pi/segments_split
			rotate <- rotate * pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind(ellipse[, 1] * cos(rotate) + ellipse[, 2] * sin(rotate), ellipse[, 2] * cos(rotate) - ellipse[, 1] * sin(rotate))
			ellipse <- cbind(center[1] + ellipse[, 1], center[2] + ellipse[, 2])
			colnames(ellipse) <- c("x", "y")
			return(as.data.frame(ellipse))
		},
		petal = function(r = 1, n= 1000, a = 4, b = 1.2, mx = 0, my = 0, rotate = 0){
			ang <- seq(0, 360, len=n+1)
			ang <- ang[1:n]
			ang <- ang*pi/180
			x <- r * cos(ang)
			y <- r* sin(ang)
			xy <- rbind(x,y)
			m <- diag(c(a,b))
			xy <- m %*% xy
			rotate <- rotate*pi/180
			m <- c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate))
			dim(m) <- c(2,2)
			xy <- m %*% xy
			xy[1, ] <- xy[1, ] + mx
			xy[2, ] <- xy[2, ] + my
			xy <- as.data.frame(t(xy))
			colnames(xy) <- c("x", "y")
			return(xy)
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


