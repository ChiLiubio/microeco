#' @title
#' Rarefaction and plotting.
#'
#' @description
#' Rarefaction based on the trans_norm class and visualization based on the ggplot2.
#'
#' @export
trans_rarefy <- R6Class(classname = "trans_rarefy",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param alphadiv default "Shannon"; alpha diversity measurement used for the rarefaction; see microtable$cal_alphadiv for all the measurement.
		#' @param depth default NULL; an integer vecotr used for the rarefying.
		#' @param ... parameters passed to \code{rarefy_samples} function of \code{microtable} class, except the sample.size parameter.
		#' @return res_rarefy stored in the object.
		#' @examples
		#' \donttest{
		#' library(microeco)
		#' data(dataset)
		#' t1 <- trans_rarefy$new(dataset = dataset, depth = c(0, 10, 50, 400, 800))
		#' }
		initialize = function(dataset = NULL, alphadiv = "Shannon", depth = NULL, ...)
			{
			res <- data.frame()
			if(length(alphadiv) > 1){
				message("The input alphadiv has multiple indexes! Only use the first one ...")
				alphadiv %<>% .[1]
			}
			for(i in depth){
				use_data <- clone(dataset)
				if(i == 0){
					suppressMessages(use_data$cal_alphadiv(measures = alphadiv))
					inter_res <- use_data$alpha_diversity[, alphadiv, drop = FALSE]
					inter_res[, 1] <- 0
					inter_res <- cbind.data.frame(SampleID = rownames(inter_res), seqnum = i, value = inter_res[, 1])
				}else{
					message("Rarefy data at depth ", i," ...")
					use_data$rarefy_samples(sample.size = i, ...)
					suppressMessages(use_data$cal_alphadiv(measures = alphadiv))
					inter_res <- use_data$alpha_diversity[, alphadiv, drop = FALSE]
					inter_res <- data.frame(SampleID = rownames(inter_res), seqnum = i, value = inter_res[, 1])
				}
				res <- rbind(res, inter_res)
			}
			colnames(res)[3] <- alphadiv
			# used for the following plotting
			self$dataset <- clone(dataset)
			self$measure <- alphadiv
			self$res_rarefy <- res
			message('The rarefied data is stored in object$res_rarefy ...')
		},
		#' @description
		#' Plotting the rarefied result.
		#'
		#' @param color_values colors used for presentation.
		#' @param color default "SampleID"; color mapping in the plot.
		#' @param show_point default TRUE; whether show the point.
		#' @param point_size default .3; point size value.
		#' @param point_alpha default .6; point alpha value.
		#' @param add_fitting default FALSE; whether add fitted line.
		#' @param x_axis_title default "Sequence number"; x axis title.
		#' @param y_axis_title default NULL; default NULL represents the measure used.
		#' @param show_legend default TRUE;	whether show the legend in the plot.
		#' @param show_samplename default FALSE; whether show the sample name in the plot.
		#' @param samplename_size default 3; the sample name text size. Only available when show_samplename is TRUE.
		#' @param samplename_color default "grey70"; sample name text color. Only available when show_samplename is TRUE.
		#' @param ... parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_rarefy(color = "Group")
		#' }
		plot_rarefy = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			color = "SampleID",
			show_point = TRUE,
			point_size = .3,
			point_alpha = .6,
			add_fitting = FALSE,
			x_axis_title = "Sequence number",
			y_axis_title = NULL,
			show_legend = TRUE,
			show_samplename = FALSE,
			samplename_size = 3,
			samplename_color = "grey30",
			...
			){
			alphadiv <- self$measure
			rarefy_data <- self$res_rarefy
			if(color != "SampleID"){
				sample_info <- self$dataset$sample_table
				if(any(colnames(sample_info) == "SampleID")){
					sample_info <- sample_info[, -which(colnames(sample_info) == "SampleID")]
				}
				sample_info <- data.frame(SampleID = rownames(sample_info), sample_info) %>% dropallfactors
				rarefy_data <- merge(rarefy_data, sample_info, by.x = 'SampleID', by.y = 'SampleID')
			}
			if(is.null(y_axis_title)){
				y_axis_title <- alphadiv
			}
			color_values <- microeco:::expand_colors(color_values, length(unique(rarefy_data[, color])))
			
			p <- ggplot(rarefy_data, aes(x = .data[["seqnum"]], y = .data[[alphadiv]], color = .data[[color]], group = .data[["SampleID"]])) +
				xlab(x_axis_title) +
				ylab(y_axis_title)
			if(show_point == T){
				p <- p + geom_point(alpha = point_alpha, size = point_size)
			}
			if(add_fitting == T){
				p <- p + geom_smooth(se = FALSE, ...)
			}else{
				p <- p + geom_line(...)
			}
			p <- p + scale_color_manual(values = color_values)

			p <- p + theme_bw()
			if(show_legend == F){
				p <- p + theme(legend.position = "none")
			}
			if(show_samplename){
				# select the max number for each sample
				num_order <- lapply(unique(rarefy_data$SampleID), function(x){
					tmp <- rarefy_data[rarefy_data$SampleID == x, "seqnum"] %>% max
					which(rarefy_data$SampleID == x & rarefy_data$seqnum == tmp)
					}) %>% unlist
				submax <- rarefy_data[num_order, ]
				
				p <- p + ggrepel::geom_text_repel(data = submax, aes(.data[["seqnum"]], .data[[alphadiv]], label = .data[["SampleID"]]), size = samplename_size, 
					color = samplename_color, parse = FALSE)
			}
			p
		},
		#' @description
		#' Print the trans_rarefy object.
		print = function() {
			cat("trans_rarefy class:\n")
			cat(paste("res_rarefy have been calculated at depths: ", paste0(unique(self$res_rarefy$seqnum), collapse = ", "), "\n"))
			invisible(self)
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
