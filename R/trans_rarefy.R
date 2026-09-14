#' @title
#' Rarefaction and plotting.
#'
#' @description
#' Rarefaction based on the trans_norm class and visualization of the rarefaction curve based on the ggplot2.
#'
#' @export
trans_rarefy <- R6Class(classname = "trans_rarefy",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param alphadiv default "Shannon"; one or more alpha diversity measurements used for the rarefaction; see microtable$cal_alphadiv for all the available measurements. 
		#'   The measurement names are case-insensitive, e.g. "shannon" is same as "Shannon". Multiple measurements are supported and share the same rarefying process.
		#' @param depth default NULL; a numeric vector used for the rarefying. 0 represents using the original data without rarefying and the diversity values are all set to 0 as the starting point of the curve. 
		#'   If NULL, a sequence of 10 depths from 0 to the maximum of sample sums (i.e. \code{max(dataset$sample_sums())}) is generated automatically.
		#'   Note that for each sample, only the depths not larger than its own sequencing depth are available, 
		#'   as the sample is removed at the larger depths by the \code{norm} function of \code{\link{trans_norm}} class, leading to an interrupted curve in the plot.
		#' @param PD default FALSE; whether add Faith's phylogenetic diversity (PD) to the alpha diversity. The calculation depends on \code{cal_alphadiv} function of \code{microtable} class, 
		#'   so the phylogenetic tree (\code{phylo_tree} in the dataset) is necessary for this option.
		#' @param ... parameters passed to \code{norm} function of \code{\link{trans_norm}} class, such as \code{method} ("rarefy" or "SRS") and \code{rngseed}, 
		#'   except the sample.size parameter which is controlled by the depth parameter.
		#' @return \code{res_rarefy} stored in the object, which is a data.frame with the columns \code{SampleID} (sample name), 
		#'   \code{seqnum} (the rarefying depth) and the alpha diversity measurements used (the column names are the standard measurement names). 
		#'   The measurements used are also stored in \code{object$measure}, and the rarefying depths in \code{object$depth}.
		#' @examples
		#' \dontrun{
		#' data(dataset)
		#' t1 <- trans_rarefy$new(dataset = dataset, depth = c(0, 10, 50, 400, 800))
		#' # the depth is generated automatically when it is not provided
		#' t1 <- trans_rarefy$new(dataset = dataset)
		#' }
		initialize = function(dataset = NULL, alphadiv = "Shannon", depth = NULL, PD = FALSE, ...)
			{
			if(is.null(dataset)){
				stop("Please provide a microtable object in the parameter dataset!")
			}
			if(! inherits(dataset, "microtable")){
				stop("The input dataset must be a microtable object! Please check it!")
			}
			max_depth <- max(dataset$sample_sums())
			if(is.null(depth)){
				if(! is.finite(max_depth) || max_depth <= 0){
					stop("No available sequencing depth is found in the dataset! Please check the otu_table in the dataset ...")
				}
				depth <- unique(round(seq(0, max_depth, length.out = 10)))
				message("The parameter depth is not provided! Use the automatically generated depths: ", paste0(depth, collapse = ", "), " ...")
			}
			depth_numeric <- suppressWarnings(as.numeric(depth))
			if(length(depth_numeric) == 0){
				stop("The parameter depth is empty! Please provide at least one depth value ...")
			}
			if(any(is.na(depth_numeric))){
				stop("All the numbers in depth should be convertible to numeric values!")
			}
			if(any(depth_numeric < 0)){
				stop("All the numbers in depth should be non-negative!")
			}
			if(! all(is.wholenumber(depth_numeric))){
				# the rarefying depth must be integers, as the non-integer value is silently truncated in the subsampling
				stop("All the numbers in depth should be integers! As the depth denotes the number of reads for the rarefying ...")
			}
			if(any(depth_numeric > max_depth)){
				stop("The following depths are larger than the maximum of sample sums (", max_depth, "): ", 
					paste0(depth_numeric[depth_numeric > max_depth], collapse = ", "), "! Please check the parameter depth ...")
			}
			if(anyDuplicated(depth_numeric)){
				message("The duplicated values in the parameter depth are removed ...")
			}
			depth <- sort(unique(depth_numeric))
			# resolve the measure names to be consistent with the output colnames of cal_alphadiv
			# the keys of measure_map should be kept consistent with the renamevec in the cal_alphadiv function of microtable class
			measure_map <- c("s.obs" = "Observed", "observed" = "Observed", "coverage" = "Coverage", "s.chao1" = "Chao1", "chao1" = "Chao1", 
				"s.ace" = "ACE", "ace" = "ACE", "shannon" = "Shannon", "simpson" = "Simpson", "invsimpson" = "InvSimpson", 
				"fisher" = "Fisher", "pielou" = "Pielou", "pd" = "PD")
			measure_key <- tolower(as.character(alphadiv))
			unmatched <- ! measure_key %in% names(measure_map)
			if(any(unmatched)){
				stop("The following measures are not supported: ", paste0(alphadiv[unmatched], collapse = ", "), 
					". The available measures: ", paste0(sort(unique(unname(measure_map))), collapse = ", "))
			}
			measures <- unique(unname(measure_map[measure_key]))
			if("PD" %in% measures && ! isTRUE(PD)){
				stop("Please set PD = TRUE when the PD measure is required!")
			}
			# parameters passed to the norm function of trans_norm class
			dots <- list(...)
			if("sample.size" %in% names(dots)){
				message("The sample.size parameter in ... is ignored! Please use the depth parameter instead ...")
				dots$sample.size <- NULL
			}
			method <- "rarefy"
			if("method" %in% names(dots)){
				method <- as.character(dots$method)[1]
				dots$method <- NULL
			}
			if(identical(method, "rarefying")){
				# keep the alias supported by the deprecated rarefy_samples function of microtable class
				method <- "rarefy"
			}
			# clone once and construct the trans_norm object once; 
			# the norm function reconstructs a new microtable object from the original data at each call, 
			# so there is no need to reset the data for each depth
			use_data <- clone(dataset)
			tmp_norm <- suppressMessages(trans_norm$new(use_data))
			# the samples in the abundance table, used for the rows of the depth 0
			sample_names_all <- use_data$sample_names()
			# measures passed to cal_alphadiv; PD is controlled by the PD parameter, so it is excluded
			cal_measures <- setdiff(measures, "PD")
			if(length(cal_measures) == 0){
				cal_measures <- "Observed"
			}
			res_list <- vector("list", length(depth))
			measures_checked <- FALSE
			for(i in seq_along(depth)){
				if(depth[i] == 0){
					# 0 denotes no rarefying; the diversity values are all set to 0
					zero_table <- matrix(0, nrow = length(sample_names_all), ncol = length(measures), dimnames = list(NULL, measures))
					res_list[[i]] <- data.frame(SampleID = sample_names_all, seqnum = 0, zero_table, check.names = FALSE, stringsAsFactors = FALSE)
				}else{
					message("Rarefy data at depth ", depth[i], " ...")
					# rarefy the data with the norm function of trans_norm class instead of the deprecated rarefy_samples function
					new_data <- do.call(tmp_norm$norm, c(list(method = method, sample.size = depth[i]), dots))
					suppressMessages(new_data$cal_alphadiv(measures = cal_measures, PD = PD))
					if(! measures_checked){
						# check the measures based on the first calculated result
						missed <- setdiff(measures, colnames(new_data$alpha_diversity))
						if(length(missed) > 0){
							stop("The following measures are not found in the calculated alpha diversity: ", paste0(missed, collapse = ", "), 
								". The available measures: ", paste0(colnames(new_data$alpha_diversity), collapse = ", "), 
								". Some measures may be unavailable at a very low rarefying depth, e.g. Fisher ...")
						}
						measures_checked <- TRUE
					}
					inter_res <- new_data$alpha_diversity[, measures, drop = FALSE]
					res_list[[i]] <- data.frame(SampleID = rownames(inter_res), seqnum = depth[i], inter_res, check.names = FALSE, stringsAsFactors = FALSE)
				}
			}
			res <- do.call(rbind, res_list)
			# check the missing values, which usually occur at a very low rarefying depth, e.g. Pielou is 0/0 when only one species is present
			na_num <- colSums(is.na(res[, measures, drop = FALSE]))
			if(any(na_num > 0)){
				message("NaN or NA is found in the rarefied result for the following measures: ", 
					paste0(names(na_num)[na_num > 0], " (", na_num[na_num > 0], " rows)", collapse = ", "), 
					". This usually happens at a very low rarefying depth. Please check the result or exclude the corresponding depths ...")
			}
			# used for the following plotting
			self$dataset <- use_data
			self$measure <- measures
			self$depth <- depth
			self$res_rarefy <- res
			message('The rarefied data is stored in object$res_rarefy ...')
		},
		#' @description
		#' Plotting the rarefied result.
		#'
		#' @param color_values colors used for presentation.
		#' @param color default "SampleID"; color mapping in the plot.
		#' @param measure default NULL; the alpha diversity measurement(s) used as the y axis; see the alphadiv parameter in \code{new} function. 
		#'   If NULL, all the measures stored in the object are used. When multiple measures are used, the plot is faceted by the measures.
		#' @param show_point default TRUE; whether show the point.
		#' @param point_size default .3; point size value.
		#' @param point_alpha default .6; point alpha value.
		#' @param add_fitting default FALSE; whether add fitted line.
		#' @param fitting_method default "lm"; the method used by \code{ggplot2::geom_smooth} function. Only available when \code{add_fitting = TRUE}. 
		#'   The default "lm" avoids the warnings generated by the default "loess" method, as there are usually only a few depths for each sample.
		#' @param fitting_formula default \code{y ~ log(x + 1)}; the formula used by \code{ggplot2::geom_smooth} function. Only available when \code{add_fitting = TRUE}. 
		#'   Note that \code{x + 1} is used instead of \code{x}, as the depth 0 makes \code{log(x)} invalid.
		#' @param x_axis_title default "Sequence number"; x axis title.
		#' @param y_axis_title default NULL; default NULL represents the measure used; when multiple measures are used, default NULL represents "Value".
		#' @param show_legend default TRUE; whether show the legend in the plot.
		#' @param show_samplename default FALSE; whether show the sample name in the plot.
		#' @param samplename_size default 3; the sample name text size. Only available when show_samplename is TRUE.
		#' @param samplename_color default "grey30"; sample name text color. Only available when show_samplename is TRUE.
		#' @param ... parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).
		#' @return ggplot.
		#' @examples
		#' \dontrun{
		#' t1$plot_rarefy(color = "Group")
		#' }
		plot_rarefy = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			color = "SampleID",
			measure = NULL,
			show_point = TRUE,
			point_size = .3,
			point_alpha = .6,
			add_fitting = FALSE,
			fitting_method = "lm",
			fitting_formula = y ~ log(x + 1),
			x_axis_title = "Sequence number",
			y_axis_title = NULL,
			show_legend = TRUE,
			show_samplename = FALSE,
			samplename_size = 3,
			samplename_color = "grey30",
			...
			){
			measures <- self$measure
			rarefy_data <- self$res_rarefy
			if(is.null(measure)){
				use_measures <- measures
			}else{
				unmatched <- ! measure %in% measures
				if(any(unmatched)){
					stop("The following measures are not found in the object: ", paste0(measure[unmatched], collapse = ", "), 
						". The available measures: ", paste0(measures, collapse = ", "))
				}
				use_measures <- measure
			}
			if(color != "SampleID"){
				sample_info <- self$dataset$sample_table
				if(any(colnames(sample_info) == "SampleID")){
					sample_info <- sample_info[, -which(colnames(sample_info) == "SampleID"), drop = FALSE]
				}
				sample_info <- data.frame(SampleID = rownames(sample_info), sample_info, stringsAsFactors = FALSE) %>% dropallfactors
				rarefy_data <- merge(rarefy_data, sample_info, by = "SampleID")
				if(! color %in% colnames(rarefy_data)){
					stop("The parameter color = '", color, "' is not found in the sample_table of the object! Available columns: ", 
						paste0(colnames(self$dataset$sample_table), collapse = ", "), " ...")
				}
			}
			# sort the data to avoid the wrong line connection when the input depth is unordered
			rarefy_data %<>% .[order(.[["SampleID"]], .[["seqnum"]]), ]
			if(length(use_measures) > 1){
				facet_plot <- TRUE
				y_var <- "value"
				rarefy_data %<>% reshape2::melt(id.vars = setdiff(colnames(rarefy_data), use_measures), variable.name = "measure", value.name = "value")
				# sort the data again as the melt operation may change the row order
				rarefy_data %<>% .[order(.[["measure"]], .[["SampleID"]], .[["seqnum"]]), ]
			}else{
				facet_plot <- FALSE
				y_var <- use_measures
			}
			if(is.null(y_axis_title)){
				y_axis_title <- if(facet_plot) "Value" else y_var
			}
			color_values <- expand_colors(color_values, length(unique(rarefy_data[[color]])))

			p <- ggplot(rarefy_data, aes(x = .data[["seqnum"]], y = .data[[y_var]], color = .data[[color]], group = .data[["SampleID"]])) +
				xlab(x_axis_title) +
				ylab(y_axis_title)
			if(isTRUE(show_point)){
				p <- p + geom_point(alpha = point_alpha, size = point_size)
			}
			if(isTRUE(add_fitting)){
				# the method and formula can be overwritten by the parameters in ... when necessary
				smooth_args <- list(...)
				if(is.null(smooth_args$method)){
					smooth_args$method <- fitting_method
				}
				if(is.null(smooth_args$formula)){
					smooth_args$formula <- fitting_formula
				}
				p <- p + do.call(geom_smooth, c(list(se = FALSE), smooth_args))
			}else{
				p <- p + geom_line(...)
			}
			p <- p + scale_color_manual(values = color_values)

			p <- p + theme_bw()
			if(facet_plot){
				p <- p + facet_wrap(~ measure, scales = "free_y")
			}
			if(! isTRUE(show_legend)){
				p <- p + theme(legend.position = "none")
			}
			if(isTRUE(show_samplename)){
				# select the max number for each sample
				num_order <- lapply(unique(rarefy_data$SampleID), function(x){
					tmp <- rarefy_data[rarefy_data$SampleID == x, "seqnum"] %>% max
					which(rarefy_data$SampleID == x & rarefy_data$seqnum == tmp)
					}) %>% unlist
				submax <- rarefy_data[num_order, ]

				p <- p + ggrepel::geom_text_repel(data = submax, aes(.data[["seqnum"]], .data[[y_var]], label = .data[["SampleID"]]), size = samplename_size,
					color = samplename_color, parse = FALSE)
			}
			p
		},
		#' @description
		#' Print the trans_rarefy object.
		print = function() {
			cat("trans_rarefy class:\n")
			cat("The measures used:", paste0(self$measure, collapse = ", "), "\n")
			if(is.null(self$depth)){
				cat("The rarefying depths are not available ...\n")
			}else{
				cat("res_rarefy have been calculated at depths:", paste0(self$depth, collapse = ", "), "\n")
			}
			invisible(self)
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
