#' @title Create \code{trans_niche} object for niche breadth and niche overlap analysis.
#'
#' @description
#' This class is a wrapper for niche breadth and niche overlap calculations,
#' including Levins' niche breadth, Pianka's niche overlap, and environment/trait-based
#' methods: OMI (Outlying Mean Index) based on Doledec et al. (2000) <doi:10.1890/0012-9658(2000)081[2914:NSICAA]2.0.CO;2>,
#' N-dimensional Hypervolume based on Blonder et al. (2014) <doi:10.1111/geb.12146> and Blonder (2018) <doi:10.1111/ecog.03187>,
#' and TPD (Trait Probability Density) based on Carmona et al. (2016) <doi:10.1016/j.tree.2016.02.003>.
#'
#' @export
trans_niche <- R6Class(classname = "trans_niche",
        public = list(
                #' @param dataset \code{\link{microtable}} object.
                #' @param filter_thres default 0; the relative abundance threshold for filtering low-abundance taxa.
                #'   Taxa with total relative abundance <= filter_thres will be removed.
                #'   This is delegated to \code{dataset$filter_taxa(rel_abund = filter_thres)}.
                #' @param env_cols default NULL; character vector of column names in \code{dataset$sample_table}
                #'   to be used as environmental variables. Only numeric columns are supported.
                #'   Passed to \code{\link{trans_env}} for processing (NA check and auto-fill).
                #' @param add_data default NULL; a data.frame of environmental variables with sample names as row names.
                #'   If provided, it takes priority over \code{env_cols}.
                #'   Passed to \code{\link{trans_env}} for processing (NA check and auto-fill).
                #' @return \code{data_rel_abund}, \code{n_samples}, \code{n_taxa}, \code{dataset}, \code{sample_table}
                #'   and optionally \code{data_env} stored in the object.
                #' @examples
                #' \donttest{
                #' data(dataset)
                #' t1 <- trans_niche$new(dataset = dataset)
                #' }
                initialize = function(dataset = NULL, filter_thres = 0, env_cols = NULL, add_data = NULL) {
                        microeco:::check_microtable(dataset)
                        use_dataset <- microeco::clone(dataset)
                        
                        if(filter_thres > 0){
                                use_dataset$filter_taxa(rel_abund = filter_thres)
                        }
                        # otu_table: rows = taxa, cols = samples
                        otu_mat <- as.matrix(use_dataset$otu_table)
                        # transpose to: rows = samples, cols = taxa
                        # this aligns with the resource state (sample) x species (taxa) convention
                        P <- t(otu_mat)
                        # compute relative abundance: each column (taxon) divided by its column sum
                        col_sums <- colSums(P)
                        col_sums[col_sums == 0] <- 1
                        rel_abund <- sweep(P, 2, col_sums, "/")
                        self$data_rel_abund <- rel_abund
                        self$n_samples <- nrow(rel_abund)
                        self$n_taxa <- ncol(rel_abund)
                        self$dataset <- use_dataset
                        self$sample_table <- use_dataset$sample_table
                        # handle environmental data via trans_env (NA check and auto-fill)
                        self$data_env <- NULL
                        if(!is.null(env_cols) || !is.null(add_data)){
                                env_obj <- trans_env$new(dataset = use_dataset, env_cols = env_cols, add_data = add_data)
                                self$data_env <- env_obj$data_env
                                message('The environmental data is stored in object$data_env ...')
                        }
                        message('The relative abundance data is stored in object$data_rel_abund ...')
                },
                #' @description
                #' Calculate niche breadth for each taxon.
                #'
                #' @param method default "levins"; the method for niche breadth calculation.
                #'   Options:
                #'   \describe{
                #'     \item{\strong{'levins'}}{Levins' niche breadth}
                #'     \item{\strong{'OMI'}}{Outlying Mean Index via ade4 package; Doledec et al. (2000) <doi:10.1890/0012-9658(2000)081[2914:NSICAA]2.0.CO;2>}
                #'     \item{\strong{'hypervolume'}}{N-dimensional hypervolume via hypervolume package; Blonder et al. (2014) <doi:10.1111/geb.12146>; Blonder (2018) <doi:10.1111/ecog.03187>}
                #'     \item{\strong{'TPD'}}{Trait Probability Density via TPD package; Carmona et al. (2016) <doi:10.1016/j.tree.2016.02.003>}
                #'   }
                #' @param ... additional parameters passed to the specific method:
                #'   \describe{
                #'     \item{OMI:}{\code{n_axis} default 2; number of niche axes to retain.}
                #'     \item{hypervolume:}{\code{dims} default 3; number of PCA dimensions for hypervolume construction.
                #'       \code{min_points} default 5; minimum number of occurrence samples required.
                #'       \code{kde.bandwidth} default "silverman"; bandwidth for kernel density estimation.
                #'       \code{quantile.requested} default 0.95; probability density quantile threshold.
                #'       \code{samples.per.point} default NULL; samples per point for hypervolume construction.}
                #'     \item{TPD:}{\code{dims} default 2; number of PCA dimensions for TPD calculation.
                #'       \code{alpha} default 0.95; probability density quantile threshold.}
                #'   }
                #' @return \code{res_niche_breadth} stored in the object. For OMI/hypervolume/TPD methods,
                #'   intermediate objects are also stored in \code{res_omi}, \code{res_hypervolume}, or \code{res_TPD}.
                cal_niche_breadth = function(method = c("levins", "OMI", "hypervolume", "TPD")[1], ...) {
                        method <- match.arg(method, c("levins", "OMI", "hypervolume", "TPD"))
                        if(method == "levins"){
                                P <- self$data_rel_abund
                                # Levins' B = 1 / sum(p_ij^2)
                                sum_p_sq <- colSums(P^2)
                                B <- ifelse(sum_p_sq == 0, NA, 1 / sum_p_sq)
                                # Standardized Levins: (B - 1) / (R - 1)
                                R <- self$n_samples
                                B_std <- if(R > 1) (B - 1) / (R - 1) else NA_real_
                                self$res_niche_breadth <- data.frame(
                                        Taxa = colnames(P),
                                        Levins = B,
                                        Levins_Standardized = B_std,
                                        stringsAsFactors = FALSE
                                )
                                message('The result is stored in object$res_niche_breadth ...')
                        }else if(method == "OMI"){
                                private$cal_niche_breadth_OMI(...)
                        }else if(method == "hypervolume"){
                                private$cal_niche_breadth_hypervolume(...)
                        }else if(method == "TPD"){
                                private$cal_niche_breadth_TPD(...)
                        }
                        invisible(self)
                },
                #' @description
                #' Calculate niche overlap matrix among taxa.
                #'
                #' @param method default "pianka"; the method for niche overlap calculation.
                #'   Options:
                #'   \describe{
                #'     \item{\strong{'pianka'}}{Pianka's niche overlap}
                #'     \item{\strong{'OMI'}}{based on OMI axes; Doledec et al. (2000) <doi:10.1890/0012-9658(2000)081[2914:NSICAA]2.0.CO;2>}
                #'     \item{\strong{'hypervolume'}}{based on hypervolume intersection; Blonder et al. (2014) <doi:10.1111/geb.12146>; Blonder (2018) <doi:10.1111/ecog.03187>}
                #'     \item{\strong{'TPD'}}{based on TPD probability density overlap; Carmona et al. (2016) <doi:10.1016/j.tree.2016.02.003>}
                #'   }
                #' @param ... additional parameters passed to the specific method:
                #'   \describe{
                #'     \item{OMI:}{\code{n_axis} default 2; number of OMI axes used for overlap calculation.}
                #'     \item{hypervolume:}{\code{overlap_metric} default "jaccard"; overlap metric,
                #'       either "jaccard" or "sorensen".}
                #'     \item{TPD:}{\code{symmetric} default TRUE; whether to compute symmetric overlap.
                #'       When TRUE, the geometric mean of the two directional overlaps is used (analogous
                #'       to Pianka's symmetric formulation). When FALSE, the raw directional overlaps
                #'       O(i,j) and O(j,i) are stored in the upper and lower triangles respectively,
                #'       producing an asymmetric matrix.}
                #'   }
                #' @return \code{res_niche_overlap} stored in the object, a matrix of dimension Taxa x Taxa.
                #'   For TPD method with \code{symmetric = FALSE}, the matrix is asymmetric where
                #'   overlap[i,j] represents the proportion of taxon i's niche overlapped by taxon j.
                cal_niche_overlap = function(method = c("pianka", "OMI", "hypervolume", "TPD")[1], ...) {
                        method <- match.arg(method, c("pianka", "OMI", "hypervolume", "TPD"))
                        if(method == "pianka"){
                                P <- self$data_rel_abund
                                # numerator: crossprod = P^T %*% P
                                # numerator[j, k] = sum_i(p_ij * p_ik)
                                numerator <- crossprod(P)
                                # denominator: sqrt(sum(p_ij^2) * sum(p_ik^2))
                                col_sq_sums <- colSums(P^2)
                                denominator <- sqrt(outer(col_sq_sums, col_sq_sums))
                                # Pianka's overlap
                                pianka_matrix <- numerator / denominator
                                # handle NaN from division by zero (taxa with all-zero abundance)
                                pianka_matrix[is.na(pianka_matrix)] <- 0
                                self$res_niche_overlap <- pianka_matrix
                                message('The result is stored in object$res_niche_overlap ...')
                        }else if(method == "OMI"){
                                private$cal_niche_overlap_OMI(...)
                        }else if(method == "hypervolume"){
                                private$cal_niche_overlap_hypervolume(...)
                        }else if(method == "TPD"){
                                private$cal_niche_overlap_TPD(...)
                        }
                        invisible(self)
                },
                #' @description
                #' Plot the niche breadth result as a horizontal bar plot.
                #'
                #' @param measure default NULL; a column name in \code{res_niche_breadth} used as the bar length or an integer to select the column name.
                #'   If NULL, the first numeric column except \code{Taxa} is used automatically.
                #' @param ntaxa default 20; number of top taxa to show, ordered by \code{measure} from high to low.
                #' @param add_prefix default TRUE; whether prepend the taxonomic name to the y-axis taxa labels.
                #' @param taxa_level default "Genus"; the taxonomic level used when \code{add_prefix} is \code{TRUE}
                #' @param sep default " : "; separator between genus name and taxa name.
                #' @param color default RColorBrewer::brewer.pal(8, "Dark2")[2]; bar fill color.
                #' @param barwidth default 0.7; bar width passed to \code{\link{geom_bar}}.
                #' @param ytext_size default 11; y-axis text size.
                #' @return ggplot2 plot.
                #' @examples
                #' \dontrun{
                #' t1$cal_niche_breadth(method = "levins")
                #' t1$plot_niche_breadth()
                #' }
                plot_niche_breadth = function(
                        measure = NULL,
                        ntaxa = 20,
                        add_prefix = TRUE,
						taxa_level = "Genus",
                        sep = " : ",
                        color = RColorBrewer::brewer.pal(8, "Dark2")[2],
                        barwidth = 0.7,
                        ytext_size = 11
                        ){
                        if(is.null(self$res_niche_breadth)){
                                stop("Please first use cal_niche_breadth() to calculate niche breadth !")
                        }
                        breadth_df <- self$res_niche_breadth
                        num_cols <- colnames(breadth_df)[unlist(lapply(breadth_df, is.numeric))]
                        num_cols <- setdiff(num_cols, "Taxa")
                        if(length(num_cols) == 0){
                                stop("No numeric columns found in res_niche_breadth!")
                        }
                        if(is.null(measure)){
							measure <- num_cols[1]
                        }else{
							if(! is.numeric(measure)){
								if(! measure %in% colnames(breadth_df)){
									stop("The measure '", measure, "' is not found in res_niche_breadth!")
								}
							}else{
								if(is.wholenumber(measure)){
									if(measure > length(num_cols)){
										message("Input measure is larger than the length of available names! Use the last one ...")
										measure <- length(num_cols)
									}
									measure <- num_cols[measure]
								}else{
									stop("Input measure should be integer when it is numeric class!")
								}
							}
                        }
                        # select top taxa
                        ordered_idx <- order(breadth_df[[measure]], decreasing = TRUE, na.last = TRUE)
                        use_n <- min(ntaxa, length(ordered_idx))
                        selected_idx <- ordered_idx[seq_len(use_n)]
                        plot_data <- breadth_df[selected_idx, c("Taxa", measure), drop = FALSE]
                        plot_data$Taxa <- as.character(plot_data$Taxa)
                        # create y-axis labels
                        if(add_prefix){
                                plot_data$Taxa_label <- private$get_taxa_label(plot_data$Taxa, taxa_level, sep = sep)
                        }else{
                                plot_data$Taxa_label <- plot_data$Taxa
                        }
                        # highest value at top
                        plot_data$Taxa_label <- factor(plot_data$Taxa_label, levels = rev(plot_data$Taxa_label))
                        p <- ggplot(plot_data, aes_meco(x = measure, y = "Taxa_label")) +
                                geom_bar(stat = "identity", fill = color, width = barwidth) +
                                xlab(measure) + ylab("") +
                                theme_bw() +
                                theme(panel.grid = element_blank(), panel.border = element_blank(),
                                        axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square"),
                                        axis.text.y = element_text(size = ytext_size))

                        p
                },
                #' @description
                #' Plot the niche overlap matrix as a heatmap.
                #'
                #' @param taxa_level default NULL; a taxonomic rank in \code{tax_table} used for filtering taxa.
                #' @param taxa_name default NULL; a character vector of taxa names at \code{taxa_level} used for filtering.
                #'   Must be provided together with \code{taxa_level}.
                #' @param color_low default "white"; color for low overlap values.
                #' @param color_high default RColorBrewer::brewer.pal(8, "Dark2")[2]; color for high overlap values.
                #' @param plot_breaks default NULL; breaks for the color scale.
                #' @param withmargin default TRUE; whether draw tile margins.
                #' @param margincolor default "white"; margin color when \code{withmargin = TRUE}.
                #' @param xtext_size default 9; x-axis text size.
                #' @param ytext_size default 9; y-axis text size.
                #' @param add_anno default FALSE; whether add taxonomic group annotation lines and labels to x and y axes.
                #'   When it is TRUE, the feature names will not be shown.
                #' @param add_anno_level default "Phylum"; the taxonomic rank used for axis annotation.
                #' @param anno_color default "grey60"; color of the annotation lines and labels.
                #' @param anno_text_size default 9; text size of the annotation labels.
                #' @param anno_line_size default 0.5; line width of the annotation lines.
                #' @return ggplot2 plot.
                #' @examples
                #' \dontrun{
                #' t1$cal_niche_overlap(method = "pianka")
                #' t1$plot_niche_overlap()
                #' }
                plot_niche_overlap = function(
                        taxa_level = NULL,
                        taxa_name = NULL,
                        color_low = "white",
                        color_high = RColorBrewer::brewer.pal(8, "Dark2")[2],
                        plot_breaks = NULL,
                        withmargin = TRUE,
                        margincolor = "white",
                        xtext_size = 9,
                        ytext_size = 9,
                        add_anno = FALSE,
                        add_anno_level = "Phylum",
                        anno_color = "grey60",
                        anno_text_size = 9,
                        anno_line_size = 0.5
                        ){
                        if(is.null(self$res_niche_overlap)){
                                stop("Please first use cal_niche_overlap() to calculate niche overlap !")
                        }
                        overlap_mat <- self$res_niche_overlap
                        if(!is.matrix(overlap_mat) && !is.data.frame(overlap_mat)){
                                stop("res_niche_overlap is not a matrix!")
                        }
                        all_taxa <- rownames(overlap_mat)
                        if(is.null(all_taxa) || length(all_taxa) == 0){
                                stop("No taxa names found in res_niche_overlap!")
                        }
                        # filtering by taxonomic level
                        if(!is.null(taxa_level) && is.null(taxa_name)){
                                stop("You provide taxa_level but no taxa_name! Please provide both or neither.")
                        }
                        if(is.null(taxa_level) && !is.null(taxa_name)){
                                stop("You provide taxa_name but no taxa_level! Please provide both or neither.")
                        }
                        if(!is.null(taxa_level)){
                                if(is.null(self$dataset$tax_table)){
                                        stop("No tax_table found in the dataset! Cannot filter by taxa_level.")
                                }
                                if(! taxa_level %in% colnames(self$dataset$tax_table)){
                                        stop("The taxa_level '", taxa_level, "' is not found in tax_table!")
                                }
                                tax_vals <- private$strip_prefix(self$dataset$tax_table[all_taxa, taxa_level, drop = TRUE])
                                selected_vals <- private$strip_prefix(taxa_name)
                                keep_taxa <- all_taxa[tax_vals %in% selected_vals]
                                if(length(keep_taxa) == 0){
                                        stop("No taxa retained after filtering by taxa_level and taxa_name!")
                                }
                                overlap_mat <- overlap_mat[keep_taxa, keep_taxa, drop = FALSE]
                                all_taxa <- keep_taxa
                        }
                        # order taxa by taxonomy
                        ordered_taxa <- private$order_taxa_by_taxonomy(all_taxa)
                        overlap_mat <- overlap_mat[ordered_taxa, ordered_taxa, drop = FALSE]
                        # convert to long format
                        plot_data <- reshape2::melt(overlap_mat)
                        colnames(plot_data) <- c("Taxa1", "Taxa2", "Overlap")
                        plot_data$Taxa1 <- factor(plot_data$Taxa1, levels = ordered_taxa)
                        plot_data$Taxa2 <- factor(plot_data$Taxa2, levels = ordered_taxa)
                        p <- ggplot(plot_data, aes_meco(x = "Taxa1", y = "Taxa2", fill = "Overlap"))
                        if(withmargin){
                                p <- p + geom_tile(colour = margincolor, linewidth = 0.5)
                        }else{
                                p <- p + geom_tile()
                        }
                        if(is.null(plot_breaks)){
                                p <- p + scale_fill_gradient(low = color_low, high = color_high, na.value = "grey90")
                        }else{
                                p <- p + scale_fill_gradient(low = color_low, high = color_high, na.value = "grey90", breaks = plot_breaks)
                        }
                        p <- p + xlab("") + ylab("") + theme_bw()
                        if(! add_anno){
                                p <- p + theme(
                                        panel.grid = element_blank(),
                                        panel.border = element_blank(),
                                        axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = xtext_size),
                                        axis.text.y = element_text(size = ytext_size)
                                )
                        }else{
                                p <- p + theme(
                                        panel.grid = element_blank(),
                                        panel.border = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank()
                                )
                        }
                        # add taxonomic group annotations to axes
                        if(add_anno){
                                if(is.null(self$dataset$tax_table)){
                                        stop("No tax_table found in the dataset! Cannot add axis annotation.")
                                }
                                if(! add_anno_level %in% colnames(self$dataset$tax_table)){
                                        stop("The add_anno_level '", add_anno_level, "' is not found in tax_table!")
                                }
                                anno_groups <- private$get_anno_groups(ordered_taxa, add_anno_level)
                                if(nrow(anno_groups) > 0){
                                        n_groups <- nrow(anno_groups)
                                        is_odd <- seq_len(n_groups) %% 2 == 1
                                        y_line <- ifelse(is_odd, -0.5, -1.5)
                                        df_anno_x <- data.frame(
                                                start = anno_groups$start,
                                                end = anno_groups$end,
                                                center = anno_groups$center,
                                                y_line = y_line,
                                                label = anno_groups$label,
                                                stringsAsFactors = FALSE
                                        )
                                        x_line <- ifelse(is_odd, -0.5, -1.5)
                                        df_anno_y <- data.frame(
                                                start = anno_groups$start,
                                                end = anno_groups$end,
                                                center = anno_groups$center,
                                                x_line = x_line,
                                                label = anno_groups$label,
                                                stringsAsFactors = FALSE
                                        )
                                        p <- p +
                                                geom_segment(
                                                        data = df_anno_x,
                                                        aes(x = start, xend = end, y = y_line, yend = y_line),
                                                        color = anno_color,
                                                        linewidth = anno_line_size,
                                                        inherit.aes = FALSE
                                                ) +
                                                geom_segment(
                                                        data = df_anno_y,
                                                        aes(x = x_line, xend = x_line, y = start, yend = end),
                                                        color = anno_color,
                                                        linewidth = anno_line_size,
                                                        inherit.aes = FALSE
                                                ) +
                                                geom_text(
                                                        data = df_anno_x,
                                                        aes(x = center, y = y_line, label = label),
                                                        size = anno_text_size / ggplot2::.pt,
                                                        color = anno_color,
                                                        vjust = 1.5,
                                                        inherit.aes = FALSE
                                                ) +
                                                geom_text(
                                                        data = df_anno_y,
                                                        aes(x = x_line, y = center, label = label),
                                                        size = anno_text_size / ggplot2::.pt,
                                                        color = anno_color,
                                                        hjust = 1,
                                                        inherit.aes = FALSE
                                                )
                                        p <- p +
                                                coord_cartesian(clip = "off") +
                                                theme(plot.margin = margin(b = 60, l = 60, t = 5, r = 5))
                                }
                        }
                        p
                },
                #' @description
                #' Print the trans_niche object.
                print = function() {
                        cat("trans_niche class object\n")
                        cat(paste0("Number of samples: ", self$n_samples, "\n"))
                        cat(paste0("Number of taxa: ", self$n_taxa, "\n"))
                        if(!is.null(self$data_env)){
                                cat(paste0("Environmental variables: ", ncol(self$data_env), " (", paste(colnames(self$data_env), collapse = ", "), ")\n"))
                        }else{
                                cat("Environmental data: Not provided\n")
                        }
                        if(!is.null(self$res_niche_breadth)){
                                cat("Niche breadth calculated: Yes\n")
                        }else{
                                cat("Niche breadth calculated: No\n")
                        }
                        if(!is.null(self$res_niche_overlap)){
                                cat("Niche overlap calculated: Yes\n")
                        }else{
                                cat("Niche overlap calculated: No\n")
                        }
                        if(!is.null(self$res_omi)){
                                cat("OMI analysis: Yes\n")
                        }
                        if(!is.null(self$res_hypervolume)){
                                cat("Hypervolume analysis: Yes\n")
                        }
                        if(!is.null(self$res_TPD)){
                                cat("TPD analysis: Yes\n")
                        }
                        invisible(self)
                }
        ),
        private = list(
                # ----- OMI method (Outlying Mean Index) -----
                cal_niche_breadth_OMI = function(n_axis = 2) {
                        if(!requireNamespace("ade4", quietly = TRUE)){
                                stop("Package 'ade4' is required for OMI method. Please install it: install.packages('ade4')")
                        }
                        if(is.null(self$data_env)){
                                stop("Environmental data is required for OMI method! Please provide env_cols or add_data when creating the object.")
                        }
                        env_data <- self$data_env
                        rel_abund <- self$data_rel_abund
                        # check that env_data and rel_abund have matching samples
                        if(!all(rownames(env_data) == rownames(rel_abund))){
                                common_samples <- intersect(rownames(env_data), rownames(rel_abund))
                                if(length(common_samples) == 0){
                                        stop("No common samples between environmental data and abundance data!")
                                }
                                env_data <- env_data[common_samples, , drop = FALSE]
                                rel_abund <- rel_abund[common_samples, , drop = FALSE]
                        }
                        # remove constant columns (zero variance)
                        env_sd <- sapply(env_data, sd, na.rm = TRUE)
                        if(any(env_sd == 0)){
                                message("Removing constant environmental variables: ", paste(names(env_sd)[env_sd == 0], collapse = ", "), " ...")
                                env_data <- env_data[, env_sd > 0, drop = FALSE]
                        }
                        if(ncol(env_data) == 0){
                                stop("No variable environmental columns left after removing constant columns!")
                        }
                        # step 1: PCA on environmental data
                        nf <- min(n_axis, ncol(env_data), nrow(env_data) - 1)
                        dudi_env <- ade4::dudi.pca(env_data, scannf = FALSE, nf = nf)
                        # step 2: prepare species abundance data.frame for niche()
                        # niche() requires Y as a data.frame: rows = samples, columns = species
                        # remove species with zero total abundance (niche() requires positive column sums)
                        col_sums <- colSums(rel_abund)
                        valid_taxa <- names(which(col_sums > 0))
                        Y <- as.data.frame(rel_abund[, valid_taxa, drop = FALSE])
                        # step 3: run niche analysis (all species at once)
                        niche_res <- ade4::niche(dudi_env, Y = Y, scannf = FALSE, nf = nf)
                        # step 4: extract niche parameters
                        params <- ade4::niche.param(niche_res)
                        # params is a matrix: rows = species, columns = inertia, OMI, Tol, Rtol, omi, tol, rtol
                        param_df <- as.data.frame(params)
                        # rename columns for clarity
                        colnames(param_df) <- c("OMI_Inertia", "OMI_Marginality", "OMI_Tol_Marginal",
                                "OMI_Residual_Tol", "OMI_Marginality_pct", "OMI_Tolerance_pct", "OMI_Residual_Tol_pct")
                        param_df$Taxa <- rownames(param_df)
                        # add back taxa with zero abundance (NA values)
                        all_taxa <- colnames(rel_abund)
                        zero_taxa <- setdiff(all_taxa, valid_taxa)
                        if(length(zero_taxa) > 0){
                                zero_df <- data.frame(
                                        OMI_Inertia = NA_real_, OMI_Marginality = NA_real_, OMI_Tol_Marginal = NA_real_,
                                        OMI_Residual_Tol = NA_real_, OMI_Marginality_pct = NA_real_,
                                        OMI_Tolerance_pct = NA_real_, OMI_Residual_Tol_pct = NA_real_,
                                        Taxa = zero_taxa, stringsAsFactors = FALSE
                                )
                                param_df <- rbind(param_df, zero_df)
                        }
                        # reorder to match original taxa order
                        param_df <- param_df[match(all_taxa, param_df$Taxa), ]
                        rownames(param_df) <- NULL
                        param_df <- param_df[, c("Taxa", setdiff(colnames(param_df), "Taxa"))]
                        self$res_niche_breadth <- param_df
                        self$res_omi <- niche_res
                        message('The OMI niche breadth result is stored in object$res_niche_breadth ...')
                        message('The OMI niche object is stored in object$res_omi ...')
                },
                cal_niche_overlap_OMI = function(n_axis = 2) {
                        if(!requireNamespace("ade4", quietly = TRUE)){
                                stop("Package 'ade4' is required for OMI method. Please install it: install.packages('ade4')")
                        }
                        if(is.null(self$data_env)){
                                stop("Environmental data is required for OMI method! Please provide env_cols or add_data when creating the object.")
                        }
                        env_data <- self$data_env
                        rel_abund <- self$data_rel_abund
                        # check sample matching
                        if(!all(rownames(env_data) == rownames(rel_abund))){
                                common_samples <- intersect(rownames(env_data), rownames(rel_abund))
                                env_data <- env_data[common_samples, , drop = FALSE]
                                rel_abund <- rel_abund[common_samples, , drop = FALSE]
                        }
                        # remove constant columns
                        env_sd <- sapply(env_data, sd, na.rm = TRUE)
                        if(any(env_sd == 0)){
                                env_data <- env_data[, env_sd > 0, drop = FALSE]
                        }
                        if(ncol(env_data) == 0){
                                stop("No variable environmental columns left!")
                        }
                        # PCA on environmental data
                        nf <- min(n_axis, ncol(env_data), nrow(env_data) - 1)
                        dudi_env <- ade4::dudi.pca(env_data, scannf = FALSE, nf = nf)
                        # prepare Y for niche analysis
                        col_sums <- colSums(rel_abund)
                        valid_taxa <- names(which(col_sums > 0))
                        Y <- as.data.frame(rel_abund[, valid_taxa, drop = FALSE])
                        niche_res <- ade4::niche(dudi_env, Y = Y, scannf = FALSE, nf = nf)
                        # use species scores on niche axes for overlap calculation
                        # niche_res$li contains species coordinates on niche axes
                        species_scores <- as.matrix(niche_res$li)
                        n_taxa <- nrow(species_scores)
                        taxa_names <- rownames(species_scores)
                        # compute overlap using Bhattacharyya coefficient on niche axes
                        # sample scores on niche axes
                        sample_scores <- as.matrix(niche_res$ls)
                        overlap_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa)
                        rownames(overlap_matrix) <- taxa_names
                        colnames(overlap_matrix) <- taxa_names
                        diag(overlap_matrix) <- 1
                        # compute weighted covariance for each species
                        cov_list <- vector("list", n_taxa)
                        for(i in seq_len(n_taxa)){
                                w_i <- Y[, taxa_names[i]]
                                w_i <- w_i / sum(w_i)
                                mu_i <- species_scores[i, ]
                                diff_i <- sweep(sample_scores, 2, mu_i)
                                cov_i <- matrix(0, nf, nf)
                                for(k in 1:nrow(sample_scores)){
                                        cov_i <- cov_i + w_i[k] * outer(diff_i[k, ], diff_i[k, ])
                                }
                                cov_list[[i]] <- cov_i
                        }
                        for(i in seq_len(n_taxa - 1)){
                                mu_i <- species_scores[i, ]
                                cov_i <- cov_list[[i]]
                                for(j in (i + 1):n_taxa){
                                        mu_j <- species_scores[j, ]
                                        cov_j <- cov_list[[j]]
                                        # Bhattacharyya distance
                                        cov_avg <- (cov_i + cov_j) / 2
                                        cov_avg_inv <- tryCatch(solve(cov_avg), error = function(e) NULL)
                                        if(is.null(cov_avg_inv)){
                                                overlap_matrix[i, j] <- 0
                                        }else{
                                                diff_mu <- mu_i - mu_j
                                                det_cov_avg <- det(cov_avg)
                                                det_cov_i <- det(cov_i)
                                                det_cov_j <- det(cov_j)
                                                D_B <- tryCatch({
                                                        as.numeric(1/8 * t(diff_mu) %*% cov_avg_inv %*% diff_mu +
                                                                1/2 * log(abs(det_cov_avg) / sqrt(abs(det_cov_i * det_cov_j))))
                                                }, error = function(e) Inf)
                                                if(is.na(D_B) || is.infinite(D_B)){
                                                        BC <- 0
                                                }else{
                                                        BC <- exp(-D_B)
                                                }
                                                BC <- max(0, min(1, BC))
                                                overlap_matrix[i, j] <- BC
                                                overlap_matrix[j, i] <- BC
                                        }
                                }
                        }
                        self$res_niche_overlap <- overlap_matrix
                        message('The OMI niche overlap result is stored in object$res_niche_overlap ...')
                },
                # ----- Hypervolume method -----
                cal_niche_breadth_hypervolume = function(dims = 3, min_points = 5,
                        kde.bandwidth = "silverman", quantile.requested = 0.95, samples.per.point = NULL) {
                        if(!requireNamespace("hypervolume", quietly = TRUE)){
                                stop("Package 'hypervolume' is required for hypervolume method. Please install it: install.packages('hypervolume')")
                        }
                        if(is.null(self$data_env)){
                                stop("Environmental data is required for hypervolume method! Please provide env_cols or add_data when creating the object.")
                        }
                        env_data <- self$data_env
                        otu_mat <- as.matrix(self$dataset$otu_table)
                        # check sample matching
                        common_samples <- intersect(rownames(env_data), colnames(otu_mat))
                        if(length(common_samples) == 0){
                                stop("No common samples between environmental data and abundance data!")
                        }
                        env_data <- env_data[common_samples, , drop = FALSE]
                        otu_mat <- otu_mat[, common_samples, drop = FALSE]
                        # remove constant columns before PCA (prcomp with scale. = TRUE fails on zero-variance columns)
                        env_sd <- sapply(env_data, sd, na.rm = TRUE)
                        if(any(env_sd == 0 | is.na(env_sd))){
                                const_cols <- names(env_sd)[env_sd == 0 | is.na(env_sd)]
                                message("Removing constant/NA-variance environmental variables for PCA: ", paste(const_cols, collapse = ", "), " ...")
                                env_data <- env_data[, env_sd > 0 & !is.na(env_sd), drop = FALSE]
                        }
                        if(ncol(env_data) == 0){
                                stop("No variable environmental columns left after removing constant columns!")
                        }
                        # PCA dimensionality reduction
                        env_pca <- prcomp(env_data, scale. = TRUE, rank. = dims)
                        actual_dims_pca <- min(dims, ncol(env_pca$x))
                        pca_scores <- as.data.frame(env_pca$x[, 1:actual_dims_pca, drop = FALSE])
                        if(actual_dims_pca < dims){
                                message("Only ", actual_dims_pca, " PCA dimensions available, using ", actual_dims_pca, " dimensions ...")
                        }
                        taxa_names <- rownames(otu_mat)
                        n_taxa <- length(taxa_names)
                        hv_list <- list()
                        volumes <- numeric(n_taxa)
                        skipped_taxa <- character(0)
                        error_taxa <- character(0)
                        for(i in seq_along(taxa_names)){
                                # get occurrence samples for this taxon
                                occ_samples <- names(which(otu_mat[i, ] > 0))
                                if(length(occ_samples) < min_points){
                                        skipped_taxa <- c(skipped_taxa, taxa_names[i])
                                        volumes[i] <- NA
                                        next
                                }
                                # extract environmental coordinates for occurrence samples
                                species_data <- pca_scores[occ_samples, , drop = FALSE]
                                # compute bandwidth externally to avoid hypervolume package
                                # internal Silverman computation dimensionality issues
                                if(is.character(kde.bandwidth) && kde.bandwidth == "silverman"){
                                        n_d <- ncol(species_data)
                                        n_r <- nrow(species_data)
                                        bw <- sapply(seq_len(n_d), function(k) sd(species_data[, k]))
                                        # replace zero/NA sd with a small positive value to avoid
                                        # degenerate kernels (species concentrated on a single PCA axis)
                                        bw[is.na(bw) | bw == 0] <- 1e-6
                                        bw <- (4 / (n_d + 2))^(1 / (n_d + 4)) *
                                              n_r^(-1 / (n_d + 4)) * bw
                                        use_bandwidth <- bw
                                }else if(is.character(kde.bandwidth) && kde.bandwidth == "plug-in"){
                                        if(!requireNamespace("ks", quietly = TRUE)){
                                                stop("Package 'ks' is required for plug-in bandwidth. Please install it: install.packages('ks')")
                                        }
                                        use_bandwidth <- ks::hpi(as.matrix(species_data))
                                }else{
                                        use_bandwidth <- kde.bandwidth
                                }
                                # build hypervolume (wrapped in tryCatch for robustness)
                                hv_args <- list(
                                        data = species_data,
                                        method = "gaussian",
                                        name = taxa_names[i],
                                        kde.bandwidth = use_bandwidth,
                                        quantile.requested = quantile.requested,
                                        verbose = FALSE
                                )
                                if(!is.null(samples.per.point)){
                                        hv_args$samples.per.point <- samples.per.point
                                }
                                hv <- tryCatch(
                                        do.call(hypervolume::hypervolume, hv_args),
                                        error = function(e) {
                                                message("Hypervolume construction failed for ", taxa_names[i], ": ", conditionMessage(e))
                                                NULL
                                        }
                                )
                                if(is.null(hv)){
                                        error_taxa <- c(error_taxa, taxa_names[i])
                                        volumes[i] <- NA
                                        next
                                }
                                hv_list[[taxa_names[i]]] <- hv
                                volumes[i] <- hypervolume::get_volume(hv)
                        }
                        if(length(skipped_taxa) > 0){
                                message("Skipped ", length(skipped_taxa), " taxa with fewer than ", min_points, " occurrence samples: ", 
									paste(head(skipped_taxa, 10), collapse = ", "), if(length(skipped_taxa) > 10) " ..." else "", " ...")
                        }
                        if(length(error_taxa) > 0){
                                message("Skipped ", length(error_taxa), " taxa due to hypervolume construction errors: ", paste(head(error_taxa, 10), collapse = ", "), 
									if(length(error_taxa) > 10) " ..." else "", " ...")
                        }
                        self$res_niche_breadth <- data.frame(
                                Taxa = taxa_names,
                                Hypervolume_Volume = volumes,
                                stringsAsFactors = FALSE
                        )
                        self$res_hypervolume <- hv_list
                        self$env_pca_scores <- pca_scores
                        message('The hypervolume niche breadth result is stored in object$res_niche_breadth ...')
                        message('The hypervolume objects are stored in object$res_hypervolume ...')
                },
                cal_niche_overlap_hypervolume = function(overlap_metric = c("jaccard", "sorensen")[1]) {
                        if(!requireNamespace("hypervolume", quietly = TRUE)){
                                stop("Package 'hypervolume' is required for hypervolume method. Please install it: install.packages('hypervolume')")
                        }
                        if(is.null(self$res_hypervolume)){
                                stop("Please first run cal_niche_breadth(method = 'hypervolume') to build hypervolumes!")
                        }
                        overlap_metric <- match.arg(overlap_metric, c("jaccard", "sorensen"))
                        hv_list <- self$res_hypervolume
                        valid_taxa <- names(hv_list)
                        n_valid <- length(valid_taxa)
                        if(n_valid == 0){
                                stop("No valid hypervolumes found! Check if cal_niche_breadth ran correctly.")
                        }
                        if(n_valid > 50){
                                message("Warning: Computing overlap for ", n_valid, " taxa (", n_valid * (n_valid - 1) / 2, " pairs) may take a long time ...")
                        }
                        # Build full all_taxa x all_taxa matrix (NA for skipped taxa)
                        all_taxa <- colnames(self$data_rel_abund)
                        n_all <- length(all_taxa)
                        overlap_matrix <- matrix(NA_real_, nrow = n_all, ncol = n_all)
                        rownames(overlap_matrix) <- all_taxa
                        colnames(overlap_matrix) <- all_taxa
                        diag(overlap_matrix) <- 1
                        # Compute pairwise overlap only for valid taxa
                        for(i in seq_len(n_valid - 1)){
                                for(j in (i + 1):n_valid){
                                        hv_set <- tryCatch({
                                                hypervolume::hypervolume_set(hv_list[[i]], hv_list[[j]], verbose = FALSE, check.memory = FALSE)
                                        }, error = function(e) NULL)
                                        ov <- 0
                                        if(!is.null(hv_set)){
                                                stats <- hypervolume::hypervolume_set_statistics(hv_set)
                                                ov <- if(overlap_metric == "jaccard") stats$jaccard else stats$sorensen
                                        }
                                        overlap_matrix[valid_taxa[i], valid_taxa[j]] <- ov
                                        overlap_matrix[valid_taxa[j], valid_taxa[i]] <- ov
                                }
                        }
                        # Set self-overlap for valid-only taxa (already 1 on diagonal)
                        self$res_niche_overlap <- overlap_matrix
                        message('The hypervolume niche overlap result is stored in object$res_niche_overlap ...')
                        if(n_valid < n_all){
                                message("Note: ", n_all - n_valid, " taxa were skipped (insufficient occurrence samples); their overlap values are NA.")
                        }
                },
                # ----- TPD method (Trait Probability Density) -----
                cal_niche_breadth_TPD = function(dims = 2, alpha = 0.95) {
                        if(!requireNamespace("TPD", quietly = TRUE)){
                                stop("Package 'TPD' is required for TPD method. Please install it: install.packages('TPD')")
                        }
                        if(is.null(self$data_env)){
                                stop("Environmental data is required for TPD method! Please provide env_cols or add_data when creating the object.")
                        }
                        env_data <- self$data_env
                        otu_mat <- as.matrix(self$dataset$otu_table)
                        # check sample matching
                        common_samples <- intersect(rownames(env_data), colnames(otu_mat))
                        if(length(common_samples) == 0){
                                stop("No common samples between environmental data and abundance data!")
                        }
                        env_data <- env_data[common_samples, , drop = FALSE]
                        otu_mat <- otu_mat[, common_samples, drop = FALSE]
                        # remove constant columns before PCA (prcomp with scale. = TRUE fails on zero-variance columns)
                        env_sd <- sapply(env_data, sd, na.rm = TRUE)
                        if(any(env_sd == 0 | is.na(env_sd))){
                                const_cols <- names(env_sd)[env_sd == 0 | is.na(env_sd)]
                                message("Removing constant/NA-variance environmental variables for PCA: ", paste(const_cols, collapse = ", "), " ...")
                                env_data <- env_data[, env_sd > 0 & !is.na(env_sd), drop = FALSE]
                        }
                        if(ncol(env_data) == 0){
                                stop("No variable environmental columns left after removing constant columns!")
                        }
                        # PCA dimensionality reduction
                        env_pca <- prcomp(env_data, scale. = TRUE, rank. = dims)
                        actual_dims_pca <- min(dims, ncol(env_pca$x))
                        pca_scores <- as.data.frame(env_pca$x[, 1:actual_dims_pca, drop = FALSE])
                        if(actual_dims_pca < dims){
                                message("Only ", actual_dims_pca, " PCA dimensions available, using ", actual_dims_pca, " dimensions ...")
                        }
                        # construct species-trait table
                        # for each taxon, extract environmental coordinates of its occurrence samples
                        taxa_names <- rownames(otu_mat)
                        species_vec <- character(0)
                        traits_list <- list()
                        for(i in seq_along(taxa_names)){
                                occ_samples <- names(which(otu_mat[i, ] > 0))
                                if(length(occ_samples) == 0) next
                                occ_env <- pca_scores[occ_samples, , drop = FALSE]
                                species_vec <- c(species_vec, rep(taxa_names[i], nrow(occ_env)))
                                traits_list[[length(traits_list) + 1]] <- occ_env
                        }
                        if(length(species_vec) == 0){
                                stop("No occurrence data found for any taxon!")
                        }
                        traits_df <- do.call(rbind, traits_list)
                        # compute TPDs
                        tpd_result <- TPD::TPDs(species = species_vec, traits = traits_df, alpha = alpha)
                        # compute REND (functional diversity indices)
                        # TPDs needs TPDc for REND; create a dummy TPDc with equal weights
                        unique_species <- unique(species_vec)
                        n_species <- length(unique_species)
                        # create community data: each species as its own "community" for individual metrics
                        sampUnit <- matrix(1, nrow = 1, ncol = n_species)
                        colnames(sampUnit) <- unique_species
                        rownames(sampUnit) <- "community"
                        tpd_comm <- TPD::TPDc(TPDs = tpd_result, sampUnit = sampUnit)
                        rend_result <- TPD::REND(TPDs = tpd_result, TPDc = tpd_comm)
                        # extract FRic, FEve, FDiv per species
                        breadth_df <- data.frame(
                                Taxa = taxa_names,
                                TPD_FRic = NA_real_,
                                TPD_FEve = NA_real_,
                                TPD_FDiv = NA_real_,
                                stringsAsFactors = FALSE
                        )
                        # rend_result$species contains per-species indices (named vectors)
                        if(!is.null(rend_result$species)){
                                fric <- rend_result$species$FRichness
                                feve <- rend_result$species$FEvenness
                                fdiv <- rend_result$species$FDivergence
                                for(i in seq_along(taxa_names)){
                                        idx <- which(names(fric) == taxa_names[i])
                                        if(length(idx) > 0){
                                                breadth_df$TPD_FRic[i] <- fric[idx[1]]
                                                breadth_df$TPD_FEve[i] <- feve[idx[1]]
                                                breadth_df$TPD_FDiv[i] <- fdiv[idx[1]]
                                        }
                                }
                        }
                        self$res_niche_breadth <- breadth_df
                        self$res_TPD <- tpd_result
                        self$env_pca_scores <- pca_scores
                        message('The TPD niche breadth result is stored in object$res_niche_breadth ...')
                        message('The TPD object is stored in object$res_TPD ...')
                },
                cal_niche_overlap_TPD = function(symmetric = TRUE) {
                        if(!requireNamespace("TPD", quietly = TRUE)){
                                stop("Package 'TPD' is required for TPD method. Please install it: install.packages('TPD')")
                        }
                        if(is.null(self$res_TPD)){
                                stop("Please first run cal_niche_breadth(method = 'TPD') to compute TPDs!")
                        }
                        tpd_result <- self$res_TPD
                        # get species names from TPDs (only taxa with enough occurrence data)
                        species_names <- names(tpd_result$TPDs)
                        n_species <- length(species_names)
                        if(n_species == 0){
                                stop("No species found in TPD result!")
                        }
                        # Build full all_taxa x all_taxa matrix (NA for skipped taxa)
                        all_taxa <- colnames(self$data_rel_abund)
                        n_all <- length(all_taxa)
                        overlap_matrix <- matrix(NA_real_, nrow = n_all, ncol = n_all)
                        rownames(overlap_matrix) <- all_taxa
                        colnames(overlap_matrix) <- all_taxa
                        diag(overlap_matrix) <- 1
                        # Pre-compute total density sums for each species.
                        total_sums <- sapply(species_names, function(sp) sum(tpd_result$TPDs[[sp]]))
                        for(i in seq_len(n_species - 1)){
                                density_i <- tpd_result$TPDs[[species_names[i]]]
                                sum_i <- total_sums[i]
                                for(j in (i + 1):n_species){
                                        density_j <- tpd_result$TPDs[[species_names[j]]]
                                        sum_j <- total_sums[j]
                                        # shared area under both density curves
                                        shared <- sum(pmin(density_i, density_j))
                                        if(symmetric){
                                                # Geometric mean of directional overlaps (Pianka-style)
                                                # O_ij = shared / sum_i ; O_ji = shared / sum_j
                                                # O_sym = sqrt(O_ij * O_ji) = shared / sqrt(sum_i * sum_j)
                                                if(sum_i > 0 && sum_j > 0){
                                                        overlap_val <- shared / sqrt(sum_i * sum_j)
                                                }else{
                                                        overlap_val <- 0
                                                }
                                                overlap_val <- max(0, min(1, overlap_val))
                                                overlap_matrix[species_names[i], species_names[j]] <- overlap_val
                                                overlap_matrix[species_names[j], species_names[i]] <- overlap_val
                                        }else{
                                                # Asymmetric: store directional overlaps separately
                                                # overlap[i,j] = proportion of taxon i overlapped by taxon j
                                                # overlap[j,i] = proportion of taxon j overlapped by taxon i
                                                o_ij <- if(sum_i > 0) shared / sum_i else 0
                                                o_ji <- if(sum_j > 0) shared / sum_j else 0
                                                overlap_matrix[species_names[i], species_names[j]] <- max(0, min(1, o_ij))
                                                overlap_matrix[species_names[j], species_names[i]] <- max(0, min(1, o_ji))
                                        }
                                }
                        }
                        self$res_niche_overlap <- overlap_matrix
                        if(symmetric){
                                message('The TPD niche overlap result (symmetric) is stored in object$res_niche_overlap ...')
                        }else{
                                message('The TPD niche overlap result (asymmetric: upper=O(i,j), lower=O(j,i)) is stored in object$res_niche_overlap ...')
                        }
                        if(n_species < n_all){
                                message("Note: ", n_all - n_species, " taxa were skipped (no occurrence data for TPD); their overlap values are NA.")
                        }
                },
                # ----- plotting helpers -----
                strip_prefix = function(x){
                        gsub("^.*__", "", as.character(x))
                },
                get_taxa_label = function(taxa_names, taxa_level, sep = " : "){
                        if(is.null(self$dataset$tax_table) || ! taxa_level %in% colnames(self$dataset$tax_table)){
                                return(taxa_names)
                        }
                        taxa_raw <- self$dataset$tax_table[taxa_names, taxa_level]
                        taxa <- private$strip_prefix(taxa_raw)
                        label <- ifelse(is.na(taxa) | taxa == "", taxa_names, paste0(taxa, sep, taxa_names))
                        label
                },
                order_taxa_by_taxonomy = function(taxa_names){
                        if(is.null(self$dataset$tax_table)){
                                return(taxa_names)
                        }
                        ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), colnames(self$dataset$tax_table))
                        if(length(ranks) == 0){
                                return(taxa_names)
                        }
                        tax_sub <- self$dataset$tax_table[taxa_names, ranks, drop = FALSE]
                        ord <- do.call(order, c(tax_sub, list(na.last = TRUE)))
                        rownames(tax_sub)[ord]
                },
                get_anno_groups = function(ordered_taxa, anno_level){
                        if(is.null(self$dataset$tax_table) || ! anno_level %in% colnames(self$dataset$tax_table)){
                                stop("No tax_table found or anno_level is not in tax_table!")
                        }
                        anno_raw <- self$dataset$tax_table[ordered_taxa, anno_level, drop = TRUE]
                        anno_vals <- private$strip_prefix(anno_raw)
                        anno_vals[is.na(anno_vals) | anno_vals == ""] <- ""
                        rle_res <- rle(anno_vals)
                        ends <- cumsum(rle_res$lengths)
                        starts <- ends - rle_res$lengths + 1
                        data.frame(
                                label = rle_res$values,
                                start = starts - 0.5,
                                end = ends + 0.5,
                                center = (starts + ends) / 2,
                                stringsAsFactors = FALSE
                        )[rle_res$values != "", , drop = FALSE]
                }
        ),
        lock_objects = FALSE,
        lock_class = FALSE
)
