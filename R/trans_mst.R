#' @title Create \code{trans_mst} object for microbial source tracking (MST) analysis.
#'
#' @description
#' This class is a wrapper for microbial source tracking based on the FEAST algorithm
#' (Fast Expectation-maximization for microbial Source Tracking; Shenhav et al. 2019, Nature Methods, <doi:10.1038/s41592-019-0431-x>).
#' It estimates the contribution proportions of potential source environments to sink samples
#' (e.g., tracking the microbial transfer from bulk soil to rhizosphere soil and then to root).
#' The class automatically constructs the metadata (Env, SourceSink, id) required by FEAST
#' from the \code{sample_table}, runs the \code{FEAST} package, and provides ggplot2-based
#' visualization. Results of multiple runs (e.g., sequential transfer steps) can be stored
#' in the same object using different labels.
#'
#' @export
trans_mst <- R6Class(classname = "trans_mst",
        public = list(
                #' @param dataset \code{\link{microtable}} object.
                #' @param env_col default NULL; a column name in \code{sample_table} used as the environment
                #'   description (the "Env" information for FEAST), e.g., "Compartment".
                #'   If NULL, the environment of each sample is represented by its sample name when
                #'   \code{cal_mst} is called with explicit sample vectors (\code{source_samples}/\code{sink_samples}).
                #' @param filter_thres default 0; the relative abundance threshold for filtering low-abundance taxa.
				#'   Taxa with total relative abundance < filter_thres will be removed.
				#'   This is delegated to \code{dataset$filter_taxa(rel_abund = filter_thres)}.
                #' @return \code{dataset}, \code{sample_table}, \code{otu_table} (taxa x samples counts),
                #'   \code{env_col}, \code{n_samples}, \code{n_taxa} stored in the object.
                #' @examples
                #' \dontrun{
                #' data(soil_microb)
                #' t1 <- trans_mst$new(dataset = soil_microb, env_col = "Compartment")
                #' }
                initialize = function(dataset = NULL, env_col = NULL, filter_thres = 0) {
                        microeco:::check_microtable(dataset)
                        use_dataset <- microeco::clone(dataset)
                        if(filter_thres > 0){
                                use_dataset$filter_taxa(rel_abund = filter_thres)
                        }
                        self$dataset <- use_dataset
                        self$sample_table <- use_dataset$sample_table
                        # otu_table: rows = taxa, cols = samples
                        self$otu_table <- as.matrix(use_dataset$otu_table)
                        self$n_samples <- ncol(self$otu_table)
                        self$n_taxa <- nrow(self$otu_table)
                        if(!is.null(env_col)){
                                if(! env_col %in% colnames(self$sample_table)){
                                        stop("The env_col '", env_col, "' is not found in sample_table!")
                                }
                                self$env_col <- env_col
                        }
                        message("The object is created. Use cal_mst() to run the source tracking analysis ...")
                },
                #' @description
                #' Run FEAST source tracking analysis.
                #'
                #' Sources and sinks can be specified either by environment values (\code{sources}/\code{sinks},
                #' matched against the \code{env_col} column of \code{sample_table}) or by explicit sample names
                #' (\code{source_samples}/\code{sink_samples}).
                #' The function automatically generates the FEAST metadata:
                #' in the same-sources mode (\code{different_sources_flag = 0}), all sinks share the same source pool;
                #' in the different-sources mode (\code{different_sources_flag = 1}), each sink uses its own source
                #' set defined by the pairing column \code{id_col} (samples sharing the same \code{id_col} value
                #' belong to the same group, e.g., the same plant individual).
                #'
                #' @param sources default NULL; character vector of environment values (in the \code{env_col} column)
                #'   used as sources, e.g., c("Bulk soil", "Rhizosphere").
                #' @param sinks default NULL; character vector of environment values used as sinks, e.g., "Endophyte".
                #' @param source_samples default NULL; character vector of sample names used as sources.
                #'   Alternative to \code{sources}; cannot be provided together with \code{sources}.
                #' @param sink_samples default NULL; character vector of sample names used as sinks.
                #'   Alternative to \code{sinks}; cannot be provided together with \code{sinks}.
                #' @param different_sources_flag default 0; 0 (FALSE) means all sinks use the same source pool;
                #'   1 (TRUE) means each sink uses its own source set paired via \code{id_col}.
                #' @param id_col default NULL; a column name in \code{sample_table} defining the pairing groups
                #'   for the different-sources mode (required when \code{different_sources_flag = 1}).
                #'   Samples sharing the same value are considered paired, e.g., "Plant_ID".
                #'   Each pairing group must contain exactly one sink sample; the sources in the group
                #'   constitute the source set of that sink.
                #' @param EM_iterations default 1000; the number of EM iterations passed to \code{\link[FEAST:FEAST]{FEAST::FEAST}}.
                #' @param COVERAGE default NULL; the rarefaction depth. If NULL, the minimum sequencing depth
                #'   across all involved source and sink samples is used, so all samples are uniformly downsampled.
                #' @param label default NULL; the label for storing the result, e.g., "bulk2rhizo".
                #'   If NULL, a label like "run1", "run2", ... is generated automatically.
                #'   Reusing an existing label overwrites the previous result with the same label.
                #' @return \code{res_mst[[label]]}: Sink x Source contribution proportions matrix
                #'   (rownames: sink sample names; colnames: source sample names + "Unknown"; each row sums to 1);
                #'   \code{res_mst_env[[label]]}: contributions aggregated by source environment;
                #'   \code{res_mst_param[[label]]}: parameters and sample information of the run.
                #' @examples
                #' \dontrun{
                #' data(soil_microb)
                #' t1 <- trans_mst$new(dataset = soil_microb, env_col = "Compartment")
                #' # same sources for all sinks
                #' t1$cal_mst(sources = "Bulk soil", sinks = "Rhizosphere", EM_iterations = 1000, label = "bulk2rhizo")
                #' # different sources for each sink, paired by Plant_ID
                #' t1$cal_mst(sources = c("Bulk soil", "Rhizosphere"), sinks = "Endophyte",
                #'   different_sources_flag = 1, id_col = "Plant_ID", label = "soil2root")
                #' }
                cal_mst = function(
                        sources = NULL,
                        sinks = NULL,
                        source_samples = NULL,
                        sink_samples = NULL,
                        different_sources_flag = 0,
                        id_col = NULL,
                        EM_iterations = 1000,
                        COVERAGE = NULL,
                        label = NULL
                        ){
                        # ----- check FEAST availability -----
                        if(!requireNamespace("FEAST", quietly = TRUE)){
                                stop("Package 'FEAST' is required for source tracking! ",
                                        "It is not on CRAN; please install it from GitHub: remotes::install_github('cozygene/FEAST').")
                        }
                        # ----- resolve samples -----
                        res <- private$resolve_samples(sources = sources, sinks = sinks,
                                source_samples = source_samples, sink_samples = sink_samples)
                        source_ids <- res$source_ids
                        sink_ids <- res$sink_ids
                        env_map <- res$env_map  # data.frame: sample_id, env, role
                        if(length(intersect(source_ids, sink_ids)) > 0){
                                stop("Some samples are assigned to both sources and sinks! Please check the inputs.")
                        }
                        different_sources_flag <- as.integer(different_sources_flag)
                        if(! different_sources_flag %in% c(0, 1)){
                                stop("different_sources_flag must be 0 (same sources) or 1 (different sources)!")
                        }
                        # ----- build FEAST metadata (Env, SourceSink, id) -----
                        # each sink gets a unique id; sources get NA (same-sources mode)
                        # or the id of their paired sink (different-sources mode)
                        metadata <- data.frame(
                                Env = env_map$env,
                                SourceSink = env_map$role,
                                id = NA_real_,
                                stringsAsFactors = FALSE
                        )
                        rownames(metadata) <- env_map$sample_id
                        sink_idx <- which(metadata$SourceSink == "Sink")
                        metadata$id[sink_idx] <- seq_along(sink_idx)
                        if(different_sources_flag == 1){
                                if(is.null(id_col)){
                                        stop("id_col is required when different_sources_flag = 1! Please provide a sample_table column defining the pairing groups.")
                                }
                                if(! id_col %in% colnames(self$sample_table)){
                                        stop("The id_col '", id_col, "' is not found in sample_table!")
                                }
                                groups <- as.character(self$sample_table[env_map$sample_id, id_col])
                                names(groups) <- env_map$sample_id
                                if(any(is.na(groups))){
                                        stop("NA values found in the id_col column '", id_col, "' for the involved samples! Please check the sample_table.")
                                }
                                # check: each group must contain at most one sink
                                sink_groups <- groups[env_map$role == "Sink"]
                                dup_groups <- names(table(sink_groups))[table(sink_groups) > 1]
                                if(length(dup_groups) > 0){
                                        stop("Each pairing group (id_col) must contain only one sink sample, but the following groups have multiple sinks: ",
                                                paste(dup_groups, collapse = ", "), " !")
                                }
                                # assign the sink's id to its paired sources
                                paired <- character(0)
                                for(sid in sink_ids){
                                        g <- groups[[sid]]
                                        pool <- source_ids[groups[source_ids] == g]
                                        if(length(pool) == 0){
                                                stop("The sink sample '", sid, "' (group '", g, "') has no source samples with the same ", id_col, " value!")
                                        }
                                        metadata[pool, "id"] <- metadata[sid, "id"]
                                        paired <- union(paired, pool)
                                }
                                # drop sources not paired with any sink
                                dropped_sources <- setdiff(source_ids, paired)
                                if(length(dropped_sources) > 0){
                                        message(length(dropped_sources), " source samples without any paired sink are excluded: ",
                                                paste(dropped_sources, collapse = ", "), " ...")
                                        source_ids <- setdiff(source_ids, dropped_sources)
                                        metadata <- metadata[! rownames(metadata) %in% dropped_sources, , drop = FALSE]
                                }
                        }
                        # ----- prepare count matrix: samples x taxa, integer -----
                        use_samples <- c(source_ids, sink_ids)
                        C <- t(self$otu_table[, use_samples, drop = FALSE])
                        if(! all(is.wholenumber(C))){
                                message("Non-integer values detected in otu_table! Values are converted by ceiling() ...")
                                C <- ceiling(C)
                        }
                        storage.mode(C) <- "integer"
                        # ----- COVERAGE: uniform rarefaction depth -----
                        zero_depth <- rownames(C)[rowSums(C) <= 0]
                        if(length(zero_depth) > 0){
                                stop("The following samples have a sequencing depth of 0, please remove them first: ",
                                        paste(zero_depth, collapse = ", "), " !")
                        }
                        if(is.null(COVERAGE)){
                                COVERAGE <- min(rowSums(C))
                                message("COVERAGE is automatically set to the minimum sequencing depth of the involved samples: ", COVERAGE, " ...")
                        }
                        if(COVERAGE <= 0){
                                stop("COVERAGE must be a positive number!")
                        }
                        # ----- run FEAST -----
                        message("Running FEAST source tracking, please wait ...")
                        props_raw <- FEAST::FEAST(
                                C = C,
                                metadata = metadata,
                                EM_iterations = EM_iterations,
                                COVERAGE = COVERAGE,
                                different_sources_flag = different_sources_flag
                        )
                        # ----- post-process: recover sample names -----
                        # FEAST output rownames: paste0(sink_id, "_", Env); colnames: paste0(source_id, "_", Env) + "Unknown"
                        sink_env <- setNames(env_map$env[env_map$role == "Sink"], sink_ids)
                        source_env <- setNames(env_map$env[env_map$role == "Source"], source_ids)
                        row_m <- match(paste0(sink_ids, "_", sink_env[sink_ids]), rownames(props_raw))
                        col_m <- match(paste0(source_ids, "_", source_env[source_ids]), colnames(props_raw))
                        if(any(is.na(row_m)) || any(is.na(col_m))){
                                stop("Unexpected FEAST output format! Cannot map the results back to sample names.")
                        }
                        props <- props_raw[row_m, c(col_m, ncol(props_raw)), drop = FALSE]
                        rownames(props) <- sink_ids
                        colnames(props) <- c(source_ids, "Unknown")
                        props_env <- private$aggregate_env(props = props, source_env = source_env[source_ids])
                        # ----- store results -----
                        if(is.null(label)){
                                label <- private$get_run_label()
                        }
                        self$res_mst[[label]] <- props
                        self$res_mst_env[[label]] <- props_env
                        self$res_mst_param[[label]] <- list(
                                sources = sources,
                                sinks = sinks,
                                source_samples = source_samples,
                                sink_samples = sink_samples,
                                source_ids = source_ids,
                                sink_ids = sink_ids,
                                source_env = source_env[source_ids],
                                sink_env = sink_env[sink_ids],
                                different_sources_flag = different_sources_flag,
                                id_col = id_col,
                                EM_iterations = EM_iterations,
                                COVERAGE = COVERAGE
                        )
                        message('The result is stored in object$res_mst[["', label, '"]] (sample level) and object$res_mst_env[["', label, '"]] (environment level) ...')
                        invisible(self)
                },
                #' @description
                #' Run sequential source tracking steps (e.g., bulk soil -> rhizosphere soil -> root).
                #'
                #' @param steps default NULL; a list where each element is a list of arguments
                #'   passed to \code{cal_mst} (without \code{label}), e.g.,
                #'   \code{list(list(sources = "Bulk soil", sinks = "Rhizosphere"),
                #'   list(sources = c("Bulk soil", "Rhizosphere"), sinks = "Endophyte"))}.
                #'   Each step is stored in \code{res_mst[["step1"]]}, \code{res_mst[["step2"]]}, ...
                #'   unless a custom \code{label} is provided inside the step.
                #' @return results stored in \code{res_mst}, \code{res_mst_env} and \code{res_mst_param}.
                #' @examples
                #' \dontrun{
                #' data(soil_microb)
                #' t1 <- trans_mst$new(dataset = soil_microb, env_col = "Compartment")
                #' t1$cal_mst_chain(list(
                #'   list(sources = "Bulk soil", sinks = "Rhizosphere"),
                #'   list(sources = c("Bulk soil", "Rhizosphere"), sinks = "Endophyte",
                #'     different_sources_flag = 1, id_col = "Plant_ID")
                #' ))
                #' }
                cal_mst_chain = function(steps = NULL){
                        if(!is.list(steps) || length(steps) == 0 || ! all(vapply(steps, is.list, logical(1)))){
                                stop("steps must be a list of lists, where each element contains the arguments for cal_mst!")
                        }
                        for(i in seq_along(steps)){
                                step <- steps[[i]]
                                label <- if(!is.null(step$label)) step$label else paste0("step", i)
                                step$label <- NULL
                                message("Running transmission step ", i, " of ", length(steps), " ...")
                                do.call(self$cal_mst, c(step, list(label = label)))
                        }
                        invisible(self)
                },
                #' @description
                #' Plot the source tracking result as a stacked barplot.
                #'
                #' @param use_run default NULL; the label (or index) of the run stored in \code{res_mst}.
                #'   If NULL, the first run is used.
                #' @param use_env default TRUE; whether to use the environment-level aggregated contributions
                #'   (\code{res_mst_env}). If FALSE, the sample-level contributions (\code{res_mst}) are shown.
                #' @param group_sink default NULL; a column name in \code{sample_table} used to group the sinks.
                #'   The plot is facetted by this column and the sinks are ordered accordingly, e.g., "Group".
                #' @param show_unknown default TRUE; whether to show the "Unknown" contribution.
                #' @param colors default NULL; a vector of fill colors for the sources.
                #'   If NULL, the Dark2 (environment level) or Paired (sample level) palette of RColorBrewer is used.
                #' @param unknown_color default "grey80"; the fill color for the "Unknown" source.
                #' @param legend_title default NULL; the legend title. If NULL, "Source" is used.
                #' @param bar_width default 0.8; the bar width passed to \code{\link{geom_bar}}.
                #' @param xtext_size default 10; x-axis text size.
                #' @param xtext_angle default 45; x-axis text angle.
                #' @param facet_scales default "free_x"; the scales parameter of \code{\link{facet_wrap}}.
                #' @return ggplot2 plot.
                #' @examples
                #' \dontrun{
                #' t1$plot_mst()
                #' t1$plot_mst(use_run = "soil2root", group_sink = "Group")
                #' }
                plot_mst = function(
                        use_run = NULL,
                        use_env = TRUE,
                        group_sink = NULL,
                        show_unknown = TRUE,
                        colors = NULL,
                        unknown_color = "grey80",
                        legend_title = NULL,
                        bar_width = 0.8,
                        xtext_size = 10,
                        xtext_angle = 45,
                        facet_scales = "free_x"
                        ){
                        run <- private$get_run(use_run = use_run)
                        label <- run$label
                        if(use_env){
                                props <- self$res_mst_env[[label]]
                                src_title <- "Source (Env)"
                        }else{
                                props <- self$res_mst[[label]]
                                src_title <- "Source"
                        }
                        if(is.null(props)){
                                stop("No result found for the run '", label, "'!")
                        }
                        # long format; melt keeps the row/column order as factor levels
                        plot_data <- reshape2::melt(
                                as.matrix(props),
                                varnames = c("Sink", "Source"),
                                value.name = "Contribution"
                        )
                        # in the different-sources mode at sample level, sources not used
                        # for a sink are NA; remove these entries for plotting
                        plot_data <- plot_data[! is.na(plot_data$Contribution), , drop = FALSE]
                        if(! show_unknown){
                                plot_data <- plot_data[plot_data$Source != "Unknown", , drop = FALSE]
                        }
                        # sink order: sample_table order; grouped by group_sink if provided
                        sink_order <- rownames(props)
                        if(!is.null(group_sink)){
                                if(! group_sink %in% colnames(self$sample_table)){
                                        stop("The group_sink '", group_sink, "' is not found in sample_table!")
                                }
                                gvals <- as.character(self$sample_table[sink_order, group_sink])
                                if(is.factor(self$sample_table[, group_sink])){
                                        glev <- intersect(levels(self$sample_table[, group_sink]), unique(gvals))
                                }else{
                                        glev <- unique(gvals)
                                }
                                plot_data$Group <- factor(gvals[match(plot_data$Sink, sink_order)], levels = glev)
                                sink_order <- sink_order[order(match(gvals, glev))]
                        }
                        plot_data$Sink <- factor(plot_data$Sink, levels = sink_order)
                        # fill colors
                        src_levels <- levels(plot_data$Source)
                        n_src <- length(src_levels)
                        has_unknown <- "Unknown" %in% src_levels
                        n_real <- n_src - ifelse(has_unknown, 1, 0)
                        if(is.null(colors)){
                                if(use_env){
                                        pal <- RColorBrewer::brewer.pal(8, "Dark2")
                                }else{
                                        pal <- RColorBrewer::brewer.pal(12, "Paired")
                                }
                                if(n_real <= length(pal)){
                                        colors <- pal[seq_len(n_real)]
                                }else{
                                        colors <- colorRampPalette(pal)(n_real)
                                }
                        }
                        if(has_unknown){
                                colors <- c(colors[seq_len(n_real)], unknown_color)
                        }
                        names(colors) <- src_levels
                        if(is.null(legend_title)){
                                legend_title <- src_title
                        }
                        p <- ggplot(plot_data, aes_meco(x = "Sink", y = "Contribution", fill = "Source")) +
                                geom_bar(stat = "identity", width = bar_width, color = "black", linewidth = 0.1) +
                                scale_fill_manual(values = colors, breaks = src_levels) +
                                xlab("Sink samples") + ylab("Contribution proportion") +
                                guides(fill = guide_legend(title = legend_title, reverse = TRUE)) +
                                theme_bw() +
                                ggplot_xtext_anglesize(xtext_angle = xtext_angle, xtext_size = xtext_size)
                        if(!is.null(group_sink)){
                                p <- p + facet_wrap(~Group, scales = facet_scales)
                        }
                        p
                },
                #' @description
                #' Print the trans_mst object.
                print = function() {
                        cat("trans_mst class object\n")
                        cat(paste0("Number of samples: ", self$n_samples, "\n"))
                        cat(paste0("Number of taxa: ", self$n_taxa, "\n"))
                        if(!is.null(self$env_col)){
                                cat(paste0("Env column: ", self$env_col, "\n"))
                        }else{
                                cat("Env column: Not provided\n")
                        }
                        if(is.null(self$res_mst)){
                                cat("Source tracking results: No\n")
                        }else{
                                cat("Source tracking runs: ", length(self$res_mst), "\n")
                                for(label in names(self$res_mst)){
                                        param <- self$res_mst_param[[label]]
                                        cat(paste0("  [", label, "] sources: ", length(param$source_ids),
                                                " samples (", paste(unique(param$source_env), collapse = ", "),
                                                "); sinks: ", length(param$sink_ids),
                                                " samples (", paste(unique(param$sink_env), collapse = ", "),
                                                "); different_sources_flag: ", param$different_sources_flag, "\n"))
                                }
                        }
                        invisible(self)
                }
        ),
        private = list(
                # resolve source/sink samples from env values or explicit sample names
                resolve_samples = function(sources = NULL, sinks = NULL, source_samples = NULL, sink_samples = NULL){
                        st <- self$sample_table
                        all_samples <- rownames(st)
                        if(!is.null(sources) && !is.null(source_samples)){
                                stop("Please provide either 'sources' (env values) or 'source_samples' (sample names), not both!")
                        }
                        if(!is.null(sinks) && !is.null(sink_samples)){
                                stop("Please provide either 'sinks' (env values) or 'sink_samples' (sample names), not both!")
                        }
                        # env value of each sample
                        if(!is.null(self$env_col)){
                                env_vals <- as.character(st[[self$env_col]])
                        }else{
                                if(!is.null(sources) || !is.null(sinks)){
                                        stop("env_col is required when sources/sinks are specified as env values! Please provide env_col when creating the object.")
                                }
                                env_vals <- all_samples
                        }
                        names(env_vals) <- all_samples
                        # resolve source samples
                        if(!is.null(sources)){
                                source_ids <- all_samples[env_vals %in% sources]
                                if(length(source_ids) == 0){
                                        stop("No samples found for the source env values: ", paste(sources, collapse = ", "), " !")
                                }
                                not_found <- setdiff(sources, unique(env_vals))
                                if(length(not_found) > 0){
                                        message("The following source env values are not found in sample_table$", self$env_col, ": ",
                                                paste(not_found, collapse = ", "), " ...")
                                }
                        }else if(!is.null(source_samples)){
                                not_found <- setdiff(source_samples, all_samples)
                                if(length(not_found) > 0){
                                        stop("The following source samples are not found in the dataset: ", paste(not_found, collapse = ", "), " !")
                                }
                                source_ids <- unique(source_samples)
                        }else{
                                stop("Please provide sources (env values) or source_samples (sample names)!")
                        }
                        # resolve sink samples
                        if(!is.null(sinks)){
                                sink_ids <- all_samples[env_vals %in% sinks]
                                if(length(sink_ids) == 0){
                                        stop("No samples found for the sink env values: ", paste(sinks, collapse = ", "), " !")
                                }
                                not_found <- setdiff(sinks, unique(env_vals))
                                if(length(not_found) > 0){
                                        message("The following sink env values are not found in sample_table$", self$env_col, ": ",
                                                paste(not_found, collapse = ", "), " ...")
                                }
                        }else if(!is.null(sink_samples)){
                                not_found <- setdiff(sink_samples, all_samples)
                                if(length(not_found) > 0){
                                        stop("The following sink samples are not found in the dataset: ", paste(not_found, collapse = ", "), " !")
                                }
                                sink_ids <- unique(sink_samples)
                        }else{
                                stop("Please provide sinks (env values) or sink_samples (sample names)!")
                        }
                        message("Source samples: ", length(source_ids), "; Sink samples: ", length(sink_ids), " ...")
                        list(
                                source_ids = source_ids,
                                sink_ids = sink_ids,
                                env_map = data.frame(
                                        sample_id = c(source_ids, sink_ids),
                                        env = c(env_vals[source_ids], env_vals[sink_ids]),
                                        role = c(rep("Source", length(source_ids)), rep("Sink", length(sink_ids))),
                                        stringsAsFactors = FALSE
                                )
                        )
                },
                # aggregate the sample-level proportions by source environment
                # NA handling: a source column is NA for sinks whose source set did not include it;
                # the aggregated value is NA only when all columns of the env are NA for that sink
                aggregate_env = function(props, source_env){
                        uenv <- unique(source_env)
                        agg <- matrix(NA_real_, nrow = nrow(props), ncol = length(uenv) + 1,
                                dimnames = list(rownames(props), c(uenv, "Unknown")))
                        for(e in uenv){
                                cols <- names(source_env)[source_env == e]
                                sub <- props[, cols, drop = FALSE]
                                all_na <- rowSums(is.na(sub)) == ncol(sub)
                                agg[! all_na, e] <- rowSums(sub[! all_na, , drop = FALSE], na.rm = TRUE)
                        }
                        agg[, "Unknown"] <- props[, "Unknown"]
                        agg
                },
                # generate the next available run label
                get_run_label = function(){
                        existing <- names(self$res_mst)
                        i <- 1
                        repeat{
                                label <- paste0("run", i)
                                if(! label %in% existing){
                                        return(label)
                                }
                                i <- i + 1
                        }
                },
                # find a run by label or index
                get_run = function(use_run = NULL){
                        if(is.null(self$res_mst)){
                                stop("No source tracking results found! Please first run cal_mst() ...")
                        }
                        if(is.null(use_run)){
                                label <- names(self$res_mst)[1]
                        }else{
                                if(is.numeric(use_run)){
                                        if(use_run < 1 || use_run > length(self$res_mst) || use_run != round(use_run)){
                                                stop("use_run index must be an integer between 1 and ", length(self$res_mst), "!")
                                        }
                                        label <- names(self$res_mst)[use_run]
                                }else{
                                        if(! use_run %in% names(self$res_mst)){
                                                stop("The run '", use_run, "' is not found in res_mst! Available runs: ",
                                                        paste(names(self$res_mst), collapse = ", "))
                                        }
                                        label <- use_run
                                }
                        }
                        list(label = label)
                }
        ),
	lock_class = FALSE,
	lock_objects = FALSE
)
