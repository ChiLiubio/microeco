#' @title trans_phylo class
#'
#' @description Publication-grade phylogenetic tree visualization built on ggtree,
#'   designed for microeco's microtable object. Supports group-based branch/tip
#'   coloring (with full clade coloring via groupOTU) and multi-layer outer ring
#'   annotations.
#'
#' @export
trans_phylo <- R6::R6Class(
        classname = "trans_phylo",

        # ===========================================================================
        # Public fields
        # ===========================================================================
        public = list(
                # =======================================================================
                # Initialize
                # =======================================================================
                #' @description Initialize a trans_phylo object
                #' @param dataset microtable object containing phylo_tree and tax_table
                #' @param group_col string, taxonomic rank column for branch/tip coloring, default "Phylum"
                #' @param ring_data data.frame, outer ring annotation data; row names or first column
                #'   should match tip labels. Supports numeric (heatmap) and categorical (color band) columns
                #' @param color_palette optional custom color palette; auto-generated if NULL
                #' @param group_prefix string, prefix to add to group names for legend display
                #' @param group_prefix_remove string, prefix to strip from group names (e.g. "p__")
                #' @param clade_coloring logical, if TRUE (default), propagate group info to internal
                #'   nodes via groupOTU so that entire clades are colored; if FALSE, groupOTU is not
                #'   applied and branch_color_by_group is automatically set to FALSE in plot_tree()
                #'   (since ggtree's stat_tree() requires groupOTU-annotated data for branch coloring).
                #'   Tip points are still colored by group via the fill aesthetic regardless.
                #' @examples
                #' \dontrun{
                #' library(microeco)
                #' library(magrittr)
                #'
                #' mt <- clone(dataset)
                #'
                #' mt$tax_table %<>% base::subset(Kingdom == "k__Bacteria")
                #' mt$tidy_dataset()
                #' mt200 <- clone(mt)
                #'
                #' # Subset to top 200 OTUs
                #' top_otus <- names(sort(mt200$taxa_sums(), decreasing = TRUE))[1:200]
                #' mt200$otu_table %<>% .[top_otus, , drop = FALSE]
                #' mt200$tidy_dataset()
                #' 
                #' # ========================================================================
                #' # Build rich ring annotation data
                #' # ========================================================================
                #' 
                #' otu <- mt200$otu_table
                #' tax <- mt200$tax_table
                #' si  <- mt200$sample_table
                #' n_otus <- nrow(otu)
                #' 
                #' # Relative abundance (mean %)
                #' total_reads <- colSums(mt$otu_table)
                #' total_rel_abund   <- sweep(mt$otu_table, 2, total_reads, "/")
                #' rel_abund   <- total_rel_abund[rownames(otu), ]
                #' mean_abund  <- rowMeans(rel_abund) * 100
                #' 
                #' # Prevalence (% of samples where detected)
                #' prevalence <- rowMeans(otu > 0) * 100
                #' 
                #' # Log10 abundance
                #' log_abund <- log10(rowMeans(otu) + 1)
                #' 
                #' # ---- Simulated ecological indicators ----
                #' # Importance index: composite of abundance and prevalence
                #' # (mimics indicators like indicator species analysis or random forest importance)
                #' importance_raw <- mean_abund * sqrt(prevalence / 100)
                #' importance <- round(importance_raw / max(importance_raw, na.rm = TRUE) * 100, 2)
                #' 
                #' # Specificity: how concentrated is an OTU in its preferred habitat
                #' # (low = generalist, high = specialist)
                #' group_cols  <- split(colnames(otu), si$Group)
                #' group_means <- sapply(group_cols, function(cols) rowMeans(rel_abund[, cols, drop = FALSE]))
                #' if (!is.matrix(group_means)) group_means <- matrix(group_means, nrow = n_otus)
                #' # Specificity = 1 - evenness (Simpson-like)
                #' specificity_raw <- 1 - apply(group_means, 1, function(p) {
                #'   p <- p / sum(p)
                #'   sum(p^2)
                #' })
                #' specificity <- round((1 - specificity_raw) / max(1 - specificity_raw, na.rm = TRUE) * 100, 2)
                #' 
                #' # Connectivity: simulated network centrality score (higher = hub OTU)
                #' # Based on co-occurrence signal approximated from abundance correlation
                #' set.seed(42)
                #' connectivity_raw <- importance_raw * runif(n_otus, 0.5, 1.5) + rnorm(n_otus, 0, 2)
                #' connectivity_raw[connectivity_raw < 0] <- 0
                #' connectivity <- round(connectivity_raw / max(connectivity_raw, na.rm = TRUE) * 100, 2)
                #' 
                #' # Dominant habitat
                #' dominant_habitat <- apply(group_means, 1, function(x) colnames(group_means)[which.max(x)])
                #' 
                #' # Salinity preference
                #' sal_cols  <- split(colnames(otu), si$Saline)
                #' sal_means <- sapply(sal_cols, function(cols) rowMeans(rel_abund[, cols, drop = FALSE]))
                #' if (!is.matrix(sal_means)) sal_means <- matrix(sal_means, nrow = n_otus)
                #' dominant_saline <- apply(sal_means, 1, function(x) colnames(sal_means)[which.max(x)])
                #' 
                #' # Class taxonomy
                #' class_col <- if ("Class" %in% colnames(tax)) tax$Class else rep("Unknown", n_otus)
                #' 
                #' # Assemble ring data
                #' ring_df <- data.frame(
                #'   label        = rownames(otu),
                #'   RelAbund     = round(mean_abund, 4),
                #'   Prevalence   = round(prevalence, 1),
                #'   LogAbund     = round(log_abund, 3),
                #'   Importance   = importance,
                #'   Specificity  = specificity,
                #'   Connectivity = connectivity,
                #'   Habitat      = dominant_habitat,
                #'   Salinity     = dominant_saline,
                #'   Class        = class_col,
                #'   stringsAsFactors = FALSE
                #' )
                #' 
                #' # Clean Class column
                #' ring_df$Class <- gsub("^c__", "", ring_df$Class)
                #' ring_df$Class[ring_df$Class == "" | is.na(ring_df$Class)] <- "Unclassified"
                #' 
                #' # ========================================================================
                #' # Custom style color palette for Phyla
                #' # ========================================================================
                #' 
                #' phyla_present <- unique(gsub("^p__", "", mt$tax_table$Phylum))
                #' n_phyla <- length(phyla_present)
                #' 
                #' phyla_colors <- c(
                #'   "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                #'   "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                #'   "#7E6148", "#B09C85", "#BB8FCE", "#5B88D4"
                #' )
                #' if (n_phyla <= length(phyla_colors)) {
                #'   phyla_colors <- phyla_colors[1:n_phyla]
                #' } else {
                #'   phyla_colors <- colorRampPalette(phyla_colors)(n_phyla)
                #' }
                #' names(phyla_colors) <- phyla_present
                #'
                #' pviz <- trans_phylo$new(
                #'   dataset             = mt,
                #'   group_col           = "Phylum",
                #'   group_prefix_remove = "p__",
                #'   ring_data           = ring_df,
                #'   color_palette       = list(group_color = phyla_colors),
                #'   clade_coloring      = TRUE
                #' )
                #' }
                initialize = function(dataset,
                                                          group_col = "Phylum",
                                                          ring_data = NULL,
                                                          color_palette = NULL,
                                                          group_prefix = NULL,
                                                          group_prefix_remove = NULL,
                                                          clade_coloring = TRUE) {

                        # --- Parameter validation ---
                        if (is.null(dataset)) {
                                stop("dataset (microtable object) cannot be NULL.")
                        }
                        if (!inherits(dataset, "microtable")) {
                                stop("dataset must be a microtable object from the microeco package.")
                        }
                        if (is.null(dataset$phylo_tree)) {
                                stop("The microtable object does not contain a phylo_tree slot.")
                        }

                        self$dataset            <- dataset
                        self$phylo_tree         <- dataset$phylo_tree
                        self$phylo_tree_original <- dataset$phylo_tree
                        self$tax_table          <- dataset$tax_table
                        self$group_col          <- group_col

                        # --- Extract tip labels ---
                        self$tip_labels <- self$phylo_tree$tip.label

                        # --- Build group mapping table ---
                        private$build_group_info(
                                group_col           = group_col,
                                group_prefix        = group_prefix,
                                group_prefix_remove = group_prefix_remove
                        )

                        # --- Store clade_coloring flag ---
                        self$clade_coloring <- clade_coloring

                        # --- Initialize ring extent tracker ---
                        self$max_ring_extent <- 0

                        # --- Propagate group info to internal nodes via groupOTU ---
                        if (clade_coloring) {
                                private$apply_groupOTU()
                        }

                        # --- Process outer ring annotation data ---
                        if (!is.null(ring_data)) {
                                self$ring_data <- private$process_ring_data(ring_data)
                        }

                        # --- Generate color palette ---
                        if (is.null(color_palette)) {
                                self$color_palette <- private$generate_palette()
                        } else {
                                self$color_palette <- color_palette
                        }

                        message("trans_phylo object initialized successfully.")
                        message(sprintf("  Tree tips: %d", length(self$tip_labels)))
                        message(sprintf("  Group levels (%s): %d", group_col, length(unique(self$group_info$Group))))
                        message(sprintf("  Clade coloring: %s", ifelse(clade_coloring, "enabled (groupOTU)", "tip-only")))
                        if (!is.null(self$ring_data)) {
                                message(sprintf("  Ring annotation columns: %d", ncol(self$ring_data) - 1))
                        }
                },

                # =======================================================================
                # plot_tree: Draw phylogenetic tree (supports multiple layouts)
                # =======================================================================
                #' @description Draw a phylogenetic tree (publication grade)
                #' @param layout tree layout: "fan", "circular", "radial", or "rectangular", default "fan"
                #' @param open_angle numeric, opening angle of the fan (0-360), default 30 (fan/circular only)
                #' @param branch_size branch line width, default 0.3
                #' @param tip_point logical, whether to show tip points, default TRUE
                #' @param tip_point_size tip point size, default 2
                #' @param tip_point_shape tip point shape, default 21 (filled circle with border)
                #' @param tip_point_stroke tip point border width, default 0.3
                #' @param tip_point_alpha tip point transparency, default 0.9
                #' @param show_tip_label logical, whether to show tip labels, default FALSE
                #' @param tip_label_size tip label font size, default 2
                #' @param tip_label_offset tip label offset, default 0.5
                #' @param show_node_label logical, whether to show internal node labels, default FALSE
                #' @param branch_color_by_group logical, color branches by group, default TRUE
                #' @param branch_color unified branch color (when branch_color_by_group = FALSE), default "grey40"
                #' @param legend_title legend title, default uses group_col
                #' @param legend_position legend position, default "right"
                #' @param legend_font_size legend font size, default 3.5
                #' @param theme_base base theme, default ggtree::theme_tree2()
                #' @return ggtree plot object
                #' @examples
                #' \dontrun{
                #' pviz$plot_tree(
                #'   layout             = "fan",
                #'   open_angle         = 18,
                #'   branch_size        = 0.15,
                #'   tip_point          = TRUE,
                #'   tip_point_size     = 1.2,
                #'   tip_point_shape    = 21,
                #'   tip_point_stroke   = 0.1,
                #'   tip_point_alpha    = 0.80,
                #'   legend_title       = "Phylum"
                #' )
                #' }
                plot_tree = function(layout               = "fan",
                                                         open_angle           = 30,
                                                         branch_size          = 0.3,
                                                         tip_point            = TRUE,
                                                         tip_point_size       = 2,
                                                         tip_point_shape      = 21,
                                                         tip_point_stroke     = 0.3,
                                                         tip_point_alpha      = 0.9,
                                                         show_tip_label       = FALSE,
                                                         tip_label_size       = 2,
                                                         tip_label_offset     = 0.5,
                                                         show_node_label      = FALSE,
                                                         branch_color_by_group = TRUE,
                                                         branch_color         = "grey40",
                                                         legend_title         = NULL,
                                                         legend_position      = "right",
                                                         legend_font_size     = 3.5,
                                                         theme_base           = NULL) {

                        # --- Build base ggtree object ---
                        # ggtree's internal stat_tree() requires lowercase 'group' for branch coloring.
                        # We use lowercase 'group' for edges (from groupOTU) and uppercase 'Group'
                        # for tip point fill (from self$group_info) to keep them isolated.
                        #
                        # IMPORTANT: When clade_coloring = FALSE, groupOTU is not applied, so the
                        # tree lacks the 'group' column. The old fallback path (using %<+% then
                        # ggtree::geom_tree(ggplot2::aes(color=Group))) is broken because ggtree's stat_tree() cannot
                        # resolve columns attached via %<+% — it only sees columns in the tree object
                        # itself (from groupOTU). Therefore, when clade_coloring = FALSE, we
                        # automatically set branch_color_by_group = FALSE. Tip points are still
                        # colored by Group via the fill aesthetic.

                        if (!self$clade_coloring && branch_color_by_group) {
                                message("Note: clade_coloring = FALSE, so branch_color_by_group is set to FALSE. ",
                                                "Tip points are still colored by group via fill.")
                                branch_color_by_group <- FALSE
                        }

                        # Unified legend title and breaks for both branch color and tip point fill
                        leg_title  <- ifelse(is.null(legend_title), self$group_col, legend_title)
                        pal        <- self$color_palette$group_color
                        # Use explicit breaks to ensure color and fill scales map identically.
                        # groupOTU's 'group' column includes "0" for unassigned internal nodes,
                        # which is NOT in pal. Setting breaks = names(pal) ensures both scales
                        # show exactly the same categories in the same order, so ggplot2 can
                        # merge them into a single legend with consistent colors.
                        leg_breaks <- names(pal)

                        if (branch_color_by_group && "group" %in% colnames(tidytree::as_tibble(self$phylo_tree))) {
                                # Tree has group column from groupOTU -> color entire clades
                                p <- ggtree::ggtree(self$phylo_tree,
                                                        layout     = layout,
                                                        open.angle = open_angle,
                                                        linewidth  = branch_size,
                                                        ggplot2::aes(color  = group))
                                p <- p + ggplot2::scale_color_manual(
                                        values   = pal,
                                        breaks   = leg_breaks,
                                        name     = leg_title,
                                        guide    = "legend",
                                        na.value = "grey70"
                                )
                        } else {
                                # No group-based branch coloring (either disabled or groupOTU not applied)
                                p <- ggtree::ggtree(self$phylo_tree,
                                                        layout     = layout,
                                                        open.angle = open_angle,
                                                        linewidth  = branch_size)
                                p <- p + ggtree::geom_tree(color = branch_color, linewidth = branch_size)
                        }

                        # --- Tip points ---
                        if (tip_point) {
                                # Use shape = 21 (fillable circle) so that 'fill' maps to Group,
                                # keeping 'color' aesthetic reserved for tree branches only.
                                # IMPORTANT: color = NA removes the border so that the perceived
                                # color of the point comes entirely from 'fill', matching the
                                # branch 'color' exactly. A visible border (e.g. "grey30") makes
                                # points appear darker/different from branches of the same group.
                                # Attach group_info for tip fill mapping
                                group_df <- self$group_info
                                p <- ggtree::`%<+%`(p, group_df)

                                p <- p + ggtree::geom_tippoint(
                                        ggplot2::aes(fill = Group),
                                        size   = tip_point_size,
                                        shape  = tip_point_shape,
                                        stroke = tip_point_stroke,
                                        alpha  = tip_point_alpha,
                                        color  = "transparent",
                                        na.rm  = TRUE
                                )
                                # Add fill scale for tip points FIRST, THEN create new fill scale boundary.
                                # ggnewscale::new_scale("fill") acts as a separator: scales before it apply
                                # to geoms before it, scales after it apply to geoms after it.
                                # The old code had new_scale BEFORE scale_fill_manual, which caused the fill
                                # scale to be assigned to ring layers instead of the tip point layer.
                                p <- p + ggplot2::scale_fill_manual(
                                        values   = pal,
                                        breaks   = leg_breaks,
                                        name     = leg_title,
                                        guide    = ggplot2::guide_legend(override.aes = list(size = 3, shape = 21, linetype = 1)),
                                        na.value = "grey70"
                                )
                                # Create new fill scale boundary for subsequent ring layers
                                p <- p + ggnewscale::new_scale("fill")
                        }

                        # --- Tip labels ---
                        if (show_tip_label) {
                                if (layout %in% c("fan", "circular", "radial")) {
                                        p <- p + ggtree::geom_tiplab2(size = tip_label_size, offset = tip_label_offset)
                                } else {
                                        p <- p + ggtree::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
                                }
                        }

                        # --- Internal node labels ---
                        if (show_node_label) {
                                p <- p + ggtree::geom_nodelab(size = 2, hjust = 0.5)
                        }

                        # --- Theme ---
                        if (is.null(theme_base)) {
                                p <- p + ggtree::theme_tree2()
                        } else {
                                p <- p + theme_base
                        }

                        # --- Legend styling ---
                        p <- p + ggplot2::theme(
                                legend.position    = legend_position,
                                legend.text        = ggplot2::element_text(size = legend_font_size),
                                legend.title       = ggplot2::element_text(size = legend_font_size + 1, face = "bold"),
                                legend.key.size    = ggplot2::unit(0.35, "cm"),
                                legend.key.spacing = ggplot2::unit(0.15, "cm"),
                                legend.box.spacing = ggplot2::unit(0.3, "cm"),
                                panel.background   = ggplot2::element_blank(),
                                plot.background    = ggplot2::element_rect(fill = "white", color = NA)
                                )

                        self$plot_obj <- p
                        message(sprintf("Tree plotted (layout: %s). Use $add_ring() to add outer annotations.", layout))
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # plot_circular: Convenience wrapper for plot_tree with fan layout
                # =======================================================================
                #' @description Convenience wrapper for plot_tree() with layout = "fan".
                #'   All parameters are passed through to plot_tree().
                #' @param ... all arguments passed to plot_tree()
                #' @return ggtree plot object
                #' @examples
                #' \dontrun{
                #' pviz$plot_circular(
                #'   open_angle       = 18,
                #'   branch_size      = 0.15,
                #'   tip_point_size   = 1.2,
                #'   tip_point_shape  = 21,
                #'   tip_point_stroke = 0.1,
                #'   tip_point_alpha  = 0.80,
                #'   legend_title     = "Phylum"
                #' )
                #' }
                plot_circular = function(...) {
                        args <- list(...)
                        if (!"layout" %in% names(args)) {
                                args$layout <- "fan"
                        }
                        do.call(self$plot_tree, args)
                },

                # =======================================================================
                # add_ring: Add one outer ring annotation layer
                # =======================================================================
                #' @description Add one outer ring annotation layer to the current plot.
                #'   Extra arguments (...) are passed to ggtreeExtra::geom_fruit() for advanced customization
                #'   (e.g. axis.params, grid.params, direction).
                #' @param col_name string, column name from ring_data to display
                #' @param ring_data optional custom data.frame for this ring (overrides self$ring_data)
                #' @param color_low color for low numeric values, default "#0D0887"
                #' @param color_high color for high numeric values, default "#F0F921"
                #' @param limits numeric vector of length 2, explicit limits for the fill scale;
                #'   default NULL means auto from data range
                #' @param categorical_colors color vector for categorical columns, default NULL (auto)
                #' @param ring_width ring width, default 0.08
                #' @param ring_offset ring offset, default 0.03
                #' @param show_legend logical, whether to show this ring's legend, default TRUE
                #' @param legend_title_custom custom legend title, default uses column name
                #' @param geom type: "bar" for bar heatmap, "point" for scatter, "star" for stars;
                #'   categorical columns automatically use color band regardless of this setting
                #' @param na_fill fill color for NA values, default "grey90"
                #' @param na_replace categorical value to replace NA/empty with, default "Unknown"
                #' @param point_size point size (only for geom = "point"), default 1.5
                #' @param ... additional arguments passed to ggtreeExtra::geom_fruit() (e.g. axis.params, grid.params)
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' # Numeric bar ring
                #' pviz$add_ring(
                #'   col_name            = "LogAbund",
                #'   geom                = "bar",
                #'   ring_width          = 0.038,
                #'   ring_offset         = 0.012,
                #'   color_low           = "#440154",
                #'   color_high          = "#FDE725",
                #'   na_fill             = "grey95",
                #'   legend_title_custom = "Log10\\nAbund."
                #' )
                #'
                #' # Categorical color band ring
                #' pviz$add_ring(
                #'   col_name            = "Habitat",
                #'   ring_width          = 0.022,
                #'   ring_offset         = 0.242,
                #'   na_fill             = "grey95",
                #'   legend_title_custom = "Habitat",
                #'   categorical_colors  = c("IW" = "#E64B35", "CW" = "#4DBBD5", "TW" = "#00A087")
                #' )
                #' }
                add_ring = function(col_name,
                                                        ring_data           = NULL,
                                                        color_low           = "#0D0887",
                                                        color_high          = "#F0F921",
                                                        limits              = NULL,
                                                        categorical_colors  = NULL,
                                                        ring_width          = 0.08,
                                                        ring_offset         = 0.03,
                                                        show_legend         = TRUE,
                                                        legend_title_custom = NULL,
                                                        geom                = "bar",
                                                        na_fill             = "grey90",
                                                        na_replace          = "Unknown",
                                                        point_size          = 1.5,
                                                        ...) {

                        if (is.null(self$plot_obj)) {
                                stop("No plot object found. Please run $plot_tree() or $plot_circular() first.")
                        }

                        # --- Determine data source ---
                        df <- if (!is.null(ring_data)) ring_data else self$ring_data
                        if (is.null(df)) {
                                stop("No ring data available. Provide ring_data in add_ring() or during initialization.")
                        }
                        if (!col_name %in% colnames(df)) {
                                stop(sprintf("Column '%s' not found in ring data.", col_name))
                        }

                        # --- Build sub data.frame ---
                        sub_df <- df[, c("label", col_name), drop = FALSE]
                        colnames(sub_df)[2] <- "value"

                        # --- Detect data type ---
                        is_numeric <- is.numeric(sub_df$value)

                        # --- Compute limits for numeric scales ---
                        if (is_numeric && is.null(limits)) {
                                limits <- range(sub_df$value, na.rm = TRUE)
                        }

                        # IMPORTANT: All ring annotations use the 'fill' aesthetic exclusively.
                        # The 'color' aesthetic is reserved for the tree's discrete Group variable.
                        # We use ggnewscale::new_scale("fill") before each ring to allow independent fill scales.
                        self$plot_obj <- self$plot_obj + ggnewscale::new_scale("fill")

                        # --- Capture and merge user-provided ... arguments ---
                        # We must handle axis.params / grid.params specially to avoid
                        # "formal argument matched by multiple actual arguments" errors.
                        # Strategy: capture ... as a list, merge defaults with user values,
                        # then use rlang::inject() to splice into ggtreeExtra::geom_fruit().
                        # rlang::inject() is used instead of do.call() because ggtreeExtra::geom_fruit()
                        # uses substitute()/deparse() on the 'geom' argument, and do.call()
                        # evaluates the function object first, making it a closure that
                        # cannot be coerced to character.
                        dot_args <- list(...)

                        # Merge axis.params: default = list(add = FALSE); user values override
                        default_axis <- list(add = FALSE)
                        if (!is.null(dot_args$axis.params)) {
                                dot_args$axis.params <- modifyList(default_axis, dot_args$axis.params)
                        } else {
                                dot_args$axis.params <- default_axis
                        }

                        # Merge grid.params: default = empty list (no grid); user values override
                        default_grid <- list()
                        if (!is.null(dot_args$grid.params)) {
                                dot_args$grid.params <- modifyList(default_grid, dot_args$grid.params)
                        } else {
                                dot_args$grid.params <- default_grid
                        }

                        if (is_numeric) {
                                # Numeric -> heatmap bar / point / star
                                if (geom == "bar") {
                                        fruit_layer <- rlang::inject(ggtreeExtra::geom_fruit(
                                                data    = sub_df,
                                                geom    = geom_col,
                                                mapping = ggplot2::aes(y = label, x = value, fill = value),
                                                offset  = ring_offset,
                                                pwidth  = ring_width,
                                                !!!dot_args
                                                )
                                        )
                                        # Build fill scale — DRY: construct once, reuse for all numeric geoms
                                        legend_name  <- ifelse(is.null(legend_title_custom), col_name, legend_title_custom)
                                        legend_guide <- if (show_legend) ggplot2::guide_legend() else ggplot2::guide_none()
                                        if (!identical(color_low, "#0D0887") || !identical(color_high, "#F0F921")) {
                                                fill_scale <- ggplot2::scale_fill_gradient(
                                                        low      = color_low,
                                                        high     = color_high,
                                                        name     = legend_name,
                                                        na.value = na_fill,
                                                        limits   = limits,
                                                        guide    = legend_guide
                                                )
                                        } else {
                                                fill_scale <- ggplot2::scale_fill_viridis_c(
                                                        option   = "D",
                                                        begin    = 0,
                                                        end      = 1,
                                                        name     = legend_name,
                                                        na.value = na_fill,
                                                        limits   = limits,
                                                        guide    = legend_guide
                                                )
                                        }
                                        self$plot_obj <- self$plot_obj +
                                                fruit_layer +
                                                fill_scale
                                } else if (geom == "point") {
                                        fruit_layer <- rlang::inject(ggtreeExtra::geom_fruit(
                                                data    = sub_df,
                                                geom    = geom_point,
                                                mapping = ggplot2::aes(y = label, x = value, fill = value),
                                                shape   = 21,
                                                size    = point_size,
                                                color   = "grey30",
                                                offset  = ring_offset,
                                                pwidth  = ring_width,
                                                !!!dot_args
                                                )
                                        )
                                        # Reuse the same fill_scale logic (DRY)
                                        self$plot_obj <- self$plot_obj +
                                                fruit_layer +
                                                fill_scale
                                } else if (geom == "star") {
                                        fruit_layer <- rlang::inject(ggtreeExtra::geom_fruit(
                                                data    = sub_df,
                                                geom    = geom_point,
                                                mapping = ggplot2::aes(y = label, x = value, fill = value),
                                                shape   = 8,
                                                size    = 2,
                                                offset  = ring_offset,
                                                pwidth  = ring_width,
                                                !!!dot_args
                                                )
                                        )
                                        # Reuse the same fill_scale logic (DRY)
                                        self$plot_obj <- self$plot_obj +
                                                fruit_layer +
                                                fill_scale
                                }
                        } else {
                                # Categorical -> color band
                                # Sanitize NA and empty strings before as.factor
                                sub_df$value[is.na(sub_df$value)] <- na_replace
                                sub_df$value[sub_df$value == ""]  <- na_replace
                                sub_df$value <- as.factor(sub_df$value)
                                n_levels     <- length(levels(sub_df$value))

                                if (is.null(categorical_colors)) {
                                        if (n_levels <= 9) {
                                                categorical_colors <- RColorBrewer::brewer.pal(max(3, n_levels), "Set1")[1:n_levels]
                                        } else if (n_levels <= 12) {
                                                categorical_colors <- RColorBrewer::brewer.pal(max(3, n_levels), "Set3")[1:n_levels]
                                        } else {
                                                categorical_colors <- colorRampPalette(
                                                  RColorBrewer::brewer.pal(9, "Set1")
                                                )(n_levels)
                                        }
                                }
                                # Ensure names match levels
                                names(categorical_colors) <- levels(sub_df$value)

                                fruit_layer <- rlang::inject(ggtreeExtra::geom_fruit(
                                        data    = sub_df,
                                        geom    = geom_tile,
                                        mapping = ggplot2::aes(y = label, x = 1, fill = value),
                                        offset  = ring_offset,
                                        pwidth  = ring_width,
                                        width   = 1,
                                        !!!dot_args
                                        )
                                )
                                self$plot_obj <- self$plot_obj +
                                        fruit_layer +
                                        ggplot2::scale_fill_manual(
                                                values   = categorical_colors,
                                                name     = ifelse(is.null(legend_title_custom), col_name, legend_title_custom),
                                                drop     = FALSE,
                                                na.value = na_fill,
                                                guide    = if (show_legend) ggplot2::guide_legend() else ggplot2::guide_none()
                                        )
                        }

                        # --- Track ring position for add_spacer() ---
                        self$max_ring_extent <- max(self$max_ring_extent, ring_offset + ring_width)

                        message(sprintf("Ring annotation '%s' added. (type: %s)", col_name,
                                                   ifelse(is_numeric, "numeric", "categorical")))
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # add_spacer: Add invisible tile spacer to prevent outermost ring stretching
                # =======================================================================
                #' @description Add an invisible spacer ring AFTER the last visible ring to
                #'   prevent ggtreeExtra from auto-stretching the outermost ring.
                #'
                #'   ggtreeExtra::geom_fruit() auto-expands the LAST (outermost) ring to
                #'   fill the remaining panel space, making it much thicker than the specified
                #'   ring_width (pwidth). This spacer uses geom_tile (same as categorical rings)
                #'   to render a uniform invisible white band as the last layer. Since geom_tile
                #'   fills a complete angular band regardless of data values, even if ggtreeExtra
                #'   auto-expands this spacer, it just adds more invisible white space — the
                #'   last VISIBLE ring stays at its correct width, and the tree layout is
                #'   completely unaffected.
                #'
                #'   IMPORTANT: This does NOT change the coordinate system, so the tree shape
                #'   and proportions remain exactly the same as without the spacer.
                #'
                #'   Call AFTER all add_ring() calls and BEFORE save_plot().
                #'   The offset is auto-calculated from the tracked ring positions.
                #' @param ring_width numeric, radial width of the spacer band, default 0.02.
                #'   This should be similar to other categorical ring widths. Even if
                #'   auto-expanded by ggtreeExtra, the white band remains invisible.
                #' @param gap numeric, small gap between the last visible ring's outer edge
                #'   and the spacer's inner edge, default 0.004
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' # Call after all add_ring() calls to prevent the outermost ring
                #' # from being auto-expanded by ggtreeExtra
                #' pviz$add_spacer()
                #' }
                add_spacer = function(ring_width = 0.02, gap = 0.004) {
                        if (is.null(self$plot_obj)) {
                                stop("No plot object found. Please run $plot_tree() or $plot_circular() first.")
                        }

                        if (self$max_ring_extent == 0) {
                                message("No rings have been added yet. Spacer is not needed.")
                                return(invisible(self$plot_obj))
                        }

                        # Auto-calculate offset: just beyond the outermost ring + small gap
                        spacer_offset <- self$max_ring_extent + gap

                        # Create spacer data: a single dummy categorical value for all tips
                        # Using the same structure as categorical rings (geom_tile)
                                spacer_df <- data.frame(
                                label  = self$tip_labels,
                                spacer = rep("_spacer", length(self$tip_labels)),
                                stringsAsFactors = FALSE
                        )

                        # Add new fill scale for this ring (required for independent scales)
                        self$plot_obj <- self$plot_obj + ggnewscale::new_scale("fill")

                        # Local reference for ggtreeExtra::geom_fruit() geom argument.
                        # .convert_to_name() requires a bare function name, not ggplot2::geom_tile.
                        geom_tile <- ggplot2::geom_tile

                        # Use geom_tile — same as categorical rings.
                        # This renders a complete uniform band regardless of data values.
                        # Even if ggtreeExtra auto-expands the last ring, the extra space
                        # is just white (invisible) — the last VISIBLE ring keeps its width.
                        # color = NA removes the tile border outlines (otherwise visible as
                        # thin radiating lines even when fill is white).
                        spacer_layer <- ggtreeExtra::geom_fruit(
                                data        = spacer_df,
                                geom        = geom_tile,
                                mapping     = ggplot2::aes(y = label, x = 1, fill = spacer),
                                color       = NA,
                                linewidth   = 0,
                                offset      = spacer_offset,
                                pwidth      = ring_width,
                                width       = 1,
                                axis.params = list(add = FALSE),
                                grid.params = list()
                        )

                        # White fill, no legend — completely invisible
                        self$plot_obj <- self$plot_obj +
                        spacer_layer +
                        ggplot2::scale_fill_manual(
                                values   = c("_spacer" = "white"),
                                guide    = "none",
                                na.value = "white"
                        )

                        message(sprintf(paste0("Spacer ring added (offset: %.3f, width: %.3f). ",
                                                                  "Uses geom_tile — invisible white band absorbs ",
                                                                  "ggtreeExtra's auto-expansion without affecting tree layout."),
                                                   spacer_offset, ring_width))
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # add_rings_batch: Batch add multiple ring annotations
                # =======================================================================
                #' @description Batch add all (or selected) columns from ring_data as outer rings
                #' @param col_names character vector of column names to add; NULL means all columns
                #' @param numeric_geom geom type for numeric columns, default "bar"
                #' @param ring_width ring width, default 0.08
                #' @param ring_offset spacing between rings, default 0.02
                #' @param color_low color for low numeric values, default "#0D0887"
                #' @param color_high color for high numeric values, default "#F0F921"
                #' @param categorical_palette categorical palette function, default NULL
                #' @param show_legend logical, whether to show legends
                #' @param ... additional arguments passed to each add_ring() call
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' pviz$add_rings_batch(
                #'   col_names    = c("RelAbund", "Prevalence", "LogAbund"),
                #'   numeric_geom = "bar",
                #'   ring_width   = 0.04,
                #'   ring_offset  = 0.02,
                #'   show_legend  = TRUE
                #' )
                #' }
                add_rings_batch = function(col_names           = NULL,
                                                                   numeric_geom        = "bar",
                                                                   ring_width          = 0.08,
                                                                   ring_offset         = 0.02,
                                                                   color_low           = "#0D0887",
                                                                   color_high          = "#F0F921",
                                                                   categorical_palette = NULL,
                                                                   show_legend         = TRUE,
                                                                   ...) {

                        if (is.null(self$ring_data) && is.null(col_names)) {
                                stop("No ring data available.")
                        }

                        df <- self$ring_data
                        all_cols <- setdiff(colnames(df), "label")

                        if (!is.null(col_names)) {
                                all_cols <- intersect(all_cols, col_names)
                        }

                        if (length(all_cols) == 0) {
                                stop("No valid columns found for ring annotation.")
                        }

                        cumulative_offset <- ring_offset

                        for (i in seq_along(all_cols)) {
                                cn <- all_cols[i]
                                self$add_ring(
                                        col_name            = cn,
                                        ring_data           = df,
                                        ring_width          = ring_width,
                                        ring_offset         = cumulative_offset,
                                        geom                = numeric_geom,
                                        color_low           = color_low,
                                        color_high          = color_high,
                                        show_legend         = show_legend,
                                        legend_title_custom = cn,
                                        ...
                                )
                                # Accumulate offset so rings are arranged from inner to outer
                                cumulative_offset <- cumulative_offset + ring_width + ring_offset
                        }

                        message(sprintf("Batch added %d ring annotation(s).", length(all_cols)))
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # add_scale_bar: Add evolutionary distance scale bar
                # =======================================================================
                #' @description Add an evolutionary distance scale bar to the plot
                #' @param x numeric, x coordinate
                #' @param y numeric, y coordinate
                #' @param label string, scale bar label
                #' @return updated plot
                #' @examples
                #' \dontrun{
                #' pviz$add_scale_bar(x = 0, y = 0, label = "0.05")
                #' }
                add_scale_bar = function(x = NULL, y = NULL, label = "0.05") {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        self$plot_obj <- self$plot_obj +
                        ggtree::geom_treescale(
                                x        = x,
                                y        = y,
                                width    = as.numeric(label),
                                label    = label,
                                fontsize = 3,
                                linesize = 0.5,
                                offset   = 0.2
                        )
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # add_highlight_clade: Highlight a specified clade
                # =======================================================================
                #' @description Add highlight background to a specified clade
                #' @param node node number (can be found in $plot_obj$data)
                #' @param fill highlight fill color, default "steelblue"
                #' @param alpha transparency, default 0.3
                #' @param extend extend amount, default 0.15
                #' @return updated plot
                #' @examples
                #' \dontrun{
                #' # Highlight clade at node 150
                #' pviz$add_highlight_clade(node = 150, fill = "steelblue", alpha = 0.3)
                #' }
                add_highlight_clade = function(node,
                                                                                fill   = "steelblue",
                                                                                alpha  = 0.3,
                                                                                extend = 0.15) {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        self$plot_obj <- self$plot_obj +
                                ggtree::geom_highlight(
                                node     = node,
                                fill     = fill,
                                alpha    = alpha,
                                extend   = extend,
                                extendto = NULL
                        )
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # add_clade_label: Add clade annotation text
                # =======================================================================
                #' @description Add annotation text to a specified clade
                #' @param node node number
                #' @param label annotation text
                #' @param color text color, default "black"
                #' @param size text size, default 4
                #' @param offset offset amount, default 0.5
                #' @param barsize bar thickness, default 0.8
                #' @return updated plot
                #' @examples
                #' \dontrun{
                #' pviz$add_clade_label(node = 150, label = "Clade A", color = "black", size = 4)
                #' }
                add_clade_label = function(node,
                                                                        label,
                                                                        color    = "black",
                                                                        size     = 4,
                                                                        offset   = 0.5,
                                                                        barsize  = 0.8) {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        self$plot_obj <- self$plot_obj +
                        ggtree::geom_cladelab(
                                node       = node,
                                label      = label,
                                color      = color,
                                size       = size,
                                offset     = offset,
                                barsize    = barsize,
                                horizontal = TRUE,
                                angle      = 0
                        )
                        invisible(self$plot_obj)
                },

                # =======================================================================
                # theme_publication: Apply publication-grade theme
                # =======================================================================
                #' @description Apply a publication-grade theme to the current plot
                #' @param base_size base font size, default 12
                #' @param base_family font family, default "sans"
                #' @param bg_color background color, default "white"
                #' @return updated plot
                #' @examples
                #' \dontrun{
                #' pviz$theme_publication(base_size = 7)
                #' }
                theme_publication = function(base_size   = 12,
                                                                          base_family = "sans",
                                                                          bg_color    = "white") {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        self$plot_obj <- self$plot_obj +
                        ggplot2::theme(
                                text                 = ggplot2::element_text(size = base_size, family = base_family),
                                plot.background      = ggplot2::element_rect(fill = bg_color, color = NA),
                                panel.background     = ggplot2::element_blank(),
                                legend.position      = "right",
                                legend.text          = ggplot2::element_text(size = base_size * 0.75),
                                legend.title         = ggplot2::element_text(size = base_size * 0.85, face = "bold"),
                                legend.key.size      = ggplot2::unit(0.4, "cm"),
                                legend.key           = ggplot2::element_rect(fill = "transparent", color = NA),
                                legend.spacing       = ggplot2::unit(0.2, "cm"),
                                legend.box           = "vertical",
                                legend.justification = c(0, 0.5),
                                plot.margin          = ggplot2::margin(t = 5, r = 15, b = 5, l = 5, unit = "pt")
                        )

                        invisible(self$plot_obj)
                },
				
                # =======================================================================
                # get_plot: Retrieve current plot object
                # =======================================================================
                #' @description Return the current ggtree plot object
                #' @return ggtree / ggplot object
                #' @examples
                #' \dontrun{
                #' p <- pviz$get_plot()
                #' print(p)
                #' }
                get_plot = function() {
                        if (is.null(self$plot_obj)) {
                                warning("No plot object available. Run $plot_tree() first.")
                        }
                        return(self$plot_obj)
                },
				
                # =======================================================================
                # save_plot: Save plot to file
                # =======================================================================
                #' @description Save the current plot to a file
                #' @param filename file path; format inferred from extension
                #' @param width plot width, default 12
                #' @param height plot height, default 12
                #' @param dpi resolution, default 300
                #' @param device graphics device, default NULL (auto)
                #' @return NULL
                #' @examples
                #' \dontrun{
                #' pviz$save_plot(
                #'   filename = "phylo_tree.png",
                #'   width    = 16,
                #'   height   = 16,
                #'   dpi      = 300
                #' )
                #' 
                #' 				
                #' }
                save_plot = function(filename,
                                                         width   = 12,
                                                         height  = 12,
                                                         dpi     = 300,
                                                         device  = NULL) {
                        if (is.null(self$plot_obj)) {
                                stop("No plot object to save. Run $plot_tree() first.")
                        }
                        ggplot2::ggsave(
                                filename = filename,
                                plot     = self$plot_obj,
                                width    = width,
                                height   = height,
                                dpi      = dpi,
                                device   = device,
                                bg       = "white"
                        )
                        message(sprintf("Plot saved to: %s", filename))
                        invisible(NULL)
                },

                # =======================================================================
                # print: Print method
                # =======================================================================
                #' @description Print a summary of the trans_phylo object
                #' @param ... ignored (for compatibility with the generic print method)
                #' @examples
                #' \dontrun{
                #' print(pviz)
                #' }
                print = function(...) {
                        cat("<trans_phylo> object\n")
                        cat(sprintf("  Phylo tree tips: %d\n", length(self$tip_labels)))
                        cat(sprintf("  Group column: %s\n", self$group_col))
                        if (!is.null(self$group_info)) {
                                cat(sprintf("  Number of groups: %d\n", length(unique(self$group_info$Group))))
                        }
                        cat(sprintf("  Clade coloring (groupOTU): %s\n", ifelse(self$clade_coloring, "yes", "no")))
                        if (!is.null(self$ring_data)) {
                                ring_cols <- setdiff(colnames(self$ring_data), "label")
                                cat(sprintf("  Ring annotation columns (%d): %s\n",
                                                        length(ring_cols), paste(ring_cols, collapse = ", ")))
                        }
                        if (!is.null(self$plot_obj)) {
                                cat("  Plot object: available\n")
                        } else {
                                cat("  Plot object: not yet created\n")
                        }
                        invisible(self)
                }
        ),

        # ===========================================================================
        # Private methods
        # ===========================================================================
        private = list(
        # Build group mapping table from tax_table
        # group_col string, taxonomic rank column name
        # group_prefix string, prefix to add to group names
        # group_prefix_remove string, prefix to strip from group names
        build_group_info = function(group_col, group_prefix, group_prefix_remove) {
                tax <- self$tax_table
                # Ensure tax_table is a data.frame
                if (!is.data.frame(tax)) {
                        tax <- as.data.frame(tax, stringsAsFactors = FALSE)
                }
                # Row names are OTU/ASV identifiers
                if (!"label" %in% colnames(tax)) {
                        tax$label <- rownames(tax)
                }

                # Check that group_col exists
                if (!group_col %in% colnames(tax)) {
                        stop(sprintf("Column '%s' not found in tax_table. Available: %s",
                                         group_col, paste(colnames(tax), collapse = ", ")))
                }

                group_df <- tax[, c("label", group_col), drop = FALSE]
                colnames(group_df)[2] <- "Group"

                # Keep only tips present in the tree
                group_df <- group_df[group_df$label %in% self$tip_labels, ]

                # Remove prefix
                if (!is.null(group_prefix_remove)) {
                        group_df$Group <- gsub(group_prefix_remove, "", group_df$Group)
                }

                # Add prefix
                if (!is.null(group_prefix)) {
                        group_df$Group <- paste0(group_prefix, group_df$Group)
                }

                # Handle NA / unknown groups
                group_df$Group[is.na(group_df$Group) | group_df$Group == ""] <- "Others"

                self$group_info <- group_df
        },
        # Propagate group info to internal nodes via groupOTU
        apply_groupOTU = function() {
                # Build a named list: group -> vector of tip labels
                group_list <- split(self$group_info$label, self$group_info$Group)

                # Apply groupOTU to annotate internal nodes with Group
                # This ensures entire clades are colored, not just tip branches
                # ggtree's stat_tree() internally looks for lowercase 'group';
                # using group_name = "Group" (uppercase) causes 'object Group not found' error.
                self$phylo_tree <- ggtree::groupOTU(self$phylo_tree, group_list, group_name = "group")
        },
        # Process and validate outer ring annotation data
        # ring_data data.frame or matrix of ring annotation data
        process_ring_data = function(ring_data) {
                if (is.data.frame(ring_data)) {
                        df <- ring_data
                } else if (is.matrix(ring_data)) {
                        df <- as.data.frame(ring_data, stringsAsFactors = FALSE)
                } else {
                        stop("ring_data must be a data.frame or matrix.")
                }

                # Identify the label column: either a column named "label",
                # row names matching tip labels, or the first column matching tip labels
                if ("label" %in% colnames(df)) {
                        # Already has label column
                } else if (all(rownames(df) %in% self$tip_labels)) {
                        df$label <- rownames(df)
                } else if (all(df[[1]] %in% self$tip_labels)) {
                        colnames(df)[1] <- "label"
                } else {
                        stop(
                          "Cannot identify tip labels in ring_data. ",
                          "Please ensure rownames or the first column matches tree tip labels."
                        )
                }

                # Keep only tips present in the tree
                df <- df[df$label %in% self$tip_labels, ]

                # Move label to first column
                df <- df[, c("label", setdiff(colnames(df), "label")), drop = FALSE]

                return(df)
        },
        # Auto-generate color palette based on group levels
        generate_palette = function() {
                groups   <- unique(self$group_info$Group)
                n_groups <- length(groups)

                # Group colors
                if (n_groups <= 9) {
                        group_color <- RColorBrewer::brewer.pal(max(3, n_groups), "Set1")[1:n_groups]
                } else if (n_groups <= 12) {
                        group_color <- RColorBrewer::brewer.pal(max(3, n_groups), "Paired")[1:n_groups]
                } else if (n_groups <= 20) {
                        group_color <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_groups)
                } else {
                        # Use viridis when more than 20 groups
                        group_color <- viridis::viridis(n_groups, option = "D")
                }
                names(group_color) <- groups

                palette <- list(
                group_color = group_color
                )
                return(palette)
        }
        ),
        lock_class = FALSE,
        lock_objects = FALSE
)
