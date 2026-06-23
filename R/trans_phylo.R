#' @title trans_phylo class
#'
#' @description Publication-grade phylogenetic tree visualization built on ggtree package (DOI: 10.1002/imt2.56),
#'   designed for the microtable object. Supports group-based branch/tip
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
                #' @param color_palette vector, optional custom color palette; 
                #'   If the color vector has names, they will be used in correspondence with the names; auto-generated if NULL
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
                #' total_rel_abund <- trans_norm$new(mt)$norm(method = "TSS")$otu_table
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
                #' niche <- trans_niche$new(mt)$cal_niche_breadth()
                #' specificity <- 1 - niche$res_niche_breadth[rownames(otu), "Levins_Standardized"]
                #' 
                #' # network metrics
                #' net <- trans_network$new(mt, cor_method = "pearson")$cal_network()$get_node_table()
                #' centrality <- net$res_node_table[rownames(otu), "closeness_centrality"]
                #' connectivity <- net$res_node_table[rownames(otu), "z"]
                #' 
                #' # Habitat
                #' mt$add_rownames2tax("OTU")
                #' mt$cal_abund()
                #' mt_diff <- trans_diff$new(dataset = mt, method = "rf", taxa_level = "OTU", group = "Group", alpha = 1)
                #' mt_diff_res <- mt_diff$res_diff 
                #' rownames(mt_diff_res) %<>% gsub(".*\\|", "", .)
                #' dominant_habitat <- mt_diff_res[rownames(otu), "Group"]
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
                #'   Centrality   = centrality,
                #'   Connectivity = connectivity,
                #'   Habitat      = dominant_habitat,
                #'   Class        = class_col,
                #'   stringsAsFactors = FALSE
                #' )
                #' 
                #' # Clean Class column
                #' ring_df$Class <- gsub("^c__", "", ring_df$Class)
                #' ring_df$Class[ring_df$Class == "" | is.na(ring_df$Class)] <- "Unclassified"
                #'
                #' pviz <- trans_phylo$new(
                #'   dataset             = mt,
                #'   group_col           = "Phylum",
                #'   group_prefix_remove = "p__",
                #'   ring_data           = ring_df,
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

                        self$params <- list(
                                open_angle   = 30,
                                ring_history = NULL
                        )

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
                                color_palette <- private$generate_palette()
                        }else{
							groups   <- unique(self$group_info$Group)
							n_groups <- length(groups)
							if(!is.vector(color_palette)) stop("Input color_palette must be a vector!")
							if(length(color_palette) < n_groups) stop("Input color_palette is not enough!", n_groups, " is need!")
							if(is.null(names(color_palette))){
								names(color_palette) <- groups
							}
						}
						self$color_palette <- color_palette

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

                        self$params$open_angle <- open_angle

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
                        pal        <- self$color_palette
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
                        invisible(self)
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
                #' @param categorical_colors color vector or palette function for categorical columns, default NULL (auto).
                #'   If a named character vector, names must correspond to factor levels present in
                #'   the data. If an unnamed character vector, it is recycled and named by the
                #'   factor levels in their alphabetical order. If a palette function (e.g.
                #'   `scales::hue_pal()` or `colorRampPalette(c(...))`), it is invoked with the
                #'   number of factor levels to produce a color vector, which is then named by
                #'   the levels.
                #' @param ring_width ring width, default NULL; see ring_width_bar and ring_width_other.
                #' @param ring_offset ring offset, default NULL; If NULL, ring_offset_first + ring_gap
                #' @param ring_gap numeric, gap between consecutive rings when ring_offset is not manually specified
                #' @param ring_offset_first numeric, starting offset for the very first ring when ring_offset is not manually specified
                #' @param ring_width_bar numeric, default ring width for numeric bar-type rings when ring_width is NULL
                #' @param ring_width_other numeric, default ring width for non-bar rings (categorical color band, point, star) when ring_width is NULL
                #' @param show_legend logical, whether to show this ring's legend, default TRUE
                #' @param legend_title_custom custom legend title, default uses column name
                #' @param ring_label label around the ring, default uses column name
                #' @param geom type: "bar" for bar heatmap, "point" for scatter, "star" for stars;
                #'   categorical columns automatically use color band regardless of this setting
                #' @param na_fill fill color for NA values, default "grey90"
                #' @param na_replace categorical value to replace NA/empty with, default "Unknown".
                #'   Only applied when categorical_colors is NULL (auto palette). Set to NULL to keep
                #'   NA values as NA, in which case they are colored with na_fill via
                #'   scale_fill_manual(na.value = na_fill). When categorical_colors is provided,
                #'   NA is always kept as NA and colored with na_fill (regardless of na_replace),
                #'   so the user's color vector does not need to cover an extra "Unknown" level.
                #' @param point_size point size (only for geom = "point"), default 1.5
                #' @param ... additional arguments passed to ggtreeExtra::geom_fruit() (e.g. axis.params, grid.params)
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' # Numeric bar ring
                #' pviz$add_ring(
                #'   col_name            = "LogAbund",
                #'   geom                = "bar",
                #'   ring_width          = 0.03,
                #'   ring_offset         = 0.02,
                #'   color_low           = "#440154",
                #'   color_high          = "#FDE725",
                #'   na_fill             = "grey95",
                #'   legend_title_custom = "Log10\\nAbund."
                #' )
                #'
                #' # Categorical color band ring
                #' pviz$add_ring(
                #'   col_name            = "Habitat",
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
									ring_width          = NULL,
									ring_offset         = NULL,
									ring_gap            = 0.03,
									ring_offset_first   = 0.02,
									ring_width_bar      = 0.03,
									ring_width_other    = 0.03,
									show_legend         = TRUE,
									legend_title_custom = NULL,
									ring_label          = NULL,
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
						if (is.null(df)) stop("No ring data available.")
						if (!col_name %in% colnames(df)) stop(sprintf("Column '%s' not found in ring data.", col_name))

						sub_df <- df[, c("label", col_name), drop = FALSE]
						colnames(sub_df)[2] <- "value"
						is_numeric <- is.numeric(sub_df$value)

						# ==========================================
						# Smart auto-layout: compute width and offset
						# ==========================================
						# 1. Smart ring_width and ring_offset assignment
						if (is.null(ring_width)) {
								ring_width <- if (is_numeric && geom == "bar") ring_width_bar else ring_width_other
						}
						if (is.null(ring_offset)) {
								if (is.null(self$params$ring_history) || nrow(self$params$ring_history) == 0) {
										ring_offset <- ring_offset_first # Default starting offset for the first ring
								} else {
										ring_offset <- ring_gap
								}
						}

						r_label <- if (!is.null(ring_label)) ring_label else if (!is.null(legend_title_custom)) legend_title_custom else col_name
						r_label <- gsub("\n", " ", r_label)
						entry <- data.frame(
							label = r_label, 
							offset = ring_offset, 
							width = ring_width, 
							is_numeric = is_numeric, # add var type for label position identification
							stringsAsFactors = FALSE)
						
						self$params$ring_history <- if (is.null(self$params$ring_history)) entry else rbind(self$params$ring_history, entry)

						if (is_numeric && is.null(limits)) {
								limits <- range(sub_df$value, na.rm = TRUE)
						}

						self$plot_obj <- self$plot_obj + ggnewscale::new_scale("fill")
						dot_args <- list(...)
						
						default_axis <- list(add = FALSE)
						dot_args$axis.params <- if (!is.null(dot_args$axis.params)) modifyList(default_axis, dot_args$axis.params) else default_axis
						
						default_grid <- list()
						dot_args$grid.params <- if (!is.null(dot_args$grid.params)) modifyList(default_grid, dot_args$grid.params) else default_grid

                        if (is_numeric) {
                                # Numeric -> heatmap bar / point / star
                                # Build fill scale once (DRY) and reuse for all numeric geoms.
                                legend_name  <- ifelse(is.null(legend_title_custom), col_name, legend_title_custom)
                                legend_guide <- if (show_legend) ggplot2::guide_colorbar() else ggplot2::guide_none()
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
                                        self$plot_obj <- self$plot_obj +
                                                fruit_layer +
                                                fill_scale
                                }
                        } else {
                                # Categorical -> color band
                                # Sanitize empty strings: always convert to NA first (safe even when
                                # na_replace is NULL — avoids deleting rows by assigning NULL to a subset).
                                sub_df$value[sub_df$value == ""] <- NA_character_

                                # NA handling strategy:
                                #   - When the user provides categorical_colors, keep NA as NA so that
                                #     scale_fill_manual(na.value = na_fill) paints it with na_fill.
                                #     This avoids creating an extra "Unknown" level that the user's
                                #     color vector does not cover, which previously caused:
                                #     "names attribute [4] must be the same length as the vector [3]".
                                #   - When colors are auto-generated, honor na_replace: if non-NULL,
                                #     NA becomes its own labeled level (default "Unknown"); if NULL,
                                #     NA stays as NA and is colored with na_fill.
                                if (is.null(categorical_colors) && !is.null(na_replace)) {
                                        sub_df$value[is.na(sub_df$value)] <- na_replace
                                }
                                sub_df$value <- as.factor(sub_df$value)
                                n_levels     <- length(levels(sub_df$value))

                                # Normalize categorical_colors input: accept a palette FUNCTION
                                # (e.g. scales::hue_pal() or colorRampPalette(c(...))) and resolve
                                # it to a named color vector before the coverage logic below.
                                # This makes add_rings_batch(categorical_palette = ...) work, since
                                # batch forwards the palette through this same parameter.
                                if (is.function(categorical_colors)) {
                                        categorical_colors <- categorical_colors(n_levels)
                                }

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
                                        # Ensure names match levels
                                        names(categorical_colors) <- levels(sub_df$value)
                                } else {
                                        # User-supplied palette (named vector, unnamed vector, or
                                        # the resolved output of a palette function). If unnamed,
                                        # assign names by factor-level order so the coverage check
                                        # below works correctly.
                                        if (is.null(names(categorical_colors))) {
                                                names(categorical_colors) <- levels(sub_df$value)[seq_along(categorical_colors)]
                                        }
                                        # User-supplied palette: validate / repair coverage for the
                                        # levels actually present in the data. NA is excluded because
                                        # it is not a factor level — it is handled by na.value downstream.
                                        data_levels <- levels(sub_df$value)
                                        # Drop user entries that don't appear in the data (avoids
                                        # spurious legend keys), preserving data level order.
                                        categorical_colors <- categorical_colors[
                                                intersect(names(categorical_colors), data_levels)
                                        ]
                                        missing_levels <- setdiff(data_levels, names(categorical_colors))
                                        if (length(missing_levels) > 0) {
                                                fallback <- if (length(missing_levels) <= 9) {
                                                        RColorBrewer::brewer.pal(max(3, length(missing_levels)), "Set1")[1:length(missing_levels)]
                                                } else if (length(missing_levels) <= 12) {
                                                        RColorBrewer::brewer.pal(max(3, length(missing_levels)), "Set3")[1:length(missing_levels)]
                                                } else {
                                                        colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(missing_levels))
                                                }
                                                categorical_colors <- c(categorical_colors, setNames(fallback, missing_levels))
                                                warning(sprintf(
                                                        "[trans_phylo$add_ring] categorical_colors missing for level(s): %s. Auto-assigned fallback colors.",
                                                        paste(missing_levels, collapse = ", ")
                                                ))
                                        }
                                }

                                # Forcefully inject a constant virtual X value to reduce the heatmap to a bar chart with values uniformly fixed at 1
                                sub_df$.dummy_x <- 1
                                
                                fruit_layer <- rlang::inject(ggtreeExtra::geom_fruit(
                                        data        = sub_df,
                                        geom        = geom_col,   # Deprecate geom_tile
                                        mapping     = ggplot2::aes(y = label, x = .dummy_x, fill = value),
                                        offset      = ring_offset,
                                        pwidth      = ring_width,
                                        orientation = "y",        # Ensure distribution along the y-axis
                                        width       = 1,          # width=1 ensures no white gaps between adjacent categorical tiles, forming a continuous color band
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

						message(sprintf("Ring annotation '%s' added. (offset: %.3f, width: %.3f)", 
										col_name, ring_offset, ring_width))
                        invisible(self)
				},

                # =======================================================================
                # add_ring_labels: Auto-label all rings at the fan opening edge
                # =======================================================================
                #' @description Automatically place text labels for all added rings at the edge of the fan opening.
                #'   Uses the ring geometry history recorded by add_ring() to compute positions.
                #'   Call AFTER all add_ring() calls and AFTER add_spacer().
                #' @param angle_offset numeric, fine-tuning angle padding inside the opening gap, default 2 (degrees)
                #' @param size text font size, default 1.8
                #' @param color text color, default "grey25"
                #' @param fontface text font style, default "italic"
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' pviz$add_ring_labels(size = 1.8, color = "grey25")
                #' }
                add_ring_labels = function(size = 1.8, color = "grey25", angle_offset = 15) {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        history <- self$params$ring_history
                        if (is.null(history) || nrow(history) == 0) return(invisible(self))

                        #
                        tree_data <- tidytree::as_tibble(self$plot_obj$data)
                        min_x <- min(tree_data$x, na.rm = TRUE)
                        max_x <- max(tree_data$x, na.rm = TRUE)
                        tree_width <- max_x - min_x
                        if (is.na(tree_width) || tree_width <= 0) tree_width <- 1.0

                        #outer_edges <- cumsum(history$offset + history$width)
                        # Calculate the relative proportional positions of each ring's center point
                        #mid_points  <- outer_edges - (history$width / 2)

						mid_points <- numeric(nrow(history))
                        current_max_x_prop <- 0
                        
                        # All layers are now uniformly stacked according to their physical left boundary
                        for (i in seq_len(nrow(history))) {
                                mid_points[i] <- current_max_x_prop + history$offset[i] + history$width[i] / 2
                                current_max_x_prop <- current_max_x_prop + history$offset[i] + history$width[i]
                        }

						# Compute absolute X coordinates: precisely align with the physical centers of the rings
						x_pos_all <- max_x + mid_points * tree_width

                        # Y coordinates: 0.5 is the valid lower boundary for the ggtree fan layout
                        y_pos_all <- rep(0.5, nrow(history))

                        # Expand the plotting area boundaries to prevent long labels on the outermost rings from being clipped
                        self$plot_obj <- self$plot_obj + ggplot2::expand_limits(x = max(x_pos_all) * 1.1)

                        OPEN_ANGLE <- if (!is.null(self$params$open_angle)) self$params$open_angle else 30
                        edge_angle_deg <- 270 + OPEN_ANGLE / 2 + angle_offset
                        tangent_angle  <- (edge_angle_deg - 90) %% 360
                        
                        # Ensure text is always readable from left to right (auto-flip orientation)
                        text_angle <- if (tangent_angle > 90 && tangent_angle < 270) tangent_angle - 180 else tangent_angle

                        for (i in seq_len(nrow(history))) {
                                self$plot_obj <- self$plot_obj +
                                        ggplot2::annotate(
                                                "text",
                                                x     = x_pos_all[i],
                                                y     = y_pos_all[i], 
                                                label = history$label[i],
                                                angle = text_angle,
                                                size  = size,
                                                # hjust = 1 Ensure that the end (tail) of the text is anchored at the computed X coordinate and radiates outward
                                                hjust = 1,  
                                                vjust = 0.5, 
                                                color = color
                                        )
                        }

                        message(sprintf("Successfully auto-labeled %d rings with dynamic cumulative layout.", nrow(history)))
                        invisible(self)
                },

                # =======================================================================
                # add_rings_batch: Batch add multiple ring annotations
                # =======================================================================
                #' @description Batch add all (or selected) columns from ring_data as outer rings.
                #'   Parameters ring_width, ring_offset, color_low, color_high, numeric_geom, show_legend
                #'   and categorical_palette all support vector input (length = number of col_names);
                #'   a single value is automatically recycled for all rings.
                #' @param col_names character vector of column names to add; NULL means all columns
                #' @param numeric_geom character or character vector; geom type(s) for numeric columns, default "bar".
                #'   Recycled if length 1; must match length of col_names otherwise.
                #' @param ring_width numeric or numeric vector; ring width(s), default 0.08.
                #'   Single value is recycled; vector must match length of col_names.
                #' @param ring_offset numeric or numeric vector; gap(s) before each ring, default 0.02.
                #'   Single value is recycled; vector must match length of col_names.
                #' @param color_low character or character vector; color(s) for low numeric values, default "#0D0887".
                #'   Single value is recycled; vector must match length of col_names.
                #' @param color_high character or character vector; color(s) for high numeric values, default "#F0F921".
                #'   Single value is recycled; vector must match length of col_names.
                #' @param categorical_palette palette function, named color vector, or list thereof, for
                #'   categorical columns. A single value is recycled for all rings; an unnamed list is
                #'   recycled positionally (length recycled to n); a NAMED list is matched to columns by
                #'   name (e.g. list(Habitat = c(...), Salinity = c(...))), and columns without a matching
                #'   entry fall back to auto-generated colors. Each entry may be either a palette function
                #'   (e.g. `scales::hue_pal()` or `colorRampPalette(c(...))`) or a named color vector.
                #'   Forwarded to each add_ring() call as `categorical_colors`. Default NULL.
                #' @param show_legend logical or logical vector; whether to show legends, default TRUE.
                #'   Single value is recycled; vector must match length of col_names.
                #' @param ... additional arguments passed to each add_ring() call
                #' @return updated ggtree plot object
                #' @examples
                #' \dontrun{
                #' pviz$add_rings_batch(
                #'   col_names    = c("RelAbund", "Prevalence", "Habitat"),
                #'   numeric_geom = "bar",
                #'   ring_width   = c(0.02, 0.03, 0.04),
                #'   ring_offset  = c(0.02, 0.01, 0.02),
                #'   color_low    = c("#440154", "#EDF8FB"),
                #'   color_high   = c("#FDE725", "#006D2C"),
                #'   categorical_palette = list(Habitat = c("IW" = "red", "CW" = "blue", "TW" = "green")),
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

                        n <- length(all_cols)
                        ring_width   <- rep_len(ring_width, n)
                        ring_offset  <- rep_len(ring_offset, n)
                        color_low    <- rep_len(color_low, n)
                        color_high   <- rep_len(color_high, n)
                        numeric_geom <- rep_len(numeric_geom, n)
                        show_legend  <- rep_len(show_legend, n)
                        cat_pal_list <- if (is.null(categorical_palette)) {
									rep_len(list(NULL), n)
							} else if (is.list(categorical_palette) && !is.null(names(categorical_palette))) {
									lapply(all_cols, function(cn) categorical_palette[[cn, exact = TRUE]])
							} else if (is.list(categorical_palette)) {
									rep_len(categorical_palette, n)
							} else {
									rep_len(list(categorical_palette), n)
							}

                        for (i in seq_along(all_cols)) {
                                cn <- all_cols[i]
                                self$add_ring(
                                        col_name            = cn,
                                        ring_data           = df,
                                        ring_width          = ring_width[i],
                                        ring_offset         = ring_offset[i],
                                        geom                = numeric_geom[i],
                                        color_low           = color_low[i],
                                        color_high          = color_high[i],
                                        categorical_colors  = cat_pal_list[[i]],
                                        show_legend         = show_legend[i],
                                        legend_title_custom = cn,
                                        ...
                                )
                        }

                        message(sprintf("Batch added %d ring annotation(s).", length(all_cols)))
                        invisible(self)
                },

                # =======================================================================
                # add_scale_bar: Add evolutionary distance scale bar
                # =======================================================================
                #' @description Add an evolutionary distance scale bar to the plot
                #' @param x numeric, x coordinate
                #' @param y numeric, y coordinate
                #' @param label string, scale bar label
                #' @param fontsize numeric, default 2, text size
                #' @param linesize numeric, default 0.5, line size
				#' @param offset numeric, default 0.3, offset of text to line
                #' @param ... additional arguments passed to ggtree::geom_treescale
                #' @return updated plot
                #' @examples
                #' \dontrun{
                #' pviz$add_scale_bar(x = 0.2, y = 0, label = "0.05")
                #' }
                add_scale_bar = function(x = NULL, y = NULL, label = "0.05", fontsize = 2, linesize = 0.5, offset = 0.3, ...) {
                        if (is.null(self$plot_obj)) stop("No plot object found.")

                        self$plot_obj <- self$plot_obj +
                        ggtree::geom_treescale(
                                x        = x,
                                y        = y,
                                width    = as.numeric(label),
                                label    = label,
                                fontsize = fontsize,
                                linesize = linesize,
								offset   = offset,
                                ...
                        )
                        invisible(self)
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
                                extend   = extend
                        )
                        invisible(self)
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
                        invisible(self)
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

                        invisible(self)
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
                        # Use viridis when more than 20
                        group_color <- viridis::viridis(n_groups, option = "D")
                }
                names(group_color) <- groups

                return(group_color)
        }
        ),
        lock_class = FALSE,
        lock_objects = FALSE
)
