#' @title
#' Multi-omics data integration analysis.
#' 
#' @description
#' This class is a wrapper for microbiome multi-omics data integration analysis based on the \code{mixOmics} package
#' <doi:10.1371/journal.pcbi.1005752>. It integrates multiple \code{\link{microtable}} objects (e.g. microbiome and metabolome data)
#' and provides methods for data preprocessing, supervised/unsupervised multi-block integration (DIABLO, sPLS),
#' parameter tuning, model evaluation, feature extraction and visualization.
#'
#' The core integration models are powered by \code{mixOmics} functions such as
#' \code{block.splsda} (DIABLO), \code{block.spls} and \code{spls}.
#'
#' @export
trans_multiomics <- R6Class(classname = "trans_multiomics",
	public = list(
		#' @description
		#' Create the \code{trans_multiomics} object.
		#'
		#' The function receives multiple \code{\link{microtable}} objects and automatically aligns samples
		#' (intersection of sample names across all microtable objects) to ensure the matrices are compatible
		#' for multi-omics integration.
		#' @param microtables default NULL; a named \code{list} of \code{\link{microtable}} objects (e.g. \code{list(microb = soil_microb, metab = soil_metab)}),
		#' 	 or a single \code{microtable} object. The names of the list are used as the omics block names.
		#' @param group default NULL; a factor or character vector indicating the group membership of each sample.
		#' 	 If NULL, the function extracts the group from \code{group_col} column of \code{sample_table} in the first microtable.
		#' @param group_col default "Group"; the column name in \code{sample_table} used to extract the group factor when \code{group} is NULL.
		#' @return \code{self} (invisible), enabling chained calls; \code{dataset_list}, \code{data_list_names} and \code{Y} are stored in the object.
		#' @examples
		#' \dontrun{
		#' data(soil_microb)
		#' data(soil_metab)
		#' t1 <- trans_multiomics$new(
		#'   microtables = list(microb = soil_microb, metab = soil_metab),
		#'   group_col = "Group"
		#' )
		#' }
		initialize = function(
			microtables = NULL,
			group = NULL,
			group_col = "Group"
			){
			if(is.null(microtables)){
				stop("Please provide microtables parameter: a list of microtable objects or a single microtable object!")
			}
			if(inherits(microtables, "microtable")){
				microtables <- list(microtables)
				names(microtables) <- "block1"
			}
			if(!inherits(microtables, "list")){
				stop("microtables must be a list of microtable objects or a single microtable object!")
			}
			lapply(microtables, function(x){
				if(!inherits(x, "microtable")){
					stop("Each element in microtables must be a microtable object!")
				}
			})
			if(is.null(names(microtables))){
				names(microtables) <- paste0("block", seq_along(microtables))
			}else{
				# fill the empty/NA names of a partially named list (e.g. list(microb = x, y))
				empty_names <- is.na(names(microtables)) | names(microtables) == ""
				if(any(empty_names)){
					names(microtables)[empty_names] <- paste0("block", which(empty_names))
				}
			}
			block_names <- names(microtables)
			# tidy each microtable (clone to protect user's original data from R6 reference semantics)
			microtables <- lapply(microtables, function(x){
				new_x <- clone(x)
				new_x$tidy_dataset()
				new_x
			})
			# align samples: intersection across all microtable objects
			sample_lists <- lapply(microtables, function(x) rownames(x$sample_table))
			common_samples <- Reduce(intersect, sample_lists)
			if(length(common_samples) == 0){
				stop("No common samples found across all microtable objects! Please check the sample names!")
			}
			removed_total <- sum(lengths(sample_lists)) - length(common_samples) * length(microtables)
			if(removed_total > 0){
				message(removed_total, " sample(s) removed to align samples across all blocks (intersection) ...")
			}
			# filter each microtable to common samples
			microtables <- lapply(microtables, function(x){
				if(nrow(x$sample_table) > length(common_samples)){
					x$sample_table <- x$sample_table[common_samples, , drop = FALSE]
					x$tidy_dataset()
				}
				x
			})
			# sort common samples consistently
			common_samples <- sort(common_samples)
			microtables <- lapply(microtables, function(x){
				x$sample_table <- x$sample_table[common_samples, , drop = FALSE]
				x$tidy_dataset()
				x
			})
			# extract Y
			if(is.null(group)){
				if(!is.null(group_col)){
					microeco:::check_table_variable(microtables[[1]]$sample_table, group_col, "group_col", "sample_table")
					group <- microtables[[1]]$sample_table[, group_col]
					# the sample_table has been sorted to common_samples above; attach names directly
					names(group) <- rownames(microtables[[1]]$sample_table)
				}
			}
			if(!is.null(group)){
				if(!is.factor(group)){
					group <- as.factor(group)
				}
				# drop the empty levels possibly left by the sample intersection (or by a named group with extra levels)
				group <- droplevels(group)
				if(is.null(names(group))){
					if(length(group) != length(common_samples)){
						stop("The length of group (", length(group), ") does not match the number of common samples (", length(common_samples), ")!")
					}
					warning("The group vector has no names; samples are matched by order. Please provide a named vector to ensure correct sample alignment!")
					names(group) <- common_samples
				}else{
					if(!setequal(names(group), common_samples)){
						stop("The names of the group vector do not match the common samples! Please check the sample names of the group parameter!")
					}
					# reorder group to match the sorted common samples
					group <- group[common_samples]
				}
				if(nlevels(group) < 2){
					stop("The group factor must have at least 2 levels for supervised analysis!")
				}
			}
			self$dataset_list <- microtables
			self$data_list_names <- block_names
			self$Y <- group
			message("Multi-omics data loaded with ", length(common_samples), " samples and ", length(microtables), " blocks: ", paste(block_names, collapse = ", "), " ...")
			invisible(self)
		},

		#' @description
		#' Preprocess the data: filter low-abundance features and apply normalization/transformation.
		#'
		#' This function reuses the \code{\link{trans_norm}} class to apply normalization methods (e.g. CLR)
		#' to each omics block. The processed matrices (rows are samples, columns are features) are stored
		#' in \code{data_list} for \code{mixOmics} input.
		#' @param method default "clr"; the normalization method passed to \code{trans_norm$norm}.
		#' 	 Available options include \code{"clr"}, \code{"rclr"}, \code{"css"}, \code{"tss"}, \code{"log"}, etc.
		#' 	 See the \code{method} parameter of \code{\link{trans_norm}} class for details.
		#' @param pseudocount default 1; add pseudocount for zero values when \code{method = "clr"}.
		#' @param filter_thres default NULL; the relative abundance threshold for filtering low-abundance features.
		#' 	 Passed to \code{rel_abund} parameter of \code{microtable$filter_taxa}. NULL means no filtering.
		#' @param filter_freq default 1; the occurrence frequency threshold for filtering features.
		#' 	 Passed to \code{freq} parameter of \code{microtable$filter_taxa}.
		#' @param ... other parameters passed to \code{norm} method of \code{\link{trans_norm}} class.
		#' @return \code{self} (invisible), enabling chained calls; \code{data_list} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$preprocess_data(method = "clr", filter_thres = 0.0001, filter_freq = 0.2)
		#' }
		preprocess_data = function(
			method = "clr",
			pseudocount = 1,
			filter_thres = NULL,
			filter_freq = 1,
			...
			){
			if(is.null(self$dataset_list)){
				stop("The dataset_list is NULL! Please create the object with microtables first!")
			}
			data_list <- list()
			for(i in seq_along(self$dataset_list)){
				block_name <- self$data_list_names[i]
				tmp_dataset <- clone(self$dataset_list[[i]])
				# filter low-abundance features
				if(!is.null(filter_thres) || filter_freq != 1){
					if(is.null(filter_thres)){
						filter_thres_use <- 0
					}else{
						filter_thres_use <- filter_thres
					}
					tmp_dataset$filter_taxa(rel_abund = filter_thres_use, freq = filter_freq)
				}
				# reuse trans_norm for normalization
				tmp_norm <- trans_norm$new(dataset = tmp_dataset)
				res_dataset <- tmp_norm$norm(method = method, pseudocount = pseudocount, ...)
				# extract transposed otu_table (rows = samples, cols = features)
				if(inherits(res_dataset, "microtable")){
					abund_matrix <- as.matrix(t(res_dataset$otu_table))
				}else{
					abund_matrix <- as.matrix(res_dataset)
				}
				# align rows to Y names if available
				if(!is.null(self$Y)){
					abund_matrix <- abund_matrix[rownames(abund_matrix) %in% names(self$Y), , drop = FALSE]
					abund_matrix <- abund_matrix[names(self$Y), , drop = FALSE]
				}
				data_list[[block_name]] <- abund_matrix
				message("Block '", block_name, "': ", nrow(abund_matrix), " samples x ", ncol(abund_matrix), " features after preprocessing ...")
			}
			self$data_list <- data_list
			message("Preprocessed data is stored in object$data_list ...")
			invisible(self)
		},

		#' @description
		#' Set the design matrix for DIABLO analysis.
		#'
		#' The design matrix defines the weights of correlation between different omics blocks.
		#' A full design (weight = 1) maximizes the covariance between blocks, while a small weight
		#' (e.g. 0.1) favors classification accuracy.
		#' @param design default NULL; a user-specified design matrix. If NULL, a matrix filled with
		#' 	 \code{weight} (with 0 on the diagonal) is created.
		#' @param weight default 0.1; the off-diagonal weight when \code{design} is NULL.
		#' @return \code{self} (invisible), enabling chained calls; \code{design} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$set_design(weight = 0.1)
		#' }
		set_design = function(
			design = NULL,
			weight = 0.1
			){
			if(is.null(self$data_list)){
				stop("The data_list is NULL! Please run preprocess_data() first!")
			}
			n_block <- length(self$data_list)
			if(is.null(design)){
				design <- matrix(weight, ncol = n_block, nrow = n_block,
					dimnames = list(self$data_list_names, self$data_list_names))
				diag(design) <- 0
			}else{
				if(!all(dim(design) == c(n_block, n_block))){
					stop("The design matrix must be ", n_block, " x ", n_block, "!")
				}
				if(is.null(rownames(design))){
					rownames(design) <- self$data_list_names
				}
				if(is.null(colnames(design))){
					colnames(design) <- self$data_list_names
				}
			}
			self$design <- design
			message("Design matrix is stored in object$design ...")
			invisible(self)
		},

		#' @description
		#' Run supervised multi-omics integration using DIABLO (block.sPLS-DA).
		#'
		#' This function wraps \code{mixOmics::block.splsda} to perform supervised multi-block sparse PLS-DA.
		#' It is suitable when the user has a clear grouping variable and wants to identify cross-omics
		#' features that discriminate between groups.
		#' @param ncomp default 2; the number of components.
		#' @param keepX default NULL; a list of variable numbers to select per block per component.
		#' 	 If NULL, all features are retained (full model). Example: \code{list(microb = c(10, 10), metab = c(15, 15))}.
		#' @param design default NULL; the design matrix. If NULL, use \code{self$design}; if still NULL, create with \code{set_design()}.
		#' @param scale default TRUE; whether to scale the data.
		#' @param ... other parameters passed to \code{mixOmics::block.splsda}.
		#' @return \code{self} (invisible), enabling chained calls; \code{model} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$run_diablo(ncomp = 2, keepX = list(microb = c(10, 10), metab = c(15, 15)))
		#' }
		run_diablo = function(
			ncomp = 2,
			keepX = NULL,
			design = NULL,
			scale = TRUE,
			...
			){
			private$check_data_list()
			if(is.null(self$Y)){
				stop("The Y (group) is NULL! Please provide group when creating the object!")
			}
			if(is.null(design)){
				if(is.null(self$design)){
					self$set_design()
				}
				design <- self$design
			}
			if(is.null(keepX)){
				keepX <- private$build_keepX_default(ncomp = ncomp)
				message("keepX is NULL; using all features (full model). For sparse model, please provide keepX or run tune_model() ...")
			}
			message("Running DIABLO (block.splsda) with ncomp = ", ncomp, " ...")
			self$model <- mixOmics::block.splsda(
				X = self$data_list,
				Y = self$Y,
				ncomp = ncomp,
				keepX = keepX,
				design = design,
				scale = scale,
				...
			)
			message("DIABLO model is stored in object$model ...")
			invisible(self)
		},

		#' @description
		#' Run unsupervised multi-omics integration using block.sPLS.
		#'
		#' This function wraps \code{mixOmics::block.spls} to perform unsupervised multi-block sparse PLS.
		#' It is suitable when there is no grouping variable and the user wants to explore the covariance
		#' structure across multiple omics data sets.
		#' Note that \code{mixOmics::block.spls} requires a numeric response matrix, so one omics block
		#' must be designated as the response block through the \code{Y_block} parameter (passed to \code{indY}).
		#' @param ncomp default 2; the number of components.
		#' @param Y_block default NULL; the index or name of the omics block used as the response (\code{indY}).
		#' 	 If NULL, the last block in \code{data_list} is used.
		#' @param keepX default NULL; a list of variable numbers to select per block per component (including the response block,
		#' 	 as \code{mixOmics} ignores keepY when \code{indY} is used). If NULL, all features are retained (full model).
		#' @param design default NULL; the design matrix. If NULL, use \code{self$design}; if still NULL, create with \code{set_design()}.
		#' @param mode default "regression"; the PLS mode. One of "regression", "canonical", "invariant", "classic".
		#' @param scale default TRUE; whether to scale the data.
		#' @param ... other parameters passed to \code{mixOmics::block.spls}.
		#' @return \code{self} (invisible), enabling chained calls; \code{model} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$run_spls(ncomp = 2, Y_block = "metab", mode = "regression")
		#' }
		run_spls = function(
			ncomp = 2,
			Y_block = NULL,
			keepX = NULL,
			design = NULL,
			mode = "regression",
			scale = TRUE,
			...
			){
			private$check_data_list()
			n_block <- length(self$data_list)
			if(n_block < 2){
				stop("run_spls requires at least 2 omics blocks! For two-omics analysis, run_spls_single() is recommended!")
			}
			# mixOmics::block.spls requires a numeric response matrix; designate one omics block as the response via indY
			if(is.null(Y_block)){
				Y_block <- n_block
				message("Y_block is NULL; using the last block '", self$data_list_names[Y_block], "' as the response block (indY) ...")
			}
			if(is.character(Y_block)){
				Y_block_name <- Y_block
				Y_block <- match(Y_block, self$data_list_names)
				if(is.na(Y_block)){
					stop("Y_block '", Y_block_name, "' not found in data_list_names: ", paste(self$data_list_names, collapse = ", "), "!")
				}
			}
			if(!is.numeric(Y_block) || length(Y_block) != 1 || Y_block < 1 || Y_block > n_block){
				stop("Y_block must be a single integer between 1 and ", n_block, " or a valid block name!")
			}
			if(is.null(design)){
				if(is.null(self$design)){
					self$set_design()
				}
				design <- self$design
			}
			if(is.null(keepX)){
				keepX <- private$build_keepX_default(ncomp = ncomp)
				message("keepX is NULL; using all features (full model). For sparse model, please provide keepX ...")
			}
			message("Running block.spls (unsupervised) with response block '", self$data_list_names[Y_block], "', ncomp = ", ncomp, " ...")
			self$model <- mixOmics::block.spls(
				X = self$data_list,
				indY = Y_block,
				ncomp = ncomp,
				keepX = keepX,
				design = design,
				mode = mode,
				scale = scale,
				...
			)
			message("block.spls model is stored in object$model ...")
			invisible(self)
		},

		#' @description
		#' Run two-omics sPLS (sparse PLS) regression for unsupervised association analysis.
		#'
		#' This function wraps \code{mixOmics::spls} to perform sparse PLS between two omics blocks.
		#' It is useful for exploring the covariation between two data sets (e.g. microbiome and metabolome).
		#' @param X_block default 1; the index or name of the X block in \code{data_list}.
		#' @param Y_block default 2; the index or name of the Y block in \code{data_list}.
		#' @param ncomp default 2; the number of components.
		#' @param mode default "regression"; the PLS mode.
		#' @param keepX default NULL; the number of variables to select in X per component.
		#' @param keepY default NULL; the number of variables to select in Y per component.
		#' @param scale default TRUE; whether to scale the data.
		#' @param ... other parameters passed to \code{mixOmics::spls}.
		#' @return \code{self} (invisible), enabling chained calls; \code{model} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$run_spls_single(X_block = "microb", Y_block = "metab", ncomp = 2)
		#' }
		run_spls_single = function(
			X_block = 1,
			Y_block = 2,
			ncomp = 2,
			mode = "regression",
			keepX = NULL,
			keepY = NULL,
			scale = TRUE,
			...
			){
			private$check_data_list()
			X_block <- private$resolve_single_block(X_block, "X_block")
			Y_block <- private$resolve_single_block(Y_block, "Y_block")
			if(X_block == Y_block){
				stop("X_block and Y_block must be different blocks!")
			}
			X_data <- self$data_list[[X_block]]
			Y_data <- self$data_list[[Y_block]]
			message("Running spls (two-omics) with X block '", X_block,
				"' and Y block '", Y_block, "', ncomp = ", ncomp, " ...")
			self$model <- mixOmics::spls(
				X = X_data,
				Y = Y_data,
				ncomp = ncomp,
				mode = mode,
				keepX = keepX,
				keepY = keepY,
				scale = scale,
				...
			)
			message("spls model is stored in object$model ...")
			invisible(self)
		},

		#' @description
		#' Tune the DIABLO model to select the optimal number of variables to select (keepX).
		#'
		#' This function wraps \code{mixOmics::tune.block.splsda} (or \code{tune.block.spls}) to perform
		#' cross-validation and identify the optimal \code{keepX} per block per component.
		#' The optimal parameters are stored in \code{best_keepX} for use in the final model.
		#' @param method default "block.splsda"; the tuning method. Currently only "block.splsda" is supported,
		#' 	 as \code{mixOmics} provides no tune function for unsupervised block.spls.
		#' @param ncomp default 2; the number of components to tune.
		#' @param test_keepX default NULL; a list of candidate keepX values per block.
		#' 	 If NULL, an adaptive grid (5\%, 10\% and 20\% of the feature number of each block) is used.
		#' @param folds default 10; the number of folds for cross-validation.
		#' 	 Must not exceed the smallest class size of \code{Y}.
		#' @param nrepeat default 1; the number of repeats. Increase to 10-50 for robust results.
		#' @param dist default "centroids.dist"; the prediction distance for DA methods.
		#' @param seed default NULL; an integer used to set the random seed before cross-validation for reproducible results.
		#' @param ... other parameters passed to \code{mixOmics::tune.block.splsda}.
		#' @return \code{self} (invisible), enabling chained calls; \code{tune_result} and \code{best_keepX} are stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$tune_model(ncomp = 2, test_keepX = list(microb = c(5, 10, 20), metab = c(5, 10, 20)), nrepeat = 1, seed = 123)
		#' }
		tune_model = function(
			method = "block.splsda",
			ncomp = 2,
			test_keepX = NULL,
			folds = 10,
			nrepeat = 1,
			dist = "centroids.dist",
			seed = NULL,
			...
			){
			private$check_data_list()
			if(!is.null(self$Y)){
				min_class <- min(table(self$Y))
				if(folds > min_class){
					stop("folds (", folds, ") exceeds the smallest class size (", min_class, ") in Y! Please reduce folds!")
				}
			}
			if(is.null(self$design)){
				self$set_design()
			}
			if(is.null(test_keepX)){
				test_keepX <- lapply(self$data_list, function(x){
					cand <- unique(round(ncol(x) * c(0.05, 0.1, 0.2)))
					cand <- cand[cand >= 1 & cand <= ncol(x)]
					if(length(cand) == 0){
						cand <- min(5, ncol(x))
					}
					cand
				})
				message("Using adaptive test_keepX grid: ", paste(unlist(lapply(test_keepX, function(x) paste(x, collapse = ","))), collapse = " | "), " ...")
			}
			if(method != "block.splsda"){
				stop("method must be 'block.splsda'! mixOmics does not provide a tune function for unsupervised block.spls!")
			}
			if(is.null(self$Y)){
				stop("The Y (group) is NULL! tune_model requires a grouping variable for block.splsda tuning!")
			}
			message("Tuning ", method, " with folds = ", folds, ", nrepeat = ", nrepeat, " ...")
			# mixOmics (>= 6.3.2) provides its own seed argument and resets the RNG internally when it is missing;
			# pass seed explicitly when supported, otherwise fall back to set.seed
			seed_supported <- "seed" %in% names(formals(mixOmics::tune.block.splsda))
			if(!is.null(seed) && !seed_supported){
				set.seed(seed)
			}
			tune_args <- list(
				X = self$data_list,
				Y = self$Y,
				ncomp = ncomp,
				test.keepX = test_keepX,
				design = self$design,
				validation = "Mfold",
				folds = folds,
				nrepeat = nrepeat,
				dist = dist,
				...
			)
			if(!is.null(seed) && seed_supported){
				tune_args$seed <- seed
			}
			tune_res <- do.call(mixOmics::tune.block.splsda, tune_args)
			self$tune_result <- tune_res
			if(!is.null(tune_res$choice.keepX)){
				self$best_keepX <- tune_res$choice.keepX
				message("Optimal keepX selected: ", paste(unlist(lapply(tune_res$choice.keepX, function(x) paste(x, collapse = ","))), collapse = " | "), " ...")
			}
			message("Tuning result is stored in object$tune_result, optimal keepX in object$best_keepX ...")
			invisible(self)
		},

		#' @description
		#' Evaluate the performance of the model using cross-validation.
		#'
		#' This function wraps \code{mixOmics::perf} to assess the classification error rate (BER) or
		#' prediction performance of the model.
		#' @param method default "Mfold"; the cross-validation method. "Mfold" or "loo".
		#' @param folds default 10; the number of folds (when method = "Mfold").
		#' 	 Must not exceed the smallest class size of \code{Y}.
		#' @param nrepeat default 1; the number of repeats.
		#' @param seed default NULL; an integer used to set the random seed before cross-validation for reproducible results.
		#' @param ... other parameters passed to \code{mixOmics::perf}.
		#' @return \code{self} (invisible), enabling chained calls; \code{perf_result} is stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$eval_performance(folds = 5, nrepeat = 1, seed = 123)
		#' }
		eval_performance = function(
			method = "Mfold",
			folds = 10,
			nrepeat = 1,
			seed = NULL,
			...
			){
			private$check_model()
			if(method == "Mfold" && !is.null(self$Y)){
				min_class <- min(table(self$Y))
				if(folds > min_class){
					stop("folds (", folds, ") exceeds the smallest class size (", min_class, ") in Y! Please reduce folds!")
				}
			}
			message("Evaluating model performance with ", method, ", folds = ", folds, ", nrepeat = ", nrepeat, " ...")
			# pass seed to the dispatched perf method when supported (mixOmics >= 6.3.2); otherwise fall back to set.seed
			seed_supported <- FALSE
			if(!is.null(seed)){
				for(cl in class(self$model)){
					f <- tryCatch(getFromNamespace(paste0("perf.", cl), "mixOmics"), error = function(e) NULL)
					if(!is.null(f) && "seed" %in% names(formals(f))){
						seed_supported <- TRUE
						break
					}
				}
				if(!seed_supported){
					set.seed(seed)
				}
			}
			perf_args <- list(object = self$model, validation = method, folds = folds, nrepeat = nrepeat, ...)
			if(!is.null(seed) && seed_supported){
				perf_args$seed <- seed
			}
			self$perf_result <- do.call(mixOmics::perf, perf_args)
			message("Performance result is stored in object$perf_result ...")
			message("Tip: run get_stability() to extract the feature selection stability from the cross-validation results ...")
			invisible(self)
		},

		#' @description
		#' Extract the selected features (variables with non-zero loadings) from the model.
		#'
		#' This function wraps \code{mixOmics::selectVar} to extract the selected variables on a given
		#' component, and returns a tidy data.frame with feature names, block, component and loading values.
		#' @param comp default 1; the component of interest.
		#' @param block default NULL; the block(s) to extract. NULL means all blocks. Can be integer index or block name.
		#' @return \code{res_features} stored in the object and returned.
		#' @examples
		#' \dontrun{
		#' features <- t1$get_features(comp = 1)
		#' }
		get_features = function(
			comp = 1,
			block = NULL
			){
			private$check_model()
			private$check_comp(comp)
			sel <- mixOmics::selectVar(object = self$model, comp = comp, block = block)
			# organize into a tidy data.frame
			res_list <- list()
			# newer mixOmics versions append a non-block element (e.g. 'comp') to the output; keep only real blocks
			blocks_available <- names(sel)[vapply(sel, is.list, logical(1))]
			for(b in blocks_available){
				block_sel <- sel[[b]]
				if(!is.null(block_sel$name)){
					feature_names <- block_sel$name
					# the column name of the loading data.frame differs across mixOmics versions (e.g. 'value.var')
					loading_values <- if(is.data.frame(block_sel$value)) block_sel$value[[1]] else as.vector(block_sel$value)
					# align the loading values with feature names
					if(length(loading_values) == length(feature_names)){
						res_list[[b]] <- data.frame(
							Feature = feature_names,
							Block = b,
							Comp = comp,
							Loading_value = loading_values,
							stringsAsFactors = FALSE
						)
					}else{
						warning("Loading values misaligned with feature names in block '", b, "' (comp = ", comp, "); Loading_value set to NA!")
						res_list[[b]] <- data.frame(
							Feature = feature_names,
							Block = b,
							Comp = comp,
							Loading_value = NA_real_,
							stringsAsFactors = FALSE
						)
					}
				}
			}
			if(length(res_list) > 0){
				res_df <- as.data.frame(dplyr::bind_rows(res_list), stringsAsFactors = FALSE)
				rownames(res_df) <- NULL
					res_df <- res_df[order(abs(res_df$Loading_value), decreasing = TRUE), ]
				}else{
					res_df <- data.frame(Feature = character(), Block = character(), Comp = integer(), Loading_value = numeric(), stringsAsFactors = FALSE)
				}
			self$res_features <- res_df
			message("Extracted ", nrow(res_df), " features. Result is stored in object$res_features ...")
			res_df
		},

		#' @description
		#' Extract the correlation/similarity matrix between selected variables from the model.
		#'
		#' This function uses the internal similarity computation of the model to derive a correlation
		#' matrix between selected variables across blocks. The result can be used for external network
		#' analysis or pathway enrichment.
		#' @param comp default 1; the component(s) used to compute the similarity (identical across all blocks).
		#' @param block default NULL; the block(s) to include. NULL means all blocks.
		#' @param cutoff default NULL; only show correlations with absolute value above this threshold.
		#' @return \code{res_correlation} stored in the object and returned.
		#' @examples
		#' \dontrun{
		#' cor_mat <- t1$get_correlation(comp = 1, cutoff = 0.7)
		#' }
		get_correlation = function(
			comp = 1,
			block = NULL,
			cutoff = NULL
			){
			private$check_model()
			private$check_comp(comp)
			# use network function to extract the similarity matrix
			# note: mixOmics >= 6.x uses plot.graph to control plotting (older versions used show),
			# and requires comp as a named list for block models
			comp_fmt <- private$format_network_comp(comp, block)
			if("plot.graph" %in% names(formals(mixOmics::network))){
				net_res <- mixOmics::network(self$model, comp = comp_fmt, blocks = block, plot.graph = FALSE)
			}else{
				net_res <- mixOmics::network(self$model, comp = comp_fmt, blocks = block, show = FALSE)
			}
			if(!is.null(net_res)){
				M <- net_res$M
				if(!is.null(cutoff)){
					M[abs(M) < cutoff] <- 0
				}
				self$res_correlation <- M
				message("Correlation matrix is stored in object$res_correlation ...")
				return(M)
			}
			message("No correlation matrix extracted ...")
			invisible(NULL)
		},

		#' @description
		#' Extract the feature selection stability from the cross-validation results.
		#'
		#' This function summarizes how often each feature was selected across the cross-validation
		#' folds and repeats (stored in \code{perf_result$features$stable} of \code{mixOmics::perf}).
		#' Features with high stability (close to 1) are robust biomarker candidates, while features
		#' with low stability may be artifacts of a particular data split.
		#' @param comp default NULL; the component(s) to extract. NULL means all components.
		#' @return \code{res_stability} stored in the object and returned; a tidy data.frame with columns
		#' 	 Feature, Block, Comp and Stability (mean selection frequency across folds and repeats).
		#' @examples
		#' \dontrun{
		#' t1$eval_performance(folds = 5, nrepeat = 10, seed = 123)
		#' stab <- t1$get_stability()
		#' }
		get_stability = function(
			comp = NULL
			){
			private$check_model()
			if(is.null(self$perf_result)){
				stop("The perf_result is NULL! Please run eval_performance() first!")
			}
			stable <- self$perf_result$features$stable
			if(is.null(stable)){
				stop("No feature stability information found in perf_result! ",
					"Feature stability requires a sparse supervised model (e.g. block.splsda with keepX) evaluated with eval_performance()!")
			}
			# structure of stable: stable[[nrep]][[block]][[compX]] = named frequency table
			res_list <- list()
			nreps <- names(stable)
			blocks <- unique(unlist(lapply(stable, names)))
			for(b in blocks){
				comps <- unique(unlist(lapply(stable, function(x){
					if(!is.null(x[[b]])) names(x[[b]]) else NULL
				})))
				for(cp in comps){
					comp_num <- as.integer(gsub("comp", "", cp))
					if(!is.null(comp) && !(comp_num %in% comp)){
						next
					}
					# collect the frequency tables across repeats
					tabs <- lapply(nreps, function(nr){
						stable[[nr]][[b]][[cp]]
					})
					tabs <- tabs[!vapply(tabs, is.null, logical(1))]
					if(length(tabs) == 0){
						next
					}
					# a feature absent in one repeat was never selected there: fill with 0 and average
					all_feats <- unique(unlist(lapply(tabs, names)))
					freq_mat <- vapply(tabs, function(tb){
						v <- rep(0, length(all_feats))
						names(v) <- all_feats
						v[names(tb)] <- as.numeric(tb)
						v
					}, numeric(length(all_feats)))
					res_list[[paste(b, cp, sep = "_")]] <- data.frame(
						Feature = all_feats,
						Block = b,
						Comp = comp_num,
						Stability = as.numeric(rowMeans(freq_mat)),
						stringsAsFactors = FALSE
					)
				}
			}
			if(length(res_list) > 0){
				res_df <- as.data.frame(dplyr::bind_rows(res_list), stringsAsFactors = FALSE)
				rownames(res_df) <- NULL
				res_df <- res_df[order(res_df$Comp, -res_df$Stability), ]
			}else{
				res_df <- data.frame(Feature = character(), Block = character(), Comp = integer(), Stability = numeric(), stringsAsFactors = FALSE)
			}
			self$res_stability <- res_df
			message("Extracted stability for ", nrow(res_df), " features. Result is stored in object$res_stability ...")
			res_df
		},

		#' @description
		#' Plot the samples projected onto the latent components.
		#'
		#' This function wraps \code{mixOmics::plotIndiv} with \code{style = "ggplot2"} to visualize
		#' the clustering of samples in the latent space. The returned object is a ggplot object that
		#' can be further customized with \code{+ theme_xx()}.
		#' @param comp default c(1, 2); the components on the x and y axes.
		#' @param group default NULL; the group factor for coloring. If NULL, use \code{self$Y}.
		#' @param blocks default "average"; which blocks to display. "average" for consensus, or block index/name.
		#' @param ellipse default FALSE; whether to add ellipses.
		#' @param legend default TRUE; whether to show the legend.
		#' @param ... other parameters passed to \code{mixOmics::plotIndiv}.
		#' @return a ggplot object.
		#' @examples
		#' \dontrun{
		#' p <- t1$plot_samples(comp = c(1, 2), ellipse = TRUE)
		#' }
		plot_samples = function(
			comp = c(1, 2),
			group = NULL,
			blocks = "average",
			ellipse = FALSE,
			legend = TRUE,
			...
			){
			private$check_model()
			private$check_comp(comp)
			if(is.null(group)){
				group <- self$Y
			}
			# note: plotIndiv treats an explicitly passed group = NULL differently from a missing group
			# (missingness check internally), so only include group when it is not NULL
			plotindiv_args <- list(object = self$model, comp = comp, blocks = blocks, style = "ggplot2",
				ellipse = ellipse, legend = legend, ...)
			if(!is.null(group)){
				plotindiv_args$group <- group
			}
			p <- tryCatch(do.call(mixOmics::plotIndiv, plotindiv_args), error = function(e) e)
			# plotIndiv with blocks = "average" fails for block.spls (indY) models in some mixOmics versions;
			# fall back to the first block in that case
			if(inherits(p, "error")){
				if(identical(blocks, "average")){
					warning("plotIndiv with blocks = 'average' failed for this model (", conditionMessage(p), "); falling back to the first block!")
					plotindiv_args$blocks <- 1
					p <- do.call(mixOmics::plotIndiv, plotindiv_args)
				}else{
					stop(p)
				}
			}
			# with style = "ggplot2", plotIndiv (mixOmics >= 6.x) returns a list with the ggplot in element 'graph'
			if(is.list(p) && !inherits(p, "ggplot") && inherits(p$graph, "ggplot")){
				return(p$graph)
			}
			# plotIndiv may return a list of ggplot objects (one per block) for multi-block models
			if(inherits(p, "list") && !inherits(p, "ggplot")){
				message("plotIndiv returned a list of plots (one per block). Use blocks = 'average' or a single block index/name to get a single ggplot object.")
			}
			p
		},

		#' @description
		#' Plot the loading weights of selected variables.
		#'
		#' This function wraps \code{mixOmics::plotLoadings} with \code{style = "ggplot2"} to visualize
		#' the contribution of selected variables on a given component.
		#' @param comp default 1; the component to plot.
		#' @param block default NULL; the block to display. NULL means the first block.
		#' @param ... other parameters passed to \code{mixOmics::plotLoadings}.
		#' @return a ggplot object.
		#' @examples
		#' \dontrun{
		#' p <- t1$plot_loadings(comp = 1, block = "microb")
		#' }
		plot_loadings = function(
			comp = 1,
			block = NULL,
			...
			){
			private$check_model()
			private$check_comp(comp)
			p <- mixOmics::plotLoadings(self$model, comp = comp, block = block, style = "ggplot2", ...)
			p
		},

		#' @description
		#' Plot the relevance network showing correlations between variables across blocks.
		#'
		#' This function wraps \code{mixOmics::network} to display the relevance network of selected
		#' variables. The network data is also stored in \code{res_network} for potential bridging with
		#' the \code{\link{trans_network}} class.
		#' @param comp default 1; the component used to compute the similarity (identical across all blocks).
		#' @param cutoff default 0.7; only show edges with absolute correlation above this threshold.
		#' @param block default NULL; the block(s) to include.
		#' @param ... other parameters passed to \code{mixOmics::network}.
		#' @return the network result (invisibly), with \code{res_network} stored.
		#' @examples
		#' \dontrun{
		#' t1$plot_cor_network(comp = 1, cutoff = 0.7)
		#' }
		plot_cor_network = function(
			comp = 1,
			cutoff = 0.7,
			block = NULL,
			...
			){
			private$check_model()
			private$check_comp(comp)
			# note: mixOmics >= 6.x uses plot.graph to control plotting (older versions used show),
			# and requires comp as a named list for block models
			comp_fmt <- private$format_network_comp(comp, block)
			if("plot.graph" %in% names(formals(mixOmics::network))){
				net_res <- mixOmics::network(self$model, comp = comp_fmt, cutoff = cutoff, blocks = block, plot.graph = TRUE, ...)
			}else{
				net_res <- mixOmics::network(self$model, comp = comp_fmt, cutoff = cutoff, blocks = block, show = TRUE, ...)
			}
			self$res_network <- net_res
			message("Network data is stored in object$res_network ...")
			invisible(net_res)
		},

		#' @description
		#' Plot the circos plot showing correlations between selected variables across blocks.
		#'
		#' This function wraps \code{mixOmics::circosPlot} to display a circos visualization of the
		#' correlations between selected variables from different blocks.
		#' @param comp default 1; the component(s) to display.
		#' @param cutoff default 0.7; only show links with absolute correlation above this threshold.
		#' @param ... other parameters passed to \code{mixOmics::circosPlot}.
		#' @return the circos plot result (invisibly).
		#' @examples
		#' \dontrun{
		#' t1$plot_circos(comp = 1, cutoff = 0.7)
		#' }
		plot_circos = function(
			comp = 1,
			cutoff = 0.7,
			...
			){
			private$check_model()
			private$check_comp(comp)
			res <- mixOmics::circosPlot(self$model, comp = comp, cutoff = cutoff, ...)
			invisible(res)
		},

		#' @description
		#' Plot the diagnostic plot for DIABLO showing the correlation between components of each block.
		#'
		#' This function wraps \code{mixOmics::plotDiablo} to check whether the correlations between
		#' components from each data set were maximized as specified in the design matrix.
		#' @param ncomp default 1; the dimension to be assessed.
		#' @param ... other parameters passed to \code{mixOmics::plotDiablo}.
		#' @return the plot (invisibly).
		#' @examples
		#' \dontrun{
		#' t1$plot_diablo(ncomp = 1)
		#' }
		plot_diablo = function(
			ncomp = 1,
			...
			){
			private$check_model()
			private$check_comp(ncomp)
			res <- mixOmics::plotDiablo(self$model, ncomp = ncomp, ...)
			invisible(res)
		},

		#' @description
		#' Plot the clustered image map (heatmap) of selected variables.
		#'
		#' This function wraps \code{mixOmics::cimDiablo} to display a clustered image map (CIM) of
		#' the selected variables across blocks.
		#' @param comp default 1; the component(s) to display.
		#' @param ... other parameters passed to \code{mixOmics::cimDiablo}.
		#' @return the CIM result (invisibly).
		#' @examples
		#' \dontrun{
		#' t1$plot_heatmap(comp = 1)
		#' }
		plot_heatmap = function(
			comp = 1,
			...
			){
			private$check_model()
			private$check_comp(comp)
			res <- mixOmics::cimDiablo(self$model, comp = comp, ...)
			invisible(res)
		},
		#' @description
		#' Build a trans-kingdom (cross-omics) association network following the TkNA framework.
		#'
		#' This function implements a lightweight Transkingdom Network Analysis (TkNA) workflow:
		#' (1) reduce nodes by differential filtering within each omics block (when \code{Y} is available),
		#' (2) build a cross-block (bipartite) association network from correlations between features of
		#' different blocks, and (3) compute the Bipartite Betweenness Centrality (BiBC) to identify
		#' cross-omics hub nodes that are most likely to mediate the interactions between omics layers.
		#' Confounders (age, BMI, batch, ...) can be regressed out from all features before computing
		#' the associations (MaAsLin2-style marginal models) to reduce spurious cross-layer edges.
		#'
		#' The resulting network provides testable mechanistic hypotheses rather than causal conclusions;
		#' prioritized hubs/edges should be validated with independent cohorts or perturbation experiments.
		#' @param corr_method default "spearman"; the association method, "spearman" or "pearson".
		#' 	 The data in \code{data_list} must be compositional-transformed (e.g. CLR via \code{preprocess_data})
		#' 	 before entering correlation-based network analysis.
		#' @param corr_thres default 0.6; the minimal absolute correlation to retain an edge.
		#' @param p_thres default 0.05; the maximal adjusted p value to retain an edge.
		#' @param p_adjust default "BH"; the p value adjustment method used for both the differential
		#' 	 filtering (within each block) and the correlation testing (across all cross-block pairs).
		#' @param diff_test default "kruskal"; the differential test for node filtering:
		#' 	 "kruskal" (Kruskal-Wallis rank sum test, any number of groups), "aov" (one-way ANOVA) or
		#' 	 "wilcox" (Wilcoxon rank sum test, exactly 2 groups). Use NULL to skip differential filtering.
		#' @param diff_p default 0.05; the adjusted p value threshold of the differential filtering.
		#' @param covariates default NULL; a data.frame of confounders (rows are samples with rownames matching
		#' 	 the sample names) to be regressed out from all features before computing correlations.
		#' 	 NULL means no confounder adjustment.
		#' @param blocks default NULL; the blocks to include (block names or indices). NULL means all blocks.
		#' @param hub_quantile default 0.9; features with a BiBC above this quantile (among features with
		#' 	 BiBC > 0) are flagged as cross-omics hubs.
		#' @return \code{self} (invisible); the network is stored in \code{res_transkingdom}
		#' 	 (a list with \code{graph}, \code{nodes} and \code{edges}).
		#' @examples
		#' \dontrun{
		#' t1$cal_transkingdom(corr_method = "spearman", corr_thres = 0.6, diff_test = "kruskal")
		#' head(t1$res_transkingdom$nodes)
		#' head(t1$res_transkingdom$edges)
		#' }
		cal_transkingdom = function(
			corr_method = "spearman",
			corr_thres = 0.6,
			p_thres = 0.05,
			p_adjust = "BH",
			diff_test = "kruskal",
			diff_p = 0.05,
			covariates = NULL,
			blocks = NULL,
			hub_quantile = 0.9
			){
			private$check_data_ready()
			if(!requireNamespace("igraph", quietly = TRUE)){
				stop("Please install igraph package first: install.packages('igraph')")
			}
			corr_method <- match.arg(corr_method, c("spearman", "pearson"))
			# ---- resolve blocks ----
			blocks_use <- private$resolve_blocks(blocks)
			if(length(blocks_use) < 2){
				stop("cal_transkingdom requires at least 2 blocks!")
			}
			samples_use <- rownames(self$data_list[[blocks_use[1]]])
			# ---- validate covariates ----
			if(!is.null(covariates)){
				if(!is.data.frame(covariates) || is.null(rownames(covariates)) || !setequal(rownames(covariates), samples_use)){
					stop("covariates must be a data.frame whose rownames match the sample names of the blocks!")
				}
				if(anyNA(covariates)){
					stop("covariates contain NA values! Please remove or impute them before confounder adjustment!")
				}
				covariates <- covariates[samples_use, , drop = FALSE]
			}
			# ---- step 1: differential filtering (TkNA node reduction) ----
			do_diff <- !is.null(diff_test) && !is.null(diff_p)
			if(do_diff && is.null(self$Y)){
				message("Y is NULL; differential filtering is skipped and all features are retained ...")
				do_diff <- FALSE
			}
			mat_list <- list()
			node_rows <- list()
			for(b in blocks_use){
				mat <- self$data_list[[b]]
				padj <- setNames(rep(NA_real_, ncol(mat)), colnames(mat))
				if(do_diff){
					pvals <- vapply(seq_len(ncol(mat)), function(j){
						private$diff_test_single(mat[, j], self$Y, diff_test)
					}, numeric(1))
					padj <- stats::p.adjust(pvals, method = p_adjust)
					names(padj) <- colnames(mat)
				}
				keep <- if(do_diff) (!is.na(padj) & padj <= diff_p) else rep(TRUE, ncol(mat))
				if(all(!keep)){
					warning("No feature passed the differential filtering (diff_p = ", diff_p,
						") in block '", b, "'; all features are retained instead!")
					keep <- rep(TRUE, ncol(mat))
				}
				mat_list[[b]] <- mat[, keep, drop = FALSE]
				node_rows[[b]] <- data.frame(Feature = colnames(mat_list[[b]]), Block = b,
					Diff_p = unname(padj[keep]), stringsAsFactors = FALSE)
			}
			nodes_df <- do.call(rbind, node_rows)
			rownames(nodes_df) <- NULL
			if(do_diff){
				message("Differential filtering (", diff_test, ", ", p_adjust, "-adjusted p <= ", diff_p, ") retained ",
					nrow(nodes_df), "/", sum(vapply(self$data_list[blocks_use], ncol, numeric(1))), " features ...")
			}
			# ---- step 2: cross-block correlations (bipartite edges) ----
			block_pairs <- utils::combn(blocks_use, 2, simplify = FALSE)
			# Candidate p values are kept as compact numeric vectors for the multiple-testing
			# correction (across all cross-block pairs); the edge table is assembled only for
			# the pairs passing the correlation threshold to reduce the memory peak.
			p_list <- list()
			keep_list <- list()
			edge_rows <- list()
			for(pr in block_pairs){
				X <- mat_list[[pr[1]]]
				Ymat <- mat_list[[pr[2]]]
				cor_res <- private$cal_block_corr(X, Ymat, corr_method, covariates)
				R <- cor_res$r
				P <- cor_res$p
				pair_id <- paste(pr, collapse = "_")
				p_list[[pair_id]] <- as.vector(P)
				# linear indices follow the column-major order of as.vector()
				keep_lin <- which(!is.na(R) & abs(R) >= corr_thres)
				keep_list[[pair_id]] <- keep_lin
				if(length(keep_lin) > 0){
					edge_rows[[pair_id]] <- data.frame(
						From = rep(rownames(R), times = ncol(R))[keep_lin],
						To = rep(colnames(R), each = nrow(R))[keep_lin],
						Corr = as.vector(R)[keep_lin],
						P = as.vector(P)[keep_lin],
						Block_from = pr[1],
						Block_to = pr[2],
						stringsAsFactors = FALSE
					)
				}
			}
			n_candidate <- sum(lengths(p_list))
			# multiple-testing correction across all cross-block pairs
			p_adj_all <- stats::p.adjust(unlist(p_list, use.names = FALSE), method = p_adjust)
			# map the adjusted p values back onto the retained edges (pairs keep their original order)
			offsets <- cumsum(lengths(p_list)) - lengths(p_list)
			keep_global <- unlist(Map(function(off, kl) off + kl, offsets, keep_list), use.names = FALSE)
			edges_df <- do.call(rbind, edge_rows)
			if(is.null(edges_df)){
				edges_df <- data.frame(From = character(), To = character(), Corr = numeric(), P = numeric(),
					Block_from = character(), Block_to = character(), P_adj = numeric(), stringsAsFactors = FALSE)
			}else{
				rownames(edges_df) <- NULL
				edges_df$P_adj <- unname(p_adj_all[keep_global])
				edges_df <- edges_df[edges_df$P_adj <= p_thres, , drop = FALSE]
			}
			if(nrow(edges_df) == 0){
				warning("No edge passed the thresholds (corr_thres = ", corr_thres, ", p_thres = ", p_thres,
					")! Consider lowering corr_thres or p_thres, or relaxing the differential filtering.")
			}
			message("Cross-block correlation (", corr_method, ", ", p_adjust, "-adjusted p <= ", p_thres,
				", |corr| >= ", corr_thres, ") retained ", nrow(edges_df), "/", n_candidate, " candidate edges ...")
			# ---- step 3: build the bipartite network and compute BiBC ----
			g <- igraph::graph_from_data_frame(
				d = edges_df[, c("From", "To", "Corr", "P", "P_adj")],
				directed = FALSE,
				vertices = nodes_df[, c("Feature", "Block", "Diff_p")]
			)
			# edge distance for path-based centrality: strongly correlated features are close
			igraph::E(g)$weight <- pmax(1 - abs(igraph::E(g)$Corr), 1e-12)
			bipc_vals <- private$cal_bipc(g, layer_vec = nodes_df$Block)
			g <- igraph::set_vertex_attr(g, "BiBC", value = unname(bipc_vals))
			# hub: top BiBC quantile among nodes lying on cross-layer paths
			bipc_pos <- bipc_vals[bipc_vals > 0]
			hub_flag <- rep(FALSE, length(bipc_vals))
			if(length(bipc_pos) > 0){
				hub_cut <- stats::quantile(bipc_pos, probs = hub_quantile)
				hub_flag <- bipc_vals >= hub_cut & bipc_vals > 0
			}
			g <- igraph::set_vertex_attr(g, "Hub", value = hub_flag)
			# ---- tidy output tables ----
			nodes_out <- igraph::as_data_frame(g, what = "vertices")
			colnames(nodes_out)[1] <- "Feature"
			nodes_out$Hub <- as.logical(nodes_out$Hub)
			edges_out <- igraph::as_data_frame(g, what = "edges")
			colnames(edges_out)[1:2] <- c("From", "To")
			rownames(nodes_out) <- NULL
			rownames(edges_out) <- NULL
			self$res_transkingdom <- list(graph = g, nodes = nodes_out, edges = edges_out)
			message("Trans-kingdom network built: ", igraph::vcount(g), " nodes, ", igraph::ecount(g),
				" edges, ", sum(hub_flag), " hub(s). Result is stored in object$res_transkingdom ...")
			invisible(self)
		},
		#' @description
		#' Plot the trans-kingdom (cross-omics) network from \code{cal_transkingdom}.
		#'
		#' Nodes are colored by omics block and sized by BiBC (bipartite betweenness centrality);
		#' edges are colored by the sign and strength of the correlation. Cross-omics hub features
		#' (flagged by \code{cal_transkingdom}) are labeled by default. The returned object is a
		#' ggplot object that can be further customized with \code{+ theme_xx()}.
		#' @param label default "hub"; which nodes to label: "hub" (only cross-omics hubs),
		#' 	 "all" (all nodes) or "none".
		#' @param layout default "fr"; the node layout: "fr" (force-directed Fruchterman-Reingold,
		#' 	 using the absolute correlations as attraction) or "circle".
		#' @param node_size_range default c(3, 8); the size range mapping BiBC to node size.
		#' @param edge_colors default c("#B2182B", "#2166AC"); colors for positive and negative correlations.
		#' @return a ggplot object.
		#' @examples
		#' \dontrun{
		#' t1$cal_transkingdom()
		#' p <- t1$plot_transkingdom(label = "hub")
		#' print(p)
		#' }
		plot_transkingdom = function(
			label = "hub",
			layout = "fr",
			node_size_range = c(3, 8),
			edge_colors = c("#B2182B", "#2166AC")
			){
			label <- match.arg(label, c("hub", "all", "none"))
			layout <- match.arg(layout, c("fr", "circle"))
			if(!requireNamespace("igraph", quietly = TRUE)){
				stop("Please install igraph package first: install.packages('igraph')")
			}
			if(!requireNamespace("ggplot2", quietly = TRUE)){
				stop("Please install ggplot2 package first: install.packages('ggplot2')")
			}
			if(is.null(self$res_transkingdom)){
				stop("The res_transkingdom is NULL! Please run cal_transkingdom() first!")
			}
			g <- self$res_transkingdom$graph
			if(igraph::ecount(g) == 0){
				warning("The trans-kingdom network has no edge; nothing to plot!")
				return(invisible(NULL))
			}
			# ---- layout ----
			if(layout == "fr"){
				coords <- igraph::layout_with_fr(g, weights = abs(igraph::E(g)$Corr))
			}else{
				coords <- igraph::layout_in_circle(g)
			}
			node_df <- self$res_transkingdom$nodes
			node_df$x <- coords[, 1]
			node_df$y <- coords[, 2]
			edge_df <- self$res_transkingdom$edges
			edge_df$x <- node_df$x[match(edge_df$From, node_df$Feature)]
			edge_df$y <- node_df$y[match(edge_df$From, node_df$Feature)]
			edge_df$xend <- node_df$x[match(edge_df$To, node_df$Feature)]
			edge_df$yend <- node_df$y[match(edge_df$To, node_df$Feature)]
			# ---- ggplot assembly ----
			p <- ggplot2::ggplot() +
				ggplot2::geom_segment(data = edge_df,
					ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = Corr),
					alpha = 0.5, show.legend = TRUE) +
				ggplot2::scale_color_gradient2(low = edge_colors[2], mid = "grey85", high = edge_colors[1],
					midpoint = 0, name = "Correlation") +
				ggplot2::geom_point(data = node_df,
					ggplot2::aes(x = x, y = y, fill = Block, size = BiBC),
					shape = 21, color = "grey30", stroke = 0.3) +
				ggplot2::scale_size_continuous(range = node_size_range, name = "BiBC")
			# RColorBrewer Set2 supports at most 8 groups; fall back to a hue palette for more blocks
			n_block_cols <- length(unique(node_df$Block))
			if(n_block_cols <= 8){
				p <- p + ggplot2::scale_fill_brewer(palette = "Set2", name = "Block")
			}else{
				p <- p + ggplot2::scale_fill_manual(values = scales::hue_pal()(n_block_cols), name = "Block")
			}
			if(label != "none"){
				label_df <- if(label == "hub") node_df[node_df$Hub, , drop = FALSE] else node_df
				if(nrow(label_df) > 0){
					p <- p + ggplot2::geom_text(data = label_df,
						ggplot2::aes(x = x, y = y, label = Feature), vjust = -0.9, size = 3)
				}
			}
			p + ggplot2::theme_classic(base_size = 12) +
				ggplot2::labs(x = NULL, y = NULL, title = "Trans-kingdom network") +
				ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
					axis.ticks = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5))
		},

		#' @description
		#' Compute the Area Under the ROC Curve (AUC) for the model.
		#'
		#' This function wraps \code{mixOmics::auroc} to perform ROC analysis for supervised models.
		#' @param comp default 1; the component to use for ROC analysis.
		#' @param ... other parameters passed to \code{mixOmics::auroc}.
		#' @return the AUC result (invisibly).
		#' @examples
		#' \dontrun{
		#' t1$cal_auroc(comp = 1)
		#' }
		cal_auroc = function(
			comp = 1,
			...
			){
			private$check_model()
			private$check_comp(comp)
			# mixOmics::auroc methods use roc.comp instead of comp; detect the dispatched method to stay version-compatible
			auroc_args <- list(object = self$model, ...)
			use_roc_comp <- FALSE
			for(cl in class(self$model)){
				f <- tryCatch(getFromNamespace(paste0("auroc.", cl), "mixOmics"), error = function(e) NULL)
				if(!is.null(f)){
					use_roc_comp <- "roc.comp" %in% names(formals(f))
					break
				}
			}
			if(use_roc_comp){
				auroc_args$roc.comp <- comp
			}else{
				auroc_args$comp <- comp
			}
			res <- do.call(mixOmics::auroc, auroc_args)
			invisible(res)
		},
		#' @description
		#' Run mediation analysis bridging to the \code{multimedia} package.
		#'
		#' This function is a thin wrapper around the \code{multimedia} package (Jiang et al. 2024,
		#' <doi:10.1101/2024.03.27.587024>) to quantify the pathway "exposure -> mediator -> outcome",
		#' e.g. microbes -> metabolites -> host phenotype. The high-dimensional blocks are taken from
		#' \code{data_list} (compositional-transformed by \code{preprocess_data}); the exposure and the
		#' outcome can alternatively be low-dimensional named vectors/factors (e.g. a treatment indicator
		#' or a continuous phenotype). The result provides direct/indirect (mediation) effect estimates
		#' that decompose the total effect, plus optional bootstrap confidence intervals.
		#'
		#' Note that the "causal" interpretation of the mediation effects relies on the assumed ordering
		#' of the chain and remains observational; validation with longitudinal or perturbation data is
		#' recommended. Variable names are replaced internally by syntactically valid names (V1, V2, ...)
		#' because multimedia builds R formulas from the variable names; the original names are restored
		#' in the returned effect tables (the mapping is stored in \code{res_mediation$name_map}).
		#' @param exposure the exposure (treatment): either the name of one omics block in
		#' 	 \code{data_list_names} (high-dimensional), or a named numeric vector/factor whose names match
		#' 	 the sample names (low-dimensional, e.g. the grouping factor).
		#' @param mediator the name of the omics block acting as the mediator (e.g. \code{"metab"}).
		#' @param outcome the outcome: either the name of one omics block in \code{data_list_names}
		#' 	 (e.g. host transcriptome block), or a named numeric vector/2-level factor whose names match
		#' 	 the sample names.
		#' @param model default "glmnet"; the mediation/outcome model engine passed to \code{multimedia}:
		#' 	 "glmnet" (regularized regression, recommended for high-dimensional blocks),
		#' 	 "lm" (linear model, only for low-dimensional exposure) or "rf" (random forest).
		#' @param features default NULL; an optional named list of feature subsets per block
		#' 	 (e.g. \code{list(microb = c("m1", "m2"))}), useful to focus on features selected by
		#' 	 \code{get_features} of a previous DIABLO analysis.
		#' @param max_features default NULL; an optional cap on the number of features per block; when a
		#' 	 block exceeds it, only the most variable features are kept.
		#' @param covariates default NULL; a data.frame of covariates (rows are samples with rownames
		#' 	 matching the sample names) passed to the \code{pretreatments} of \code{multimedia}.
		#' @param n_boot default 0; the number of bootstrap resamples for confidence intervals.
		#' 	 0 skips the bootstrap (recommended for a first exploration; use 1000 for final results).
		#' @return \code{self} (invisible); the result is stored in \code{res_mediation}.
		#' @examples
		#' \dontrun{
		#' # microbes -> metabolites -> host phenotype
		#' t1$run_mediation(exposure = "microb", mediator = "metab", outcome = pheno_vector,
		#'   model = "glmnet", n_boot = 200)
		#' t1$res_mediation$indirect_pathwise
		#' }
		run_mediation = function(
			exposure,
			mediator,
			outcome,
			model = "glmnet",
			features = NULL,
			max_features = NULL,
			covariates = NULL,
			n_boot = 0
			){
			if(!requireNamespace("multimedia", quietly = TRUE)){
				stop("Please install the multimedia package first: install.packages('multimedia')")
			}
			private$check_data_ready()
			samples_use <- rownames(self$data_list[[1]])
			# ---- resolve mediator ----
			if(!is.character(mediator) || length(mediator) != 1 || !(mediator %in% self$data_list_names)){
				stop("mediator must be the name of one omics block in data_list_names: ",
					paste(self$data_list_names, collapse = ", "), "!")
			}
			mediator_df <- as.data.frame(private$mediation_block_matrix(mediator, features, max_features))
			# ---- resolve exposure ----
			if(is.character(exposure) && length(exposure) == 1){
				if(!(exposure %in% self$data_list_names)){
					stop("exposure block '", exposure, "' not found in data_list_names: ",
						paste(self$data_list_names, collapse = ", "), "!")
				}
				if(exposure == mediator){
					stop("exposure and mediator must be different blocks!")
				}
				exposure_df <- as.data.frame(private$mediation_block_matrix(exposure, features, max_features))
				if(model == "lm" && ncol(exposure_df) > 1){
					warning("The exposure block '", exposure, "' is high-dimensional; consider model = 'glmnet' ",
						"because 'lm' may overfit or fail with more predictors than samples!")
				}
			}else{
				exposure_df <- private$mediation_vector(exposure, samples_use, "exposure")
				if(is.factor(exposure_df[[1]])){
					if(nlevels(exposure_df[[1]]) != 2){
						stop("The exposure factor must have exactly 2 levels (treatment vs control)! Please encode a multi-level exposure as numeric indicators first!")
					}
					exposure_df[[1]] <- as.numeric(exposure_df[[1]]) - 1
					message("The exposure factor is converted to 0/1 numeric for the regression models ...")
				}
			}
			# ---- resolve outcome ----
			if(is.character(outcome) && length(outcome) == 1 && outcome %in% self$data_list_names){
				if(outcome == mediator){
					stop("outcome and mediator must be different blocks!")
				}
				outcome_df <- as.data.frame(private$mediation_block_matrix(outcome, features, max_features))
			}else{
				outcome_df <- private$mediation_vector(outcome, samples_use, "outcome")
				if(is.factor(outcome_df[[1]])){
					if(nlevels(outcome_df[[1]]) != 2){
						stop("The outcome factor must have exactly 2 levels! Please encode a multi-level outcome as numeric variables first!")
					}
					outcome_df[[1]] <- as.numeric(outcome_df[[1]]) - 1
					message("The outcome factor is converted to 0/1 numeric for the regression models ...")
				}
			}
			# ---- covariates -> pretreatments ----
			pretreat_names <- NULL
			if(!is.null(covariates)){
				if(!is.data.frame(covariates) || is.null(rownames(covariates)) || !setequal(rownames(covariates), samples_use)){
					stop("covariates must be a data.frame whose rownames match the sample names of the blocks!")
				}
				if(anyNA(covariates)){
					stop("covariates contain NA values! Please remove or impute them before mediation analysis!")
				}
				cov_df <- covariates[samples_use, , drop = FALSE]
				pretreat_names <- colnames(cov_df)
			}else{
				cov_df <- NULL
			}
			# ---- assemble the mediation data.frame ----
			med_df <- cbind(exposure_df, mediator_df, outcome_df)
			if(!is.null(cov_df)){
				med_df <- cbind(med_df, cov_df)
			}
			if(any(duplicated(colnames(med_df)))){
				stop("Duplicated variable names across blocks! Please rename the features before mediation analysis!")
			}
			# ---- sanitize variable names: multimedia engines build R formulas from the variable
			# names, and metabolite names like "trehalose+d-Glucoheptose 1" are not syntactically
			# valid. To be robust against any internal name mangling, plain synthetic names
			# (V1, V2, ...) are used internally and restored in the returned effect tables ----
			n_exp <- ncol(exposure_df)
			n_med <- ncol(mediator_df)
			n_out <- ncol(outcome_df)
			n_cov <- if(is.null(cov_df)) 0 else ncol(cov_df)
			safe_names <- paste0("V", seq_len(ncol(med_df)))
			name_map <- setNames(colnames(med_df), safe_names)  # synthetic -> original
			if(!identical(safe_names, colnames(med_df))){
				message("Variable names are replaced internally by syntactically valid names for the ",
					"regression formulas; the results use the original names (see res_mediation$name_map) ...")
				colnames(med_df) <- safe_names
			}
			exper <- multimedia::mediation_data(med_df,
				outcomes = safe_names[n_exp + n_med + seq_len(n_out)],
				treatments = safe_names[seq_len(n_exp)],
				mediators = safe_names[n_exp + seq_len(n_med)],
				pretreatments = if(n_cov > 0) safe_names[n_exp + n_med + n_out + seq_len(n_cov)] else NULL)
			model_spec <- switch(model,
				lm = multimedia::lm_model(),
				glmnet = multimedia::glmnet_model(),
				rf = multimedia::rf_model(),
				stop("model must be one of 'lm', 'glmnet', 'rf'!")
			)
			message("Fitting the multimedia mediation model (", model, ") with ", nrow(med_df), " samples, ",
				ncol(exposure_df), " exposure, ", ncol(mediator_df), " mediator and ", ncol(outcome_df), " outcome variable(s) ...")
			mm <- multimedia::multimedia(exper, model_spec)
			fit <- multimedia::estimate(mm, exper)
			res <- list(
				fit = fit,
				direct_effect = private$restore_med_names(multimedia::direct_effect(fit), name_map),
				indirect_overall = private$restore_med_names(multimedia::indirect_overall(fit), name_map),
				indirect_pathwise = private$restore_med_names(multimedia::indirect_pathwise(fit), name_map),
				name_map = name_map
			)
			if(n_boot > 0){
				message("Running bootstrap with B = ", n_boot, " resamples (may take a while) ...")
				res$bootstrap <- private$restore_med_names(multimedia::bootstrap(mm, exper, B = n_boot,
					fs = list(multimedia::direct_effect, multimedia::indirect_overall, multimedia::indirect_pathwise)), name_map)
			}
			self$res_mediation <- res
			message("Mediation result is stored in object$res_mediation (direct_effect, indirect_overall, indirect_pathwise) ...")
			invisible(self)
		},

		#' @description
		#' Export feature tables formatted for knowledge-based mechanism tools (MIMOSA2 / AMON / anansi).
		#'
		#' This function writes the abundance tables of the selected blocks into the input formats expected
		#' by external knowledge-constrained tools, bridging layer 1 (data-driven association) with layer 2
		#' (biochemistry-constrained mechanism inference):
		#' \itemize{
		#' \item MIMOSA2: KO (or other microbial feature) abundance table + metabolite abundance table,
		#' 	 both tab-separated with features as rows and samples as columns, used to link community
		#' 	 metabolic potential to measured metabolites.
		#' \item AMON: metabolite KEGG ID list + KO ID list (one ID per line), used to attribute metabolite
		#' 	 origins (microbial vs host) on the KEGG reaction network.
		#' \item anansi: KO abundance table + metabolite abundance table (features as rows, samples as columns)
		#' 	 for constrained KO-metabolite association testing.
		#' }
		#' The exported tables use the original (non-transformed) abundances stored in \code{dataset_list},
		#' because these tools expect raw counts or relative abundances rather than CLR-transformed values.
		#' @param target default "mimosa2"; the target tool: "mimosa2", "amon" or "anansi".
		#' @param file default NULL; the path prefix of the exported files (without extension), e.g.
		#' 	 \code{"./bridge/input"} produces \code{"./bridge/input_ko.txt"} and \code{"./bridge/input_metabolites.txt"}.
		#' @param ko_block default NULL; the block holding the microbial features (KO/gene/taxon abundances).
		#' 	 NULL uses the first block of \code{data_list_names}.
		#' @param metab_block default NULL; the block holding the metabolite abundances. NULL auto-detects
		#' 	 a block whose name contains "metab" (case-insensitive) and falls back to the second block.
		#' @param metab_id_table default NULL; an optional 2-column data.frame mapping the metabolite
		#' 	 names in the data (column 1) to KEGG compound IDs (column 2). If NULL, the original
		#' 	 metabolite names are used and a reminder is printed that MIMOSA2/AMON/anansi expect KEGG IDs.
		#' @param features default NULL; an optional named list of feature subsets per block
		#' 	 (e.g. \code{list(microb = c("K00001", ...), metab = c("C00095", ...))}); useful to export only
		#' 	 features prioritized by a previous DIABLO analysis (see \code{get_features}).
		#' @return a named character vector of the written file paths (invisibly).
		#' @examples
		#' \dontrun{
		#' t1$export_bridge(target = "mimosa2", file = "./bridge/mimosa_input",
		#'   ko_block = "microb", metab_block = "metab")
		#' }
		export_bridge = function(
			target = "mimosa2",
			file = NULL,
			ko_block = NULL,
			metab_block = NULL,
			metab_id_table = NULL,
			features = NULL
			){
			if(is.null(self$dataset_list)){
				stop("The dataset_list is NULL! Please create the object with microtables first!")
			}
			target <- match.arg(target, c("mimosa2", "amon", "anansi"))
			if(is.null(file)){
				stop("Please provide the file path prefix for the exported files!")
			}
			if(!dir.exists(dirname(file)) && dirname(file) != ""){
				dir.create(dirname(file), recursive = TRUE)
			}
			# ---- resolve blocks ----
			if(is.null(ko_block)){
				ko_block <- self$data_list_names[1]
				message("ko_block is NULL; using the first block '", ko_block, "' ...")
			}
			if(is.null(metab_block)){
				cand <- grep("metab", self$data_list_names, value = TRUE, ignore.case = TRUE)
				if(length(cand) > 0){
					metab_block <- cand[1]
				}else{
					metab_block <- self$data_list_names[2]
				}
				message("metab_block is NULL; using block '", metab_block, "' ...")
			}
			ko_block <- private$resolve_single_block(ko_block, "ko_block")
			metab_block <- private$resolve_single_block(metab_block, "metab_block")
			if(ko_block == metab_block){
				stop("ko_block and metab_block must be different blocks!")
			}
			# ---- extract original abundance tables (features x samples) ----
			ko_tab <- private$bridge_block_table(ko_block, features)
			metab_tab <- private$bridge_block_table(metab_block, features)
			if(!is.null(metab_id_table)){
				if(!is.data.frame(metab_id_table) || ncol(metab_id_table) != 2){
					stop("metab_id_table must be a 2-column data.frame: original metabolite names and KEGG compound IDs!")
				}
				id_map <- setNames(as.character(metab_id_table[[2]]), as.character(metab_id_table[[1]]))
				matched <- rownames(metab_tab) %in% names(id_map)
				if(any(matched)){
					rownames(metab_tab)[matched] <- id_map[rownames(metab_tab)[matched]]
				}
				if(sum(matched) < nrow(metab_tab)){
					warning(sum(matched), "/", nrow(metab_tab),
						" metabolites mapped to KEGG IDs; unmapped metabolites keep their original names ...")
				}
			}else{
				message("Note: MIMOSA2/AMON/anansi work best with KEGG compound IDs; provide metab_id_table to map the metabolite names ...")
			}
			files_written <- c()
			if(target %in% c("mimosa2", "anansi")){
				# anansi consumes the same two abundance tables as MIMOSA2
				ko_file <- paste0(file, "_ko.txt")
				metab_file <- paste0(file, "_metabolites.txt")
				utils::write.table(ko_tab, ko_file, sep = "\t", quote = FALSE, col.names = NA)
				utils::write.table(metab_tab, metab_file, sep = "\t", quote = FALSE, col.names = NA)
				files_written <- c(files_written, ko_file, metab_file)
			}
			if(target == "amon"){
				ko_id_file <- paste0(file, "_ko_ids.txt")
				metab_id_file <- paste0(file, "_metabolite_ids.txt")
				writeLines(rownames(ko_tab), ko_id_file)
				writeLines(rownames(metab_tab), metab_id_file)
				files_written <- c(files_written, ko_id_file, metab_id_file)
			}
			names(files_written) <- NULL
			message("Exported files for '", target, "': ", paste(files_written, collapse = ", "),
				"\nPlease check the format requirements of the target tool before uploading ...")
			invisible(files_written)
		}
	),
	private = list(
		# check mixOmics package availability
		check_mixomics = function(){
			if(!requireNamespace("mixOmics", quietly = TRUE)){
				stop("Please install mixOmics package: if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('mixOmics')")
			}
		},
		# check mixOmics availability and that data_list is ready
		check_data_list = function(){
			private$check_mixomics()
			if(is.null(self$data_list)){
				stop("The data_list is NULL! Please run preprocess_data() first!")
			}
		},
		# check mixOmics availability and that a model has been fitted
		check_model = function(){
			private$check_mixomics()
			if(is.null(self$model)){
				stop("The model is NULL! Please run run_diablo(), run_spls() or run_spls_single() first!")
			}
		},
		# check that the component index does not exceed the fitted model
		check_comp = function(comp){
			ncomp_model <- self$model$ncomp
			if(is.null(ncomp_model)){
				return(invisible(NULL))
			}
			if(any(comp > max(ncomp_model)) || any(comp < 1)){
				stop("comp (", paste(comp, collapse = ", "), ") is out of range! The model has ", max(ncomp_model), " component(s)!")
			}
		},
		# convert a numeric component vector to the format required by mixOmics::network:
		# block models (block.spls/block.pls) need a named list with identical components per block,
		# while two-block models (mixo_spls) take a plain numeric vector
		format_network_comp = function(comp, block = NULL){
			if(inherits(self$model, c("block.spls", "block.pls"))){
				block_names <- names(self$model$X)
				if(!is.null(block)){
					if(is.character(block)){
						block_use <- intersect(block, block_names)
					}else{
						block_use <- block_names[block]
					}
					if(length(block_use) > 0){
						block_names <- block_use
					}
				}
				comp_list <- lapply(block_names, function(x) comp)
				names(comp_list) <- block_names
				comp_list
			}else{
				comp
			}
		},
		# construct default keepX list (all features retained)
		build_keepX_default = function(ncomp){
			keepX_list <- lapply(self$data_list, function(x){
				rep(ncol(x), ncomp)
			})
			keepX_list
		},
		# check that data_list is ready (without requiring mixOmics, e.g. for cal_transkingdom)
		check_data_ready = function(){
			if(is.null(self$data_list)){
				stop("The data_list is NULL! Please run preprocess_data() first!")
			}
		},
		# resolve a vector of block names/indices into block names
		resolve_blocks = function(blocks){
			if(is.null(blocks)){
				return(self$data_list_names)
			}
			if(is.character(blocks)){
				blocks_idx <- match(blocks, self$data_list_names)
				if(any(is.na(blocks_idx))){
					stop("blocks not found in data_list_names: ",
						paste(blocks[is.na(blocks_idx)], collapse = ", "), "!")
				}
			}else{
				blocks_idx <- blocks
				if(any(is.na(match(blocks_idx, seq_along(self$data_list_names))))){
					stop("blocks must be valid block indices or names!")
				}
			}
			self$data_list_names[blocks_idx]
		},
		# resolve a single block name/index into a block name
		resolve_single_block = function(block, arg_name){
			if(is.character(block)){
				if(!(block %in% self$data_list_names)){
					stop(arg_name, " '", block, "' not found in data_list_names: ",
						paste(self$data_list_names, collapse = ", "), "!")
				}
				return(block)
			}
			if(!is.numeric(block) || length(block) != 1 || block < 1 || block > length(self$data_list_names)){
				stop(arg_name, " must be a single block name or index!")
			}
			self$data_list_names[block]
		},
		# differential test for one feature (TkNA node reduction)
		diff_test_single = function(x, Y, method){
			if(stats::sd(x) == 0){
				return(NA_real_)
			}
			switch(method,
				kruskal = stats::kruskal.test(x, Y)$p.value,
				aov = summary(stats::aov(x ~ Y))[[1]][1, "Pr(>F)"],
				wilcox = {
					if(nlevels(droplevels(Y)) != 2){
						stop("diff_test = 'wilcox' requires exactly 2 groups in Y! Use 'kruskal' or 'aov' for multi-group data!")
					}
					stats::wilcox.test(x ~ droplevels(Y))$p.value
				},
				stop("diff_test must be one of 'kruskal', 'aov', 'wilcox' or NULL!")
			)
		},
		# MaAsLin2-style confounder removal: regress the covariate design matrix out of each feature
		residualize_matrix = function(mat, cov_mat){
			res <- apply(mat, 2, function(y){
				stats::lm.fit(cov_mat, y)$residuals
			})
			# apply(., 2, f) returns a matrix with one column per feature
			rownames(res) <- rownames(mat)
			res
		},
		# batch correlation (with p values) between the features of two blocks;
		# spearman is computed as pearson on ranks; covariates are regressed out first
		cal_block_corr = function(X, Ymat, corr_method, covariates = NULL){
			n_cov <- 0
			if(!is.null(covariates)){
				# build the design matrix once: factors expand to dummy columns, so the degrees of
				# freedom must count the non-intercept columns of the design matrix, not ncol(covariates)
				cov_mat <- stats::model.matrix(stats::as.formula("~ ."), data = covariates)
				n_cov <- ncol(cov_mat) - 1
				X <- private$residualize_matrix(X, cov_mat)
				Ymat <- private$residualize_matrix(Ymat, cov_mat)
			}
			if(corr_method == "spearman"){
				X <- apply(X, 2, rank)
				Ymat <- apply(Ymat, 2, rank)
			}
			R <- suppressWarnings(stats::cor(X, Ymat))
			n <- nrow(X)
			df <- n - 2 - n_cov
			if(df <= 0){
				stop("Insufficient samples (", n, ") to compute correlations with ", n_cov, " covariate(s)!")
			}
			tstat <- abs(R) * sqrt(df / pmax(1 - R^2, 0))
			P <- 2 * stats::pt(-tstat, df)
			P[is.na(R)] <- NA_real_
			list(r = R, p = P)
		},
		# Bipartite betweenness centrality (BiBC) following the TkNA framework:
		# a modified Brandes algorithm (Dijkstra-based) that only counts shortest paths whose two
		# endpoints belong to different blocks, so that nodes bridging omics layers are highlighted.
		# Edge weights are interpreted as distances (strongly correlated features are close).
		cal_bipc = function(g, layer_vec){
			nv <- igraph::vcount(g)
			bipc <- numeric(nv)
			if(nv <= 2){
				return(bipc)
			}
			el <- igraph::as_edgelist(g, names = FALSE)
			if(nrow(el) == 0){
				return(bipc)
			}
			w <- igraph::E(g)$weight
			if(is.null(w)){
				w <- rep(1, nrow(el))
			}
			layer <- as.character(layer_vec)
			# adjacency list (each entry: matrix of [neighbor index, edge weight] rows)
			adj_list <- vector("list", nv)
			for(e in seq_len(nrow(el))){
				a <- el[e, 1]
				b <- el[e, 2]
				adj_list[[a]] <- rbind(adj_list[[a]], c(b, w[e]))
				adj_list[[b]] <- rbind(adj_list[[b]], c(a, w[e]))
			}
			tol <- 1e-9
			for(s in seq_len(nv)){
				dist <- rep(Inf, nv)
				sigma <- numeric(nv)
				pred <- vector("list", nv)
				dist[s] <- 0
				sigma[s] <- 1
				done <- rep(FALSE, nv)
				# ---- Dijkstra from s ----
				repeat{
					cand <- which(!done & is.finite(dist))
					if(length(cand) == 0){
						break
					}
					u <- cand[which.min(dist[cand])]
					done[u] <- TRUE
					nb <- adj_list[[u]]
					if(is.null(nb)){
						next
					}
					for(k in seq_len(nrow(nb))){
						v <- nb[k, 1]
						alt <- dist[u] + nb[k, 2]
						if(alt < dist[v] - tol){
							dist[v] <- alt
							sigma[v] <- sigma[u]
							pred[[v]] <- u
						}else if(abs(alt - dist[v]) <= tol){
							if(!(u %in% pred[[v]])){
								sigma[v] <- sigma[v] + sigma[u]
								pred[[v]] <- c(pred[[v]], u)
							}
						}
					}
				}
				# ---- back-propagate dependencies, counting only cross-layer endpoints ----
				delta <- numeric(nv)
				ord <- order(dist, decreasing = TRUE)
				ord <- ord[is.finite(dist[ord])]
				for(v in ord){
					if(v == s){
						next
					}
					# a cross-layer endpoint v contributes 1 to each of its predecessors;
					# delta[v] holds v's contribution as an intermediate node to further endpoints
					coeff <- delta[v] + (layer[v] != layer[s])
					pv <- pred[[v]]
					if(length(pv) > 0 && sigma[v] > 0){
						delta[pv] <- delta[pv] + (sigma[pv] / sigma[v]) * coeff
					}
				}
				bipc <- bipc + delta
			}
			names(bipc) <- igraph::vertex_attr(g, "name")
			bipc
		},
		# extract (and optionally subset/limit) a block matrix for mediation analysis
		mediation_block_matrix = function(block, features, max_features){
			mat <- self$data_list[[block]]
			if(!is.null(features) && !is.null(features[[block]])){
				keep <- intersect(features[[block]], colnames(mat))
				if(length(keep) == 0){
					stop("No feature of the '", block, "' block matches the provided features subset!")
				}
				mat <- mat[, keep, drop = FALSE]
			}
			if(!is.null(max_features) && ncol(mat) > max_features){
				vars <- apply(mat, 2, var)
				keep <- order(vars, decreasing = TRUE)[seq_len(max_features)]
				mat <- mat[, keep, drop = FALSE]
				message("Block '", block, "' reduced to the ", max_features, " most variable features for the mediation analysis ...")
			}
			mat
		},
		# rename the entries of a character vector through a map (names(x) -> x)
		rename_chr = function(v, map){
			if(is.null(v)){
				return(v)
			}
			idx <- match(v, names(map))
			v[!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
			v
		},
		# recursively restore the original variable names in the effect tables returned by
		# multimedia: names, colnames, rownames, dimnames, character cell values and factor levels
		restore_med_names = function(x, map){
			if(is.matrix(x)){
				dn <- dimnames(x)
				if(!is.null(dn)){
					dimnames(x) <- lapply(dn, private$rename_chr, map)
				}
				return(x)
			}
			if(is.data.frame(x)){
				# variable references may also appear as cell values (e.g. the outcome/mediator
				# columns of the effect tables)
				for(j in seq_len(ncol(x))){
					if(is.character(x[[j]])){
						idx <- match(x[[j]], names(map))
						if(any(!is.na(idx))){
							x[[j]][!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
						}
					}else if(is.factor(x[[j]])){
						lv <- levels(x[[j]])
						idx <- match(lv, names(map))
						if(any(!is.na(idx))){
							levels(x[[j]])[!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
						}
					}
				}
				colnames(x) <- private$rename_chr(colnames(x), map)
				rownames(x) <- private$rename_chr(rownames(x), map)
				return(x)
			}
			if(is.list(x)){
				names(x) <- private$rename_chr(names(x), map)
				for(i in seq_along(x)){
					x[[i]] <- private$restore_med_names(x[[i]], map)
				}
				return(x)
			}
			if(is.atomic(x) && !is.null(names(x))){
				names(x) <- private$rename_chr(names(x), map)
			}
			x
		},
		# align a named exposure/outcome vector to the samples and wrap it into a 1-column data.frame
		mediation_vector = function(vec, samples_use, var_name){
			if(is.null(names(vec))){
				if(length(vec) != length(samples_use)){
					stop("The length of ", var_name, " (", length(vec),
						") does not match the number of samples (", length(samples_use), ")!")
				}
				warning(var_name, " has no names; samples are matched by order. Please provide a named vector to ensure correct alignment!")
				out <- vec
			}else{
				if(!setequal(names(vec), samples_use)){
					stop("The names of ", var_name, " do not match the sample names of the blocks!")
				}
				out <- vec[samples_use]
			}
			res <- data.frame(out, stringsAsFactors = TRUE)
			colnames(res) <- var_name
			rownames(res) <- samples_use
			res
		},
		# extract an original (non-transformed) abundance table (features x samples) for export_bridge
		bridge_block_table = function(block, features){
			tab <- as.matrix(self$dataset_list[[block]]$otu_table)
			if(!is.null(features) && !is.null(features[[block]])){
				keep <- intersect(features[[block]], rownames(tab))
				if(length(keep) == 0){
					stop("No feature of the '", block, "' block matches the provided features subset!")
				}
				tab <- tab[keep, , drop = FALSE]
			}
			tab
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
