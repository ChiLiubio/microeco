#' @title 
#' Create \code{trans_classifier} object for machine-learning-based model prediction.
#'
#' @description
#' This class is a wrapper for methods of machine-learning-based classification or regression models, including data pre-processing, feature selection, 
#' data split, model training, prediction, confusionMatrix and ROC (Receiver Operator Characteristic) or PR (Precision-Recall) curve.
#'
#' Author(s): Felipe Mansoldo and Chi Liu
#'
#' @export
trans_classifier <- R6::R6Class(classname = "trans_classifier",
	public = list(
		#' @description
		#' Create a trans_classifier object.
		#' 
		#' @param dataset an object of \code{\link{microtable}} class.
		#' @param x.predictors default "Genus"; character vector or data.frame; a character vector represents selecting the corresponding data from \code{microtable$taxa_abund}; 
		#'   data.frame denotes other customized input. See the following available options:
		#'   \describe{
		#'     \item{\strong{'Genus'}}{use Genus level table in \code{microtable$taxa_abund}, or other specific taxonomic rank, e.g., 'Phylum'.
		#'        At most one input level (e.g., ASV) is allowed to be absent in the names of taxa_abund list; 
		#'        in that case the function will use \code{otu_table} to calculate relative abundance of features.}
		#'     \item{\strong{'all'}}{use all the levels stored in \code{microtable$taxa_abund}. Only available when it is the sole input.}
		#'     \item{\strong{multiple levels}}{use several levels at the same time, e.g., \code{c("Genus", "Species")}; 
		#'       the corresponding tables in \code{microtable$taxa_abund} are combined by rows as the predictors.}
		#'     \item{\strong{other input}}{must be a data.frame object. It should have the same format with the tables in microtable$taxa_abund, i.e. rows are features; 
		#'       columns are samples with same names in sample_table.}
		#'   }
		#' @param y.response default NULL; the response variable in \code{sample_table} of input \code{microtable} object.
		#' @param n.cores default 1; the number of CPU threads used in the model training and the feature selection. 
		#'   The parallel backend of the foreach package is registered only temporarily during the computation and the 
		#'   original backend is restored afterwards, so that the global parallel setting is not permanently changed.
		#' @return \code{data_feature} and \code{data_response} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_classifier$new(
		#' 		dataset = dataset, 
		#' 		x.predictors = "Genus",
		#' 		y.response = "Group")
		#' # use multiple taxonomic levels as the predictors
		#' t1 <- trans_classifier$new(
		#' 		dataset = dataset, 
		#' 		x.predictors = c("Genus", "Species"),
		#' 		y.response = "Group")
		#' }
		initialize = function(dataset,
				x.predictors = "Genus",
				y.response = NULL,
				n.cores = 1
			){
			check_microtable(dataset)
			sampleinfo <- dataset$sample_table
			if(is.null(y.response)){
				stop("No y.response provided!")
			}
			if(!y.response %in% colnames(sampleinfo)){
				stop("Input y.response must be a column name of dataset$sample_table!")
			}
			
			response_data <- sampleinfo[, y.response]
			if(is.numeric(response_data)){
				self$type <- "Regression"
				self$data_response <- response_data
				message("Regression for ", y.response)
			}else{
				self$type <- "Classification"
				if(nlevels(as.factor(response_data)) < 2) {
					stop("Response variable must have at least 2 factors!")
				}
				if(nlevels(as.factor(response_data)) == 2){
					ClassificationCase = "2 groups"
				}else{
					ClassificationCase = "2+ groups"
				}
				# self$ClassificationCase <- ClassificationCase
				message("Classification type = ", ClassificationCase)

				ClassNames <- make.names(as.character(response_data), unique = FALSE)
				MapNames <- data.frame(OriginalNames = as.character(response_data), 
					ClassNames = ClassNames, stringsAsFactors = FALSE)
				# as.character is necessary here, as identical() on a factor and a character vector is always FALSE
				if(!identical(MapNames$OriginalNames, MapNames$ClassNames)){
					message("Factor names are non-standard. A correction was made and the change map was saved in object$data_MapNames ...")
				}
				self$data_MapNames <- MapNames
				self$data_response <- ClassNames
			}
			# x.predictors must be character or data.frame
			if(is.character(x.predictors)){
				if(length(x.predictors) < 1 || any(is.na(x.predictors))){
					stop("Provided x.predictors must be a character vector without NA or a data.frame !")
				}
				if(any(duplicated(x.predictors))){
					stop("Duplicated names are found in the provided x.predictors: ", 
						paste(unique(x.predictors[duplicated(x.predictors)]), collapse = ", "), " ! Please check the input!")
				}
				if(is.null(dataset$taxa_abund)){
					message("No taxa_abund found in the dataset. Calculate the relative abundance ...")
					dataset$cal_abund()
				}
				# exact matching instead of grepl, otherwise any level name containing 'all' would be mis-matched
				if(length(x.predictors) == 1 && tolower(trimws(x.predictors)) == "all"){
					abund_table <- do.call(rbind, unname(dataset$taxa_abund))
				}else{
					if(any(tolower(trimws(x.predictors)) == "all")){
						stop("'all' can only be used as the sole x.predictors, e.g. x.predictors = 'all' ! Please check the input!")
					}
					# the levels that are not in taxa_abund are calculated from the otu_table; only one of them is allowed, 
					# as all of them are generated from the same row names of otu_table and would be duplicated features otherwise
					missing_levels <- x.predictors[! x.predictors %in% names(dataset$taxa_abund)]
					if(length(missing_levels) > 1){
						stop("Multiple unknown levels are found in the provided x.predictors: ", paste(missing_levels, collapse = ", "), 
							" ! Only one unknown level can be calculated from the otu_table at a time, as the other ones would ", 
							"generate duplicated features. Please check the input!")
					}
					if(length(missing_levels) == 1){
						message("x.predictors: ", missing_levels, " is not found in the names of taxa_abund list. Use the features in otu_table ...")
						if("add_rownames2tax" %in% names(dataset)){
							dataset$add_rownames2tax(use_name = missing_levels)
						}else{
							dataset$add_rownames2taxonomy(use_name = missing_levels)
						}
						message("Calculate relative abundance of features ...")
						suppressMessages(dataset$cal_abund())
					}
					if(length(x.predictors) == 1){
						abund_table <- dataset$taxa_abund[[x.predictors]]
					}else{
						abund_table <- do.call(rbind, unname(dataset$taxa_abund[x.predictors]))
						message("Multiple levels are used as the predictors: ", paste(x.predictors, collapse = ", "), 
							". Their tables are combined by rows ...")
					}
				}
			}else{
				# first check the data.frame
				if(!is.data.frame(x.predictors)){
					stop("Provided x.predictors is neither character nor data.frame !")
				}
				abund_table <- x.predictors
				if(! any(colnames(abund_table) %in% rownames(sampleinfo))){
					stop("Please make sure the column names in provided x.predictors are sample names!")
				}
			}
			# remove meaningless things
			abund_table %<>% {
				.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]
			}
			DataX <- abund_table %>% t() %>% as.data.frame(check.names = FALSE)
			message("Total feature numbers: ", ncol(DataX))
			
			# the parallel backend is not registered here permanently; it is registered temporarily in the functions that 
			# need it (cal_train, cal_feature_sel and cal_caretList) and the original backend is restored afterwards
			if(is.null(n.cores) || length(n.cores) != 1 || is.na(n.cores) || n.cores < 1){
				stop("Provided n.cores must be a positive integer !")
			}
			self$n.cores <- n.cores
			# use data_feature to make it easily remember and search
			self$data_feature <- DataX
			message("The feature table is stored in object$data_feature ...")
			message("The response variable is stored in object$data_response ...")
		},
		#' @description
		#' Split data for training and testing.
		#' 
		#' @param prop.train default 3/4; the ratio of the data used for the training.
		#' @return \code{data_train} and \code{data_test} in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_split(prop.train = 3/4)
		#' }
		cal_split = function(prop.train = 3/4){
			######################    DATA SPLIT: TRAIN and TEST
			message("Creating training set and testing set ...")
			data_response <- self$data_response
			if(self$type == "Classification"){
				data_response %<>% factor
			}
			data_feature <- self$data_feature

			data_all <- data.frame(Response = data_response, data_feature, check.names = FALSE)
			SplitData <- rsample::initial_split(data_all, prop = prop.train, strata = "Response")
			train_data <- rsample::training(SplitData)
			test_data <- rsample::testing(SplitData)
			message("Stratified sampling with the proportion of ", prop.train*100 ,"% for the training set ...")

			######################    DATA SPLIT end
			self$data_train <- train_data
			self$data_test <- test_data
			message("Training and testing data are stored in object$data_train and object$data_test respectively ...")
			invisible(self)
		},
		#' @description
		#' Pre-process (centering, scaling etc.) of features based on the caret::preProcess function. 
		#' 	 See \href{https://topepo.github.io/caret/pre-processing.html}{https://topepo.github.io/caret/pre-processing.html} for more details.
		#' 	 Some methods of \code{preProcess}, such as 'nzv', 'zv', 'pca' and 'ica', change the number of features. 
		#' 	 In that case the \code{data_train} and \code{data_test} are rebuilt with the new features and the 
		#' 	 corresponding information is printed, so that the feature names always match the data. 
		#' 	 The existing \code{data_preProcess} is reused when this function is called more than once; please assign 
		#' 	 \code{NULL} to \code{object$data_preProcess} if a new preprocess model is needed.
		#' 
		#' @param ... parameters pass to \code{preProcess} function of caret package.
		#' @return \code{data_preProcess}, \code{data_train} and \code{data_test} in the object.
		#' 	 \code{data_preProcess} is the return data generated by the \code{preProcess} function of caret package based on the training data.
		#' 	 \code{data_train} and \code{data_test} are preprocessed training and testing data based on the \code{data_preProcess}.
		#' @examples
		#' \dontrun{
		#' # "nzv" removes near zero variance predictors
		#' t1$cal_preProcess(method = c("center", "scale", "nzv"))
		#' }
		cal_preProcess = function(...){
			self <- private$check_training_data(self)
			data_train <- self$data_train
			
			if(is.null(self$data_preProcess)){
				preProcess_res <- caret::preProcess(data_train[, -1, drop = FALSE], ...)
				self$data_preProcess <- preProcess_res
				message("Preprocess model is stored in object$data_preProcess ...")
			}else{
				preProcess_res <- self$data_preProcess
				message("The existing preprocess model in object$data_preProcess is reused and the arguments provided in ", 
					"this call are ignored. Assign NULL to object$data_preProcess first if a new preprocess model is needed ...")
			}
			
			self$data_train <- private$apply_preProcess(preProcess_res, data_train, "training")
			message("Training data is preprocessed and reassigned to object$data_train ...")
			if(!is.null(self$data_test)){
				self$data_test <- private$apply_preProcess(preProcess_res, self$data_test, "testing")
				message("Testing data is preprocessed and reassigned to object$data_test ...")
			}
			invisible(self)
		},
		#' @description
		#' Perform feature selection.
		#' 	 See \href{https://topepo.github.io/caret/feature-selection-overview.html}{https://topepo.github.io/caret/feature-selection-overview.html} for more details.
		#' 	 An error is raised when no feature is selected, as the following model training can not be performed 
		#' 	 in that case. Please consider increasing \code{boruta.repetitions}, decreasing \code{boruta.pValue} or 
		#' 	 using another \code{x.predictors}.
		#' 
		#' @param boruta.maxRuns default 300; maximal number of importance source runs; passed to the \code{maxRuns} parameter in \code{Boruta} function of Boruta package.
		#'   Note that the Boruta function requires a value greater than 10.
		#' @param boruta.pValue default 0.01; p value passed to the pValue parameter in \code{Boruta} function of Boruta package.
		#' @param boruta.repetitions default 4; repetition runs for the feature selection.
		#' @param boruta.min.repetitions default 2; a feature is selected only when it is confirmed as important in at least 
		#'   this number of the \code{boruta.repetitions} independent runs. It must not be greater than \code{boruta.repetitions}.
		#' @param seed default NULL; the random seed used to make the repeated Boruta runs reproducible. When it is NULL, 
		#'   the random seeds are generated automatically and the result may change between the runs.
		#' @param ... parameters pass to \code{Boruta} function of Boruta package.
		#' @return optimized \code{data_train} and \code{data_test} in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
		#' }
		cal_feature_sel = function(
			boruta.maxRuns = 300,
			boruta.pValue = 0.01,
			boruta.repetitions = 4,
			boruta.min.repetitions = 2,
			seed = NULL,
			...
			){
			self <- private$check_training_data(self)
			
			data_input <- self$data_train
			data_x <- data_input[, -1, drop = FALSE]
			data_y <- data_input[, 1]

			if(self$type == "Classification"){
				data_y %<>% factor
			}
			
			###################### ----------------
			# checks of the input parameters
			if(!is.numeric(boruta.maxRuns) || length(boruta.maxRuns) != 1 || is.na(boruta.maxRuns) || boruta.maxRuns <= 10){
				stop("Provided boruta.maxRuns must be a single number greater than 10 (requirement of the Boruta function) !")
			}
			if(!is.numeric(boruta.repetitions) || length(boruta.repetitions) != 1 || is.na(boruta.repetitions) || boruta.repetitions < 1){
				stop("Provided boruta.repetitions must be a positive integer !")
			}
			boruta.repetitions <- as.integer(boruta.repetitions)
			if(!is.numeric(boruta.min.repetitions) || length(boruta.min.repetitions) != 1 || is.na(boruta.min.repetitions) || 
				boruta.min.repetitions < 1){
				stop("Provided boruta.min.repetitions must be a positive integer !")
			}
			boruta.min.repetitions <- as.integer(boruta.min.repetitions)
			if(boruta.repetitions < boruta.min.repetitions){
				stop("Provided boruta.repetitions (", boruta.repetitions, ") is smaller than boruta.min.repetitions (", 
					boruta.min.repetitions, ") ! A feature should be confirmed as important in at least boruta.min.repetitions runs, ", 
					"so no feature can be selected in this situation. Please increase boruta.repetitions or decrease ", 
					"boruta.min.repetitions ...")
			}
			###################### ----------------
			######################    BORUTA
			# one Boruta run; the environment is set to baseenv() so that the whole trans_classifier object is not 
			# serialized to the parallel workers
			boruta_run_once <- function(run_seed, x_data, y_data, maxRuns_value, pValue_value, extra_args){
				set.seed(run_seed)
				boruta.res <- do.call(Boruta::Boruta, c(list(x = x_data, y = y_data, 
					maxRuns = maxRuns_value, pValue = pValue_value), extra_args))
				boruta.stats <- data.frame(Boruta::attStats(boruta.res))
				rownames(boruta.stats[boruta.stats$decision == 'Confirmed', , drop = FALSE])
			}
			environment(boruta_run_once) <- baseenv()
			# a fixed seed makes the feature selection reproducible
			if(is.null(seed)){
				run_seeds <- sample.int(.Machine$integer.max, boruta.repetitions)
			}else{
				if(!is.numeric(seed) || length(seed) != 1 || is.na(seed)){
					stop("Provided seed must be a single number !")
				}
				run_seeds <- as.integer(seed) + seq_len(boruta.repetitions) - 1L
			}
			extra_args <- list(...)
			use_cores <- private$get_usable_cores(boruta.repetitions)
			message("Running Feature Selection (Boruta) based on the training data ...")
			if(use_cores > 1){
				# a PSOCK cluster is used instead of parallel::mclapply so that the parallel computation also works on Windows
				message("Performing ", boruta.repetitions, " repetitions on ", use_cores, " cores ...")
				boruta_cluster <- parallel::makeCluster(use_cores)
				on.exit(parallel::stopCluster(boruta_cluster), add = TRUE)
				boruta.list <- parallel::parLapply(boruta_cluster, run_seeds, boruta_run_once, 
					x_data = data_x, y_data = data_y, maxRuns_value = boruta.maxRuns, 
					pValue_value = boruta.pValue, extra_args = extra_args)
			}else{
				boruta.list <- lapply(run_seeds, boruta_run_once, x_data = data_x, y_data = data_y, 
					maxRuns_value = boruta.maxRuns, pValue_value = boruta.pValue, extra_args = extra_args)
			}

			boruta.final <- as.data.frame(table(unlist(boruta.list)))
			boruta.list.top <- as.character(boruta.final[which(boruta.final$Freq >= boruta.min.repetitions), 1])
			boruta.n.features <- length(unique(boruta.list.top))
			message("End of Feature Selection - Total of selected features = ", boruta.n.features)
			if(boruta.n.features < 1){
				stop("No feature was selected by Boruta ! A feature is selected only when it is confirmed as important in ", 
					"at least ", boruta.min.repetitions, " of the ", boruta.repetitions, " repetition run(s). Please consider ", 
					"increasing boruta.repetitions, decreasing boruta.pValue, increasing boruta.maxRuns, or using another ", 
					"x.predictors ...")
			}
			######################    BORUTA end
			###################### ----------------
			data_output <- data_input[, c(colnames(data_input)[1], boruta.list.top), drop = FALSE]
			self$data_train <- data_output
			
			if(is.null(self$data_test)){
				message("Selected features are reassigned to object$data_train ...")
			}else{
				data_input <- self$data_test
				data_output <- data_input[, c(colnames(data_input)[1], boruta.list.top), drop = FALSE]
				self$data_test <- data_output
				message("Selected features are reassigned to object$data_train and object$data_test ...")
			}
			invisible(self)
		},
		#' @description
		#' Control parameters for the following training. Please see \code{trainControl} function of caret package for details.
		#' 
		#' @param method default 'repeatedcv'; 'repeatedcv': Repeated k-Fold cross validation; 
		#' 	 see method parameter in \code{trainControl} function of \code{caret} package for available options.
		#' @param classProbs default TRUE; should class probabilities be computed for classification models?;
		#' 	 see classProbs parameter in \code{caret::trainControl} function.
		#' @param savePredictions default 'final'; see \code{savePredictions} parameter in \code{caret::trainControl} function.
		#'   'final' means that only the predictions of the final model with the best tuning parameters are saved, which is 
		#'   required by the \code{cal_ROC} function when \code{input = "train"}. TRUE (or "all") saves the predictions of 
		#'   all the tuning parameter combinations and is not recommended, as the repeated predictions of the same sample 
		#'   are not independent observations and would distort the ROC analysis.
		#' @param ... parameters pass to \code{trainControl} function of caret package.
		#' @return \code{trainControl} in the object.
		#' @examples
		#' \dontrun{
		#' t1$set_trainControl(method = 'repeatedcv')
		#' }
		set_trainControl = function(
			method = 'repeatedcv',
			classProbs = TRUE,
			savePredictions = "final",
			...
			){
			if(classProbs){
				if(self$type == "Regression"){
					message("classProbs is set to FALSE as the response variable is numeric (regression) ...")
					classProbs <- FALSE
				}
			}
			trainControl <- caret::trainControl(method = method,
								   classProbs = classProbs,
								   savePredictions = savePredictions,
								   ...)
			message('Generating trainControl setting stored in object$trainControl ...')
			self$trainControl <- trainControl
			invisible(self)
		},
		#' @description
		#' Run the model training. Please see \href{https://topepo.github.io/caret/available-models.html}{https://topepo.github.io/caret/available-models.html} for available models.
		#' 	 The mtry and ntree optimization is performed for the "rf" method in both the classification and the 
		#' 	 regression. The mean Accuracy of the resamples is used to select the best ntree for the classification and 
		#' 	 the mean RMSE for the regression. The selected mtry and ntree are stored in \code{object$train_params}.
		#' 
		#' @param method default "rf"; "rf": random forest; see method in \code{train} function of caret package for other options.
		#' 	  For method = "rf", the \code{tuneGrid} is set: \code{expand.grid(mtry = seq(from = 1, to = max.mtry))}
		#' @param max.mtry default 2; for method = "rf"; maximum mtry used in the \code{tuneGrid} to do hyperparameter tuning to optimize the model.
		#' @param ntree default 500; for method = "rf"; Number of trees to grow. 
		#' 	  The default 500 is same with the \code{ntree} parameter in \code{randomForest} function in randomForest package.
		#' 	  When it is a vector with more than one element, the function will try to optimize the model to select a best one, such as \code{c(100, 500, 1000)}.
		#' @param ... parameters pass to \code{caret::train} function.
		#' @return \code{res_train} in the object.
		#' @examples
		#' \dontrun{
		#' # random forest
		#' t1$cal_train(method = "rf")
		#' # Support Vector Machines with Radial Basis Function Kernel
		#' t1$cal_train(method = "svmRadial", tuneLength = 15)
		#' }
		cal_train = function(
			method = "rf",
			max.mtry = 2,
			ntree = 500,
			...
			){
			self <- private$check_training_data(self)
			train_data <- self$data_train
			
			trControl <- self$trainControl
			if(is.null(trControl)){
				trControl <- caret::trainControl()
			}
			# temporarily register the parallel backend; the original backend is restored when the function exits
			old_backend <- private$register_cores()
			on.exit(private$restore_cores(old_backend), add = TRUE)
			
			###################### ----------------
			if(method == "rf"){
				# Optimization of RF parameters
				message("Optimization of Random Forest parameters ...")
				modellist <- list()
				
				tuneGrid <- expand.grid(mtry = seq(from = 1, to = max.mtry))
				
				if(length(ntree) > 1){
					for (test_ntree in ntree){
						fit <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = test_ntree, ...)
						key <- toString(test_ntree)
						modellist[[key]] <- fit
					}
					# compare results; the mean accuracy is used for the classification and the mean RMSE for the regression
					results.tune1 <- caret::resamples(modellist)
					res.tune1 <- summary(results.tune1)
					stat_list <- res.tune1$statistics
					if(self$type == "Classification" && "Accuracy" %in% names(stat_list)){
						use_metric <- "Accuracy"
					}else{
						if("RMSE" %in% names(stat_list)){
							use_metric <- "RMSE"
						}else{
							use_metric <- names(stat_list)[1]
						}
					}
					metric_table <- as.data.frame(stat_list[[use_metric]])
					if(use_metric %in% c("RMSE", "MAE")){
						best_key <- rownames(metric_table)[which.min(metric_table$Mean)]
					}else{
						best_key <- rownames(metric_table)[which.max(metric_table$Mean)]
					}
					best_key <- best_key[1]
					fit <- modellist[[best_key]]
					ntree <- as.numeric(best_key)
					message("ntree used: ", ntree, " (selected by the mean ", use_metric, " of the resamples)")
				}else{
					fit <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = ntree, ...)
				}
				message("best mtry: ", fit$bestTune$mtry)
				######################Optimization of RF parameters end				

				tuneGrid <- expand.grid(.mtry=fit$bestTune$mtry)
				res_train <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = ntree, ...)
				self$train_params <- list(mtry = fit$bestTune$mtry, ntree = ntree)
			}else{
				res_train <- caret::train(Response ~ ., data = train_data, method = method, trControl = trControl, ...)
			}
			self$res_train <- res_train
			message('The training result is stored in object$res_train ...')
			self$train_method <- method
			invisible(self)
		},
		#' @description
		#' Get feature importance from the training model.
		#' @param rf_feature_sig default FALSE; whether calculate feature significance in 'rf' model using \code{rfPermute} package; 
		#'    only available for \code{method = "rf"} in \code{cal_train} function.
		#' @param ... parameters pass to \code{varImp} function of caret package. 
		#'    If \code{rf_feature_sig} is TRUE and \code{train_method} is "rf", the parameters will be passed to \code{rfPermute} function of rfPermute package.
		#' @return \code{res_feature_imp} in the object. One row for each predictor variable. The column(s) are different importance measures.
		#'   When \code{rf_feature_sig = FALSE}, the importance is calculated with the \code{varImp} function of caret package and 
		#'   the column is usually 'Overall' for both the classification and the regression. When \code{rf_feature_sig = TRUE}, 
		#'   the importance values generated by the \code{rfPermute} package are used, such as 'MeanDecreaseAccuracy' and 
		#'   'MeanDecreaseGini' for the classification, and the corresponding p values are stored in the columns with the 
		#'   '.pval' suffix.
		#' @examples
		#' \dontrun{
		#' t1$cal_feature_imp()
		#' }
		cal_feature_imp = function(rf_feature_sig = FALSE, ...){
			if(is.null(self$res_train)){
				stop("Please first run cal_train to train the model !")
			}
			if(rf_feature_sig && self$train_method != "rf"){
				message("The rf_feature_sig parameter is only available for the model trained with method = 'rf' ! ", 
					"The feature significance is not calculated and the importance is obtained with varImp ...")
				rf_feature_sig <- FALSE
			}
			if(self$train_method == "rf"){
				if(rf_feature_sig){
					train_data <- self$data_train
					# replace feature names with simplified character
					match_table <- data.frame(rawname = colnames(train_data)[2:ncol(train_data)], 
						replacename = paste0("r", 1:(ncol(train_data) - 1)), stringsAsFactors = FALSE)
					colnames(train_data)[2:ncol(train_data)] <- match_table$replacename
					if(is.null(self$train_params)){
						rfp_res <- rfPermute::rfPermute(Response ~ ., data = train_data, ...)
					}else{
						rfp_res <- rfPermute::rfPermute(Response ~ ., data = train_data, ntree = self$train_params[["ntree"]], mtry = self$train_params[["mtry"]], ...)
					}
					res_feature_imp <- rfPermute::importance(rfp_res, scale = TRUE) %>% as.data.frame(check.names = FALSE)
					rownames(res_feature_imp) <- match_table[match(rownames(res_feature_imp), match_table[, 2]), 1]
				}else{
					res_feature_imp <- private$get_varImp(...)
				}
			}else{
				res_feature_imp <- private$get_varImp(...)
			}
			self$res_feature_imp <- res_feature_imp
			message('The feature importance is stored in object$res_feature_imp ...')
			invisible(self)
		},
		#' @description
		#' Bar plot for feature importance.
		#' @param rf_sig_show default NULL; "MeanDecreaseAccuracy" (Default) or "MeanDecreaseGini" for random forest classification;
		#' 	  "\%IncMSE" (Default) or "IncNodePurity" for random forest regression;
		#' 	  Only available when \code{rf_feature_sig = TRUE} in function \code{cal_feature_imp}, 
		#' 	  which generate "MeanDecreaseGini" (and "MeanDecreaseAccuracy") or "\%IncMSE" (and "IncNodePurity") in the column names of \code{res_feature_imp};
		#' 	  Function can also generate "Significance" according to the p value.
		#' @param show_sig_group default FALSE; whether show the features with different significant groups;
		#' 	  Only available when "Significance" is found in the data.
		#' @param ... parameters pass to \code{plot_diff_bar} function of \code{trans_diff} package.
		#' @return \code{ggplot2} object.
		#' @examples
		#' \dontrun{
		#' t1$plot_feature_imp(use_number = 1:20, coord_flip = FALSE)
		#' }
		plot_feature_imp = function(rf_sig_show = NULL, show_sig_group = FALSE, ...){
			if(is.null(self$res_feature_imp)){
				message("The res_feature_imp is not found! It is necessary for the visualization! Call the cal_feature_imp function automatically with default settings ... ")
				self$cal_feature_imp()
			}
			tmp <- data.frame(Taxa = rownames(self$res_feature_imp), self$res_feature_imp, check.names = FALSE)
			tmp$Taxa %<>% gsub("`", "", ., fixed = TRUE) %>% gsub("\\.(.__)", "\\|\\1", .)
			if(! "Value" %in% colnames(tmp)){
				if("Overall" %in% colnames(tmp)){
					colnames(tmp)[colnames(tmp) == "Overall"] <- "Value"
				}else{
					if(any(c("MeanDecreaseAccuracy", "MeanDecreaseGini", "%IncMSE", "IncNodePurity") %in% colnames(tmp))){
						if(is.null(rf_sig_show)){
							if(self$type == "Classification"){
								rf_sig_show = "MeanDecreaseAccuracy"
							}else{
								rf_sig_show = "%IncMSE"
							}
						}else{
							if(self$type == "Classification"){
								if(! rf_sig_show %in% c("MeanDecreaseAccuracy", "MeanDecreaseGini")){
									stop("Provided rf_sig_show must be one of 'MeanDecreaseAccuracy' and 'MeanDecreaseGini'!")
								}
							}else{
								if(! rf_sig_show %in% c("%IncMSE", "IncNodePurity")){
									stop("Provided rf_sig_show must be one of '%IncMSE' and 'IncNodePurity'!")
								}
							}
						}
						colnames(tmp)[colnames(tmp) == rf_sig_show] <- "Value"
						colnames(tmp)[colnames(tmp) == paste0(rf_sig_show, ".pval")] <- "pvalue"
						tmp$Significance <- generate_p_siglabel(tmp$pvalue, nonsig = "ns")
						tmp %<>% .[, c("Taxa", "Value", "pvalue", "Significance")]
					}else{
						if(is.numeric(tmp[, 2])){
							colnames(tmp)[2] <- "Value"
						}else{
							stop("The res_feature_imp format can not be correctly recognized!")
						}
					}
				}
			}
			if(show_sig_group){
				if("Significance" %in% colnames(tmp)){
					tmp$Group <- tmp$Significance %>% factor(., levels = c("***", "**", "*", "ns"))
				}
			}
			
			suppressMessages(trans_diff_tmp <- trans_diff$new(dataset = NULL))
			trans_diff_tmp$res_diff <- tmp
			trans_diff_tmp$method <- "rf"
			g1 <- trans_diff_tmp$plot_diff_bar(...)
			if(is.null(rf_sig_show)){
				g1 <- g1 + ylab("Value")
			}else{
				g1 <- g1 + ylab(rf_sig_show)
			}
			g1
		},
		#' @description
		#' Run the prediction.
		#' 
		#' @param positive_class default NULL; see positive parameter in \code{confusionMatrix} function of caret package;
		#' If positive_class is NULL, use the first group in data as the positive class automatically.
		#' @return \code{res_predict} and the performance data stored in the object.
		#' 	  The \code{res_predict} is the predicted result for \code{data_test}.
		#' 	  For the classification, \code{res_confusion_fit} and \code{res_confusion_stats} are stored in the object.
		#' 	  For the regression, \code{res_regression_stats} with RMSE, Rsquared and MAE is stored in the object.
		#' 	  Several evaluation metrics in \code{res_confusion_fit} are defined as follows:
		#' 	 \deqn{Accuracy = \frac{TP + TN}{TP + TN + FP + FN}}
		#'   \deqn{Sensitivity = Recall = TPR = \frac{TP}{TP + FN}}
		#'   \deqn{Specificity = TNR = 1 - FPR = \frac{TN}{TN + FP}}
		#'   \deqn{Precision = \frac{TP}{TP + FP}}
		#' 	 \deqn{Prevalence = \frac{TP + FN}{TP + TN + FP + FN}}
		#' 	 \deqn{F1-Score = \frac{2 * Precision * Recall}{Precision + Recall}}
		#' 	 \deqn{Kappa = \frac{Accuracy - Pe}{1 - Pe}}
		#'   where TP is true positive; TN is ture negative; FP is false positive; and FN is false negative;
		#'   FPR is False Positive Rate; TPR is True Positive Rate; TNR is True Negative Rate;
		#'   Pe is the hypothetical probability of chance agreement on the classes for reference and prediction in the confusion matrix.
		#'   Accuracy represents the ratio of correct predictions.
		#'   Precision identifies how the model accurately predicted the positive classes.
		#'   Recall (sensitivity) measures the ratio of actual positives that are correctly identified by the model.
		#'   F1-score is the weighted average score of recall and precision. The value at 1 is the best performance and at 0 is the worst.
		#'   Prevalence represents how often positive events occurred.
		#'   Kappa identifies how well the model is predicting.
		#' @examples
		#' \dontrun{
		#' t1$cal_predict()
		#' }
		cal_predict = function(positive_class = NULL){
			###################### ----------------
			######################    Evaluation for the test set
			if(is.null(self$res_train)){
				stop("Please first run cal_train to train the model !")
			}
			fit.best <- self$res_train
			test_data <- self$data_test
			if(is.null(test_data)){
				stop("No testing data is found! Please first run cal_split function!")
			}

			fit.best.predict <- predict(fit.best, test_data[, -1, drop = FALSE])
			self$res_predict <- fit.best.predict
			message('The result of model prediction is stored in object$res_predict ...')

			######################    end: Evaluation for the test set
			###################### ----------------
			if(self$type == "Classification"){
				if (is.null(positive_class)){
					positive_class <- levels(test_data[, 1])[1]
				}
				positive_class.display <- self$data_MapNames %>% dplyr::filter(ClassNames %in% positive_class) %>% 
						dplyr::select(OriginalNames) %>% unique() %>% dplyr::pull()
				if(length(positive_class.display) < 1){
					positive_class.display <- positive_class
				}

				message('Calculating confusionMatrix with positive class = ', positive_class.display, " ...")

				confusion.fit.best <- caret::confusionMatrix(as.factor(fit.best.predict), 
											as.factor(test_data[,1]), 
											positive = positive_class)

				self$res_confusion_fit <- confusion.fit.best
				message('The result of confusionMatrix is stored in object$res_confusion_fit ...')
				# only the ratio statistics are displayed as percentages; Kappa and the p values are not ratios
				confusion.data.sts <- confusion.fit.best$overall
				Confusion.Sts <- data.frame("Overall Statistics" = private$format_confusion_stats(confusion.data.sts))
				rownames(Confusion.Sts) <- names(confusion.data.sts)
				self$res_confusion_stats <- Confusion.Sts
				message('The statistics of confusionMatrix is stored in object$res_confusion_stats ...')
				message('Model prediction Accuracy = ',Confusion.Sts$Overall.Statistics[1])
			}else{
				# the regression performance is evaluated with RMSE, Rsquared and MAE
				regression.stats <- caret::postResample(pred = fit.best.predict, obs = test_data[, 1])
				Confusion.Sts <- data.frame("Overall Statistics" = round(as.numeric(regression.stats), 4))
				rownames(Confusion.Sts) <- names(regression.stats)
				self$res_regression_stats <- Confusion.Sts
				message('The regression performance is stored in object$res_regression_stats ...')
				message('Model prediction RMSE = ', Confusion.Sts$Overall.Statistics[1], ", Rsquared = ", 
					Confusion.Sts$Overall.Statistics[2], ", MAE = ", Confusion.Sts$Overall.Statistics[3])
			}
			invisible(self)
		},
		#' @description
		#' Plot the cross-tabulation of observed and predicted classes with associated statistics based on the results of function \code{cal_predict}.
		#' 
		#' @param plot_confusion default TRUE; whether plot the confusion matrix.
		#' @param plot_statistics default TRUE; whether plot the statistics.
		#' @param as_ggplot default FALSE; whether return a \code{ggplot} object instead of the \code{gtable} object 
		#'   assembled by the gridExtra package. As the statistics table is plotted with the \code{tableGrob} function, 
		#'   only the confusion matrix is returned when \code{as_ggplot = TRUE}; the statistics table is available in 
		#'   \code{object$res_confusion_stats}.
		#' @return \code{gtable} object assembled by the gridExtra package when \code{as_ggplot = FALSE} (default), 
		#'   or \code{ggplot} object when \code{as_ggplot = TRUE}.
		#' @examples
		#' \dontrun{
		#' t1$plot_confusionMatrix()
		#' }
		plot_confusionMatrix = function(
			plot_confusion = TRUE, 
			plot_statistics = TRUE,
			as_ggplot = FALSE
			){
			if(self$type == "Regression"){
				stop("The function can only be available for the Classification !")
			}
			if(is.null(self$res_confusion_fit)){
				stop("Please first run cal_predict to get the prediction performance !")
			}
			if(!isTRUE(plot_confusion) && !isTRUE(plot_statistics)){
				stop("At least one of plot_confusion and plot_statistics must be TRUE !")
			}
			
			p1 <- ggplot(data = as.data.frame(self$res_confusion_fit$table) ,
					aes(x = Reference, y = Prediction)) +
				geom_tile(aes(fill = log(Freq)), colour = "white") +
				scale_fill_gradient(low = "white", high = "steelblue") +
				geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
				theme(legend.position = "none")

			Confusion.Sts <- self$res_confusion_stats

			p2 <- gridExtra::tableGrob(Confusion.Sts)
			if(as_ggplot){
				if(plot_statistics){
					message("as_ggplot = TRUE: only the confusion matrix is returned as a ggplot object; the statistics ", 
						"table is available in object$res_confusion_stats ...")
				}
				return(p1 + ggtitle("Confusion Matrix"))
			}
			if(plot_confusion == TRUE & plot_statistics == TRUE){
				p3 <- gridExtra::grid.arrange(p1, p2,nrow = 1, ncol = 2, 
					top=grid::textGrob("Confusion Matrix and Statistics",gp=grid::gpar(fontsize=15,font=0.5)))
			}
			if(plot_confusion == TRUE & plot_statistics == FALSE){
				p3 <- gridExtra::grid.arrange(p1,nrow = 1, ncol = 1, 
					top=grid::textGrob("Confusion Matrix",gp=grid::gpar(fontsize=15,font=0.5)))
			}
			if(plot_confusion == FALSE & plot_statistics == TRUE){
				p3 <- gridExtra::grid.arrange(p2,nrow = 1, ncol = 1, 
					top=grid::textGrob("Statistics",gp=grid::gpar(fontsize=15,font=0.5)))
			}
			p3
		},
		#' @description
		#' Get ROC (Receiver Operator Characteristic) curve data and the performance data.
		#' 	 When \code{input = "train"}, the predictions are obtained from \code{object$res_train$pred}, which requires 
		#' 	 \code{savePredictions} in \code{set_trainControl} to be 'final' (default). When the predictions of several 
		#' 	 tuning parameter combinations are saved (e.g. \code{savePredictions = TRUE}), only the first prediction of 
		#' 	 each sample is used, as the repeated predictions of a sample are not independent observations.
		#' 
		#' @param input default "pred"; 'pred' or 'train'; 'pred' represents using prediction results;
		#'   'train' represents using training results.
		#' @return a list \code{res_ROC} stored in the object. It has two tables: \code{res_roc} and \code{res_pr}. AUC: Area Under the ROC Curve.
		#'   For the definition of metrics, please refer to the return part of function \code{cal_predict}.
		#' @examples
		#' \dontrun{
		#' t1$cal_ROC()
		#' }
		cal_ROC = function(input = "pred"){
			if(self$type == "Regression"){
				stop("The function can only be available for the Classification !")
			}
			input <- match.arg(input, c("pred", "train"))
			if(is.null(self$res_train)){
				stop("Please first run cal_train to train the model !")
			}
			fit.best <- self$res_train
			train_method <- fit.best$method
			
			if(input == "pred"){
				test_data <- self$data_test
				if(is.null(test_data)){
					stop("No testing data is found! Please first run cal_split function!")
				}
				prediction_prob <- predict(fit.best, test_data[, -1, drop = FALSE] , type="prob")
				class_names <- levels(droplevels(test_data[, 1])) #drop because sometimes there is empty classes
				true_label <- test_data[, 1]
				prediction_prob <- prediction_prob[, class_names, drop = FALSE]
			}else{
				# use the prediction data in the training part
				if(is.null(fit.best$pred)){
					stop("No prediction of the training data is found! Please set savePredictions = 'final' in the ", 
						"set_trainControl function before cal_train ...")
				}
				class_names <- fit.best$levels
				if(!all(class_names %in% colnames(fit.best$pred))){
					stop("The class names are not found in the columns of object$res_train$pred ! Please check the training ", 
						"result and the classProbs parameter in set_trainControl ...")
				}
				pred_table <- fit.best$pred
				repeated_rows <- duplicated(pred_table$rowIndex)
				if(any(repeated_rows)){
					message("Multiple predictions are found for the same sample (usually caused by savePredictions = TRUE ", 
						"or 'all' with several tuning parameter combinations). Only the first prediction of each sample is ", 
						"used in the ROC analysis. Please consider set_trainControl(savePredictions = 'final') before ", 
						"cal_train for a more rigorous result ...")
					pred_table <- pred_table[!repeated_rows, , drop = FALSE]
				}
				prediction_prob <- pred_table[, class_names, drop = FALSE]
				true_label <- pred_table$obs
			}
			# use multiROC package
			label_df <- lapply(class_names, function(x){ifelse(true_label == x, 1, 0)}) %>% 
				do.call(cbind, .) %>%
				as.data.frame %>%
				`colnames<-`(paste0(class_names, "_true"))
			prob_df <- prediction_prob %>% `colnames<-`(paste0(colnames(.), "_pred_", train_method))
			use_df <- cbind(label_df, prob_df)

			roc_res <- multiROC::multi_roc(use_df, force_diag = TRUE)
			pr_res <- multiROC::multi_pr(use_df, force_diag = TRUE)

			plot_roc_df <- multiROC::plot_roc_data(roc_res)
			plot_pr_df <- multiROC::plot_pr_data(pr_res)

			# store the results
			res_ROC <- list()
			res_ROC$res_roc <- plot_roc_df
			res_ROC$res_pr <- plot_pr_df
			self$res_ROC <- res_ROC
			message('Specificity-sensitivity data is stored in object$res_ROC$res_roc ...')
			message('Recall-Precision is stored in object$res_ROC$res_pr ...')
			invisible(self)
		},
		#' @description
		#' Plot ROC curve.
		#' 
		#' @param plot_type default c("ROC", "PR")[1]; 'ROC' represents ROC (Receiver Operator Characteristic) curve; 
		#'   'PR' represents PR (Precision-Recall) curve.
		#' @param plot_group default "all"; 'all' represents all the classes in the model;
		#' 	 'add' represents all adding micro-average and macro-average results, see 
		#' 	 \href{https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html}{https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html};
		#' 	 other options should be one or more class names, same with the names in Group column of res_ROC$res_roc from cal_ROC function.
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors used in the plot.
		#' @param add_AUC default TRUE; whether add AUC in the legend.
		#' @param plot_method default FALSE; If TRUE, show the method in the legend though only one method is found.
		#' @param ... parameters pass to \code{geom_path} function of ggplot2 package.
		#' @return \code{ggplot2} object.
		#' @examples
		#' \dontrun{
		#' t1$plot_ROC(size = 1, alpha = 0.7)
		#' }
		plot_ROC = function(
			plot_type = c("ROC", "PR")[1],
			plot_group = "all",
			color_values = RColorBrewer::brewer.pal(8, "Dark2"), 			
			add_AUC = TRUE,
			plot_method = FALSE,
			...
			){
			if(is.null(self$res_ROC)){
				message("The res_ROC is not found! It is necessary for the visualization! Call the cal_ROC function automatically with default settings ... ")
				self$cal_ROC()
			}
			plot_type <- match.arg(plot_type, c("ROC", "PR"))
			if(plot_type == "ROC"){
				plot_data <- self$res_ROC$res_roc
			}else{
				plot_data <- self$res_ROC$res_pr
			}
			if(plot_group != "add"){
				if(plot_group == "all"){
					plot_data %<>% .[! .$Group %in% c("Micro", "Macro"), ]
				}else{
					if(!any(plot_data$Group %in% plot_group)){
						stop("Please input the correct plot_group !")
					}else{
						plot_data %<>% .[.$Group %in% plot_group, ]
					}
				}
			}
			if(add_AUC){
				plot_data$Group <- paste0(plot_data$Group, "\n AUC = ", round(plot_data$AUC, 2))
			}
		
			if(plot_type == "ROC"){
				# annotate is used instead of geom_segment, otherwise all the aesthetics have length 1 but the 
				# data has multiple rows, which triggers a warning of ggplot2 in each plot
				p <- ggplot(plot_data, aes(x = 1-Specificity, y = Sensitivity)) + 
					annotate("segment", x = 0, y = 0, xend = 1, yend = 1, colour = "grey", linetype = "dashed")
			}else{
				p <- ggplot(plot_data, aes(x = Recall, y = Precision))
			}
			if(length(unique(plot_data$Method)) > 1){
				p <- p + geom_path(aes(color = Group, linetype = Method), ...)
			}else{
				if(plot_method){
					p <- p + geom_path(aes(color = Group, linetype = Method), ...)
				}else{
					p <- p + geom_path(aes(color = Group), ...)
				}
			}
			p <- p + theme_bw() + 
				coord_equal() +
				xlim(0, 1) +
				ylim(0, 1) +
				scale_color_manual(values = color_values) +
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
				theme(legend.title = element_blank())
			p
		},
		#' @description
		#' Use \code{caretList} function of caretEnsemble package to run multiple models. For the available models, please run \code{names(getModelInfo())}.
		#' 
		#' @param ... parameters pass to \code{caretList} function of \code{caretEnsemble} package.
		#' @return \code{res_caretList_models} in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_caretList(methodList = c('rf', 'svmRadial'))
		#' }
		cal_caretList = function(...){
			if(!requireNamespace("caretEnsemble", quietly = TRUE)){
				stop("Please first install caretEnsemble package from CRAN!")
			}
			use_trainControl <- self$trainControl
			train_data <- self$data_train
			if(is.null(train_data)){
				message("No training data is found! The reason is function cal_split is not performed! Use all the samples for the training ...")
				train_data <- data.frame(Response = self$data_response, self$data_feature, check.names = FALSE)
				if(self$type == "Classification"){
					train_data$Response %<>% as.factor
				}
				self$data_train <- train_data
			}
			
			# temporarily register the parallel backend; the original backend is restored when the function exits
			old_backend <- private$register_cores()
			on.exit(private$restore_cores(old_backend), add = TRUE)
			models <- caretEnsemble::caretList(Response ~ ., data = train_data, trControl = use_trainControl, ...)
			self$res_caretList_models <- models
			message('Models are stored in object$res_caretList_models ...')
			invisible(self)
		},
		#' @description
		#' Use \code{resamples} function of caret package to collect the metric values based on the \code{res_caretList_models} data.
		#' 
		#' @param ... parameters pass to \code{resamples} function of \code{caret} package.
		#' @return \code{res_caretList_resamples} list and \code{res_caretList_resamples_reshaped} table in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_caretList_resamples()
		#' }
		cal_caretList_resamples = function(...){
			if(is.null(self$res_caretList_models)){
				stop("Please first run cal_caretList function!")
			}
			res_caretList_models <- self$res_caretList_models
			results <- caret::resamples(res_caretList_models, ...)
			self$res_caretList_resamples <- results
			message('Raw resamples results are stored in object$res_caretList_resamples ...')

			all_data <- results$values
			model_names <- colnames(all_data) %>% .[. != "Resample"] %>% gsub("~.*", "", .) %>% unique
			metric_names <- colnames(all_data) %>% .[. != "Resample"] %>% gsub(".*~", "", .) %>% unique
			reshaped_data <- lapply(model_names, function(x){
				extract_value <- all_data[, paste0(x, "~", metric_names), drop = FALSE]
				colnames(extract_value) <- metric_names
				data.frame(Model = x, extract_value)
			}) %>% do.call(rbind, .) %>% as.data.frame(check.names = FALSE)
			reshaped_data <- reshape2::melt(reshaped_data, id.vars = "Model")
			colnames(reshaped_data) <- c("Model", "Metric", "Value")
			self$res_caretList_resamples_reshaped <- reshaped_data
			message('Reshaped metric values are stored in object$res_caretList_resamples_reshaped ...')
			invisible(self)
		},
		#' @description
		#' Visualize the metric values based on the \code{res_caretList_resamples_reshaped} data.
		#' 
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the box.
		#' @param ... parameters pass to \code{geom_boxplot} function of \code{ggplot2} package.
		#' @return ggplot object.
		#' @examples
		#' \dontrun{
		#' t1$plot_caretList_resamples()
		#' }
		plot_caretList_resamples = function(color_values = RColorBrewer::brewer.pal(8, "Dark2"), ...){
			if(is.null(self$res_caretList_resamples_reshaped)){
				stop("Please first run cal_caretList_resamples function!")
			}
			reshaped_data <- self$res_caretList_resamples_reshaped
			
			p <- ggplot(reshaped_data, aes(x = Model, y = Value, colour = Model)) + 
				geom_boxplot(...) +
				scale_color_manual(values = color_values) +
				facet_grid(Metric ~ ., drop = TRUE, scale = "free", space = "fixed")
			p
		},
		#' @description
		#' Print the trans_classifier object.
		print = function() {
			cat("trans_classifier class:\n")
			cat("The type of the response variable:", self$type, "\n")
			if(is.null(self$data_feature)){
				cat("The feature table is not available ...\n")
			}else{
				cat("The feature number of the original input:", ncol(self$data_feature), "\n")
			}
			if(is.null(self$data_train)){
				cat("The training data is not available ...\n")
			}else{
				cat("The dimension of data_train:", nrow(self$data_train), "samples x ", 
					ncol(self$data_train) - 1, "features\n")
			}
			if(!is.null(self$data_test)){
				cat("The dimension of data_test:", nrow(self$data_test), "samples x ", 
					ncol(self$data_test) - 1, "features\n")
			}
			if(!is.null(self$train_method)){
				cat("The model is trained with the method:", self$train_method, "\n")
			}
			invisible(self)
		}
	),
	private = list(
		check_training_data = function(self){
			train_data <- self$data_train
			if(is.null(train_data)){
				message("No training data is found! The reason is that the cal_split function is not performed! Use all the samples for the training ...")
				train_data <- data.frame(Response = self$data_response, self$data_feature, check.names = FALSE)
				if(self$type == "Classification"){
					train_data$Response %<>% as.factor
				}
				self$data_train <- train_data
			}
			self
		},
		# apply the preprocess model and rebuild the whole table, as the number of the features may be changed by 
		# some methods (e.g. nzv and pca). A positional assignment such as data[, -1] <- predict(...) would silently 
		# recycle the columns and leave the old feature names in this case.
		apply_preProcess = function(preProcess_res, input_data, data_type){
			response_name <- colnames(input_data)[1]
			old_features <- colnames(input_data)[-1]
			pred_data <- predict(preProcess_res, newdata = input_data[, -1, drop = FALSE])
			pred_data <- as.data.frame(pred_data, check.names = FALSE)
			if(ncol(pred_data) < 1){
				stop("No feature is left after the preprocessing ! Please check the method parameter in the ", 
					"cal_preProcess function ...")
			}
			output_data <- data.frame(input_data[, 1, drop = FALSE], pred_data, check.names = FALSE)
			colnames(output_data)[1] <- response_name
			if(ncol(output_data) != ncol(input_data)){
				new_features <- colnames(output_data)[-1]
				message("The feature number in the ", data_type, " data is changed from ", ncol(input_data) - 1, " to ", 
					ncol(output_data) - 1, " after the preprocessing ...")
				if(all(grepl("^(PC|IC)[0-9]+$", new_features))){
					message("The new features are the principal or independent components: ", 
						paste(new_features, collapse = ", "), " ...")
				}else{
					removed_features <- setdiff(old_features, new_features)
					if(length(removed_features) > 0){
						message("Features removed by the preprocessing: ", paste(removed_features, collapse = ", "), " ...")
					}
				}
			}
			output_data
		},
		# register the parallel backend temporarily; the previous backend is returned so that it can be restored
		register_cores = function(){
			cores <- self$n.cores
			if(is.null(cores) || length(cores) != 1 || is.na(cores) || cores <= 1){
				return(NULL)
			}
			old_backend <- list(
				registered = tryCatch(foreach::getDoParRegistered(), error = function(e) FALSE),
				name = tryCatch(foreach::getDoParName(), error = function(e) NA_character_),
				workers = tryCatch(foreach::getDoParWorkers(), error = function(e) 1L)
			)
			doParallel::registerDoParallel(cores)
			message("Registering cores = ", cores, " for the parallel computation ...")
			old_backend
		},
		# restore the previous parallel backend
		restore_cores = function(old_backend){
			if(is.null(old_backend)){
				return(invisible(NULL))
			}
			if(isTRUE(old_backend$registered) && !is.na(old_backend$name) && 
				grepl("doParallel|doSnow|doSNOW|doMC", old_backend$name)){
				doParallel::registerDoParallel(old_backend$workers)
			}else{
				if(requireNamespace("foreach", quietly = TRUE)){
					foreach::registerDoSEQ()
				}
			}
			invisible(NULL)
		},
		# the number of the cores that can be used for a task with n_tasks independent runs
		get_usable_cores = function(n_tasks){
			cores <- self$n.cores
			if(is.null(cores) || length(cores) != 1 || is.na(cores) || cores < 2){
				return(1L)
			}
			as.integer(min(cores, n_tasks))
		},
		# the feature importance of some models (e.g. svmRadial) is only available for the train object, 
		# not for the finalModel; fall back to the train object in that case
		get_varImp = function(...){
			res_feature_imp <- tryCatch(caret::varImp(self$res_train$finalModel, ...), error = function(e) e)
			if(inherits(res_feature_imp, "error")){
				res_feature_imp <- tryCatch(caret::varImp(self$res_train, ...), error = function(e2){
					stop("The feature importance can not be calculated for the model trained with method = '", 
						self$train_method, "' ! Original error: ", conditionMessage(e))
				})
			}
			res_feature_imp
		},
		# only the ratio statistics are formatted as percentages; Kappa and the p values are not ratios
		format_confusion_stats = function(overall_stats){
			stat_names <- names(overall_stats)
			stat_values <- as.numeric(overall_stats)
			ratio_stats <- c("Accuracy", "AccuracyLower", "AccuracyUpper", "AccuracyNull")
			pvalue_stats <- c("AccuracyPValue", "McnemarPValue")
			vapply(seq_along(stat_names), function(i){
				use_value <- stat_values[i]
				if(is.na(use_value)){
					return("NA")
				}
				if(stat_names[i] %in% ratio_stats){
					return(paste0(round(use_value * 100, 2), "%"))
				}
				if(stat_names[i] %in% pvalue_stats){
					if(use_value < 0.001){
						return("< 0.001")
					}
					return(as.character(round(use_value, 3)))
				}
				as.character(round(use_value, 3))
			}, character(1))
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)

