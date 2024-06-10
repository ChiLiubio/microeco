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
		#' @param x.predictors default "Genus"; character string or data.frame; a character string represents selecting the corresponding data from microtable$taxa_abund; 
		#'   data.frame denotes other customized input. See the following available options:
		#'   \describe{
		#'     \item{\strong{'Genus'}}{use Genus level table in microtable$taxa_abund, or other specific taxonomic rank, e.g. 'Phylum'}
		#'     \item{\strong{'all'}}{use all the taxa stored in microtable$taxa_abund}
		#'     \item{\strong{other input}}{must be a data.frame; It should have the same format with the data.frame in microtable$taxa_abund, i.e. rows are features; 
		#'       cols are samples with same names in sample_table}
		#'   }
		#' @param y.response default NULL; the response variable in \code{sample_table} of input \code{microtable} object.
		#' @param n.cores default 1; the CPU thread used.
		#' @return \code{data_feature} and \code{data_response} stored in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_classifier$new(
		#' 		dataset = dataset, 
		#' 		x.predictors = "Genus",
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

				ClassNames <- make.names(response_data, unique = F)
				MapNames <- data.frame(OriginalNames = response_data, ClassNames = ClassNames)
				if(!identical(MapNames$OriginalNames, MapNames$ClassNames)){
					message("Factor names are non-standard. A correction was made and the change map was saved in object$data_MapNames ...")
				}
				self$data_MapNames <- MapNames
				self$data_response <- ClassNames
			}
			# x.predictors must be character or data.frame
			if(is.character(x.predictors)){
				if(is.null(dataset$taxa_abund)){
					message("No taxa_abund found in the dataset. Calculate the relative abundance ...")
					dataset$cal_abund()
				}
				if (grepl("all", x.predictors, ignore.case = TRUE)) {
					abund_table <- do.call(rbind, unname(dataset$taxa_abund))
				}else{
					if(! x.predictors %in% names(dataset$taxa_abund)){
						message("x.predictors is not found in the taxa_abund list. Use the features in otu_table ...")
						dataset$add_rownames2taxonomy(use_name = x.predictors)
						suppressMessages(dataset$cal_abund())
					}
					abund_table <- dataset$taxa_abund[[x.predictors]]
				}
			}else{
				# first check the data.frame
				if(!is.data.frame(x.predictors)){
					stop("Provided x.predictors is neither character nor data.frame !")
				}
				abund_table <- x.predictors
				# maybe more checking later
			}
			# remove meaningless things
			abund_table %<>% {
				.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]
			}
			DataX <- abund_table %>% t() %>% as.data.frame()			
			message("Total feature numbers: ", ncol(DataX))
			
			if(n.cores > 1){
				message("Registering cores = ", n.cores)
				doParallel::registerDoParallel(n.cores)
			}
			# use data_feature to make it easily remember and search
			self$data_feature <- DataX
			message("The feature table is stored in object$data_feature ...")
			message("The response variable is stored in object$data_response ...")
		},
		#' @description
		#' Pre-process (centering, scaling etc.) of the feature data based on the caret::preProcess function. 
		#' 	 See \href{https://topepo.github.io/caret/pre-processing.html}{https://topepo.github.io/caret/pre-processing.html} for more details.
		#' 
		#' @param ... parameters pass to \code{preProcess} function of caret package.
		#' @return preprocessed \code{data_feature} in the object.
		#' @examples
		#' \dontrun{
		#' # "nzv" removes near zero variance predictors
		#' t1$cal_preProcess(method = c("center", "scale", "nzv"))
		#' }
		cal_preProcess = function(...){
			raw_feature <- self$data_feature
			preProcess_res <- caret::preProcess(raw_feature, ...)
			new_data <- predict(preProcess_res, newdata = raw_feature)
			self$data_feature <- new_data
			message("The converted feature table is stored in object$data_feature ...")
		},
		#' @description
		#' Perform feature selection.
		#' 	 See \href{https://topepo.github.io/caret/feature-selection-overview.html}{https://topepo.github.io/caret/feature-selection-overview.html} for more details.
		#' 
		#' @param boruta.maxRuns default 300; maximal number of importance source runs; passed to the \code{maxRuns} parameter in \code{Boruta} function of Boruta package.
		#' @param boruta.pValue default 0.01; p value passed to the pValue parameter in \code{Boruta} function of Boruta package.
		#' @param boruta.repetitions default 4; repetition runs for the feature selection.
		#' @param ... parameters pass to \code{Boruta} function of Boruta package.
		#' @return optimized \code{data_feature} in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
		#' }
		cal_feature_sel = function(
			boruta.maxRuns = 300,
			boruta.pValue = 0.01,
			boruta.repetitions = 4,
			...
			){
			# ClassNames
			data_response <- self$data_response
			if(self$type == "Classification"){
				data_response <- factor(data_response)
			}
			DataX <- self$data_feature
			
			###################### ----------------
			######################    BORUTA
			boruta.list <- list()
			boura.fs <- function(i){
				boruta.res <- Boruta::Boruta(x = DataX, y = data_response, 
					maxRuns = boruta.maxRuns, pValue = boruta.pValue, ...)
				boruta.stats <- data.frame(Boruta::attStats(boruta.res))
				boruta.list[[i]] <- rownames(boruta.stats[boruta.stats$decision =='Confirmed', ])
			}
			message("Running Feature Selection (Boruta) ...")
			boruta.list <- parallel::mclapply(1:boruta.repetitions, boura.fs)

			boruta.final <- as.data.frame(table(unlist(boruta.list)))
			#boruta.store.top <- as.character(boruta.store[which(boruta.store$Freq>10),1])
			boruta.list.top <- as.character(boruta.final[which(boruta.final$Freq >= 2), 1])
			boruta.n.features <- length(unique(boruta.list.top))
			message("End of Feature Selection - Total of selected features = ", boruta.n.features)
			######################    BORUTA end
			###################### ----------------
			self$data_feature <- DataX[, boruta.list.top]
			message("The selected features is reassigned to object$data_feature ...")
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
		},
		#' @description
		#' Control parameters for the following training. Please see \code{trainControl} function of caret package for details.
		#' 
		#' @param method default 'repeatedcv'; 'repeatedcv': Repeated k-Fold cross validation; 
		#' 	 see method parameter in \code{trainControl} function of \code{caret} package for available options.
		#' @param classProbs default TRUE; should class probabilities be computed for classification models?;
		#' 	 see classProbs parameter in \code{caret::trainControl} function.
		#' @param savePredictions default TRUE; see \code{savePredictions} parameter in \code{caret::trainControl} function.
		#' @param ... parameters pass to \code{trainControl} function of caret package.
		#' @return \code{trainControl} in the object.
		#' @examples
		#' \dontrun{
		#' t1$set_trainControl(method = 'repeatedcv')
		#' }
		set_trainControl = function(
			method = 'repeatedcv',
			classProbs = TRUE,
			savePredictions = TRUE,
			...
			){
			if(classProbs){
				if(self$type == "Regression"){
					classProbs <- FALSE
				}
			}
			trainControl <- caret::trainControl(method = method,
								   classProbs = classProbs,
								   savePredictions = savePredictions,
								   ...)
			message('Generating trainControl setting stored in object$trainControl ...')
			self$trainControl <- trainControl
		},
		#' @description
		#' Run the model training.
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
			train_data <- self$data_train
			trControl <- self$trainControl
			
			###################### ----------------
			if(method == "rf" & self$type == "Classification"){
				# Optimization of RF parameters
				message("Optimization of Random Forest parameters ...")
				modellist <- list()
				# capture the parameters
				all_parameters <- c(as.list(environment()), list(...))
				
				tuneGrid <- expand.grid(mtry = seq(from = 1, to = max.mtry))
				
				if(length(ntree) > 1){
					for (test_ntree in ntree){
						fit <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = test_ntree, ...)
						key <- toString(test_ntree)
						modellist[[key]] <- fit
					}
					# compare results
					results.tune1 <- caret::resamples(modellist)
					res.tune1 <- summary(results.tune1)
					res.tune1 <- as.data.frame(res.tune1$statistics$Accuracy)
					ntree <- as.numeric(rownames(res.tune1)[which(res.tune1$Mean == max(res.tune1$Mean))])[1]
					fit <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = ntree, ...)
					message("ntree used:", ntree)
				}else{
					fit <- caret::train(Response ~ ., data = train_data, method = method, tuneGrid = tuneGrid, trControl = trControl, ntree = ntree, ...)
				}
				message("best mtry:", fit$bestTune$mtry)
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
		},
		#' @description
		#' Get feature importance from the training model.
		#' @param rf_feature_sig default FALSE; whether calculate feature significance; only available for \code{method = "rf"} in \code{cal_train} function.
		#' @param ... parameters pass to \code{varImp} function of caret package. 
		#'    If \code{rf_feature_sig} is TURE and \code{train_method} is "rf", the parameters will be passed to \code{rfPermute} function of rfPermute package.
		#' @return \code{res_feature_imp} in the object. One row for each predictor variable. The column(s) are different importance measures.
		#'   For the method 'rf', it is MeanDecreaseGini (classification) or IncNodePurity (regression).
		#' @examples
		#' \dontrun{
		#' t1$cal_feature_imp()
		#' }
		cal_feature_imp = function(rf_feature_sig = FALSE, ...){
			if(is.null(self$res_train)){
				stop("Please first run cal_train to train the model !")
			}
			if(self$train_method == "rf"){
				if(rf_feature_sig){
					train_data <- self$data_train
					# replace feature names with simplified character
					match_table <- data.frame(rawname = colnames(train_data)[2:ncol(train_data)], replacename = paste0("r", 1:(ncol(train_data) - 1)))
					colnames(train_data)[2:ncol(train_data)] <- match_table$replacename
					rfp_res <- rfPermute::rfPermute(Response ~ ., data = train_data, ntree = self$train_params[["ntree"]], mtry = self$train_params[["mtry"]], ...)
					res_feature_imp <- rfPermute::importance(rfp_res, scale = TRUE) %>% as.data.frame(check.names = FALSE)
					rownames(res_feature_imp) <- match_table[match(rownames(res_feature_imp), match_table[, 2]), 1]
				}else{
					res_feature_imp <- caret::varImp(self$res_train$finalModel, ...)
				}
			}else{
				res_feature_imp <- caret::varImp(self$res_train$finalModel, ...)
			}
			self$res_feature_imp <- res_feature_imp
			message('The feature importance is stored in object$res_feature_imp ...')
		},
		#' @description
		#' Bar plot for feature importance.
		#' @param rf_sig_show default "MeanDecreaseGini"; "MeanDecreaseGini" or "MeanDecreaseAccuracy"; 
		#' 	  Only available when \code{rf_feature_sig = TRUE} in function \code{cal_feature_imp}, 
		#' 	  which generate "MeanDecreaseGini" and "MeanDecreaseAccuracy" in the column names of \code{res_feature_imp};
		#' 	  Function can also generate "Significance" according to the p value.
		#' @param show_sig_group default FALSE; whether show the features with different significant groups;
		#' 	  Only available when "Significance" is found in the data.
		#' @param ... parameters pass to \code{plot_diff_bar} function of \code{trans_diff} package.
		#' @return \code{ggplot2} object.
		#' @examples
		#' \dontrun{
		#' t1$plot_feature_imp(use_number = 1:20, coord_flip = FALSE)
		#' }
		plot_feature_imp = function(rf_sig_show = "MeanDecreaseGini", show_sig_group = FALSE, ...){
			if(is.null(self$res_feature_imp)){
				stop("Please first run function cal_feature_imp !")
			}
			tmp <- data.frame(Taxa = rownames(self$res_feature_imp), self$res_feature_imp)
			tmp$Taxa %<>% gsub("`", "", ., fixed = TRUE) %>% gsub("\\.(.__)", "\\|\\1", .)
			if(! "Value" %in% colnames(tmp)){
				if("Overall" %in% colnames(tmp)){
					colnames(tmp)[colnames(tmp) == "Overall"] <- "Value"
				}else{
					if(any(c("MeanDecreaseAccuracy", "MeanDecreaseGini") %in% colnames(tmp))){
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
			g1 <- trans_diff_tmp$plot_diff_bar(...) + ylab("Value")
			g1
		},
		#' @description
		#' Run the prediction.
		#' 
		#' @param positive_class default NULL; see positive parameter in \code{confusionMatrix} function of caret package;
		#' If positive_class is NULL, use the first group in data as the positive class automatically.
		#' @return \code{res_predict}, \code{res_confusion_fit} and \code{res_confusion_stats} stored in the object.
		#' 	  The \code{res_predict} is the predicted result for \code{data_test}.
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

			fit.best.predict <- predict(fit.best, test_data[, 2:ncol(test_data)])
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

				message('Calculating confusionMatrix with positive class = ', positive_class.display, " ...")

				confusion.fit.best <- caret::confusionMatrix(as.factor(fit.best.predict), 
											as.factor(test_data[,1]), 
											positive = positive_class)

				self$res_confusion_fit <- confusion.fit.best
				message('The result of confusionMatrix is stored in object$res_confusion_fit ...')
				confusion.data.sts <- data.frame(confusion.fit.best$overall)
				Confusion.Sts <- data.frame("Overall Statistics" = paste0(round(confusion.data.sts[,1],2) * 100,"%")  )
				rownames(Confusion.Sts) <- rownames(confusion.data.sts)
				self$res_confusion_stats <- Confusion.Sts
				message('The statistics of confusionMatrix is stored in object$res_confusion_stats ...')
				message('Model prediction Accuracy = ',Confusion.Sts$Overall.Statistics[1])
			}
		},
		#' @description
		#' Plot the cross-tabulation of observed and predicted classes with associated statistics based on the results of function \code{cal_predict}.
		#' 
		#' @param plot_confusion default TRUE; whether plot the confusion matrix.
		#' @param plot_statistics default TRUE; whether plot the statistics.
		#' @return \code{ggplot} object.
		#' @examples
		#' \dontrun{
		#' t1$plot_confusionMatrix()
		#' }
		plot_confusionMatrix = function(
			plot_confusion = TRUE, 
			plot_statistics = TRUE
			){
			if(self$type == "Regression"){
				stop("The function can only be available for the Classification !")
			}
			if(is.null(self$res_confusion_fit)){
				stop("Please first run cal_predict to get the prediction performance !")
			}
			
			p1 <- ggplot(data = as.data.frame(self$res_confusion_fit$table) ,
					aes(x = Reference, y = Prediction)) +
				geom_tile(aes(fill = log(Freq)), colour = "white") +
				scale_fill_gradient(low = "white", high = "steelblue") +
				geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
				theme(legend.position = "none")

			Confusion.Sts <- self$res_confusion_stats

			p2 <- gridExtra::tableGrob(Confusion.Sts)
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
				prediction_prob <- predict(fit.best, test_data[, 2:ncol(test_data)] , type="prob")
				class_names <- levels(droplevels(test_data[, 1])) #drop because sometimes there is empty classes
				true_label <- test_data[, 1]
			}else{
				# use the prediction data in the training part
				class_names <- fit.best$levels
				prediction_prob <- fit.best$pred[, class_names]
				true_label <- fit.best$pred$obs
			}
			# use multiROC package
			label_df <- lapply(class_names, function(x){ifelse(true_label == x, 1, 0)}) %>% 
				do.call(cbind, .) %>%
				as.data.frame %>%
				`colnames<-`(paste0(class_names, "_true"))
			prob_df <- prediction_prob %>% `colnames<-`(paste0(colnames(.), "_pred_", train_method))
			use_df <- cbind(label_df, prob_df)

			roc_res <- multiROC::multi_roc(use_df, force_diag = T)
			pr_res <- multiROC::multi_pr(use_df, force_diag = T)

			plot_roc_df <- multiROC::plot_roc_data(roc_res)
			plot_pr_df <- multiROC::plot_pr_data(pr_res)

			# store the results
			res_ROC <- list()
			res_ROC$res_roc <- plot_roc_df
			res_ROC$res_pr <- plot_pr_df
			self$res_ROC <- res_ROC
			message('Specificity-sensitivity data is stored in object$res_ROC$res_roc ...')
			message('Recall-Precision is stored in object$res_ROC$res_pr ...')
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
				stop("Please first run cal_ROC to get the data for ROC curve !")
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
				p <- ggplot(plot_data, aes(x = 1-Specificity, y = Sensitivity)) + 
					geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour = 'grey', linetype = 'dashed')
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
		#' @return \code{res_caretList_models} object.
		#' @examples
		#' \dontrun{
		#' t1$cal_caretList(methodList = c('rf', 'svmRadial'))
		#' }
		cal_caretList = function(
			...
			){
			if(!require(caretEnsemble)){
				stop("Please first install caretEnsemble package from CRAN!")
			}
			if(is.null(self$trainControl)){
				stop("Please first run set_trainControl function!")
			}
			use_trainControl <- self$trainControl
			train_data <- self$data_train
			
			models <- caretList(Response ~ ., data = train_data, trControl = use_trainControl, ...)
			self$res_caretList_models <- models
			message('Models are stored in object$res_caretList_models ...')
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)

