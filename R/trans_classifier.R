#' @title 
#' Create trans_classifier object for machine-learning-based model prediction.
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
		#' Create the trans_classifier object.
		#' 
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param x.predictors default "all"; character string or data.frame; a character string represents selecting the corresponding data from microtable$taxa_abund; 
		#'   data.frame represents other customized input. See the following available options:
		#'   \describe{
		#'     \item{\strong{'all'}}{use all the taxa stored in microtable$taxa_abund}
		#'     \item{\strong{'Genus'}}{use Genus level table in microtable$taxa_abund, or other specific taxonomic rank, e.g. 'Phylum'}
		#'     \item{\strong{other input}}{must be a data.frame; It should have the same format with the data.frame in microtable$taxa_abund, i.e. rows are features; 
		#'       cols are samples with same names in sample_table}
		#'   }
		#' @param y.response default NULL; the response variable in sample_table.
		#' @param n.cores default 1; the CPU thread used.
		#' @return data_feature and data_response in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_classifier$new(
		#' 		dataset = dataset, 
		#' 		x.predictors = "Genus",
		#' 		y.response = "Group")
		#' }
		initialize = function(dataset = NULL,
				x.predictors = "all",
				y.response = NULL,
				n.cores = 1
			){
			check_microtable(dataset)
			sampleinfo <- dataset$sample_table
			if(is.null(y.response)){
				stop("No y.response provided!")
			}
			if(!y.response %in% colnames(sampleinfo)){
				stop("Input y.response must be dataset$sample_table!")
			}
			# parse y.response
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
				if (grepl("all", x.predictors, ignore.case = TRUE)) {
					abund_table <- do.call(rbind, unname(dataset$taxa_abund))
				}else{
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
		#' @param ... parameters pass to preProcess function of caret package.
		#' @return converted data_feature in the object.
		#' @examples
		#' \dontrun{
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
		#' @param boruta.maxRuns default 300; maximal number of importance source runs; passed to the maxRuns parameter in Boruta function of Boruta package.
		#' @param boruta.pValue default 0.01; p value passed to the pValue parameter in Boruta function of Boruta package.
		#' @param boruta.repetitions default 4; repetition runs for the feature selection.
		#' @param ... parameters pass to Boruta function of Boruta package.
		#' @return optimized data_feature in the object.
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
		#' @param prop.train default 3/4; the ratio of the dataset used for the training.
		#' @return data_train and data_test in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_split(prop.train = 3/4)
		#' }
		cal_split = function(prop.train = 3/4){
			###################### ----------------
			######################    DATA SPLIT: TRAIN and TEST
			######################
			message("Creating training set and testing set ...")
			data_response <- self$data_response
			if(self$type == "Classification"){
				data_response %<>% factor
			}
			DataX <- self$data_feature

			data_all <- data.frame(Response = data_response, DataX)
			SplitData <- rsample::initial_split(data_all, prop = prop.train, strata = "Response")
			train_data <- rsample::training(SplitData)
			test_data <- rsample::testing(SplitData)
			message("Stratified sampling with the proportion of ", prop.train*100 ,"% for the training set ...")

			###################### 
			######################    DATA SPLIT end
			###################### ----------------
			self$data_train <- train_data
			self$data_test <- test_data
			message("Training and testing data are stored in object$data_train and object$data_test respectively ...")
		},
		#' @description
		#' Control parameters for the following training. See trainControl function of caret package for details.
		#' 
		#' @param method default 'repeatedcv'; 'repeatedcv': Repeated k-Fold cross validation; 
		#' 	 see method parameter in \code{trainControl} function of \code{caret} package for available options.
		#' @param classProbs default TRUE; should class probabilities be computed for classification models?;
		#' 	 see classProbs parameter in \code{caret::trainControl} function.
		#' @param savePredictions default TRUE; see \code{savePredictions} parameter in \code{caret::trainControl} function.
		#' @param ... parameters pass to trainControl function of caret package.
		#' @return trainControl in the object.
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
		#' @param method default "rf"; "rf": random forest; see method in caret::train function for other options.
		#' @param max.mtry default 2; for method = "rf"; maximum mtry used for the tunegrid to do hyperparameter tuning to optimize the model.
		#' @param max.ntree default 200; for method = "rf"; maximum number of trees used to optimize the model.
		#' @param ... parameters pass to \code{caret::train} function.
		#' @return res_train in the object.
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
			max.ntree = 200,
			...
			){
			train_data <- self$data_train
			trControl <- self$trainControl
			
			###################### ----------------
			if(method == "rf" & self$type == "Classification"){
				# Optimization of RF parameters
				message("Optimization of Random Forest parameters ...")

				tunegrid <- expand.grid(.mtry=seq(from =1, to = max.mtry) )
				modellist<- list()
				
				for (ntree in c(100, max.ntree)) {
					fit <- caret::train(Response ~ ., data = train_data, method = method,
									  tuneGrid = tunegrid, trControl = trControl, ntree = ntree, ...)
					key <- toString(ntree)
					modellist[[key]] <- fit
				}
				# compare results
				results.tune1 <- caret::resamples(modellist)
				res.tune1 <- summary(results.tune1)
				res.tune1 <- as.data.frame(res.tune1$statistics$Accuracy)
				#summary(results)

				ntree <- as.numeric(rownames(res.tune1)[which(res.tune1$Mean == max(res.tune1$Mean))])[1]
				#tunegrid <- expand.grid(.mtry=seq(from = 1, to=4, by = 0.5))
				modellist <- list()

				fit <- caret::train(Response ~ ., data = train_data, method = method, 
					tuneGrid = tunegrid, 
					trControl = trControl, ntree = ntree)

				message("ntree used:", ntree)
				message("best mtry:", fit$bestTune$mtry)

				tunegrid <- expand.grid(.mtry=fit$bestTune$mtry)
				res_train <- caret::train(x=train_data[,2:ncol(train_data)], y = train_data[,1], method = method, 
										 tuneGrid = tunegrid, trControl = trControl, ntree = ntree, ...)

				######################Optimization of RF parameters end				
			}else{
				res_train <- caret::train(Response ~ ., data = train_data, method = method, trControl = trControl, ...)
			}
			self$res_train <- res_train
			message('The training result is stored in object$res_train ...')
			self$train_method <- method
		},
		#' @description
		#' Get feature importance from the training model.
		#' @param ... parameters pass to varImp function of caret package.
		#' @return res_feature_imp in the object. One row for each predictor variable. The column(s) are different importance measures.
		#'   For the method 'rf', it is MeanDecreaseGini (classification) or IncNodePurity (regression).
		#' @examples
		#' \dontrun{
		#' t1$cal_feature_imp()
		#' }
		cal_feature_imp = function(...){
			if(is.null(self$res_train)){
				stop("Please first run cal_train to train the model !")
			}
			res_feature_imp <- caret::varImp(self$res_train$finalModel, ...)
			
			self$res_feature_imp <- res_feature_imp
			message('The feature importance is stored in object$res_feature_imp ...')
		},
		#' @description
		#' Bar plot for feature importance.
		#' @param ... parameters pass to \code{plot_diff_bar} function of \code{trans_diff} package.
		#' @return ggplot2 object.
		#' @examples
		#' \dontrun{
		#' t1$plot_feature_imp(use_number = 1:20, coord_flip = FALSE)
		#' }
		plot_feature_imp = function(...){
			if(is.null(self$res_feature_imp)){
				stop("Please first run cal_feature_imp !")
			}
			tmp <- data.frame(Taxa = rownames(self$res_feature_imp), Value = self$res_feature_imp[, 1])
			tmp$Taxa %<>% gsub("\\.(.__)", "\\|\\1", .)
			
			suppressMessages(trans_diff_tmp <- trans_diff$new(dataset = NULL))
			trans_diff_tmp$res_diff <- tmp
			trans_diff_tmp$method <- "rf"
			trans_diff_tmp$plot_diff_bar(coord_flip = FALSE, ...)
		},
		#' @description
		#' Run the prediction.
		#' 
		#' @param positive_class default NULL; see positive parameter in confusionMatrix function of caret package;
		#' If positive_class is NULL, use the first group in data as the positive class automatically.
		#' @return res_predict, res_confusion_fit and res_confusion_stats stored in the object.
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
		#' Plot the cross-tabulation of observed and predicted classes with associated statistics.
		#' 
		#' @param plot_confusion default TRUE; whether plot the confusion matrix.
		#' @param plot_statistics default TRUE; whether plot the statistics.
		#' @return ggplot object.
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
		#' @return a list res_ROC stored in the object.
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
		#' @param ... parameters pass to geom_path function of ggplot2 package.
		#' @return ggplot2 object.
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
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)


