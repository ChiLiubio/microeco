#' @title 
#' Create trans_classifier object for machine-learning-based model prediction.
#'
#' @description
#' This class is a wrapper for methods of machine-learning-based classification models.
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
		#' @param x.predictors default "all"; character or data.frame; character represents "all" or other specific taxonomic rank (such as "Genus" or "level_1");
		#'   "all" represent all the taxa stored in microtable$taxa_abund; data.frame represents customized data. If data.frame is provided, it should be have same
		#'   format with the data.frame in microtable$taxa_abund, i.e. rows are features; cols are samples with same names in sample_table.
		#' @param y.response default NULL; the response variable in sample_table.
		#' @param n.cores default 1; the thread used.
		#' @return data_feature and data_response in the object.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' t1 <- trans_classifier$new(
		#' 		dataset = dataset, 
		#' 		x.predictors = "all",
		#' 		y.response = "Group")
		#' }
		initialize = function(dataset = NULL,
				x.predictors = "all",
				y.response = NULL,
				n.cores = 1
			) {
			if (is.null(dataset)) {
				stop("No dataset provided!")
			}
			sampleinfo <- dataset$sample_table
			if(is.null(y.response)){
				stop("No y.response provided!")
			}
			if ( nlevels(as.factor(sampleinfo[, y.response])) < 2 ) {
				stop("Response variable must have at least 2 factors!")
			}
			# first judge whether x.predictors is character
			if(is.character(class(x.predictors))){
				if (grepl("all", x.predictors, ignore.case = TRUE)) {
					abund_table <- do.call(rbind, unname(dataset$taxa_abund))
				}else{
					abund_table <- dataset$taxa_abund[[x.predictors]]
				}
			}else{
				# first check the data.frame
				if(! is.data.frame(x.predictors)){
					stop("Provided x.predictors is neither character nor data.frame !")
				}
				abund_table <- x.predictors
				# maybe more checking later
			}
			# remove meaningless things
			abund_table %<>% {
				.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.)),]
			}
			
			if (nlevels(as.factor(sampleinfo[, y.response]))==2){
				ClassificationCase = "2 groups"
			}else{
				ClassificationCase = "2+ groups"
			}
			# self$ClassificationCase <- ClassificationCase
			message("Classification type = ", ClassificationCase)

			ClassNames <- make.names(sampleinfo[, y.response], unique = F)
			MapNames <- data.frame(OriginalNames = sampleinfo[, y.response], ClassNames = ClassNames)
			self$MapNames <- MapNames
			DataX <- abund_table %>% t() %>% as.data.frame()
			
			message("Total feature numbers: ", ncol(DataX))
			
			if(all.equal(MapNames$OriginalNames, MapNames$ClassNames) != TRUE){
				message("Factor names are non-standard. A correction was made and the change map was saved in object$MapNames")
			}
			# message("Start clasification for: ", y.response, " ...")
			if(n.cores > 1){
				message("Registering cores = ", n.cores)
				doParallel::registerDoParallel(n.cores)
			}
			# use data_feature to make it easily remember and search
			self$data_feature <- DataX
			message("The feature table is stored in object$data_feature ...")
			self$data_response <- y.response
			message("The response variable is stored in object$data_response ...")
		},
		#' @description
		#' Perform feature selection.
		#' 
		#' @param boruta.maxRuns default 300; maximal number of importance source runs; passed to the maxRuns parameter in Boruta function of Boruta package.
		#' @param boruta.pValue default 0.01; p value passed to the pValue parameter in Boruta function of Boruta package.
		#' @param boruta.repetitions default 4; repetition runs for the feature selection.
		#' @param ... parameters pass to Boruta function of Boruta package.
		#' @return optimized data_feature in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.01)
		#' }
		cal_feature_sel = function(
			boruta.maxRuns = 300,
			boruta.pValue = 0.01,
			boruta.repetitions = 4,
			...
			){
			# ClassNames
			ClassNames <- self$MapNames$ClassNames
			DataX <- self$data_feature
			
			###################### ----------------
			######################    BORUTA
			boruta.list <- list()
			boura.fs <- function(i){
				boruta.res <- Boruta::Boruta(x = DataX, y = factor(ClassNames), 
					maxRuns = boruta.maxRuns, pValue = boruta.pValue, ...)
				boruta.stats <- data.frame(Boruta::attStats(boruta.res))
				boruta.list[[i]] <- rownames(boruta.stats[boruta.stats$decision =='Confirmed',])
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
			
			# reassign
			self$data_feature <- DataX[,boruta.list.top]
			message("The selected features is reassigned to object$data_feature ...")
		},
		#' @description
		#' Split data for training and testing.
		#' 
		#' @param prop.train default 3/4; the ratio of the dataset used for the training.
		#' @return optimized data_feature in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_split(prop.train = 3/4)
		#' }
		cal_split = function(prop.train = 3/4){
			###################### ----------------
			######################    DATA SPLIT: TRAIN and TEST
			######################
			message("Creating training set and testing set ...")
			ClassNames <- self$MapNames$ClassNames
			DataX <- self$data_feature
			y.response <- self$data_response

			DataX.boruta <- data.frame(Class = factor(ClassNames), DataX)
			SplitData <- rsample::initial_split(DataX.boruta, prop = prop.train, strata = "Class")
			train_data <- rsample::training(SplitData)
			test_data <- rsample::testing(SplitData)
			message("Stratified sampling using the variable ", y.response,
			", with the proportion of ", prop.train*100 ,"% for the training set ...")

			###################### 
			######################    DATA SPLIT end
			###################### ----------------
			self$data_train <- train_data
			self$data_test <- test_data
			message("Training and testing data are stored in object$data_train and object$data_test respectively ...")
		},
		#' @description
		#' Set trainControl for the following training.
		#' 
		#' @param method default 'repeatedcv'; the method used in trainControl function of caret package.
		#' @param ... parameters pass to trainControl function of caret package.
		#' @return trainControl in the object.
		#' @examples
		#' \donttest{
		#' t1$set_trainControl()
		#' }
		set_trainControl = function(
			method = 'repeatedcv',
			...
			){
			trainControl <- caret::trainControl(method = method,
								   classProbs = TRUE,
								   savePredictions = TRUE,
								   ...)
			message('Generating trainControl setting stored in object$trainControl ...')
			self$trainControl <- trainControl
		},
		#' @description
		#' Run the training.
		#' 
		#' @param method default "rf"; representing the random forest method; see method in train function of caret package.
		#' @param metric default "Accuracy"; see metric in train function of caret package.
		#' @param max.mtry default 2; maximum mtry.
		#' @param max.ntree default 200; Number of trees to grow; pass to the ntree parameter of randomForest function in randomForest package.
		#' @param ... parameters pass to train function of caret package.
		#' @return res_train in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_train()
		#' }
		cal_train = function(
			method = "rf",
			metric = "Accuracy",
			max.mtry = 2,
			max.ntree = 200,
			...
			){
			train_data <- self$data_train
			control <- self$trainControl
			
			###################### ----------------
			######################
			if(method == "rf"){
				# Optimization of RF parameters
			  set.seed(12345)
				message("Optimization of Random Forest parameters ...")

				tunegrid <- expand.grid(.mtry=seq(from =1, to = max.mtry) )
				modellist<- list()
				#
				for (ntree in c(100, max.ntree)) {
					fit <- caret::train(Class~., data=train_data, method = method, metric=metric, 
									  tuneGrid=tunegrid, trControl=control, ntree=ntree, ...)
					key <- toString(ntree)
					modellist[[key]] <- fit
				}
				# compare results
				results.tune1 <- caret::resamples(modellist)
				res.tune1 <- summary(results.tune1)
				res.tune1 <- as.data.frame(res.tune1$statistics$Accuracy)
				#summary(results)

				ntree = as.numeric(rownames(res.tune1)[which(res.tune1$Mean == max(res.tune1$Mean))])[1]
				#tunegrid <- expand.grid(.mtry=seq(from = 1, to=4, by = 0.5))
				modellist <- list()

				fit <- caret::train(Class~., data=train_data, method = method, 
					metric=metric, tuneGrid=tunegrid, 
					trControl=control, ntree = ntree)

				message("ntree used:", ntree)
				message("best mtry:", fit$bestTune$mtry)

				tunegrid <- expand.grid(.mtry=fit$bestTune$mtry)
				fit.best <- caret::train(x=train_data[,2:ncol(train_data)], y=train_data[,1], method= method, 
										 metric=metric, tuneGrid=tunegrid, trControl=control, ntree=ntree, ...)
				self$res_train <- fit.best
				self$train_method <- method
				message('The training result is stored in object$res_train ...')
				
				######################Optimization of RF parameters end				
			}
			###################### ----------------
		},
		#' @description
		#' Get feature importance from the training model.
		#' @param ... parameters pass to the evaluating function; If "rf" used, pass to randomForest::importance.
		#' @return res_feature_imp in the object. One row for each predictor variable. The column(s) are different importance measures.
		#' @examples
		#' \donttest{
		#' t1$cal_feature_imp()
		#' }
		cal_feature_imp = function(...){
			if(self$train_method == "rf"){
				res_feature_imp <- randomForest::importance(self$res_train$finalModel, ...)
			}
			self$res_feature_imp <- res_feature_imp
			message('The feature importance evaluating result is stored in object$res_feature_imp ...')
		},
		#' @description
		#' Run the prediction.
		#' 
		#' @param positive_class default NULL; see positive parameter in confusionMatrix function of caret package.
		#' @return res_confusion_fit in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_predict()
		#' }
		cal_predict = function(positive_class = NULL){
			###################### ----------------
			######################    Evaluation for the test set
			######################
			fit.best <- self$res_train
			test_data <- self$data_test
			MapNames <- self$MapNames

			set.seed(12345)
			fit.best.predict <- predict(fit.best, test_data[, 2:ncol(test_data)])
			self$fit.best.predict <- fit.best.predict
			message('The result of model prediction is stored in object$fit.best.predict ...')

			###################### 
			######################    end: Evaluation for the test set
			###################### ----------------			
			if (is.null(positive_class)){
				positive_class <- levels(test_data[, 1])[1]
			}
			positive_class.display <- MapNames %>% dplyr::filter(ClassNames %in% positive_class) %>% 
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
			self$confusion_stats <- Confusion.Sts
			message('The statistics result of confusionMatrix is stored in object$confusion_stats ...')
			message('Model prediction Accuracy = ',Confusion.Sts$Overall.Statistics[1])
		},
		#' @description
		#' Plot the cross-tabulation of observed and predicted classes with associated statistics.
		#' 
		#' @param plot_confusion default TRUE; whether plot the confusion matrix.
		#' @param plot_statistics default TRUE; whether plot the statistics.
		#' @return ggplot object in the object.
		#' @examples
		#' \donttest{
		#' t1$plot_confusionMatrix()
		#' }
		plot_confusionMatrix = function(
			plot_confusion = TRUE, 
			plot_statistics = TRUE
			){

# 			color_values = RColorBrewer::brewer.pal(8, "Dark2")
			p1 <- ggplot(data = as.data.frame(self$res_confusion_fit$table) ,
					aes(x = Reference, y = Prediction)) +
				geom_tile(aes(fill = log(Freq)), colour = "white") +
				scale_fill_gradient(low = "white", high = "steelblue") +
				geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
				theme(legend.position = "none")# +

			confusion.data.sts <-data.frame(self$res_confusion_fit$overall)
			Confusion.Sts <- data.frame("Overall Statistics" = paste0(round(confusion.data.sts[,1],2) * 100,"%")  )
			rownames(Confusion.Sts) <- rownames(confusion.data.sts)

			p2 <- gridExtra::tableGrob(Confusion.Sts)

			if(plot_confusion == TRUE & plot_statistics == TRUE){
				gridExtra::grid.arrange(p1, p2,nrow = 1, ncol = 2, 
					top=grid::textGrob("Confusion Matrix and Statistics",gp=grid::gpar(fontsize=15,font=0.5)))
			}
			if(plot_confusion == TRUE & plot_statistics == FALSE){
				gridExtra::grid.arrange(p1,nrow = 1, ncol = 1, 
					top=grid::textGrob("Confusion Matrix",gp=grid::gpar(fontsize=15,font=0.5)))
			}
			if(plot_confusion == FALSE & plot_statistics == TRUE){
				gridExtra::grid.arrange(p2,nrow = 1, ncol = 1, 
					top=grid::textGrob("Statistics",gp=grid::gpar(fontsize=15,font=0.5)))
			}
		},
		#' @description
		#' Get ROC curve data and the performance data.
		#' 
		#' @param ... parameters pass to plot.performance function of ROCR package.
		#' @return a list including res_perf, res_auc_perf and all_perf_table stored in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_ROC()
		#' }
		cal_ROC = function(...){
			test_data <- self$data_test
			fit.best <- self$res_train

			prediction_for_roc_curve <- predict(fit.best, test_data[, 2:ncol(test_data)] , type="prob")
			classes <- levels(droplevels(test_data[, 1])) #drop because sometimes there is empty classes
			# store the results
			all_perf <- list()
			all_auc_perf <- list()
			all_perf_table <- data.frame()
			res_ROC <- list()
			
			for(i in 1:length(classes)){
			  # select observations of class[i]
			  true_values <- ifelse(test_data[, 1] == classes[i], 1, 0)
			  # performance for class[i]
			  pred <- ROCR::prediction(prediction_for_roc_curve[, i], true_values)
			  perf <- ROCR::performance(pred, "tpr", "fpr")
			  all_perf[[i]] <- perf
			  
			  gg_df <- data.frame(x = attributes(perf)$x.values[[1]], y = attributes(perf)$y.values[[1]], Group = classes[i])
			  all_perf_table <- rbind(all_perf_table, gg_df)
			  #print(classes[i])
			  auc.perf <- ROCR::performance(pred, measure = "auc")
			  all_auc_perf[[i]] <- auc.perf
			  #print(auc.perf@y.values)
			  message('AUC of ', classes[i], " = ", auc.perf@y.values)
			}
			res_ROC$res_perf <- all_perf
			res_ROC$all_auc_perf <- all_auc_perf
			res_ROC$all_perf_table <- all_perf_table
			self$res_ROC <- res_ROC
			message('Class performance of TPR-FPR is stored in object$res_ROC$res_perf ...')
			message('Class performance of AUC is stored in object$res_ROC$res_auc_perf ...')
			message('Class performance of AUC is stored in object$res_ROC$all_perf_table ...')
		},
		#' @description
		#' Plot ROC curve.
		#' 
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors used in the plot.
		#' @param ... parameters pass to geom_line function of ggplot2 package.
		#' @return ggplot2 object; res_perf and res_auc_perf stored in the object.
		#' @examples
		#' \donttest{
		#' t1$plot_ROC(size = 1, alpha = 0.7)
		#' }
		plot_ROC = function(color_values = RColorBrewer::brewer.pal(8, "Dark2"), ...){
			all_perf_table <- self$res_ROC$all_perf_table
			p <- ggplot(data = all_perf_table, aes(x = x, y = y, color = Group, group = Group)) +
				geom_line(...) + 
#				geom_ribbon(aes(x, ymin = 0, ymax = y)) +
				scale_color_manual(values = color_values) +
				xlab("False positive rate") +
				ylab("True positive rate") +
				theme_bw()
			p
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)


