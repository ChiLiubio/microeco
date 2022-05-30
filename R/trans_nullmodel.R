#' @title
#' Create trans_nullmodel object for phylogeny- and taxonomy-based null model analysis.
#'
#' @description
#' This class is a wrapper for a series of null model related approaches, 
#' including the mantel correlogram analysis of phylogenetic signal, beta nearest taxon index (betaNTI), 
#' beta net relatedness index (betaNRI), NTI, NRI and RCbray calculations;
#' See Stegen et al. (2013) <10.1038/ismej.2013.93> and Liu et al. (2017) <doi:10.1038/s41598-017-17736-w> for the algorithms and applications.
#'
#' @export
trans_nullmodel <- R6Class(classname = "trans_nullmodel",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param filter_thres default 0; the relative abundance threshold. 
		#' @param taxa_number default NULL; how many taxa the user want to keep, if provided, filter_thres parameter will be forcible invalid.
		#' @param group default NULL; which group column name in sample_table is selected.
		#' @param select_group default NULL; the group name, used following the group to filter samples.
		#' @param env_cols default NULL; number or name vector to select the environmental data in dataset$sample_table. 
		#' @param add_data default NULL; provide environmental data table additionally.
		#' @param complete_na default FALSE; whether fill the NA in environmental data based on the method in mice package.
		#' @return data_comm and data_tree in object.
		#' @examples
		#' data(dataset)
		#' data(env_data_16S)
		#' t1 <- trans_nullmodel$new(dataset, filter_thres = 0.0005, add_data = env_data_16S)
		initialize = function(
			dataset = NULL,
			filter_thres = 0,
			taxa_number = NULL,
			group = NULL,
			select_group = NULL,
			env_cols = NULL,
			add_data = NULL,
			complete_na = FALSE
			){
			use_set <- clone(dataset)
			if(!is.null(group)){
				use_set$sample_table <- base::subset(use_set$sample_table, use_set$sample_table[, group] %in% select_group)
			}
			use_set$tidy_dataset()

			if(!is.null(taxa_number)){
				use_set$otu_table %<>% {.[names(sort(apply(., 1, sum), decreasing = TRUE)[1:taxa_number]), ]}
			}else{
				if(filter_thres != 0){
					use_abund <- use_set$otu_table
					use_set$otu_table <- use_abund[apply(use_abund, 1, sum)/sum(use_abund) > filter_thres, ]
				}
			}
			use_set$tidy_dataset()
			comm <- t(use_set$otu_table)
			if(!is.null(use_set$phylo_tree)){
				tre <- use_set$phylo_tree
			}else{
				tre <- NULL
			}
			self$data_comm <- comm
			self$data_tree <- tre
			self$sample_table <- use_set$sample_table
			if(!is.null(env_cols) | !is.null(add_data)){
				if(is.null(add_data)){
					env_data <- use_set$sample_table[, env_cols, drop = FALSE]
				}else{
					env_data <- add_data[rownames(add_data) %in% rownames(use_set$sample_table), , drop = FALSE]
				}
				env_data[env_data == ""] <- NA
				if(complete_na == T){
					env_data <- dropallfactors(env_data, unfac2num = TRUE)
					env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
					env_data %<>% mice::mice(print = FALSE) %>% mice::complete(., 1)
				}
				env_data <- dropallfactors(env_data, unfac2num = TRUE)
				env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.numeric(as.factor(x)) else x})
				self$env_data <- env_data
			}
		},
		#' @description
		#' Calculate mantel correlogram.
		#'
		#' @param use_env default NULL; numeric or character vector to select env_data; if provide multiple variables or NULL, 
		#' 	 use PCA (principal component analysis) to reduce dimensionality.
		#' @param break.pts default seq(0, 1, 0.02); see break.pts parameter in \code{\link{mantel.correlog}} of vegan package.
		#' @param cutoff default FALSE; see cutoff parameter in \code{\link{mantel.correlog}}.
		#' @param ... parameters pass to \code{\link{mantel.correlog}}
		#' @return res_mantel_corr in object.
		#' @examples
		#' \donttest{
		#' t1$cal_mantel_corr(use_env = "pH")
		#' }
		cal_mantel_corr = function(
			use_env = NULL, 
			break.pts = seq(0, 1, 0.02), 
			cutoff=FALSE, 
			...
			){
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}else{
				dis <- cophenetic(self$data_tree)
			}
			comm <- self$data_comm
			dis %<>% .[colnames(comm), colnames(comm)]
			env_data <- self$env_data
			comm %<>% .[rownames(env_data), ]
			comm <- apply(comm, 2, function(x) x/sum(x))
			# if use_env select multiple variables, use PCA to reduce dimensionality
			if(is.null(use_env)){
				if(ncol(env_data) > 1){
					value_use <- scores(rda(env_data), choices = 1)$sites
				}else{
					value_use <- as.matrix(env_data)
				}
			}else{
				if(length(use_env) > 1){
					value_use <- scores(rda(env_data[, use_env]), choices = 1)$sites
				}else{
					value_use <- as.matrix(env_data[, use_env, drop = FALSE])
				}
			}
			
			message("---------------- ", Sys.time()," : Start ----------------")
			value_use <- decostand(value_use, method = "range", MARGIN = 2)
			niche_matrix <- t(comm) %*% value_use
			niche_matrix <- as.matrix(dist(niche_matrix))
			trenic_matrix <- as.matrix(dis)[rownames(niche_matrix), rownames(niche_matrix)]
			res_mantel_corr <- mantel.correlog(
				niche_matrix, 
				trenic_matrix, 
				break.pts = break.pts, 
				cutoff = cutoff, 
				...
				)
			message("---------------- ", Sys.time()," : Finish ----------------")
			self$res_mantel_corr <- res_mantel_corr
			message('The result is stored in object$res_mantel_corr ...')
		},
		#' @description
		#' Plot mantel correlogram.
		#'
		#' @param point_shape default 22; the number for selecting point shape type; see ggplot2 manual for the number meaning.
		#' @param point_size default 3; the point size.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_mantel_corr()
		#' }
		plot_mantel_corr = function(point_shape = 22, point_size = 3){
			if(is.null(self$res_mantel_corr)){
				stop("Please first run cal_mantel_corr function to get data !")
			}
			plot_data <- self$res_mantel_corr$mantel.res %>% as.data.frame
			plot_data <- plot_data[, -4]
			colnames(plot_data) <- c("index", "n.dist", "correlation", "Pr")
			plot_data %<>% {.[!is.na(.$Pr), ]}
			plot_data$significance <- cut(plot_data$Pr, breaks=c(-Inf, 0.05, Inf), right=FALSE, label=c("P.adjust < 0.05", "P.adjust >= 0.05"))
			if(sum(plot_data$Pr < 0.05) > 0){
				color_values <- c("black", "white")
			}else{
				color_values <- c("white")
			}
			g <- ggplot(plot_data, aes(x = index, y = correlation, group = 1, fill = significance)) +
				theme_bw() +
				theme(panel.grid = element_blank()) +
				geom_line(linetype = "dashed") +
				geom_point(shape = point_shape, size = point_size) +
				geom_hline(aes(yintercept = 0), linetype = "dotted") +
				scale_fill_manual(values = color_values) +
				scale_x_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
				ylab("Mantel correlation") +
				xlab("Phylogenetic distance") +
				theme(axis.text = element_text(size = 13), axis.title = element_text(size = 17)) +
				theme(legend.position = "right", legend.text = element_text(size = rel(1))) +
				theme(legend.background = element_rect(fill="white", colour="grey60")) +
				guides(fill = guide_legend(title = NULL, reverse = FALSE))
			
			g
		},
		#' @description
		#' Calculate betaMPD (mean pairwise distance). Same with comdist in picante package, but faster.
		#'
		#' @param abundance.weighted default TRUE; whether use abundance-weighted method.
		#' @return res_betampd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_betampd(abundance.weighted = TRUE)
		#' }
		cal_betampd = function(abundance.weighted = TRUE){
			comm <- self$data_comm
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}else{
				dis <- cophenetic(self$data_tree) %>% .[colnames(comm), colnames(comm)]
			}
			self$res_betampd <- private$betampd(comm = comm, dis = dis, abundance.weighted = abundance.weighted)
			message('The result is stored in object$res_betampd ...')
		},
		#' @description
		#' Calculate betaMNTD (mean nearest taxon distance). Same with comdistnt in picante package, but faster.
		#'
		#' @param abundance.weighted default TRUE; whether use abundance-weighted method.
		#' @param exclude.conspecifics default FALSE; see exclude.conspecifics parameter in comdistnt function of picante package.
		#' @param use_iCAMP default FALSE; whether use bmntd.big function of iCAMP package to calculate betaMNTD. 
		#' 	  This method can store the phylogenetic distance matrix on the disk to lower the memory spending and perform the calculation parallelly.
		#' @param use_iCAMP_force default FALSE; whether use bmntd.big function of iCAMP package automatically when the feature number is large.
		#' @param iCAMP_tempdir default NULL; the temporary directory used to place the large tree file; If NULL; use the system user tempdir.
		#' @param ... paremeters pass to iCAMP::pdist.big function.
		#' @return res_betamntd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_betamntd(abundance.weighted = TRUE)
		#' }
		cal_betamntd = function(abundance.weighted = TRUE, exclude.conspecifics = FALSE, use_iCAMP = FALSE, 
			use_iCAMP_force = TRUE, iCAMP_tempdir = NULL, ...){
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}
			comm <- self$data_comm
			if(! use_iCAMP){
				if(ncol(comm) > 5000){
					if(use_iCAMP_force == T){
						use_iCAMP <- TRUE
						message("The feature number is larger than 5000. Automatically change use_iCAMP parameter to be TRUE and ",
							"use iCAMP package for large matrix and parallel computing. Change use_iCAMP_force = FALSE to skip this method ...")
					}
				}
			}
			if(use_iCAMP){
				tree <- self$data_tree
				if(is.null(iCAMP_tempdir)){
					iCAMP_tempdir <- tempdir()
				}
				message("The temporary directory used for big tree is in ", iCAMP_tempdir, " ...")
				if(!require("iCAMP")){
					stop("iCAMP package is not installed! Please first install it !")
				}
				pd.big <- iCAMP::pdist.big(tree = tree, wd = iCAMP_tempdir, ...)
				res_betamntd <- iCAMP::bmntd.big(comm = comm, pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
					abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics) %>%
					as.matrix
			}else{
				dis <- cophenetic(self$data_tree) %>% .[colnames(comm), colnames(comm)]
				res_betamntd <- private$betamntd(
					comm = comm, 
					dis = dis,
					abundance.weighted = abundance.weighted, 
					exclude.conspecifics = exclude.conspecifics
				)
			}
			self$res_betamntd <- res_betamntd
			message('The result is stored in object$res_betamntd ...')
		},
		#' @description
		#' Calculate standardized effect size of betaMPD, i.e. beta net relatedness index (betaNRI).
		#'
		#' @param runs default 1000; simulation runs.
		#' @param null.model default "taxa.labels"; The available options include "taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", 
		#' 	  "independentswap"and "trialswap"; see null.model parameter of ses.mntd function in picante package for the algorithm details.
		#' @param abundance.weighted default TRUE; whether use weighted abundance.
		#' @param iterations default 1000; iteration number for part null models to perform; see iterations parameter of picante::randomizeMatrix function.
		#' @return res_ses_betampd in object.
		#' @examples
		#' \donttest{
		#' # only run 50 times for the example; default 1000
		#' t1$cal_ses_betampd(runs = 50, abundance.weighted = TRUE)
		#' }
		cal_ses_betampd = function(
			runs = 1000, 
			null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")[1],
			abundance.weighted = TRUE,
			iterations = 1000
			){
			comm <- self$data_comm
			null.model <- match.arg(null.model, c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"))
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}else{
				dis <- cophenetic(self$data_tree) %>% .[colnames(comm), colnames(comm)]
			}
			message("---------------- ", Sys.time()," : Start ----------------")
			cat("Calculate observed betaMPD ...\n")
			betaobs <- private$betampd(comm = comm, dis = dis, abundance.weighted = abundance.weighted) %>% as.dist
			all_samples <- rownames(comm)
			betaobs_vec <- as.vector(betaobs)
			cat("Simulate betaMPD ...\n")
			beta_rand <- sapply(seq_len(runs), function(x){
				private$show_run(x = x, runs = runs)
				rand_data <- private$null_model(null.model = null.model, comm = comm, dis = dis, tip.label = NULL, iterations = iterations)
				as.dist(private$betampd(comm = rand_data$comm, dis = rand_data$dis, abundance.weighted = abundance.weighted))
			}, simplify = "array")
			message("---------------- ", Sys.time()," : End ----------------")
			
			beta_rand_mean <- apply(X = beta_rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
			beta_rand_sd <- apply(X = beta_rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - beta_rand_mean)/beta_rand_sd
			self$res_ses_betampd <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message('The result is stored in object$res_ses_betampd ...')
		},
		#' @description
		#' Calculate standardized effect size of betaMNTD, i.e. beta nearest taxon index (betaNTI).
		#'
		#' @param runs default 1000; simulation number of null model.
		#' @param null.model default "taxa.labels"; The available options include "taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", 
		#' 	  "independentswap"and "trialswap"; see null.model parameter of ses.mntd function in picante package for the algorithm details.
		#' @param abundance.weighted default TRUE; whether use abundance-weighted method.
		#' @param exclude.conspecifics default FALSE; see comdistnt in picante package.
		#' @param use_iCAMP default FALSE; whether use bmntd.big function of iCAMP package to calculate betaMNTD. 
		#' 	  This method can store the phylogenetic distance matrix on the disk to lower the memory spending and perform the calculation parallelly.
		#' @param use_iCAMP_force default FALSE; whether to make use_iCAMP to be TRUE when the feature number is large.
		#' @param iCAMP_tempdir default NULL; the temporary directory used to place the large tree file; If NULL; use the system user tempdir.
		#' @param nworker default 2; the CPU thread number.
		#' @param iterations default 1000; iteration number for part null models to perform; see iterations parameter of picante::randomizeMatrix function.
		#' @return res_ses_betamntd in object.
		#' @examples
		#' \donttest{
		#' # only run 50 times for the example; default 1000
		#' t1$cal_ses_betamntd(runs = 50, abundance.weighted = TRUE, exclude.conspecifics = FALSE)
		#' }
		cal_ses_betamntd = function(
			runs = 1000, 
			null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")[1],
			abundance.weighted = TRUE, 
			exclude.conspecifics = FALSE,
			use_iCAMP = FALSE, 
			use_iCAMP_force = TRUE, 
			iCAMP_tempdir = NULL, 
			nworker = 2,
			iterations = 1000
			){
			comm <- self$data_comm
			null.model <- match.arg(null.model, c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"))
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}
			if(! use_iCAMP){
				if(ncol(comm) > 5000){
					if(use_iCAMP_force == T){
						use_iCAMP <- TRUE
						message("The feature number is larger than 5000. Automatically change use_iCAMP parameter to be TRUE and ",
							"use iCAMP package for large matrix and parallel computing. Change use_iCAMP_force = FALSE to skip this method ...")
					}
				}
			}
			all_samples <- rownames(comm)
			
			message("---------------- ", Sys.time()," : Start ----------------")
			if(use_iCAMP){
				tree <- self$data_tree
				if(is.null(iCAMP_tempdir)){
					iCAMP_tempdir <- tempdir()
				}
				message("The temporary directory used for big tree is in ", iCAMP_tempdir, " ...")
				if(!require("iCAMP")){
					stop("iCAMP package is not installed! Please first install it !")
				}
				pd.big <- iCAMP::pdist.big(tree = tree, wd = iCAMP_tempdir, nworker = nworker)
				cat("Calculate observed betaMNTD ...\n")
				betaobs <- iCAMP::bmntd.big(comm=comm, pd.desc = pd.big$pd.file,
										 pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
										 abundance.weighted = abundance.weighted, exclude.consp = exclude.conspecifics)
				betaobs_vec <- as.vector(betaobs)
				cat("Simulate betaMNTD ...\n")
				beta_rand <- sapply(seq_len(runs), function(x){
					private$show_run(x = x, runs = runs)
					rand_data <- private$null_model(null.model = null.model, comm = comm, dis = NULL, tip.label = pd.big$tip.label, iterations = iterations)
					iCAMP::bmntd.big(comm = rand_data$comm, pd.desc = pd.big$pd.file,
						pd.spname = rand_data$tip.label, pd.wd = pd.big$pd.wd,
						abundance.weighted = abundance.weighted, exclude.consp = exclude.conspecifics)
				}, simplify = "array")
			}else{
				dis <- cophenetic(self$data_tree) %>% .[colnames(comm), colnames(comm)]
				cat("Calculate observed betaMNTD ...\n")
				betaobs <- private$betamntd(
					comm = comm, 
					dis = dis, 
					abundance.weighted = abundance.weighted, 
					exclude.conspecifics = exclude.conspecifics
					) %>% as.dist
				betaobs_vec <- as.vector(betaobs)
				cat("Simulate betaMNTD ...\n")
				beta_rand <- sapply(seq_len(runs), function(x){
					private$show_run(x = x, runs = runs)
					rand_data <- private$null_model(null.model = null.model, comm = comm, dis = dis, tip.label = NULL, iterations = iterations)
					as.dist(private$betamntd(
						comm = rand_data$comm, 
						dis = rand_data$dis, 
						abundance.weighted = abundance.weighted, 
						exclude.conspecifics = exclude.conspecifics
						)
					)
				}, simplify = "array")
			}
			beta_rand_mean <- apply(X = beta_rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
			beta_rand_sd <- apply(X = beta_rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - beta_rand_mean)/beta_rand_sd
			res_ses_betamntd <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message("---------------- ", Sys.time()," : Finish ----------------")
			self$res_ses_betamntd <- res_ses_betamntd
			message('The result is stored in object$res_ses_betamntd ...')
		},
		#' @description
		#' Calculate Bray–Curtis-based Raup–Crick (RCbray).
		#'
		#' @param runs default 1000; simulation runs.
		#' @param verbose default TRUE; whether show the calculation process message.
		#' @param null.model default "independentswap"; see more available options in randomizeMatrix function of picante package.
		#' @return res_rcbray in object.
		#' @examples
		#' \donttest{
		#' # only run 50 times for the example; default 1000
		#' t1$cal_rcbray(runs = 50)
		#' }
		cal_rcbray = function(runs = 1000, verbose = TRUE, null.model = "independentswap") {
			comm <- self$data_comm
			betaobs_vec <- as.vector(vegdist(comm, method="bray"))
			all_samples <- rownames(comm)
			beta_rand <- sapply(seq_len(runs), function(x){
				if(verbose){
					private$show_run(x = x, runs = runs)
				}
				vegdist(picante::randomizeMatrix(comm, null.model = null.model), "bray")
			}, simplify = "array") %>% as.data.frame
			beta_rand[, (runs + 1)] <- betaobs_vec
			beta_obs_z <- apply(X = beta_rand, MARGIN = 1, FUN = function(x){sum(x > x[length(x)])/length(x)})
			beta_obs_z <- (beta_obs_z - 0.5) * 2
			self$res_rcbray <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message('The result is stored in object$res_rcbray ...')
		},
		#' @description
		#' Infer the ecological processes according to ses.betaMNTD ses.betaMPD and rcbray.
		#'
		#' @param use_betamntd default TRUE; whether use ses.betaMNTD; if false, use ses.betaMPD.
		#' @return res_rcbray in object.
		#' @examples
		#' \donttest{
		#' t1$cal_process(use_betamntd = TRUE)
		#' }
		cal_process = function(use_betamntd = TRUE){
			if(use_betamntd == T){
				ses_phylo_beta <- self$res_ses_betamntd
				if(is.null(ses_phylo_beta)){
					stop("ses_betamntd not calculated! Please first run cal_ses_betamntd function!")
				}
			}else{
				ses_phylo_beta <- self$res_ses_betampd
				if(is.null(ses_phylo_beta)){
					stop("ses_betampd not calculated! Please first run cal_ses_betampd function!")
				}
			}
			ses_comm <- self$res_rcbray
			if(is.null(ses_comm)){
				stop("RCbray not calculated! Please first run cal_rcbray function!")
			}
			self$res_process <- private$percen_proc(ses_phylo_beta = ses_phylo_beta, ses_comm = ses_comm)
			message('The result is stored in object$res_process ...')
		},
		#' @description
		#' Calculates Nearest Relative Index (NRI), equivalent to -1 times the standardized effect size of MPD.
		#'
		#' @param null.model default "taxa.labels"; Null model to use; see null.model parameter in ses.mpd function of picante package for available options.
		#' @param abundance.weighted default FALSE; Should mean nearest relative distances for each species be weighted by species abundance?
		#' @param runs default 999; Number of randomizations.
		#' @param ... paremeters pass to ses.mpd function in picante package.
		#' @return res_NRI in object, equivalent to -1 times ses.mpd.
		#' @examples
		#' \donttest{
		#' # only run 50 times for the example; default 999
		#' t1$cal_NRI(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 50)
		#' }
		cal_NRI = function(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, ...){
			samp <- self$data_comm
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}
			dis <- cophenetic(self$data_tree) %>% .[colnames(samp), colnames(samp)]
			res <- picante::ses.mpd(samp, dis, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, ...)
			res$NRI <- res$mpd.obs.z * (-1)
			self$res_NRI <- res
			message('The result is stored in object$res_NRI ...')
		},
		#' @description
		#' Calculates Nearest Taxon Index (NTI), equivalent to -1 times the standardized effect size of MNTD.
		#'
		#' @param null.model default "taxa.labels"; Null model to use; see null.model parameter in ses.mntd function of picante package for available options.
		#' @param abundance.weighted default FALSE; Should mean nearest taxon distances for each species be weighted by species abundance?
		#' @param runs default 999; Number of randomizations.
		#' @param ... paremeters pass to ses.mntd function in picante package.
		#' @return res_NTI in object, equivalent to -1 times ses.mntd.
		#' @examples
		#' \donttest{
		#' # only run 50 times for the example; default 999
		#' t1$cal_NTI(null.model = "taxa.labels", abundance.weighted = TRUE, runs = 50)
		#' }
		cal_NTI = function(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, ...){
			samp <- self$data_comm
			if(is.null(self$data_tree)){
				stop("Phylogenetic tree is required! Please see the phylo_tree parameter of microtable class!")
			}
			dis <- cophenetic(self$data_tree) %>% .[colnames(samp), colnames(samp)]
			res <- picante::ses.mntd(samp, dis, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, ...)
			res$NTI <- res$mntd.obs.z * (-1)
			self$res_NTI <- res
			message('The result is stored in object$res_NTI ...')
		},
		#' @description
		#' Calculates the (normalised) mean number of checkerboard combinations (C-score) using C.score function in bipartite package.
		#'
		#' @param by_group default NULL; one column name or number in sample_table; calculate C-score for different groups separately.
		#' @param ... paremeters pass to C.score function in bipartite package.
		#' @return results directly.
		#' @examples
		#' \dontrun{
		#' t1$cal_Cscore()
		#' }
		cal_Cscore = function(by_group = NULL, ...){
			comm <- self$data_comm
			if(is.null(by_group)){
				bipartite::C.score(comm, ...)
			}else{
				sample_table <- self$sample_table
				lapply(unique(sample_table[, by_group]), function(x){
					use_comm <- comm[sample_table[, by_group] %in% x, ]
					use_comm %<>% .[, apply(., 2, sum) >0, drop = FALSE]
					bipartite::C.score(use_comm, ...)
				})
			}
		},
		#' @description
		#' Calculate normalized stochasticity ratio (NST) based on the tNST function of NST package.
		#'
		#' @param group a colname of sample_table; 
		#' 	  the function can select the data from sample_table to generate a one-column (n x 1) matrix and provide it to the group parameter of tNST function. 
		#' @param ... paremeters pass to tNST function of NST package; see the documents of tNST function for more details.
		#' @return .
		#' @examples
		#' \dontrun{
		#' t1$cal_tNST(group = "Group", dist.method = "bray", output.rand = TRUE, SES = TRUE)
		#' }
		cal_tNST = function(group, ...){
			comm <- self$data_comm
			group <- self$sample_table[, group, drop = FALSE]
			res <- NST::tNST(comm = comm, group = group, ...)
			self$res_tNST <- res
			message('The result is stored in object$res_tNST ...')
		},
		#' @description
		#' Test the significance of NST difference between each pair of groups.
		#'
		#' @param method default "nst.boot"; "nst.boot" or "nst.panova"; see NST::nst.boot function or NST::nst.panova function for the details.
		#' @param ... paremeters pass to NST::nst.boot when method = "nst.boot" or NST::nst.panova when method = "nst.panova"
		#' @return .
		#' @examples
		#' \dontrun{
		#' t1$cal_tNST_test()
		#' }
		cal_tNST_test = function(method = "nst.boot", ...){
			if(is.null(self$res_tNST)){
				stop("Please first run cal_tNST function!")
			}else{
				if(is.null(self$res_tNST$details)){
					stop("Please first run cal_tNST function with the parameter: output.rand = TRUE ")
				}
			}
			if(method == "nst.boot"){
				NST::nst.boot(nst.result = self$res_tNST, ...)
			}else{
				NST::nst.panova(nst.result = self$res_tNST, ...)
			}
		}
	),
	private = list(
		show_run = function(x, runs){
			if(x %% 10 == 0){
				cat(paste0("Runs: ", x, " of ", runs,"\n"))
			}
		},
		fin_matrix = function(all_samples, beta_obs_z){
			res <- data.frame(t(combn(all_samples, 2)), beta_obs_z)
			colnames(res) <- c("S1", "S2", "distance")
			res1 <- rbind.data.frame(
				res, 
				data.frame(S1 = res$S2, S2 = res$S1, distance = res$distance), 
				data.frame(S1 = all_samples, S2 = all_samples, distance = 0)
				)
			res1 <- reshape2::dcast(res1, S1~S2, value.var = "distance") %>% 
				`row.names<-`(.[,1]) %>% 
				.[, -1, drop = FALSE] %>%
				.[all_samples, all_samples] %>% 
				as.matrix
			
			res1
		},
		betampd = function(comm = NULL, dis = NULL, abundance.weighted = FALSE){
			dis %<>% .[colnames(comm), colnames(comm)]
			if (abundance.weighted == F) {
				comm <- decostand(comm, method = "pa")
			}
			comm <- decostand(comm, method="total", MARGIN=1)
			all_samples <- rownames(comm)
			# use cpp instead of base
			# matrix_multi <- function(comm_use, dis_use, ag_vector){eigenMapMatMult(eigenMapMatMult(comm_use, dis_use), ag_vector)}
			matrix_multi <- function(comm_use, dis_use, ag_vector){
				(comm_use %*% dis_use) %*% ag_vector
			}
			res <- data.frame()
			rm_samples <- c()
			for(sample_name in all_samples[-length(all_samples)]){
				rm_samples <- c(rm_samples, sample_name)
				ag_vector <- comm[sample_name, , drop = FALSE] %>% t
				rownames(ag_vector) <- colnames(comm)
				ag_vector %<>% .[.[, 1] != 0, , drop = FALSE]
				dis_use <- dis[, rownames(ag_vector), drop = FALSE]
				comm_use <- comm[!rownames(comm) %in% rm_samples, , drop = FALSE]
				wd <- matrix_multi(comm_use = comm_use, dis_use = dis_use, ag_vector = ag_vector)
				inter_res <- data.frame(S1 = rownames(comm_use), S2 = sample_name, distance = wd[, 1])
				res <- rbind.data.frame(res, inter_res)
			}
			res1 <- rbind.data.frame(
				res, 
				data.frame(S1 = res$S2, S2 = res$S1, distance = res$distance), 
				data.frame(S1 = all_samples, S2 = all_samples, distance = 0)
				)
			res1 <- reshape2::dcast(res1, S1 ~ S2, value.var = "distance") %>% 
				`row.names<-`(.[,1]) %>% 
				.[, -1, drop = FALSE]
			as.matrix(res1[all_samples, all_samples])
		},
		# from v0.7.5, use the method of iCAMP to calculate betamntd
		betamntd = function(
			comm = NULL, 
			dis = NULL, 
			abundance.weighted = FALSE, 
			exclude.conspecifics = FALSE
			){
			dis %<>% .[colnames(comm), colnames(comm)]
			if(!abundance.weighted){
				comm[comm > 0] <- 1
			}
			if(! inherits(dis, "matrix")){
				dis %<>% as.matrix
			}
			nd_matrix <- comm
			if(exclude.conspecifics){
				diag(dis) <- NA
				for(i in seq_len(nrow(comm))){
					id <- comm[i, ] == 0
					nd_matrix[i, ] <- apply(dis[!id, , drop=FALSE], 2, min, na.rm = TRUE)
				}
			}else{
				for(i in seq_len(nrow(comm))){
					id <- comm[i, ] == 0
					nd_matrix[i, !id] <- 0
					nd_matrix[i, id] <- apply(dis[!id, id, drop = FALSE], 2, min)
				}
			}
			if(abundance.weighted){
				comm_stand <- comm/rowSums(comm)
				res <- as.matrix(nd_matrix) %*% (t(comm_stand)) %>% 
					{(. + t(.))/2}
			}else{
				res_inter <- as.matrix(nd_matrix) %*% (t(comm))
				res <- rowSums(comm) %>%
					matrix(., nrow = nrow(comm), ncol = nrow(comm)) %>% 
					{. + t(.)} %>%
					{(res_inter + t(res_inter))/.}
			}
			diag(res) <- 0
			res
		},
		null_model = function(null.model, comm = NULL, dis = NULL, tip.label = NULL, iterations = 1000){
				if(is.null(comm)){
					stop("comm should not be NULL!")
				}
			if(null.model %in% c("taxa.labels", "phylogeny.pool")){
				if(is.null(dis)){
					if(is.null(tip.label)){
						stop("tip.label should not be NULL when null.model is 'taxa.labels' or 'phylogeny.pool'!")
					}
					tip.label <- sample(tip.label)
				}else{
					dis <- picante::taxaShuffle(dis)
				}
				if(null.model == "phylogeny.pool"){
					comm <- picante::randomizeMatrix(comm, null.model = "richness")
				}
			}else{
				if(null.model == "sample.pool"){
					comm <- picante::randomizeMatrix(comm, null.model = "richness", iterations = iterations)
				}else{
					comm <- picante::randomizeMatrix(comm, null.model = null.model, iterations = iterations)
				}
			}
			list(comm = comm, dis = dis, tip.label = tip.label)
		},
		percen_proc = function(ses_phylo_beta, ses_comm){
			phylo_vec <- as.vector(as.dist(ses_phylo_beta))
			com_vec <- as.vector(as.dist(ses_comm))
			type_use <- c("variable selection", "homogeneous selection", "dispersal limitation", "homogeneous dispersal", "drift")
			type_value <- c(sum(phylo_vec > 2), 
							sum(phylo_vec < -2),
							sum(com_vec > 0.95 & abs(phylo_vec) <= 2),
							sum(com_vec < -0.95 & abs(phylo_vec) <= 2),
							sum(abs(com_vec) <= 0.95 & abs(phylo_vec) <= 2))
			type_value %<>% {. * 100 / sum(.)}
			res <- data.frame(process = type_use, percentage = type_value)
			res
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
