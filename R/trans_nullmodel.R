#' @title
#' Create trans_nullmodel object.
#'
#' @description
#' This class is a wrapper for a series of null model and phylogeny related approaches, including the mantel correlogram analysis of phylogenetic signal, betaNTI, betaNRI and RCbray calculations;
#' see Stegen et al. (2013) <10.1038/ismej.2013.93> and Liu et al. (2017) <doi:10.1038/s41598-017-17736-w>. 
#'
#' @export
trans_nullmodel <- R6Class(classname = "trans_nullmodel",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param filter_thres default 0; the relative abundance threshold. 
		#' @param taxa_number default NULL; how many taxa you want to use, if set, filter_thres parameter invalid.
		#' @param group default NULL; which group column name in sample_table is selected.
		#' @param select_group default NULL; the group name, used following the group to filter samples.
		#' @param env_cols default NULL; number or name vector to select the environmental data in dataset$sample_table. 
		#' @param add_data default NULL; provide environmental data table additionally.
		#' @param complete_na default FALSE; whether fill the NA in environmental data.
		#' @return intermediate files in object.
		#' @examples
		#' data(dataset)
		#' data(env_data_16S)
		#' t1 <- trans_nullmodel$new(dataset, taxa_number = 100, add_data = env_data_16S)
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
				use_set$tidy_dataset()
			}
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
				dis <- cophenetic(tre)
				dis <- dis[colnames(comm), colnames(comm)]
			}else{
				dis <- NULL
			}
			self$comm <- comm
			self$dis <- dis
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
		#' @param use_env default NULL; numeric or character vector to select env_data; if provide multiple variables or NULL, use PCA to reduce dimensionality.
		#' @param break.pts default seq(0, 1, 0.02); see \code{\link{mantel.correlog}}
		#' @param cutoff default FALSE; see cutoff in \code{\link{mantel.correlog}}
		#' @param ... parameters pass to \code{\link{mantel.correlog}}
		#' @return res_mantel_corr in object.
		#' @examples
		#' \donttest{
		#' t1$cal_mantel_corr(use_env = "pH")
		#' }
		cal_mantel_corr = function(use_env = NULL, break.pts = seq(0, 1, 0.02), cutoff=FALSE, ...){
			dis <- self$dis
			comm <- self$comm
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
			value_use <- decostand(value_use, method = "range", MARGIN = 2)
			niche_matrix <- t(comm) %*% value_use
			niche_matrix <- as.matrix(dist(niche_matrix))
			trenic_matrix <- as.matrix(dis)[rownames(niche_matrix), rownames(niche_matrix)]
			res_mantel_corr <- mantel.correlog(niche_matrix, trenic_matrix, break.pts = break.pts, cutoff = cutoff, ...)
			self$res_mantel_corr <- res_mantel_corr
			message('The result is stored in object$res_mantel_corr !')
		},
		#' @description
		#' Plot mantel correlogram.
		#'
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_mantel_corr()
		#' }
		plot_mantel_corr = function(){
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
			g <- ggplot(plot_data, aes(x=index, y=correlation, group = 1, fill=significance)) +
				theme_bw() +
				theme(panel.grid=element_blank()) +
				geom_line(linetype="dashed") +
				geom_point(shape=22, size=3) +
				geom_hline(aes(yintercept= 0), linetype="dotted") +
				scale_fill_manual(values=color_values) +
				guides(fill=FALSE) +
				scale_x_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
				ylab("Mantel correlation") +
				xlab("Phylogenetic distance") +
				theme(axis.text=element_text(size=13), axis.title=element_text(size=17)) +
				theme(legend.position="right", legend.text=element_text(size=rel(1))) +
				theme(legend.background=element_rect(fill="white", colour="grey60")) +
				guides(fill=guide_legend(title = NULL, reverse=FALSE))
			g
		},
		#' @description
		#' Calculate betaMPD. Faster than comdist in picante package.
		#'
		#' @param abundance.weighted default FALSE; whether use weighted abundance
		#' @return res_betampd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_betampd(abundance.weighted=FALSE)
		#' }
		cal_betampd = function(abundance.weighted=FALSE){
			dis <- self$dis
			if(is.null(dis)){
				stop("Phylogenetic tree is required!")
			}
			comm <- self$comm
			if (abundance.weighted == F) {
				comm <- decostand(comm, method="pa")
			}
			comm <- decostand(comm, method="total", MARGIN=1)
			self$res_betampd <- private$betampd(comm = comm, dis = dis)
			message('The result is stored in object$res_betampd !')
		},
		#' @description
		#' Calculate betaMNTD. Faster than comdistnt in picante package.
		#'
		#' @param abundance.weighted default FALSE; whether use weighted abundance
		#' @param exclude.conspecifics default FALSE; see comdistnt in picante package.
		#' @return res_betamntd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_betamntd(abundance.weighted=FALSE)
		#' }
		cal_betamntd = function(abundance.weighted = FALSE, exclude.conspecifics = FALSE){
			dis <- self$dis
			if(is.null(dis)){
				stop("Phylogenetic tree is required!")
			}
			comm <- self$comm
			self$res_betamntd <- private$betamntd(comm = comm, dis = dis, abundance.weighted=abundance.weighted, exclude.conspecifics = exclude.conspecifics)
			message('The result is stored in object$res_betamntd !')
		},
		#' @description
		#' Calculate ses.betaMPD (betaNRI).
		#'
		#' @param runs default 1000; simulation runs.
		#' @param abundance.weighted default FALSE; whether use weighted abundance.
		#' @param verbose default TRUE; whether show the calculation process message.
		#' @return res_ses_betampd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ses_betampd(runs = 100, abundance.weighted = FALSE)
		#' }
		cal_ses_betampd = function(runs=1000, abundance.weighted = FALSE, verbose = TRUE) {

			comm <- self$comm
			dis <- self$dis
			if(is.null(dis)){
				stop("Phylogenetic tree is required!")
			}
			if (abundance.weighted == F) {
				comm <- decostand(comm, method="pa")
			}
			comm <- decostand(comm, method="total", MARGIN=1)
			if(verbose){
				cat("Calculate observed betaMPD.\n")
			}
			betaobs <- private$betampd(comm = comm, dis = dis) %>% as.dist
			all_samples <- rownames(comm)
			betaobs_vec <- as.vector(betaobs)
			if(verbose){
				cat("Simulate betaMPD.\n")
			}
			beta_rand <- sapply(seq_len(runs), function(x){
				if(verbose){
					private$show_run(x = x, runs = runs)
				}
				as.dist(private$betampd(comm = comm, dis = picante::taxaShuffle(dis)))
			}, simplify = "array")
			beta_rand_mean <- apply(X = beta_rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
			beta_rand_sd <- apply(X = beta_rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - beta_rand_mean)/beta_rand_sd
			self$res_ses_betampd <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message('The result is stored in object$res_ses_betampd !')
		},
		#' @description
		#' Calculate ses.betaMNTD (betaNTI).
		#'
		#' @param runs default 1000; simulation runs.
		#' @param abundance.weighted default FALSE; whether use weighted abundance
		#' @param exclude.conspecifics default FALSE; see comdistnt in picante package.
		#' @param verbose default TRUE; whether show the calculation process message.
		#' @return res_ses_betamntd in object.
		#' @examples
		#' \donttest{
		#' t1$cal_ses_betamntd(runs = 100, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
		#' }
		cal_ses_betamntd = function(runs=1000, abundance.weighted = FALSE, exclude.conspecifics = FALSE, verbose = TRUE) {
			comm <- self$comm
			dis <- self$dis
			if(is.null(dis)){
				stop("Phylogenetic tree is required!")
			}
			all_samples <- rownames(comm)
			if(verbose){
				cat("Calculate observed betaMNTD.\n")
			}
			betaobs <- private$betamntd(comm = comm, dis = dis, abundance.weighted = abundance.weighted, exclude.conspecifics = exclude.conspecifics) %>% as.dist
			betaobs_vec <- as.vector(betaobs)
			if(verbose){
				cat("Simulate betaMNTD.\n")
			}
			beta_rand <- sapply(seq_len(runs), function(x){
				if(verbose){
					private$show_run(x = x, runs = runs)
				}
				as.dist(private$betamntd(comm = comm, dis = picante::taxaShuffle(dis), abundance.weighted = abundance.weighted, 
					exclude.conspecifics = exclude.conspecifics))
			}, simplify = "array")
			beta_rand_mean <- apply(X = beta_rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
			beta_rand_sd <- apply(X = beta_rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - beta_rand_mean)/beta_rand_sd
			self$res_ses_betamntd <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message('The result is stored in object$res_ses_betamntd !')
		},
		#' @description
		#' Calculate rcbray.
		#'
		#' @param runs default 1000; simulation runs.
		#' @param verbose default TRUE; whether show the calculation process message.
		#' @return res_rcbray in object.
		#' @examples
		#' \donttest{
		#' t1$cal_rcbray(runs=200)
		#' }
		cal_rcbray = function(runs=1000, verbose = TRUE) {
			comm <- self$comm
			betaobs_vec <- as.vector(vegdist(comm, method="bray"))
			all_samples <- rownames(comm)
			beta_rand <- sapply(seq_len(runs), function(x){
				if(verbose){
					private$show_run(x = x, runs = runs)
				}
				vegdist(picante::randomizeMatrix(comm, "independentswap"), "bray")
			}, simplify = "array") %>% as.data.frame
			beta_rand[, (runs + 1)] <- betaobs_vec
			beta_obs_z <- apply(X = beta_rand, MARGIN = 1, FUN = function(x){sum(x > x[length(x)])/length(x)})
			beta_obs_z <- (beta_obs_z - 0.5) * 2
			self$res_rcbray <- private$fin_matrix(all_samples = all_samples, beta_obs_z = beta_obs_z)
			message('The result is stored in object$res_rcbray !')
		},
		#' @description
		#' Infer the processes according to ses.betaMNTD ses.betaMPD and rcbray.
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
					stop("ses_betamntd not calculated!")
				}
			}else{
				ses_phylo_beta <- self$res_ses_betampd
				if(is.null(ses_phylo_beta)){
					stop("ses_betampd not calculated!")
				}
			}
			ses_comm <- self$res_rcbray
			if(is.null(ses_comm)){
				stop("RCbray not calculated!")
			}
			self$res_process <- private$percen_proc(ses_phylo_beta = ses_phylo_beta, ses_comm = ses_comm)
			message('The result is stored in object$res_process !')
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
			res1 <- rbind.data.frame(res, data.frame(S1 = res$S2, S2 = res$S1, distance = res$distance), 
				data.frame(S1 = all_samples, S2 = all_samples, distance = 0))
			res1 <- reshape2::dcast(res1, S1~S2, value.var = "distance") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE] %>%
				.[all_samples, all_samples] %>% as.matrix
			res1
		},
		betampd = function(comm = NULL, dis = NULL){
			all_samples <- rownames(comm)
			# use cpp instead of base
			matrix_multi <- function(comm_use, dis_use, ag_vector){(comm_use %*% dis_use) %*% ag_vector}
			# matrix_multi <- function(comm_use, dis_use, ag_vector){eigenMapMatMult(eigenMapMatMult(comm_use, dis_use), ag_vector)}
			
			res <- data.frame()
			rm_samples <- c()
			for(sample_name in all_samples[-length(all_samples)]){
				rm_samples <- c(rm_samples, sample_name)
				ag_vector <- comm[sample_name, , drop = FALSE] %>% t
				rownames(ag_vector) <- colnames(comm)
				ag_vector %<>% .[.[,1] != 0, , drop = FALSE]
				dis_use <- dis[, rownames(ag_vector), drop = FALSE]
				comm_use <- comm[!rownames(comm) %in% rm_samples, , drop = FALSE]
				wd <- matrix_multi(comm_use = comm_use, dis_use = dis_use, ag_vector = ag_vector)
				inter_res <- data.frame(S1 = rownames(comm_use), S2 = sample_name, distance = wd[, 1])
				res <- rbind.data.frame(res, inter_res)
			}
			res1 <- rbind.data.frame(res, data.frame(S1 = res$S2, S2 = res$S1, distance = res$distance), 
				data.frame(S1 = all_samples, S2 = all_samples, distance = 0))
			res1 <- reshape2::dcast(res1, S1~S2, value.var = "distance") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
			as.matrix(res1[all_samples, all_samples])
		},
		betamntd = function(comm = NULL, dis = NULL, abundance.weighted = FALSE, exclude.conspecifics = FALSE
			){
			all_samples <- rownames(comm)
			comm <- decostand(comm, method="total", MARGIN=1)
			comm <- rownames_to_column(as.data.frame(comm, stringsAsFactors = FALSE)) %>% reshape2::melt(id.vars = "rowname")
			colnames(comm) <- c("Sample", "Taxa", "Abund")
			comm %<>% .[.$Abund != 0, ]
			com_group <- combn(all_samples, 2)
			colnames(com_group) <- unlist(lapply(seq_len(ncol(com_group)), function(x) paste0(sort(com_group[, x]), collapse = "_betamntd_")))
			res <- data.frame()
			if(exclude.conspecifics == T){
				dis[dis == 0] <- NA
			}
			for(sample_name in all_samples){
				comm_use <- comm[comm$Sample != sample_name, ]
				dis_use <- dis[unique(as.character(comm_use$Taxa)), comm[comm$Sample == sample_name, "Taxa"]]
				each_taxa <- apply(dis_use, 1, min, na.rm=TRUE)
				comm_use$mindis <- each_taxa[as.character(comm_use$Taxa)]
				if(abundance.weighted == T){
					inter_res <- comm_use %>% dplyr::group_by(Sample) %>% dplyr::summarise(sum_dis = weighted.mean(mindis, Abund)) %>% as.data.frame		
				}else{
					inter_res <- comm_use %>% dplyr::group_by(Sample) %>% dplyr::summarise(sum_dis = sum(mindis), num = dplyr::n()) %>% as.data.frame		
				}
				inter_res$ag <- sample_name
				res <- rbind.data.frame(res, inter_res)
			}
			res$com_name <- unlist(lapply(seq_len(nrow(res)), function(x) {
				paste0(sort(unlist(res[x, c("Sample", "ag")])), collapse = "_betamntd_")
			}))
			if(abundance.weighted == T){
				res1 <- res %>% dplyr::group_by(com_name) %>% dplyr::summarise(distance = mean(sum_dis)) %>% as.data.frame
			}else{
				res1 <- res %>% dplyr::group_by(com_name) %>% dplyr::summarise(distance = sum(sum_dis)/sum(num)) %>% as.data.frame
			}
			res1 <- data.frame(t(com_group[, res1$com_name]), res1$distance)
			colnames(res1) <- c("S1", "S2", "distance")
			res1 <- rbind.data.frame(res1, data.frame(S1 = res1$S2, S2 = res1$S1, distance = res1$distance), 
				data.frame(S1 = all_samples, S2 = all_samples, distance = 0))
			res1 <- reshape2::dcast(res1, S1~S2, value.var = "distance") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
			as.matrix(res1[all_samples, all_samples])
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


