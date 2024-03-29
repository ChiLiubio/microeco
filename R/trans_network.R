#' @title
#' Create \code{trans_network} object for network analysis.
#'
#' @description
#' This class is a wrapper for a series of network analysis methods, 
#' including the network construction, network attributes analysis,
#' eigengene analysis, network subsetting, node and edge properties, network visualization and other operations.
#'
#' @export
trans_network <- R6Class(classname = "trans_network",
	public = list(
		#' @description
		#' Create the \code{trans_network} object, store the important intermediate data 
		#'   and calculate correlations if \code{cor_method} parameter is not NULL.
		#' 
		#' @param dataset default NULL; the object of \code{\link{microtable}} class. Default NULL means customized analysis.
		#' @param cor_method default NULL; NULL or one of "bray", "pearson", "spearman", "sparcc", "bicor", "cclasso" and "ccrepe";
		#'   All the methods refered to \code{NetCoMi} package are performed based on \code{netConstruct} function of \code{NetCoMi} package and require
		#'   \code{NetCoMi} to be installed from Github (\href{https://github.com/stefpeschel/NetCoMi}{https://github.com/stefpeschel/NetCoMi});
		#'   For the algorithm details, please see Peschel et al. 2020 Brief. Bioinform <doi: 10.1093/bib/bbaa290>;
		#'   \describe{
		#'     \item{\strong{NULL}}{NULL denotes non-correlation network, i.e. do not use correlation-based network. 
		#'       If so, the return res_cor_p list will be NULL.}
		#'     \item{\strong{'bray'}}{1-B, where B is Bray-Curtis dissimilarity; based on \code{vegan::vegdist} function}
		#'     \item{\strong{'pearson'}}{Pearson correlation; If \code{use_WGCNA_pearson_spearman} and \code{use_NetCoMi_pearson_spearman} are both FALSE, 
		#'       use the function \code{cor.test} in R; \code{use_WGCNA_pearson_spearman = TRUE} invoke \code{corAndPvalue} function of \code{WGCNA} package; 
		#'       \code{use_NetCoMi_pearson_spearman = TRUE} invoke \code{netConstruct} function of \code{NetCoMi} package}
		#'     \item{\strong{'spearman'}}{Spearman correlation; other details are same with the 'pearson' option}
		#'     \item{\strong{'sparcc'}}{SparCC algorithm (Friedman & Alm, PLoS Comp Biol, 2012, <doi:10.1371/journal.pcbi.1002687>);
		#'     	 use NetCoMi package when \code{use_sparcc_method = "NetCoMi"}; use \code{SpiecEasi} package when \code{use_sparcc_method = "SpiecEasi"} 
		#'     	 and require \code{SpiecEasi} to be installed from Github
		#'     	 (\href{https://github.com/zdk123/SpiecEasi}{https://github.com/zdk123/SpiecEasi})}
		#'     \item{\strong{'bicor'}}{Calculate biweight midcorrelation efficiently for matrices based on \code{WGCNA::bicor} function; 
		#'       This option can invoke \code{netConstruct} function of \code{NetCoMi} package;
		#'       Make sure \code{WGCNA} and \code{NetCoMi} packages are both installed}
		#'     \item{\strong{'cclasso'}}{Correlation inference of Composition data through Lasso method based on \code{netConstruct} function of \code{NetCoMi} package; 
		#'     	 for details, see \code{NetCoMi::cclasso} function}
		#'     \item{\strong{'ccrepe'}}{Calculates compositionality-corrected p-values and q-values for compositional data 
		#'     	 using an arbitrary distance metric based on \code{NetCoMi::netConstruct} function; also see \code{NetCoMi::ccrepe} function}
		#'   }
		#' @param use_WGCNA_pearson_spearman default FALSE; whether use WGCNA package to calculate correlation when \code{cor_method} = "pearson" or "spearman".
		#' @param use_NetCoMi_pearson_spearman default FALSE; whether use NetCoMi package to calculate correlation when \code{cor_method} = "pearson" or "spearman".
		#'   The important difference between NetCoMi and others is the features of zero handling and data normalization; See <doi: 10.1093/bib/bbaa290>.
		#' @param use_sparcc_method default \code{c("NetCoMi", "SpiecEasi")[1]}; 
		#'   use \code{NetCoMi} package or \code{SpiecEasi} package to perform SparCC when \code{cor_method = "sparcc"}.
		#' @param taxa_level default "OTU"; taxonomic rank; 'OTU' denotes using feature abundance table; 
		#' 	  other available options should be one of the colnames of \code{tax_table} of input dataset.
		#' @param filter_thres default 0; the relative abundance threshold.
		#' @param nThreads default 1; the CPU thread number; available when \code{use_WGCNA_pearson_spearman = TRUE} or \code{use_sparcc_method = "SpiecEasi"}.
		#' @param SparCC_simu_num default 100; SparCC simulation number for bootstrap when \code{use_sparcc_method = "SpiecEasi"}.
		#' @param env_cols default NULL; numeric or character vector to select the column names of environmental data in dataset$sample_table;
		#'   the environmental data can be used in the correlation network (as the nodes) or \code{FlashWeave} network.
		#' @param add_data default NULL; provide environmental variable table additionally instead of \code{env_cols} parameter; rownames must be sample names.
		#' @param ... parameters pass to \code{NetCoMi::netConstruct} for other operations, such as zero handling and/or data normalization 
		#' 	 when cor_method and other parameters refer to \code{NetCoMi} package. 
		#' @return \code{res_cor_p} list with the correlation (association) matrix and p value matrix. Note that when \code{cor_method} and other parameters
		#'    refer to \code{NetCoMi} package, the p value table are all zero as the significant associations have been selected.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' # for correlation network
		#' t1 <- trans_network$new(dataset = dataset, cor_method = "pearson", 
		#' 		taxa_level = "OTU", filter_thres = 0.0002)
		#' # for non-correlation network
		#' t1 <- trans_network$new(dataset = dataset, cor_method = NULL)
		#' }
		initialize = function(
			dataset = NULL,
			cor_method = NULL,
			use_WGCNA_pearson_spearman = FALSE,
			use_NetCoMi_pearson_spearman = FALSE,
			use_sparcc_method = c("NetCoMi", "SpiecEasi")[1],
			taxa_level = "OTU",
			filter_thres = 0,
			nThreads = 1,
			SparCC_simu_num = 100,
			env_cols = NULL,
			add_data = NULL,
			...
			){
			if(is.null(dataset)){
				message("Input dataset not provided. Please run the functions with your other customized data!")
			}else{
				use_dataset <- clone(dataset)
				if(!is.null(env_cols)){
					env_data <- use_dataset$sample_table[, env_cols, drop = FALSE]
				}
				if(!is.null(add_data)){
					env_data <- add_data[rownames(add_data) %in% rownames(use_dataset$sample_table), ]
				}
				if(!is.null(env_cols) | !is.null(add_data)){
					use_dataset$sample_table %<>% .[rownames(.) %in% rownames(env_data), ]
					use_dataset$tidy_dataset(main_data = TRUE)
					env_data %<>% .[rownames(use_dataset$sample_table), ] %>%
						dropallfactors(unfac2num = TRUE)
					env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
					env_data[] <- lapply(env_data, as.numeric)
					self$data_env <- env_data
				}
				if(taxa_level != "OTU"){
					use_dataset <- use_dataset$merge_taxa(taxa = taxa_level)
				}
				rel_abund <- use_dataset$otu_table %>% {apply(., 1, sum)/sum(.)}
				use_abund <- use_dataset$otu_table %>% {.[rel_abund > filter_thres, ]}
				rel_abund %<>% {. * 100} %>% .[rownames(use_abund)]
				
				private$check_filter_number(use_abund, param = "filter_thres")
				use_abund %<>% t %>% as.data.frame
				
				if((!is.null(cor_method)) & (!is.null(env_cols) | !is.null(add_data))){
					use_abund <- cbind.data.frame(use_abund, env_data)
				}
				if(!is.null(cor_method)){
					cor_method <- match.arg(cor_method, c("bray", "pearson", "spearman", "sparcc", "bicor", "cclasso", "ccrepe"))
					if(cor_method == "bray"){
						tmp <- vegan::vegdist(t(use_abund), method = "bray") %>% as.matrix
						tmp <- 1 - tmp
						cor_result <- private$get_cor_p_list(tmp)
					}
					if(cor_method %in% c("pearson", "spearman")){
						if(use_NetCoMi_pearson_spearman){
							private$check_NetCoMi()
							netConstruct_raw <- netConstruct(data = use_abund, measure = cor_method, ...)
							cor_result <- private$get_cor_p_list(netConstruct_raw$assoMat1)
						}else{
							if(use_WGCNA_pearson_spearman){
								cor_result <- WGCNA::corAndPvalue(x = use_abund, method = cor_method, nThreads = nThreads)
							}else{
								cor_result <- private$cal_corr(inputtable = use_abund, cor_method = cor_method)
							}
						}
					}
					if(cor_method == "sparcc"){
						use_sparcc_method <- match.arg(use_sparcc_method, c("NetCoMi", "SpiecEasi"))
						if(use_sparcc_method == "NetCoMi"){
							private$check_NetCoMi()
							netConstruct_raw <- netConstruct(data = use_abund, measure = cor_method, ...)
							cor_result <- private$get_cor_p_list(netConstruct_raw$assoMat1)
						}else{
							try_find <- try(find.package("SpiecEasi"), silent = TRUE)
							if(inherits(try_find, "try-error")){
								stop("SpiecEasi package is used for the SparCC calculation, but it is not installed! See https://github.com/zdk123/SpiecEasi for the installation")
							}
							bootres <- SpiecEasi::sparccboot(use_abund, ncpus = nThreads, R = SparCC_simu_num)
							cor_result <- SpiecEasi::pval.sparccboot(bootres)
							use_names <- colnames(bootres$data)
							com_res <- t(combn(use_names, 2))
							res <- cbind.data.frame(com_res, cor = cor_result$cors, p = cor_result$pvals, stringsAsFactors = FALSE)
							res_cor <- private$vec2mat(datatable = res, use_names = use_names, value_var = "cor", rep_value = 1)
							res_p <- private$vec2mat(datatable = res, use_names = use_names, value_var = "p", rep_value = 0)
							cor_result <- list(cor = res_cor, p = res_p)
						}
					}
					if(cor_method %in% c("bicor", "cclasso", "ccrepe")){
						private$check_NetCoMi()
						netConstruct_raw <- netConstruct(data = use_abund, measure = cor_method, ...)
						cor_result <- private$get_cor_p_list(netConstruct_raw$assoMat1)
					}
					self$res_cor_p <- cor_result
					message('The correlation result list is stored in object$res_cor_p ...')
				}else{
					self$res_cor_p <- NULL
				}
				self$sample_table <- use_dataset$sample_table
				# store taxonomic table for the following analysis
				self$tax_table <- use_dataset$tax_table
				self$data_abund <- use_abund
				self$data_relabund <- rel_abund
			}
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Construct network based on the \code{igraph} package or \code{SpiecEasi} package or \code{julia FlashWeave} package or \code{beemStatic} package.
		#'
		#' @param network_method default "COR"; "COR", "SpiecEasi", "gcoda", "FlashWeave" or "beemStatic"; 
		#'   \code{network_method = NULL} means skipping the network construction for the customized use.
		#'   The option details: 
		#'   \describe{
		#'     \item{\strong{'COR'}}{correlation-based network; use the correlation and p value matrices in \code{res_cor_p} list stored in the object; 
		#'     	  See Deng et al. (2012) <doi:10.1186/1471-2105-13-113> for other details}
		#'     \item{\strong{'SpiecEasi'}}{\code{SpiecEasi} network; relies on algorithms of sparse neighborhood and inverse covariance selection;
		#'     	  belong to the category of conditional dependence and graphical models;
		#'     	  see \href{https://github.com/zdk123/SpiecEasi}{https://github.com/zdk123/SpiecEasi} for installing the R package; 
		#'     	  see Kurtz et al. (2015) <doi:10.1371/journal.pcbi.1004226> for the algorithm details}
		#'     \item{\strong{'gcoda'}}{hypothesize the logistic normal distribution of microbiome data; use penalized maximum likelihood method to estimate
		#'     	  the sparse structure of inverse covariance for latent normal variables to address the high dimensionality of the microbiome data;
		#'     	  belong to the category of conditional dependence and graphical models;
		#'     	  depend on the R \code{NetCoMi} package \href{https://github.com/stefpeschel/NetCoMi}{https://github.com/stefpeschel/NetCoMi}; 
		#'     	  see FANG et al. (2017) <doi:10.1089/cmb.2017.0054> for the algorithm details}
		#'     \item{\strong{'FlashWeave'}}{\code{FlashWeave} network; Local-to-global learning framework; belong to the category of conditional dependence and graphical models;
		#'        good performance on heterogenous datasets to find direct associations among taxa;
		#'        see \href{https://github.com/meringlab/FlashWeave.jl}{https://github.com/meringlab/FlashWeave.jl} for installing \code{julia} language and 
		#'        \code{FlashWeave} package; julia must be in the computer system env path, otherwise the program can not find it;
		#'        see Tackmann et al. (2019) <doi:10.1016/j.cels.2019.08.002> for the algorithm details}
		#'     \item{\strong{'beemStatic'}}{\code{beemStatic} network;
		#'        extend generalized Lotka-Volterra model to cases of cross-sectional datasets to infer interaction among taxa based on expectation-maximization algorithm;
		#'        see \href{https://github.com/CSB5/BEEM-static}{https://github.com/CSB5/BEEM-static} for installing the R package;
		#'        see Li et al. (2021) <doi:10.1371/journal.pcbi.1009343> for the algorithm details}
		#'   }
		#' @param COR_p_thres default 0.01; the p value threshold for the correlation-based network.
		#' @param COR_p_adjust default "fdr"; p value adjustment method, see \code{method} parameter of \code{p.adjust} function for available options,
		#' 	  in which \code{COR_p_adjust = "none"} means giving up the p value adjustment.
		#' @param COR_weight default TRUE; whether use correlation coefficient as the weight of edges; FALSE represents weight = 1 for all edges.
		#' @param COR_cut default 0.6; correlation coefficient threshold for the correlation network.
		#' @param COR_optimization default FALSE; whether use random matrix theory (RMT) based method to determine the correlation coefficient; 
		#' 	  see https://doi.org/10.1186/1471-2105-13-113
		#' @param COR_optimization_low_high default \code{c(0.01, 0.8)}; the low and high value threshold used for the RMT optimization; only useful when COR_optimization = TRUE.
		#' @param COR_optimization_seq default 0.01; the interval of correlation coefficient used for RMT optimization; only useful when COR_optimization = TRUE.
		#' @param SpiecEasi_method default "mb"; either 'glasso' or 'mb';see spiec.easi function in package SpiecEasi and https://github.com/zdk123/SpiecEasi.
		#' @param FlashWeave_tempdir default NULL; The temporary directory used to save the temporary files for running FlashWeave; If not assigned, use the system user temp.
		#' @param FlashWeave_meta_data default FALSE; whether use env data for the optimization, If TRUE, the function automatically find the \code{env_data} in the object and
		#'   generate a file for meta_data_path parameter of FlashWeave package.
		#' @param FlashWeave_other_para default \code{"alpha=0.01,sensitive=true,heterogeneous=true"}; the parameters passed to julia FlashWeave package;
		#'   user can change the parameters or add more according to FlashWeave help document;
		#'   An exception is meta_data_path parameter as it is generated based on the data inside the object, see FlashWeave_meta_data parameter for the description.
		#' @param beemStatic_t_strength default 0.001; for network_method = "beemStatic"; the threshold used to limit the number of interactions (strength);
		#'   same with the t.strength parameter in showInteraction function of beemStatic package.
		#' @param beemStatic_t_stab default 0.8; for network_method = "beemStatic"; 
		#'   the threshold used to limit the number of interactions (stability); same with the t.stab parameter in showInteraction function of beemStatic package.
		#' @param add_taxa_name default "Phylum"; one or more taxonomic rank name; used to add taxonomic rank name to network node properties.
		#' @param delete_unlinked_nodes default TRUE; whether delete the nodes without any link.
		#' @param usename_rawtaxa_when_taxalevel_notOTU default FALSE; whether use OTU name as representatives of taxa when \code{taxa_level != "OTU"}.
		#'   Default \code{FALSE} means using taxonomic information of \code{taxa_level} instead of OTU name.
		#' @param ... parameters pass to \code{SpiecEasi::spiec.easi} when \code{network_method = "SpiecEasi"};
		#'   pass to \code{NetCoMi::netConstruct} when \code{network_method = "gcoda"}; 
		#'   pass to \code{beemStatic::func.EM} when \code{network_method = "beemStatic"}.
		#' @return \code{res_network} stored in object.
		#' @examples
		#' \dontrun{
		#' # for correlation network
		#' t1 <- trans_network$new(dataset = dataset, cor_method = "pearson", 
		#' 		taxa_level = "OTU", filter_thres = 0.001)
		#' t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
		#' t1 <- trans_network$new(dataset = dataset, cor_method = NULL, filter_thres = 0.003)
		#' t1$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb")
		#' t1 <- trans_network$new(dataset = dataset, cor_method = NULL, filter_thres = 0.005)
		#' t1$cal_network(network_method = "beemStatic")
		#' t1 <- trans_network$new(dataset = dataset, cor_method = NULL, filter_thres = 0.001)
		#' t1$cal_network(network_method = "FlashWeave")
		#' }
		cal_network = function(
			network_method = c("COR", "SpiecEasi", "gcoda", "FlashWeave", "beemStatic")[1],
			COR_p_thres = 0.01,
			COR_p_adjust = "fdr",
			COR_weight = TRUE,
			COR_cut = 0.6,
			COR_optimization = FALSE,
			COR_optimization_low_high = c(0.01, 0.8),
			COR_optimization_seq = 0.01,
			SpiecEasi_method = "mb",
			FlashWeave_tempdir = NULL,
			FlashWeave_meta_data = FALSE,
			FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=true",
			beemStatic_t_strength = 0.001,
			beemStatic_t_stab = 0.8,
			add_taxa_name = "Phylum",
			delete_unlinked_nodes = TRUE,
			usename_rawtaxa_when_taxalevel_notOTU = FALSE,
			...
			){
			private$check_igraph()
			if(!is.null(network_method)){
				sampleinfo <- self$sample_table
				taxa_level <- self$taxa_level
				taxa_table <- self$tax_table
				
				network_method <- match.arg(network_method, c("COR", "SpiecEasi", "gcoda", "FlashWeave", "beemStatic"))
				
				message("---------------- ", Sys.time()," : Start ----------------")
				if(network_method == "COR"){
					if(is.null(self$res_cor_p)){
						stop("The res_cor_p list in the object is NULL! Please check the created object!")
					}
					cortable <- self$res_cor_p$cor
					raw_p <- self$res_cor_p$p
					if(ncol(cortable) != ncol(raw_p)){
						stop("Correlation table and p value table have different column numbers !")
					}
					raw_vector_p <- raw_p %>% as.dist %>% as.numeric
					message("Perform p value adjustment with ", COR_p_adjust, " method ...")
					adp_raw <- p.adjust(raw_vector_p, method = COR_p_adjust)
					use_names <- colnames(raw_p)
					names_combn <- t(combn(use_names, 2))
					table_convert <- cbind.data.frame(names_combn, adjust.p = adp_raw, stringsAsFactors = FALSE)
					adp <- private$vec2mat(datatable = table_convert, use_names = use_names, value_var = "adjust.p", rep_value = 0)
					if(! identical(colnames(cortable), colnames(adp))){
						adp <- adp[colnames(cortable), colnames(cortable)]
					}
					if(COR_optimization == T){
						message("Start COR optimizing ...")
						tc1 <- private$rmt(cortable, low_thres = COR_optimization_low_high[1], high_thres = COR_optimization_low_high[2], seq_by = COR_optimization_seq)
						message("The optimized COR threshold: ", tc1, "...\n")
					}
					else {
						tc1 <- COR_cut
					}
					diag(cortable) <- 0
					cor_matrix <- as.matrix(cortable)
					if(!any(abs(cortable) >= tc1)){
						stop("All the correlation coefficients are smaller than the threshold! Please lower the COR_cut parameter!")
					}
					cor_matrix[abs(cortable) >= tc1] <- 1
					cor_matrix[adp > COR_p_thres] <- 0
					if(!any(cor_matrix == 1)){
						stop("All the correlation coefficients larger than COR_cut parameter are not significant under current COR_p_thres parameter!")
					}
					cor_matrix[cor_matrix != 1] <- 0
					network <- graph.adjacency(cor_matrix, mode = "undirected")
					edges <- t(sapply(1:ecount(network), function(x) ends(network, x)))
					E(network)$label <- unlist(lapply(seq_len(nrow(edges)), function(x) ifelse(cortable[edges[x, 1], edges[x, 2]] > 0, "+", "-")))
					if(COR_weight == T){
						E(network)$weight <- unlist(lapply(seq_len(nrow(edges)), function(x) abs(cortable[edges[x, 1], edges[x, 2]])))
					}else{
						E(network)$weight <- rep.int(1, ecount(network))
					}
				}
				if(network_method == "SpiecEasi"){
					if(!require("SpiecEasi")){
						stop("SpiecEasi package is not installed! See https://github.com/zdk123/SpiecEasi ")
					}
					SpiecEasi_method <- match.arg(SpiecEasi_method, c("glasso", "mb"))
					use_abund <- self$data_abund %>% as.matrix
					spieceasi_fit <- SpiecEasi::spiec.easi(use_abund, method = SpiecEasi_method, ...)
					if(SpiecEasi_method == "glasso"){
						assoMat <- stats::cov2cor(as.matrix(getOptCov(spieceasi_fit))) * SpiecEasi::getRefit(spieceasi_fit)
					}else{
						assoMat <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spieceasi_fit), mode = "ave")
					}
					assoMat %<>% as.matrix
					network <- adj2igraph(assoMat)
					V(network)$name <- colnames(use_abund)
				}
				if(network_method == "gcoda"){
					if(!require("NetCoMi")){
						stop("NetCoMi package is not installed! Please see https://github.com/stefpeschel/NetCoMi ")
					}
					netConstruct_raw <- netConstruct(self$data_abund, measure = network_method, ...)
					self$res_gcoda_raw <- netConstruct_raw
					message("The raw file is stored in object$res_gcoda_raw ...")
					asso_matrix <- netConstruct_raw$assoMat1
					network <- SpiecEasi::adj2igraph(asso_matrix)
					V(network)$name <- colnames(asso_matrix)
				}
				if(network_method == "FlashWeave"){
					use_abund <- self$data_abund
					oldwd <- getwd()
					# working directory is not changed when quit
					on.exit(setwd(oldwd))

					if(is.null(FlashWeave_tempdir)){
						tem_dir <- tempdir()
					}else{
						tem_dir <- FlashWeave_tempdir
						if(!dir.exists(tem_dir)){
							stop("The input temporary directory: ", tem_dir, " does not exist!")
						}
					}
					setwd(tem_dir)
					write.table(use_abund, "taxa_table_FlashWeave.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
					L1 <- "using FlashWeave\n"
					L2 <- 'data_path = "taxa_table_FlashWeave.tsv"\n'
					if(FlashWeave_meta_data){
						if(is.null(self$data_env)){
							stop("FlashWeave_meta_data is TRUE, but object$data_env not found! 
								Please use env_cols or add_data parameter of trans_network$new to provide the metadata when creating the object!")
						}
						meta_data <- self$data_env
						write.table(meta_data, "meta_table_FlashWeave.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
						L3 <- 'meta_data_path = "meta_table_FlashWeave.tsv"\n'
					}else{
						L3 <- "\n"
					}
					if(FlashWeave_meta_data){
						L4 <- paste0(gsub(",$|,\\s+$", "", paste0("netw_results = learn_network(data_path, meta_data_path, ", FlashWeave_other_para)), ")\n")
					}else{
						L4 <- paste0(gsub(",$|,\\s+$", "", paste0("netw_results = learn_network(data_path, ", FlashWeave_other_para)), ")\n")
					}
					L5 <- 'save_network("network_FlashWeave.gml", netw_results)'
					L <- paste0(L1, L2, L3, L4, L5)
					openfile <- file("calculate_network.jl", "wb")
					write(L, file = openfile)
					close(openfile)
					message("The temporary files are in ", tem_dir, " ...")
					message("Check Julia ...")
					check_out <- system(command = "julia -help", intern = TRUE)
					if(length(check_out) == 0){
						stop("Julia language not found in the system path! Install it from https://julialang.org/downloads/ and add bin directory to the path!")
					}
					message("Run the FlashWeave ...")
					system("julia calculate_network.jl")
					network <- read_graph("network_FlashWeave.gml", format = "gml")
					network <- set_vertex_attr(network, "name", value = V(network)$label)
				}
				if(network_method == "beemStatic"){
					if(!require("beemStatic")){
						stop("beemStatic package is not installed! See https://github.com/CSB5/BEEM-static ")
					}
					use_abund <- self$data_abund %>% t %>% as.data.frame
					taxa_low <- apply(use_abund, 1, function(x) sum(x != 0)) %>% .[which.min(.)]
					message("The feature table has ", nrow(use_abund), " taxa. The taxon with the lowest occurrence frequency was ", names(taxa_low)[1], 
						", found in ", taxa_low[1], " samples of total ", ncol(use_abund), " samples. If an error occurs because of low frequency, ",
						"please filter more taxa with low abundance using filter_thres parameter when creating the trans_network object ...")
					beem.out <- func.EM(use_abund, ...)
					self$res_beemStatic_raw <- beem.out
					message('beemStatic result is stored in object$res_beemStatic_raw ...')
					# based on the showInteraction function
					b <- t(beem2param(beem.out)$b.est)
					diag(b) <- 0
					if(!is.null(beem.out$resample)){
						b[t(beem.out$resample$b.stab < beemStatic_t_stab)] <- 0
					}
					b[abs(b) < beemStatic_t_strength] <- 0
					network <- graph.adjacency(b, mode = 'directed', weighted = 'weight')
					V(network)$name <- rownames(use_abund)
				}
				if(ecount(network) == 0){
					stop("No edge found in the network! Please check the input parameters!")
				}
				if(network_method != "COR"){
					E(network)$label <- ifelse(E(network)$weight > 0, '+', '-')
				}
				E(network)$weight <- abs(E(network)$weight)
				message("---------------- ", Sys.time()," : Finish ----------------")
				
				if(taxa_level != "OTU"){
					delete_nodes <- taxa_table %>% 
						.[grepl("__$|uncultured", .[, taxa_level]), ] %>% 
						rownames %>% 
						.[. %in% V(network)$name]
					network %<>% delete_vertices(delete_nodes)
				}
				if(delete_unlinked_nodes){
					network <- private$rm_unlinked_node(network)	
				}
				V(network)$taxa <- V(network)$name
				if(!is.null(add_taxa_name)){
					if(!is.null(taxa_table)){
						for(i in add_taxa_name){
							if(i %in% colnames(taxa_table)){
								network <- set_vertex_attr(network, i, value = V(network)$name %>% 
									taxa_table[., i] %>% 
									gsub("^.__", "", .))
							}else{
								message("Skip adding taxonomy: ", i, " to node as it is not in colnames of tax_table ...")
							}
						}
					}else{
						message('Skip adding taxonomy to node as tax_table is not found ...')
					}
				}
				V(network)$RelativeAbundance <- self$data_relabund[V(network)$name]

				if(taxa_level != "OTU"){
					if(usename_rawtaxa_when_taxalevel_notOTU){
						network <- set_vertex_attr(network, taxa_level, value = V(network)$name %>% 
							taxa_table[., taxa_level] %>% 
							gsub("^.__", "", .))
					}else{
						network <- set_vertex_attr(network, "name", value = V(network)$name %>% 
							taxa_table[., taxa_level] %>% 
							gsub("^.__", "", .))
					}
				}
				
				self$res_network <- network
				message('The result network is stored in object$res_network ...')
			}else{
				message('No network_method selected! Please manually assign object$res_network!')
			}
		},
		#' @description
		#' Calculate network modules and add module names to the network node properties.
		#'
		#' @param method default "cluster_fast_greedy"; the method used to find the optimal community structure of a graph;
		#' 	 the following are available functions (options) from igraph package: \cr
		#' 	 \code{"cluster_fast_greedy"}, \code{"cluster_walktrap"}, \code{"cluster_edge_betweenness"}, \cr
		#' 	 \code{"cluster_infomap"}, \code{"cluster_label_prop"}, \code{"cluster_leading_eigen"}, \cr
		#' 	 \code{"cluster_louvain"}, \code{"cluster_spinglass"}, \code{"cluster_optimal"}. \cr
		#' 	 For the details of these functions, please see the help document, such as \code{help(cluster_fast_greedy)};
		#' 	 Note that the default \code{"cluster_fast_greedy"} method can not be applied to directed network. 
		#' 	 If directed network is provided, the function can automatically switch the default method from \code{"cluster_fast_greedy"} to \code{"cluster_walktrap"}.
		#' @param module_name_prefix default "M"; the prefix of module names; module names are made of the module_name_prefix and numbers;
		#'   numbers are assigned according to the sorting result of node numbers in modules with decreasing trend.
		#' @return \code{res_network} with modules, stored in object.
		#' @examples
		#' \donttest{
		#' t1 <- trans_network$new(dataset = dataset, cor_method = "pearson", 
		#' 		taxa_level = "OTU", filter_thres = 0.0002)
		#' t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
		#' t1$cal_module(method = "cluster_fast_greedy")
		#' }
		cal_module = function(method = "cluster_fast_greedy", module_name_prefix = "M"){
			private$check_igraph()
			private$check_network()
			# add modules
			network <- self$res_network
			if(!is.character(method)){
				stop("The parameter method must be character!")
			}
			if(method == "cluster_fast_greedy" & is_directed(network)){
				message('The default method "cluster_fast_greedy" can not be applied to directed network! ',
					'Automatically switch to method "cluster_walktrap" ...')
				message('Invoke cluster_walktrap function to find densely connected subgraphs ...')
				method <- "cluster_walktrap"
			}else{
				message('Use ', method, ' function to partition modules ...')
			}
			# use NSE
			res_member <- parse(text = paste0(method, "(network)")) %>% eval
			
			mod1 <- as.character(res_member$membership)
			mod2 <- sort(table(mod1), decreasing = TRUE)
			for(i in seq_along(mod2)){
				mod1[mod1 == names(mod2)[i]] <- paste0(module_name_prefix, i)
			}
			message('Totally, ', length(mod2), ' modules are idenfified ...')
			network <- set_vertex_attr(network, "module", value = mod1)
			message('Modules are assigned in network with attribute name -- module ...')
			self$res_network <- network
		},
		#' @description
		#' Save network as gexf style, which can be opened by Gephi (\href{https://gephi.org/}{https://gephi.org/}).
		#'
		#' @param filepath default "network.gexf"; file path to save the network.
		#' @return None
		#' @examples
		#' \dontrun{
		#' t1$save_network(filepath = "network.gexf")
		#' }
		save_network = function(filepath = "network.gexf"){
			if(!require("rgexf")){
				stop("Please first install rgexf package with command: install.packages('rgexf') !")
			}
			private$check_igraph()
			private$check_network()
			private$saveAsGEXF(network = self$res_network, filepath = filepath)
		},
		#' @description
		#' Calculate network properties.
		#'
		#' @return \code{res_network_attr} stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_network_attr()
		#' }
		cal_network_attr = function(){
			private$check_igraph()
			private$check_network()
			self$res_network_attr <- private$network_attribute(self$res_network)
			message('Result is stored in object$res_network_attr ...')
		},
		#' @description
		#' Get the node property table. The properties may include the node names, modules allocation, degree, betweenness, abundance, 
		#'   taxonomy, within-module connectivity and among-module connectivity <doi:10.1016/j.geoderma.2022.115866>.
		#'
		#' @param node_roles default TRUE; whether calculate node roles, i.e. Module hubs, Network hubs, Connectors and Peripherals <doi:10.1016/j.geoderma.2022.115866>.
		#' @return \code{res_node_table} in object; Abundance expressed as a percentage; 
		#'   betweenness_centrality: betweenness centrality; betweenness_centrality: closeness centrality; eigenvector_centrality: eigenvector centrality; 
		#'  z: within-module connectivity; p: among-module connectivity.
		#' @examples
		#' \donttest{
		#' t1$get_node_table(node_roles = TRUE)
		#' }
		get_node_table = function(node_roles = TRUE){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			# Add abundance info
			sum_abund <- self$data_relabund
			node_table <- data.frame(name = V(network)$name) %>% `rownames<-`(.[, 1])
			node_table$degree <- igraph::degree(network)[rownames(node_table)]
			node_table$betweenness_centrality <- igraph::betweenness(network)[rownames(node_table)]
			node_table$closeness_centrality <- igraph::closeness(network)[rownames(node_table)]
			eigenvec_centralraw <- igraph::eigen_centrality(network)
			node_table$eigenvector_centrality <- eigenvec_centralraw$vector[rownames(node_table)]
			# Same with the above operation to make the names corresponded
			if(self$taxa_level != "OTU"){
				# create a replace_table to match the taxa name and marker name when taxa_level is not "OTU"
				replace_table <- data.frame(V(network)$name, V(network)$taxa, stringsAsFactors = FALSE) %>% `row.names<-`(.[,1])
				node_table$Abundance <- sum_abund[replace_table[rownames(node_table), 2]]
			}else{
				node_table$Abundance <- sum_abund[rownames(node_table)]
			}
			if(!is.null(V(network)$module)){
				node_table$module <- V(network)$module
				unique_level <- node_table$module %>% unique
				order_level <- unique_level %>% gsub("M", "", .) %>% as.numeric %>% order
				node_table$module %<>% factor(levels = unique_level[order_level])
			}
			if(node_roles){
				res_module_roles <- private$module_roles(network) %>% cbind.data.frame(name = rownames(.), .) %>% .[, c("name", "z", "p", "taxa_roles")]
				node_table %<>% dplyr::left_join(., res_module_roles, by = c("name" = "name")) %>% `rownames<-`(.$name)
			}
			if(self$taxa_level != "OTU"){
				node_table %<>% cbind.data.frame(., 
					self$tax_table[replace_table[rownames(.), 2], 1:which(colnames(self$tax_table) %in% self$taxa_level), drop = FALSE])
			}else{
				if(!is.null(self$tax_table)){
					node_table %<>% cbind.data.frame(., self$tax_table[rownames(.), , drop = FALSE])
				}
			}
			self$res_node_table <- node_table
			message('Result is stored in object$res_node_table ...')
		},
		#' @description
		#' Get the edge property table, including connected nodes, label and weight.
		#'
		#' @return \code{res_edge_table} in object.
		#' @examples
		#' \donttest{
		#' t1$get_edge_table()
		#' }
		get_edge_table = function(){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			res_edge_table <- igraph::as_data_frame(network, what = "edges")
			colnames(res_edge_table)[1:2] <- c("node1", "node2")
			if(! "weight" %in% colnames(res_edge_table)){
				res_edge_table$weight <- NA
			}
			self$res_edge_table <- res_edge_table
			message('Result is stored in object$res_edge_table ...')
		},
		#' @description
		#' Get the adjacency matrix from the network graph.
		#'
		#' @param ... parameters passed to as_adjacency_matrix function of \code{igraph} package.
		#' @return \code{res_adjacency_matrix} in object.
		#' @examples
		#' \donttest{
		#' t1$get_adjacency_matrix(attr = "weight")
		#' }
		get_adjacency_matrix = function(...){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			self$res_adjacency_matrix <- as_adjacency_matrix(network, ...) %>% as.matrix
			message('Result is stored in object$res_adjacency_matrix ...')
		},
		#' @description
		#' Plot the network based on a series of methods from other packages, such as \code{igraph}, \code{ggraph} and \code{networkD3}. 
		#' The networkD3 package provides dynamic network. It is especially useful for a glimpse of the whole network structure and finding 
		#' the interested nodes and edges in a large network. In contrast, the igraph and ggraph methods are suitable for relatively small network.
		#'
		#' @param method default "igraph"; The available options:
		#'   \describe{
		#'     \item{\strong{'igraph'}}{call \code{plot.igraph} function in \code{igraph} package for a static network; see plot.igraph for the parameters}
		#'     \item{\strong{'ggraph'}}{call \code{ggraph} function in \code{ggraph} package for a static network}
		#'     \item{\strong{'networkD3'}}{use forceNetwork function in \code{networkD3} package for a dynamic network; see forceNetwork function for the parameters}
		#'   }
		#' @param node_label default "name"; node label shown in the plot for \code{method = "ggraph"} or \code{method = "networkD3"}; 
		#'   Please see the column names of object$res_node_table, which is the returned table of function object$get_node_table;
		#'   User can select other column names in res_node_table.
		#' @param node_color default NULL; node color assignment for \code{method = "ggraph"} or \code{method = "networkD3"}; 
		#'   Select a column name of \code{object$res_node_table}, such as "module".
		#' @param ggraph_layout default "fr"; for \code{method = "ggraph"}; see \code{layout} parameter of \code{create_layout} function in \code{ggraph} package.
		#' @param ggraph_node_size default 2; for \code{method = "ggraph"}; the node size.
		#' @param ggraph_node_text default TRUE; for \code{method = "ggraph"}; whether show the label text of nodes.
		#' @param ggraph_text_color default NULL; for \code{method = "ggraph"}; a column name of object$res_node_table used to assign label text colors.
		#' @param ggraph_text_size default 3; for \code{method = "ggraph"}; the node label text size.
		#' @param networkD3_node_legend default TRUE; used for \code{method = "networkD3"}; logical value to enable node colour legends;
		#'   Please see the legend parameter in networkD3::forceNetwork function.
		#' @param networkD3_zoom default TRUE; used for \code{method = "networkD3"}; logical value to enable (TRUE) or disable (FALSE) zooming;
		#'   Please see the zoom parameter in networkD3::forceNetwork function.
		#' @param ... parameters passed to \code{plot.igraph} function when \code{method = "igraph"} or forceNetwork function when \code{method = "networkD3"}.
		#' @return network plot.
		#' @examples
		#' \donttest{
		#' t1$plot_network(method = "igraph", layout = layout_with_kk)
		#' t1$plot_network(method = "ggraph", node_color = "module")
		#' t1$plot_network(method = "networkD3", node_color = "module")
		#' }
		plot_network = function(method = c("igraph", "ggraph", "networkD3")[1], 
			node_label = "name", 
			node_color = NULL, 
			ggraph_layout = "fr",
			ggraph_node_size = 2,
			ggraph_node_text = TRUE,
			ggraph_text_color = NULL,
			ggraph_text_size = 3,
			networkD3_node_legend = TRUE, 
			networkD3_zoom = TRUE, 
			...
			){
			private$check_network()
			method <- match.arg(method, c("igraph", "ggraph", "networkD3"))
			if(method == "igraph"){
				network <- self$res_network
				g <- igraph::plot.igraph(network, ...)
			}else{
				if(is.null(self$res_node_table)){
					message("Run get_node_table function to obtain the node property table ...")
					self$get_node_table()
				}
				node_table <- self$res_node_table
				if(!is.null(node_color)){
					if(node_color == "module"){
						if(!any(colnames(node_table) %in% "module")){
							stop("Please first run cal_module function to get modules !")
						}
					}else{
						if(!any(colnames(node_table) %in% node_color)){
							stop("The node_color provided is not found in the object$res_node_table !")
						}
					}
				}
			}
			if(method == "ggraph"){
				if(!require("ggraph")){
					stop("Please first install ggraph package with the command: install.packages('ggraph') !")
				}
				network <- self$res_network
				network_layout <- create_layout(network, layout = ggraph_layout, ...)
				# add more node properties
				node_table %<>% {.[, c("name", colnames(.)[!colnames(.) %in% colnames(network_layout)])]}
				use_network_layout <- dplyr::left_join(network_layout, node_table, by = c("name" = "name"))

				g <- ggraph(use_network_layout)
				if(is_directed(network)){
					g <- g + geom_edge_arc(aes(col = label, width = weight), arrow = arrow(length = unit(2, 'mm')), strength = 0.2, alpha = 0.5)
				}else{
					g <- g + geom_edge_link(aes(col = label, width = weight), alpha = 0.8)
				}
				g <- g + geom_node_point(aes_meco(colour = node_color), size = ggraph_node_size, alpha = 0.5)
				if(ggraph_node_text){
					g <- g + geom_node_text(aes_meco(colour = ggraph_text_color, label = node_label), size = ggraph_text_size, repel = TRUE)
				}
				g <- g + scale_edge_width(range = c(0.5, 2)) +
					theme_void()
			}
			if(method == "networkD3"){
				try_find <- try(find.package("networkD3"), silent = TRUE)
				if(inherits(try_find, "try-error")){
					stop("Please first install networkD3 package with command: install.packages('networkD3') !")
				}
				if(is.null(self$res_edge_table)){
					message("Run get_edge_table function to obtain the edge property table ...")
					self$get_edge_table()
				}
				edge_table <- self$res_edge_table
				replace_table <- data.frame(number = 1:nrow(node_table), name = rownames(node_table))
				rownames(replace_table) <- replace_table$name
				rownames(node_table) <- replace_table$number
				edge_table$node1 <- replace_table[edge_table$node1, 1] - 1
				edge_table$node2 <- replace_table[edge_table$node2, 1] - 1
				
				if(is.null(node_color)){
					stop("networkD3 require the node_color input. Please select one column name of object$res_node_table !")
				}else{
					g <- networkD3::forceNetwork(Links = edge_table, Nodes = node_table, Source = "node1", Target = "node2",
						NodeID = node_label, zoom = networkD3_zoom, Group = node_color, legend = networkD3_node_legend, ...)
				}
			}
			g
		},
		#' @description
		#' Calculate eigengenes of modules, i.e. the first principal component based on PCA analysis, and the percentage of variance <doi:10.1186/1471-2105-13-113>.
		#'
		#' @return \code{res_eigen} and \code{res_eigen_expla} in object.
		#' @examples
		#' \donttest{
		#' t1$cal_eigen()
		#' }
		cal_eigen = function(){
			private$check_igraph()
			private$check_network()
			use_abund <- self$data_abund
			if(is.null(self$res_node_table)){
				message("Run get_node_table function to get the node property table ...")
				self$get_node_table()
			}
			node_table <- self$res_node_table
			# calculate eigengene for each module
			res_eigen <- list()
			res_eigen_expla <- c()
			for(i in unique(as.character(node_table$module))){
				tax_names <- rownames(node_table[as.character(node_table$module) == i, ])
				if(length(tax_names) < 3){
					next
				}
				if(self$taxa_level != "OTU"){
					network <- self$res_network
					replace_table <- data.frame(V(network)$name, V(network)$taxa, stringsAsFactors = FALSE) %>% `row.names<-`(.[,1])
					tax_names <- replace_table[tax_names, 2]
				}
				sel_abund <- use_abund[, tax_names]
				pca_model <- rda(sel_abund)
				sel_scores <- scores(pca_model, choices = 1)$sites %>% as.data.frame
				colnames(sel_scores)[1] <- i
				res_eigen[[i]] <- sel_scores
				expla <- paste0(round(pca_model$CA$eig/pca_model$CA$tot.chi*100, 1)[1], "%")
				names(expla) <- i
				res_eigen_expla <- c(res_eigen_expla, expla)
			}
			res_eigen <- do.call(cbind, res_eigen)
			self$res_eigen <- res_eigen
			self$res_eigen_expla <- res_eigen_expla
			message('Result is stored in object$res_eigen and object$res_eigen_expla ...')
		},
		#' @description
		#' Plot the classification and importance of nodes, see object$res_node_table for the variable names used in the parameters.
		#'
		#' @param use_type default 1; 1 or 2; 1 represents taxa roles area plot; 2 represents the layered plot with taxa as x axis.
		#' @param roles_color_background default FALSE; for use_type=1; TRUE: use background colors for each area; FALSE: use classic point colors.
		#' @param roles_color_values default NULL; for use_type=1; color palette for background or points.
		#' @param add_label default FALSE; for use_type = 1; whether add labels for the points.
		#' @param add_label_group default "Network hubs"; If add_label = TRUE; which part of tax_roles is used to show labels; character vectors.
		#' @param add_label_text default "name"; If add_label = TRUE; which column of object$res_node_table is used to label the text.
		#' @param label_text_size default 4; The text size of the label.
		#' @param label_text_color default "grey50"; The text color of the label.
		#' @param label_text_italic default FALSE; whether use italic style for the label text.
		#' @param label_text_parse default FALSE; whether parse the label text. See the parse parameter in \code{ggrepel::geom_text_repel} function.
		#' @param plot_module default FALSE; for use_type=1; whether plot the modules information.
		#' @param x_lim default c(0, 1); for use_type=1; x axis range when roles_color_background = FALSE.
		#' @param use_level default "Phylum"; for use_type=2; used taxonomic level in x axis.
		#' @param show_value default c("z", "p"); for use_type=2; used variable in y axis.
		#' @param show_number default 1:10; for use_type=2; showed number in x axis, sorting according to the nodes number.
		#' @param plot_color default "Phylum"; for use_type=2; used variable for color.
		#' @param plot_shape default "taxa_roles"; for use_type=2; used variable for shape.
		#' @param plot_size default "Abundance"; for use_type=2; used for point size; a fixed number (e.g. 5) is also available.
		#' @param color_values default RColorBrewer::brewer.pal(12, "Paired"); for use_type=2; color vector
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); for use_type=2; shape vector, see ggplot2 tutorial for the shape meaning.
		#' @param ... parameters pass to geom_point.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_taxa_roles(roles_color_background = FALSE)
		#' }
		plot_taxa_roles = function(
			use_type = c(1, 2)[1],
			roles_color_background = FALSE,
			roles_color_values = NULL,
			add_label = FALSE, 
			add_label_group = "Network hubs", 
			add_label_text = "name", 
			label_text_size = 4, 
			label_text_color = "grey50", 
			label_text_italic = FALSE,
			label_text_parse = FALSE,
			plot_module = FALSE,
			x_lim = c(0, 1),
			use_level = "Phylum",
			show_value = c("z", "p"),
			show_number = 1:10,
			plot_color = "Phylum",
			plot_shape = "taxa_roles",
			plot_size = "Abundance",
			color_values = RColorBrewer::brewer.pal(12, "Paired"),
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			...
			){
			if(is.null(self$res_node_table)){
				message("Run get_node_table function to get the node property table ...")
				self$get_node_table()
			}
			if(use_type == 1){
				res <- private$plot_roles_1(node_roles = self$res_node_table, 
					roles_color_background = roles_color_background,
					roles_color_values = roles_color_values, 
					add_label = add_label, 
					add_label_group = add_label_group, 
					add_label_text = add_label_text, 
					label_text_size = label_text_size, 
					label_text_color = label_text_color, 
					label_text_italic = label_text_italic,
					label_text_parse = label_text_parse,
					module = plot_module,
					x_lim = x_lim,
					...
				)
			}
			if(use_type == 2){
				res <- private$plot_roles_2(node_roles = self$res_node_table, 
					plot_color = plot_color,
					plot_shape = plot_shape,
					use_level = use_level, 
					show_value = show_value, 
					show_number = show_number,
					plot_size = plot_size,
					color_values = color_values, 
					shape_values = shape_values,
					...
				)
			}
			res
		},
		#' @description
		#' Subset of the network.
		#'
		#' @param node default NULL; provide the node names that you want to use in the sub-network.
		#' @param edge default NULL; provide the edge name needed; must be one of "+" or "-".
		#' @param rm_single default TRUE; whether remove the nodes without any edge in the sub-network.
		#'   So this function can also be used to remove the nodes withou any edge when node and edge are both NULL.
		#' @return a new network
		#' @examples
		#' \donttest{
		#' t1$subset_network(node = t1$res_node_table %>% base::subset(module == "M1") %>% 
		#'   rownames, rm_single = TRUE)
		#' # return a sub network that contains all nodes of module M1
		#' }
		subset_network = function(node = NULL, edge = NULL, rm_single = TRUE){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			if(!is.null(node)){
				nodes_raw <- V(network)$name
				delete_nodes <- nodes_raw %>% .[! . %in% node]
				sub_network <- delete_vertices(network, delete_nodes)
			}
			if(!is.null(edge)){
				label_raw <- E(network)$label
				sub_network <- delete_edges(network, which(label_raw != edge))
			}
			if(is.null(node) & is.null(edge)){
				sub_network <- network
			}
			# whether remove the single node without edges
			if(rm_single == T){
				sub_network <- private$rm_unlinked_node(sub_network)
			}
			sub_network
		},
		#' @description
		#' Fit degrees to a power law distribution. First, perform a bootstrapping hypothesis test to determine whether degrees follow a power law distribution.
		#' If the distribution follows power law, then fit degrees to power law distribution and return the parameters.
		#'
		#' @param ... parameters pass to bootstrap_p function in poweRlaw package.
		#' @return \code{res_powerlaw_p} and \code{res_powerlaw_fit}; see \code{poweRlaw::bootstrap_p} function for the bootstrapping p value details;
		#'   see \code{igraph::fit_power_law} function for the power law fit return details.
		#' @examples
		#' \donttest{
		#' t1$cal_powerlaw()
		#' }
		cal_powerlaw = function(...){
			private$check_network()
			network <- self$res_network
			degree_dis <- igraph::degree(network)
			if(!require("poweRlaw")){
				stop("Please first install poweRlaw package from CRAN !")
			}
			resdispl <- poweRlaw::displ$new(degree_dis + 1)
			est_xmin <- poweRlaw::estimate_xmin(resdispl)
			message('Estimated lower bound of degree: ', est_xmin$xmin)
			resdispl$setXmin(est_xmin)
			message('Perform bootstrapping ...')
			bootstrap_res <- poweRlaw::bootstrap_p(resdispl, ...)
			self$res_powerlaw_p <- bootstrap_res
			message('Bootstrap result is stored in object$res_powerlaw_p ...')
			message('Bootstrap p value: ', bootstrap_res$p)
			if(bootstrap_res$p < 0.05){
				message("The p value < 0.05; Degrees do not follow power law distribution ...")
			}else{
				message("Degrees follow power law distribution ...")
				res_powerlaw_fit <- fit_power_law(degree_dis + 1, xmin = est_xmin$xmin)
				message("The estimated alpha: ", res_powerlaw_fit$alpha)
				self$res_powerlaw_fit <- res_powerlaw_fit
				message('Powerlaw fitting result is stored in object$res_powerlaw_fit ...')
			}
		},
		#' @description
		#' This function is used to sum the links number from one taxa to another or in the same taxa, for example, at Phylum level.
		#' This is very useful to fast see how many nodes are connected between different taxa or within the taxa.
		#'
		#' @param taxa_level default "Phylum"; taxonomic rank.
		#' @return \code{res_sum_links_pos} and \code{res_sum_links_neg} in object.
		#' @examples
		#' \donttest{
		#' t1$cal_sum_links(taxa_level = "Phylum")
		#' }
		cal_sum_links = function(taxa_level = "Phylum"){
			if(is.null(self$tax_table)){
				stop("The tax_table is required! Please check your trans_network object when creating it!")
			}else{
				taxa_table <- self$tax_table
			}
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			
			if(self$taxa_level != "OTU"){
				replace_table <- data.frame(V(network)$name, V(network)$taxa, stringsAsFactors = FALSE) %>% `row.names<-`(.[, 2])
				taxa_table %<>% .[rownames(.) %in% replace_table[, 2], ]
				rownames(taxa_table) <- replace_table[rownames(taxa_table), 1]
			}
			if(is.null(self$res_edge_table)){
				self$get_edge_table()
			}
			link_table <- self$res_edge_table
			
			if(is.null(E(network)$label)){
				message('No edge label found. All edges are viewed as positive links ...')
				self$res_sum_links_pos <- private$sum_link(taxa_table = taxa_table, link_table = link_table, taxa_level = taxa_level)
				message('Results are stored in object$res_sum_links_pos ...')
			}else{				
				if(! any(c("+", "-") %in% link_table[, "label"])){
					stop("Please check the edge labels! The labels should be + or - !")
				}
				if("+" %in% link_table[, "label"]){
					link_table_use <- link_table[link_table[, "label"] %in% "+", ]
					self$res_sum_links_pos <- private$sum_link(taxa_table = taxa_table, link_table = link_table_use, taxa_level = taxa_level)
					message('The positive results are stored in object$res_sum_links_pos ...')
				}else{
					message('No positive edges found ...')
				}
				if("-" %in% link_table[, "label"]){
					link_table_use <- link_table[link_table[, "label"] %in% "-", ]
					self$res_sum_links_neg <- private$sum_link(taxa_table = taxa_table, link_table = link_table_use, taxa_level = taxa_level)
					message('The negative results are stored in object$res_sum_links_neg ...')
				}else{
					message('No negative edges found ...')
				}
			}
		},
		#' @description
		#' Plot the summed linkages among taxa.
		#'
		#' @param plot_pos default TRUE; If TRUE, plot the summed positive linkages; If FALSE, plot the summed negative linkages.
		#' @param plot_num default NULL; number of taxa presented in the plot.
		#' @param color_values default RColorBrewer::brewer.pal(8, "Dark2"); colors palette for taxa.
		#' @param method default c("chorddiag", "circlize")[1]; chorddiag package <https://github.com/mattflor/chorddiag> or circlize package.
		#' @param ... pass to \code{chorddiag::chorddiag} function when \code{method = "chorddiag"} or 
		#'	 \code{circlize::chordDiagram} function when \code{method = "circlize"}.
		#'	 Note that for \code{circlize::chordDiagram} function, \code{keep.diagonal}, \code{symmetric} and \code{self.link} parameters have been fixed to fit the input data.
		#' @return please see the invoked function.
		#' @examples
		#' \dontrun{
		#' test1$plot_sum_links(method = "chorddiag", plot_pos = TRUE, plot_num = 10)
		#' test1$plot_sum_links(method = "circlize", transparency = 0.2, 
		#' 	  annotationTrackHeight = circlize::mm_h(c(5, 5)))
		#' }
		plot_sum_links = function(plot_pos = TRUE, plot_num = NULL, color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			method = c("chorddiag", "circlize")[1],
			...){
			method <- match.arg(method, c("chorddiag", "circlize"))
			
			if(is.null(self$res_sum_links_pos) & is.null(self$res_sum_links_neg)){
				stop("Please first run cal_sum_links function!")
			}
			if(plot_pos == T){
				message("Extract the positive link information ...")
				if(is.null(self$res_sum_links_pos)){
					stop("res_sum_links_pos in object is NULL! No positive links can be used!")
				}else{
					use_data <- self$res_sum_links_pos
				}
			}else{
				message("Extract the negative link information ...")
				if(is.null(self$res_sum_links_neg)){
					stop("res_sum_links_neg in object is NULL! No negative links can be used!")
				}else{
					use_data <- self$res_sum_links_neg
				}
			}
			if(!is.null(plot_num)){
				if(plot_num > ncol(use_data)){
					message("The plot_num provided is larger than the total taxa number. Use the taxa number instead of it ...")
					plot_num <- ncol(use_data)
				}
				use_data %<>% .[1:plot_num, 1:plot_num]
			}
			if(is.null(color_values)){
				stop("Please provide the color_values parameter!")
			}else{
				if(nrow(use_data) > length(color_values)){
					message("The taxa number ", nrow(use_data), " is larger than the length of input color_values ", length(color_values), 
						". Only select ", length(color_values), " taxa ...")
					use_data %<>% .[1:length(color_values), 1:length(color_values)]
				}
				color_values %<>% .[1:nrow(use_data)]
			}
			if(method == "chorddiag"){
				chorddiag::chorddiag(use_data, groupColors = color_values, ...)
			}else{
				circlize::chordDiagram(use_data, grid.col = color_values, keep.diagonal = TRUE, symmetric = TRUE, 
					self.link = 1, ...)
			}
		},
		#' @description
		#' Generate random networks, compare them with the empirical network and get the p value of topological properties.
		#' The generation of random graph is based on the \code{erdos.renyi.game} function of igraph package.
		#' The numbers of vertices and edges in the random graph are same with the empirical network stored in the object.
		#'
		#' @param runs default 100; simulation number of random network.
		#' @param output_sim default FALSE; whether output each simulated network result.
		#' @return a data.frame with the following components:
		#' \describe{
		#'   \item{\code{Observed}}{Topological properties of empirical network}
		#'   \item{\code{Mean_sim}}{Mean of properties of simulated networks}
		#'   \item{\code{SD_sim}}{SD of properties of simulated networks}
		#'   \item{\code{p_value}}{Significance, i.e. p values}
		#' }
		#' When \code{output_sim = TRUE}, the columns from the five to the last are each simulated result.
		#' @examples
		#' \dontrun{
		#' t1$random_network(runs = 100)
		#' }
		random_network = function(runs = 100, output_sim = FALSE){
			res_network <- self$res_network
			if(is.null(res_network)){
				stop("Please first run cal_network function to obtain a network!")
			}
			if(is.null(self$res_network_attr)){
				message("Topological properties have not been calculated! First run cal_network_attr function ...")
				self$cal_network_attr()
			}
			res <- data.frame(emp = self$res_network_attr[, 1])
			print("Start simulations:")
			for(i in seq_len(runs)){
				if(i %% 10 == 0){
					print(i)
				}
				rand <- igraph::erdos.renyi.game(vcount(res_network), ecount(res_network), type = 'gnm')
				suppressMessages(tmp <- trans_network$new(dataset = NULL))
				suppressMessages(tmp$cal_network(network_method = NULL))
				tmp$res_network <- rand
				suppressMessages(tmp$cal_network_attr())
				res %<>% cbind(., tmp$res_network_attr)
			}
			colnames(res) <- c("emp", paste0("sim", seq_len(runs)))
			sig_right <- apply(res, 1, function(x){(sum(x >= x[1]))/length(x)})
			sig_left <- apply(res, 1, function(x){(sum(x <= x[1]))/length(x)})
			p_value <- sapply(seq_along(sig_right), function(x){ifelse(sig_right[x] <= sig_left[x], sig_right[x], sig_left[x])})
			output <- data.frame(Observed = self$res_network_attr[, 1], 
				Mean_sim = apply(res[, 2:ncol(res)], 1, mean),
				SD_sim = apply(res[, 2:ncol(res)], 1, sd),
				p_value = p_value
			)
			if(output_sim){
				output %<>% cbind.data.frame(., res[, 2:ncol(res)])
			}
			output
		},
		#' @description
		#' Transform classifed features to community-like microtable object for further analysis, such as module-taxa table.
		#'
		#' @param use_col default "module"; which column to use as the 'community'; must be one of the name of res_node_table from function \code{get_node_table}.
		#' @param abundance default TRUE; whether sum abundance of taxa. TRUE: sum the abundance for a taxon across all samples; 
		#' 	  FALSE: sum the frequency for a taxon across all samples.
		#' @return a new \code{\link{microtable}} class.
		#' @examples
		#' \donttest{
		#' t2 <- t1$trans_comm(use_col = "module")
		#' }
		trans_comm = function(use_col = "module", abundance = TRUE){
			private$check_igraph()
			if(use_col == "module"){
				if(is.null(V(self$res_network)$module)){
					stop("Please first run cal_module function to get node modules!")
				}
			}
			if(is.null(self$res_node_table)){
				message("Run get_node_table function to get the node property table ...")
				self$get_node_table()
			}
			if(!use_col %in% colnames(self$res_node_table)){
				stop("Provided use_col must be one of the colnames of object$res_node_table !")
			}
			if(inherits(self$res_node_table[, use_col], "numeric")){
				stop("The selected column-", use_col, " must not be numeric! Please check it!")
			}
			res_node_table <- self$res_node_table
			if(any(is.na(res_node_table[, use_col]))){
				message("Filter the taxa with NA in ", use_col, " ...")
				res_node_table %<>% .[!is.na(.[, use_col]), ]
			}
			abund_table <- self$data_abund
			if(abundance == F){
				abund_table[abund_table > 1] <- 1
			}
			tax_table <- self$tax_table
			feature_abund <- apply(abund_table, 2, sum)
			tm1 <- cbind.data.frame(res_node_table[, c("name", use_col)], abund = feature_abund[res_node_table$name])
			tm2 <- reshape2::dcast(tm1, reformulate(use_col, "name"), value.var = "abund")
			tm2[is.na(tm2)] <- 0
			rownames(tm2) <- tm2[, 1]
			tm2 %<>% .[, -1]
			microtable$new(otu_table = tm2, tax_table = tax_table, auto_tidy = TRUE)
		},
		#' @description
		#' Print the trans_network object.
		print = function() {
			cat("trans_network object:\n")
			if(!is.null(self$res_network)){
				cat("res_network: \n")
				igraph::print.igraph(self$res_network)
				cat("\n")
			}else{
				cat("res_network: NULL\n")
			}
			invisible(self)
		}
		),
	private = list(
		check_filter_number = function(input, param = "filter_thres"){
			if(nrow(input) == 0){
				stop("After filtering, no feature is remained! Please try to lower ", param, "!")
			}
			if(nrow(input) == 1){
				stop("After filtering, only one feature is remained! Please try to lower ", param, "!")
			}
			message("After filtering, ", nrow(input), " features are remained ...")
		},
		check_NetCoMi = function(){
			if(!require("NetCoMi")){
				stop("NetCoMi package is not installed! Please see https://github.com/stefpeschel/NetCoMi ")
			}
		},
		check_igraph = function(){
			if(!require("igraph")){
				stop("Please first install igraph package with the command: install.packages('igraph') !")
			}
		},
		check_network = function(){
			if(is.null(self$res_network)){
				stop("No network found! Please first run cal_network function!")
			}
		},
		# x must be symmetrical matrix
		get_cor_p_list = function(x){
			res_p <- x
			res_p[res_p != 0] <- 0
			list(cor = x, p = res_p)
		},
		# convert long format to symmetrical matrix
		# The first and second columns must be names
		vec2mat = function(datatable, use_names, value_var, rep_value){
			if(!inherits(datatable[, value_var], "numeric")){
				datatable[, value_var] %<>% as.numeric
			}
			use_table <- datatable[, c(1:2)]
			use_table[, value_var] <- datatable[, value_var]
			colnames(use_table) <- c("t1", "t2", "value")
			res_table <- rbind.data.frame(
				use_table, 
				data.frame(t1 = use_table[, 2], t2 = use_table[, 1], value = use_table[, 3]), 
				data.frame(t1 = use_names, t2 = use_names, value = rep(rep_value, length(use_names)))
			)
			res_table <- reshape2::dcast(res_table, t1~t2, value.var = "value") %>% 
				`row.names<-`(.[,1]) %>% 
				.[, -1] %>% 
				.[use_names, use_names] %>% 
				as.matrix
			res_table
		},
		cal_corr = function(inputtable, cor_method){
			use_names <- colnames(inputtable)
			com_res <- combn(use_names, 2)
			res <- lapply(seq_len(ncol(com_res)), function(x){
				cor_res <- suppressWarnings(cor.test(inputtable[, com_res[1, x]], inputtable[, com_res[2, x]], method = cor_method));
				c(t1 = com_res[1, x], t2 = com_res[2, x], cor = unname(cor_res$estimate), p = unname(cor_res$p.value))
			})
			res <- do.call(rbind, res) %>% 
				as.data.frame(stringsAsFactors = FALSE)
			
			res_cor <- private$vec2mat(datatable = res, use_names = use_names, value_var = "cor", rep_value = 1)
			res_p <- private$vec2mat(datatable = res, use_names = use_names, value_var = "p", rep_value = 0)
			res <- list(cor = res_cor, p = res_p)
			res
		},
		# RMT optimization
		rmt = function(cormat, low_thres = 0.4, high_thres = 0.8, seq_by = 0.01){
			s <- seq(0, 3, 0.1)
			pois <- exp(-s)
			ps <- c()
			seq_value <- c()
			for(i in seq(low_thres, high_thres, seq_by)){
				cormat1 <- abs(cormat)
				cormat1[cormat1 < i] <- 0  
				eigen_res <- sort(eigen(cormat1)$value)
				check_res <- tryCatch(ssp <- smooth.spline(eigen_res, control.spar = list(low = 0, high = 3)), error = function(e) { skip_to_next <- TRUE})
				if(rlang::is_true(check_res)) {
					next
				}
				nnsd1 <- density(private$nnsd(ssp$y))
				nnsdpois <- density(private$nnsd(pois))
				chival1 <- sum((nnsd1$y - nnsdpois$y)^2/nnsdpois$y/512)
				ps <- c(ps, chival1)
				seq_value <- c(seq_value, i)
				if((i*100) %% 5 == 0){
					print(i)
				}
			}
			res <- data.frame(ps, seq_value)
			tc <- res[which.min(res[, 1]), 2]
			tc
		},
		nnsd = function(x){
			abs(diff(x))
		},
		saveAsGEXF = function(network, filepath = "network.gexf"){
			require("rgexf")
			nodes <- data.frame(cbind(V(network), V(network)$name))
			edges <- get.edges(network, 1:ecount(network))
			node_attr_name <- setdiff(list.vertex.attributes(network), "name")
			node_attr <- data.frame(sapply(node_attr_name, function(attr) sub("&", "&", get.vertex.attribute(network, attr))))
			node_attr$RelativeAbundance %<>% as.numeric
			edge_attr_name <- setdiff(list.edge.attributes(network), "weight")
			edge_attr <- data.frame(sapply(edge_attr_name, function(attr) sub("&", "&", get.edge.attribute(network, attr))))
			# combine all graph attributes into a meta-data
			graphAtt <- sapply(list.graph.attributes(network), function(attr) sub("&", "&",get.graph.attribute(network, attr)))
			output_gexf <- write.gexf(nodes, edges,
				edgesLabel = as.data.frame(E(network)$label),
				edgesWeight = E(network)$weight,
				nodesAtt = node_attr,
				edgesAtt = edge_attr,
				meta=c(list(creator="trans_network class", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
			cat(output_gexf$graph, file = filepath)
		},
		network_attribute = function(x){
			if(is_directed(x)){
				ms <- cluster_walktrap(x)
			}else{
				ms <- cluster_fast_greedy(x)
			}
			res <- data.frame(
				Vertex = round(vcount(x), 0), 
				Edge = round(ecount(x), 0), 
				Average_degree = sum(igraph::degree(x))/length(igraph::degree(x)), 
				Average_path_length = average.path.length(x), 
				Network_diameter = round(diameter(x, directed = FALSE), 0), 
				Clustering_coefficient = transitivity(x), 
				Density = graph.density(x), 
				Heterogeneity = sd(igraph::degree(x))/mean(igraph::degree(x)), 
				Centralization = centr_degree(x)$centralization,
				Modularity = modularity(ms)
				)
			res <- base::as.data.frame(t(res))
			colnames(res) <- NULL
			res
		},
		# modified based on microbiomeSeq (http://www.github.com/umerijaz/microbiomeSeq) 
		module_roles = function(comm_graph){
			td <- igraph::degree(comm_graph) %>% data.frame(taxa = names(.), total_links = ., stringsAsFactors = FALSE)
			wmd <- private$within_module_degree(comm_graph)
			z <- private$zscore(wmd)
			# NaN may generate in Zi for modules with very few nodes
			if(any(is.nan(z$z))){
				message('The nodes (', sum(is.nan(z$z)),') with NaN in z will be filtered ...')
			}
			amc <- private$among_module_connectivity(comm_graph)
			pc <- private$participation_coeffiecient(amc, td)
			zp <- data.frame(z, pc[rownames(z), , drop = FALSE])
			nod_roles <- private$assign_module_roles(zp)
			nod_roles
		},
		#compute within-module degree
		within_module_degree = function(comm_graph){
			mods <- get.vertex.attribute(comm_graph, "module")
			if(is.null(mods)){
				stop("No modules found! Please first calculate network modules using function cal_module !")
			}
			modvs <- data.frame("taxon" = V(comm_graph)$name, "mod" = mods, stringsAsFactors = FALSE)
			res <- data.frame()
			for(mod in unique(modvs$mod)){
				mod_nodes <- subset(modvs$taxon, modvs$mod == mod)
				mod_subgraph <- induced.subgraph(graph = comm_graph, vids = mod_nodes)
				mod_degree <- igraph::degree(mod_subgraph)
				for(i in mod_nodes){
					ki <- mod_degree[names(mod_degree) == i]
					tmp <- data.frame(module = mod, taxa = names(ki), mod_links = ki)
					res <- rbind(res, tmp)
				}
			}
			res
		},
		#compute within-module degree z-score which
		#measures how well-connected a node is to other nodes in the module.
		zscore = function(within_md){
			ksi_bar <- aggregate(mod_links ~ module, data = within_md, FUN = mean)
			ksi_sigma <- aggregate(mod_links ~ module, data = within_md, FUN = sd)
			z <- NULL
			for(i in 1:dim(within_md)[1]){
				mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == within_md$module[i])]
				mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == within_md$module[i])]
				z[i] <- (within_md$mod_links[i] - mod_mean)/mod_sig
			}
			z <- data.frame(row.names = rownames(within_md), z, module = within_md$module)
			z
		},
		#calculate the degree (links) of each node to nodes in other modules
		among_module_connectivity = function(comm_graph){
			modvs <- data.frame(taxon = V(comm_graph)$name, mod = get.vertex.attribute(comm_graph, "module"), stringsAsFactors = FALSE)
			edges <- t(sapply(1:ecount(comm_graph), function(x) ends(comm_graph, x)))
			res <- lapply(modvs$taxon, function(x){
						sapply(unique(c(edges[edges[, 1] == x, 2], edges[edges[, 2] == x, 1])), function(y){
							c(taxa = x, module = modvs[modvs$taxon == y, "mod"])
						})
					})
			res <- do.call(cbind, res) %>% 
				t %>% 
				as.data.frame(stringsAsFactors = FALSE) %>% 
				dplyr::group_by(taxa) %>% 
				dplyr::count(module) %>% 
				as.data.frame(stringsAsFactors = FALSE)
			colnames(res)[colnames(res) == "n"] <- c("mod_links")
			res
		},
		#The participation coefficient of a node measures how well a node is distributed
		# in the entire network. It is close to 1 if its links are uniformly
		#distributed among all the modules and 0 if all its links are within its own module.
		participation_coeffiecient = function(among_mc, total_degree){
			p <- c()
			for(i in total_degree$taxa){
				ki <- subset(total_degree$total_links, total_degree$taxa == i)
				taxa_each_mod_degree <- subset(among_mc$mod_links, among_mc$taxa == i)
				p[i] <- 1 - (sum((taxa_each_mod_degree)^2)/ki^2)
			}
			as.data.frame(p)
		},
		assign_module_roles = function(zp){
			zp <- na.omit(zp)
			if(nrow(zp) == 0){
				stop("No valid result obtained for this network!")
			}
			zp$taxa_roles <- rep("", nrow(zp))
			for(i in 1:nrow(zp)){
				if(zp[i, "z"] <= 2.5){
					if(zp[i, "p"] <= 0.62){
						zp[i, "taxa_roles"] <- "Peripheral nodes"
					}else{
						zp[i, "taxa_roles"] <- "Connectors"
					}
				}
				else{
					if(zp[i, "p"] <= 0.62){
						zp[i, "taxa_roles"] <- "Module hubs"
					}else{
						zp[i, "taxa_roles"] <- "Network hubs"
					}
				}
			}
			zp
		},
		plot_roles_1 = function(node_roles, roles_color_background, add_label, add_label_group, add_label_text, 
			label_text_size, label_text_color, label_text_italic, label_text_parse,
			roles_color_values, module, x_lim, ...
			){
			if(module == T){
				all_modules <- unique(as.character(node_roles$module))
				node_roles$module <- factor(as.character(node_roles$module), levels = stringr::str_sort(all_modules, numeric = TRUE))
			}
			x1 <- c(0, 0.62, 0, 0.62)
			x2 <- c(0.62, 1, 0.62, 1)
			y1 <- c(-Inf,-Inf, 2.5, 2.5)
			y2 <- c(2.5,2.5, Inf, Inf)
			lab <- c("Peripheral nodes", "Connectors", "Module hubs", "Network hubs")
			lab <- factor(lab, levels = rev(lab))
			if(is.null(roles_color_values)){roles_color_values <- rev(c("grey80", RColorBrewer::brewer.pal(3, "Dark2")))}
			
			p <- ggplot() + theme_bw()
			if(roles_color_background){
				p <- p + geom_rect(data=NULL, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab, alpha = .4))
				p <- p + guides(fill = guide_legend(title = "Roles"), alpha = "none")
				p <- p + scale_fill_manual(values = roles_color_values)
				if(module == T){
					p <- p + geom_point(data = node_roles, aes(x = p, y = z, shape = module), ...) + 
						guides(shape = guide_legend(title = "Module"))
				}else{
					p <- p + geom_point(data = node_roles, aes(x = p, y = z), ...)
				}
				p <- p + theme(strip.background = element_rect(fill = "white"))
			}else{
				node_roles$taxa_roles %<>% factor(., levels = rev(lab))
				if(module == T){
					p <- p + geom_point(data = node_roles, aes(x = p, y = z, color = taxa_roles, shape = module), ...) + 
						guides(shape = guide_legend(title = "Module"))
				}else{
					p <- p + geom_point(data = node_roles, aes(x = p, y = z, color = taxa_roles), ...)
				}
				p <- p + xlim(x_lim[1], x_lim[2])
				p <- p + geom_hline(yintercept = 2.5, linetype = "dashed") +
					geom_vline(xintercept = 0.62, linetype = "dashed") +
					scale_color_manual(values = roles_color_values)
				p <- p + guides(color = guide_legend(title = "Roles"), alpha = "none")
			}
			if(add_label){
				if(!add_label_text %in% colnames(node_roles)){
					stop("The provided add_label_text parameter must be one of colnames of object$res_node_table!")
				}else{
					label_data <- node_roles[node_roles$taxa_roles %in% add_label_group, ]
					if(nrow(label_data) == 0){
						message("No label need to be added!")
					}else{
						if(label_text_italic == T){
							label_data[, add_label_text] %<>% paste0("italic('", .,"')")
							label_text_parse <- TRUE
						}
						p <- p + ggrepel::geom_text_repel(
							data = label_data, 
							aes(.data[["p"]], .data[["z"]], label = .data[[add_label_text]]), 
							size = label_text_size, 
							color = label_text_color, 
							segment.alpha = .01, 
							parse = label_text_parse
						)
					}
				}
			}
			p <- p + xlab("Among-module connectivity") + 
				ylab("Within-module connectivity") +
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

			p
		},
		# plot z and p according to the taxonomic levels
		plot_roles_2 = function(
			node_roles, 
			use_level = "Phylum", 
			show_value = c("z", "p"), 
			show_number = 1:10, 
			plot_color = "Phylum", 
			plot_shape = "taxa_roles",
			plot_size = "Abundance",
			color_values = RColorBrewer::brewer.pal(12, "Paired"),
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14),
			...
			){
			node_roles <- cbind.data.frame(ID = rownames(node_roles), node_roles)
			use_data <- reshape2::melt(node_roles, id.vars = c("ID", use_level), measure.vars = show_value, variable.name = "variable")
			node_roles[, use_level] <- NULL
			use_data <- suppressWarnings(dplyr::left_join(use_data, node_roles, by=c("ID" = "ID")))
			use_data %<>% dropallfactors
			# replace the names in strip text
			if(any(show_value %in% "z")){
				use_data[use_data$variable == "z", "variable"] <- "Within-module connectivity"
			}
			if(any(show_value %in% "p")){
				use_data[use_data$variable == "p", "variable"] <- "Among-module connectivity"
			}
			# filter, subset and sort the shown taxa
			use_data <- use_data[!grepl("__$", as.character(use_data[, use_level])), ]
			use_data[, use_level] %<>% gsub(".__", "", .)
			use_taxa <- table(as.character(use_data[, use_level])) %>% 
				sort(decreasing = TRUE) %>% 
				names %>% 
				.[show_number]
			use_data %<>% .[as.character(.[, use_level]) %in% use_taxa, ]
			use_data[, use_level] %<>% factor(., levels = use_taxa)

			p <- ggplot(use_data, aes_meco(x = use_level, y = "value", colour = plot_color, shape = plot_shape, size = plot_size)) + 
				geom_point(position = "jitter", ...) +
				scale_color_manual(values = color_values) +
				scale_shape_manual(values = shape_values) +				
				facet_grid(variable ~ ., drop = TRUE, scale = "free", space = "fixed") +
				theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 10))
			
			if(plot_color == use_level){
				p <- p + guides(color = "none")
			}
			if(plot_shape == use_level){
				p <- p + guides(shape = "none")
			}
			if(!is.null(plot_size)){
				p <- p + guides(size = guide_legend(title = "Abundance(%)"))
			}
			p
		},
		sum_link = function(taxa_table, link_table, taxa_level){
			# first obtain the taxa names
			all_names <- taxa_table[rownames(taxa_table) %in% unique(c(link_table[,1], link_table[,2])), ] %>%
				{table(.[, taxa_level])} %>%
				sort(., decreasing = TRUE) %>% 
				names
			com_group <- expand.grid(all_names, all_names)
			colnames(com_group) <- c("C1", "C2")
			# assign rownames irrespective of the order
			rownames(com_group) <- apply(com_group, 1, function(x) paste0(x, collapse = "-"))
			# get the unifrom combined name without regard to the order
			com_group$uni_name <- apply(com_group, 1, function(x) paste0(sort(x), collapse = "-"))
			com_group1 <- com_group[, -c(1, 2), drop = FALSE]
			res <- link_table
			# use taxa name to replace the species name
			res[, 1] <- taxa_table[res[, 1], taxa_level]
			res[, 2] <- taxa_table[res[, 2], taxa_level]
			res$pname <- paste(res[, 1], res[, 2], sep = "-")
			res %<>% dplyr::group_by(pname) %>% 
				dplyr::summarise(count = dplyr::n()) %>%
				as.data.frame(stringsAsFactors = FALSE)
			res <- dplyr::left_join(res, rownames_to_column(com_group1), by = c("pname" = "rowname")) %>%
				dplyr::group_by(uni_name) %>% 
				dplyr::summarise(sum_count = sum(count)) %>%
				as.data.frame(stringsAsFactors = FALSE)
			res <- dplyr::left_join(res, com_group, by = c("uni_name" = "uni_name"))
			res <- reshape2::dcast(res, C1~C2, value.var = "sum_count") %>%
				`row.names<-`(.[,1]) %>%
				.[, -1, drop = FALSE] %>%
				.[all_names, all_names, drop = FALSE] %>%
				as.matrix
			res[is.na(res)] <- 0			
			res
		},
		rm_unlinked_node = function(input_network){
			nodes_raw <- V(input_network)$name
			edges <- t(sapply(1:ecount(input_network), function(x) ends(input_network, x)))
			delete_nodes <- nodes_raw %>% .[! . %in% as.character(c(edges[,1], edges[,2]))]
			if(length(delete_nodes) > 0){
				input_network %<>% delete_vertices(delete_nodes)
			}
			input_network
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
