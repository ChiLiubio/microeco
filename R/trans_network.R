#' @title
#' Create trans_network object for co-occurrence network analysis.
#'
#' @description
#' This class is a wrapper for a series of network analysis methods, 
#' including the network construction approaches, network attributes analysis,
#' eigengene analysis, network subsetting, node and edge properties extraction, network plotting, and other network operations.
#'
#' @export
trans_network <- R6Class(classname = "trans_network",
	public = list(
		#' @description
		#' This function is used to create the trans_network object, store the important intermediate data 
		#'   and calculate correlations if cal_cor parameter is selected.
		#' 
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param cor_method default "pearson"; "pearson", "spearman" or "kendall"; correlation algorithm, only use for correlation-based network.
		#' @param cal_cor default "base"; "base", "WGCNA", "SparCC" or NA; correlation method; NA represents no correlation calculation, 
		#' 	  used for non-correlation based network, such as SpiecEasi and FlashWeave methods.
		#' @param taxa_level default "OTU"; taxonomic rank; 'OTU' represents using feature table directly; 
		#' 	  other available options should be one of the colnames of microtable$tax_table.
		#' @param filter_thres default 0; the relative abundance threshold.
		#' @param nThreads default 1; the thread number used for "WGCNA" and SparCC.
		#' @param SparCC_simu_num default 100; SparCC simulation number for bootstrap.
		#' @param env_cols default NULL; numeric or character vector to select the column names of environmental data in dataset$sample_table;
		#'   the environmental data can be used in the correlation network (as the nodes) or FlashWeave network.
		#' @param add_data default NULL; provide environmental table additionally instead of env_cols parameter; rownames must be sample names.
		#' @return res_cor_p list; include the correlation matrix and p value matrix.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' # for correlation network
		#' t1 <- trans_network$new(dataset = dataset, cal_cor = "base", 
		#' 		taxa_level = "OTU", filter_thres = 0.0001)
		#' # for other network
		#' t1 <- trans_network$new(dataset = dataset, cal_cor = NA)
		#' }
		initialize = function(
			dataset = NULL,
			cor_method = c("pearson", "spearman", "kendall")[1],
			cal_cor = c("base", "WGCNA", "SparCC", NA)[1],
			taxa_level = "OTU",
			filter_thres = 0,
			nThreads = 1,
			SparCC_simu_num = 100,
			env_cols = NULL,
			add_data = NULL
			){
			#cor_method <- match.arg(cor_method)
			dataset1 <- clone(dataset)
			if(!is.null(env_cols)){
				env_data <- dataset1$sample_table[, env_cols, drop = FALSE]
			}
			if(!is.null(add_data)){
				env_data <- add_data[rownames(add_data) %in% rownames(dataset1$sample_table), ]
			}
			if(!is.null(env_cols) | !is.null(add_data)){
				dataset1$sample_table %<>% .[rownames(.) %in% rownames(env_data), ]
				dataset1$tidy_dataset(main_data = TRUE)
				env_data %<>% .[rownames(dataset1$sample_table), ] %>%
					dropallfactors(unfac2num = TRUE)
				env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
				env_data[] <- lapply(env_data, as.numeric)
				self$env_data <- env_data
			}
			if(taxa_level != "OTU"){
				dataset1 <- dataset1$merge_taxa(taxa = taxa_level)
			}
			# transform each object
			self$use_sampleinfo <- dataset1$sample_table
			# store taxonomic table for the following analysis
			self$use_tax <- dataset1$tax_table
			use_abund <- dataset1$otu_table %>% 
				{.[apply(., 1, sum)/sum(.) > filter_thres, ]} %>%
				t %>%
				as.data.frame
			
			if( (!is.na(cal_cor)) & (!is.null(env_cols) | !is.null(add_data))){
				use_abund <- cbind.data.frame(use_abund, env_data)
			}
			if(!is.na(cal_cor)){
				if(cal_cor == "base"){
					cor_result <- private$cal_corr(inputtable = use_abund, cor_method = cor_method)
				}
				if(cal_cor == "WGCNA"){
					cor_result <- WGCNA::corAndPvalue(x = use_abund, method = cor_method, nThreads = nThreads)
				}
				if(cal_cor == "SparCC"){
					try_find <- try(find.package("SpiecEasi"), silent = TRUE)
					if(inherits(try_find, "try-error")){
						stop("SpiecEasi package is used for the SparCC calculation, but it is not installed! See https://github.com/zdk123/SpiecEasi for the installation")
					}
					bootres <- SpiecEasi::sparccboot(use_abund, ncpus = nThreads, R = SparCC_simu_num)
					cor_result <- SpiecEasi::pval.sparccboot(bootres)
					# reshape the results
					use_names <- colnames(bootres$data)
					com_res <- t(combn(use_names, 2))
					res <- cbind.data.frame(com_res, cor = cor_result$cors, p = cor_result$pvals, stringsAsFactors = FALSE)
					res_cor <- private$vec2mat(datatable = res, use_names = use_names, value_var = "cor", rep_value = 1)
					res_p <- private$vec2mat(datatable = res, use_names = use_names, value_var = "p", rep_value = 0)
					cor_result <- list(cor = res_cor, p = res_p)
				}
				self$res_cor_p <- cor_result
				message('The correlation result list is stored in object$res_cor_p ...')
			}else{
				self$res_cor_p <- NULL
			}
			
			self$use_abund <- use_abund
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Calculate network based on the correlation method or SpiecEasi package or julia FlashWeave package or beemStatic package.
		#'
		#' @param network_method default "COR"; "COR", "SpiecEasi", "FlashWeave" or "beemStatic"; The option details: 
		#'   \describe{
		#'     \item{\strong{'COR'}}{correlation-based network; use the correlation and p value matrixes in object$res_cor_p returned from trans_network$new; 
		#'     	  See Deng et al. (2012) <doi:10.1186/1471-2105-13-113> for other details}
		#'     \item{\strong{'SpiecEasi'}}{SpiecEasi network; relies on algorithms for sparse neighborhood and inverse covariance selection;
		#'     	  belong to the category of conditional dependence and graphical models;
		#'     	  see \href{https://github.com/zdk123/SpiecEasi}{https://github.com/zdk123/SpiecEasi} for installing the R package; 
		#'     	  see Kurtz et al. (2015) <doi:10.1371/journal.pcbi.1004226> for the algorithm details}
		#'     \item{\strong{'FlashWeave'}}{FlashWeave network; Local-to-global learning framework; belong to the category of conditional dependence and graphical models;
		#'        good performance on heterogenous datasets to find direct associations among taxa;
		#'        see \href{https://github.com/meringlab/FlashWeave.jl}{https://github.com/meringlab/FlashWeave.jl} for installing julia language and FlashWeave package;
		#'        julia must be in the computer system env path, otherwise the program can not find julia;
		#'        see Tackmann et al. (2019) <doi:10.1016/j.cels.2019.08.002> for the algorithm details}
		#'     \item{\strong{'beemStatic'}}{beemStatic network;
		#'        extend generalized Lotka-Volterra model to cases of cross-sectional datasets to infer interaction among taxa based on expectation-maximization algorithm;
		#'        see \href{https://github.com/CSB5/BEEM-static}{https://github.com/CSB5/BEEM-static} for installing the R package;
		#'        see Li et al. (2021) <doi:10.1371/journal.pcbi.1009343> for algorithm details}
		#'   }
		#' @param COR_p_thres default 0.01; the p value threshold for the correlation-based network.
		#' @param COR_p_adjust default "fdr"; p value adjustment method, see method of p.adjust function for available options.
		#' @param COR_weight default TRUE; whether use correlation coefficient as the weight of edges; FALSE represents weight = 1 for all edges.
		#' @param COR_cut default 0.6; correlation coefficient threshold for the correlation network.
		#' @param COR_optimization default FALSE; whether use random matrix theory to optimize the choice of correlation coefficient, see https://doi.org/10.1186/1471-2105-13-113
		#' @param COR_low_threshold default 0.4; the lowest correlation coefficient threshold, only useful when COR_optimization = TRUE.
		#' @param SpiecEasi_method default "mb"; either 'glasso' or 'mb';see spiec.easi function in package SpiecEasi and https://github.com/zdk123/SpiecEasi.
		#' @param FlashWeave_tempdir default NULL; The temporary directory used to save the temporary files for running FlashWeave; If not assigned, use the system user temp.
		#' @param FlashWeave_meta_data default FALSE; whether use env data for the optimization, If TRUE, the function automatically find the object$env_data in the object and
		#'   generate a file for meta_data_path parameter of FlashWeave.
		#' @param FlashWeave_other_para default "alpha=0.01,sensitive=true,heterogeneous=true"; the parameters used for FlashWeave;
		#'   user can change the parameters or add more according to FlashWeave help document;
		#'   An exception is meta_data_path parameter as it is generated based on the data inside the object, see FlashWeave_meta_data parameter for the description.
		#' @param beemStatic_t_strength default 0.001; for network_method = "beemStatic"; the threshold used to limit the number of interactions (strength);
		#'   same with the t.strength parameter in showInteraction function of beemStatic package.
		#' @param beemStatic_t_stab default 0.8; for network_method = "beemStatic"; 
		#'   the threshold used to limit the number of interactions (stability); same with the t.stab parameter in showInteraction function of beemStatic package.
		#' @param add_taxa_name default "Phylum"; NULL or a taxonomic rank name; used to add taxonomic rank name to network node properties.
		#' @param usename_rawtaxa_when_taxalevel_notOTU default FALSE; whether replace the name of nodes using the taxonomic information.
		#' @param ... paremeters pass to spiec.easi function of SpiecEasi package when network_method = "SpiecEasi" or 
		#'   func.EM function of beemStatic package when network_method = "beemStatic".
		#' @return res_network stored in object.
		#' @examples
		#' \donttest{
		#' # for correlation network
		#' t1 <- trans_network$new(dataset = dataset, cal_cor = "base", 
		#' 		taxa_level = "OTU", filter_thres = 0.0001)
		#' t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
		#' t1 <- trans_network$new(dataset = dataset, cal_cor = NA)
		#' t1$cal_network(network_method = "SpiecEasi")
		#' t1$cal_network(network_method = "beemStatic")
		#' t1$cal_network(network_method = "FlashWeave")
		#' }
		cal_network = function(
			network_method = c("COR", "SpiecEasi", "FlashWeave", "beemStatic")[1],
			COR_p_thres = 0.01,
			COR_p_adjust = "fdr",
			COR_weight = TRUE,
			COR_cut = 0.6,
			COR_optimization = FALSE,
			COR_low_threshold = 0.4,
			SpiecEasi_method = "mb",
			FlashWeave_tempdir = NULL,
			FlashWeave_meta_data = FALSE,
			FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=true",
			beemStatic_t_strength = 0.001,
			beemStatic_t_stab = 0.8,
			add_taxa_name = "Phylum",
			usename_rawtaxa_when_taxalevel_notOTU = FALSE,
			...
			){
			private$check_igraph()
			sampleinfo <- self$use_sampleinfo
			taxa_level <- self$taxa_level
			taxa_table <- self$use_tax
			
			network_method <- match.arg(network_method, c("COR", "SpiecEasi", "FlashWeave", "beemStatic"))
			
			message("---------------- ", Sys.time()," : Start ----------------")
			if(network_method == "COR"){
				cortable <- self$res_cor_p$cor
				# p adjustment for the converted vector
				raw_p <- self$res_cor_p$p
				if(ncol(cortable) != ncol(raw_p)){
					stop("Correlation table and p value table have different column numbers !")
				}
				raw_vector_p <- raw_p %>% as.dist %>% as.numeric
				adp_raw <- p.adjust(raw_vector_p, method = COR_p_adjust)
				# to matrix
				use_names <- colnames(raw_p)
				names_combn <- t(combn(use_names, 2))
				table_convert <- cbind.data.frame(names_combn, adjust.p = adp_raw, stringsAsFactors = FALSE)
				adp <- private$vec2mat(datatable = table_convert, use_names = use_names, value_var = "adjust.p", rep_value = 0)
				# make sure same names between cortable and adp
				if(! identical(colnames(cortable), colnames(adp))){
					adp <- adp[colnames(cortable), colnames(cortable)]
				}
				if(COR_optimization == T) {
					#find out threshold of correlation 
					tc1 <- private$rmt(cortable)
					tc1 <- ifelse(tc1 > COR_low_threshold, tc1, COR_low_threshold)
					message("The optimized COR threshold: ", tc1, "...\n")
				}
				else {
					tc1 <- COR_cut
				}
				diag(cortable) <- 0
				cor_matrix <- as.matrix(cortable)
				cor_matrix[abs(cortable) >= tc1] <- 1
				cor_matrix[adp >= COR_p_thres] <- 0
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
				use_abund <- self$use_abund %>% as.matrix
				# calculate SpiecEasi network, reference https://github.com/zdk123/SpiecEasi
				network <- spiec.easi(use_abund, method = SpiecEasi_method, ...)
				network <- adj2igraph(getRefit(network))
				V(network)$name <- colnames(use_abund)
				E(network)$label <- unlist(lapply(E(network)$weight, function(x) ifelse(x > 0, "+", "-")))
			}
			if(grepl("FlashWeave", network_method, ignore.case = TRUE)){
				use_abund <- self$use_abund
				oldwd <- getwd()
				# make sure working directory can not be changed by the function when quit.
				on.exit(setwd(oldwd))

				if(is.null(FlashWeave_tempdir)){
					tem_dir <- tempdir()
				}else{
					# check the directory
					tem_dir <- FlashWeave_tempdir
					if(!dir.exists(tem_dir)){
						stop("The input temporary directory: ", tem_dir, " does not exist!")
					}
				}
				setwd(tem_dir)
				write.table(use_abund, "taxa_table_FlashWeave.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
				L1 <- "using FlashWeave\n"
				L2 <- 'data_path = "taxa_table_FlashWeave.tsv"\n'
				if(FlashWeave_meta_data == T){
					if(is.null(self$env_data)){
						stop("FlashWeave_meta_data is TRUE, but object$env_data not found! 
							Please use env_cols or add_data parameter of trans_network$new to provide the metadata when creating the object!")
					}
					meta_data <- self$env_data
					write.table(meta_data, "meta_table_FlashWeave.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
					L3 <- 'meta_data_path = "meta_table_FlashWeave.tsv"\n'
				}else{
					L3 <- "\n"
				}
				if(FlashWeave_meta_data == T){
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
				message("Run the FlashWeave ...")
				system("julia calculate_network.jl")
				network <- read_graph("network_FlashWeave.gml", format = "gml")
				network <- set_vertex_attr(network, "name", value = V(network)$label)
				E(network)$label <- unlist(lapply(E(network)$weight, function(x) ifelse(x > 0, "+", "-")))
				E(network)$weight <- abs(E(network)$weight)
			}
			if(network_method == "beemStatic"){
				if(!require("beemStatic")){
					stop("beemStatic package is not installed! See https://github.com/CSB5/BEEM-static ")
				}
				use_abund <- self$use_abund %>% t %>% as.data.frame
				taxa_low <- apply(use_abund, 1, function(x) sum(x != 0)) %>% .[which.min(.)]
				message("The feature table have ", nrow(use_abund), " taxa. The taxa with the lowest occurrence frequency was ", names(taxa_low)[1], 
					", found in ", taxa_low[1], " samples of total ", ncol(use_abund), " samples. If an error occurs because of low frequency, ",
					"please filter more taxa with low abundance using filter_thres parameter when creating the trans_network object ...")
				beem.out <- func.EM(use_abund, ...)
				self$res_beemStatic_raw <- beem.out
				message('beemStatic result is stored in object$res_beemStatic_raw ...')
				# modified based on the showInteraction function
				b <- t(beem2param(beem.out)$b.est)
				diag(b) <- 0
				if(!is.null(beem.out$resample)){
					b[t(beem.out$resample$b.stab < beemStatic_t_stab)] <- 0
				}
				b[abs(b) < beemStatic_t_strength] <- 0
				network <- graph.adjacency(b, mode='directed', weighted='weight')
				V(network)$name <- rownames(use_abund)
				V(network)$RelativeAbundance <- rowMeans(beemStatic:::tss(use_abund))
				E(network)$label <- ifelse(E(network)$weight > 0, '+', '-')
				E(network)$weight <- abs(E(network)$weight)
			}
			message("---------------- ", Sys.time()," : Finish ----------------")
			
			nodes_raw <- data.frame(cbind(V(network), V(network)$name))
			# delete uncultured taxa when the taxa level is not OTU
			if(taxa_level != "OTU"){
				delete_nodes <- taxa_table %>% 
					.[grepl("__$|uncultured", .[, taxa_level]), ] %>% 
					rownames %>% 
					.[. %in% rownames(nodes_raw)]
				network %<>% delete_vertices(delete_nodes)
				nodes_raw <- data.frame(cbind(V(network), V(network)$name))
			}
			edges <- t(sapply(1:ecount(network), function(x) ends(network, x)))
			delete_nodes <- rownames(nodes_raw) %>% 
				.[! . %in% as.character(c(edges[,1], edges[,2]))]
			network %<>% delete_vertices(delete_nodes)
			V(network)$taxa <- V(network)$name
			if(!is.null(add_taxa_name)){
				network <- set_vertex_attr(network, add_taxa_name, value = V(network)$name %>% 
					taxa_table[., add_taxa_name] %>% 
					gsub("^.__", "", .))
			}
			if(taxa_level != "OTU"){
				if(usename_rawtaxa_when_taxalevel_notOTU == T){
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
		},
		#' @description
		#' Calculate network modules and add module names to the network node properties.
		#'
		#' @param method default "cluster_fast_greedy"; the method used to find the optimal community structure of a graph;
		#' 	 the following are available functions (options) from igraph package: "cluster_fast_greedy", "cluster_optimal",
		#' 	 "cluster_edge_betweenness", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen",
		#' 	 "cluster_louvain", "cluster_spinglass", "cluster_walktrap". 
		#' 	 For the details of these functions, see the help document, such as help(cluster_fast_greedy);
		#' 	 Note that the default "cluster_fast_greedy" method can only be used for undirected network. 
		#' 	 If the user selects network_method = "beemStatic" in cal_network function or provides other directed network, 
		#' 	 please use cluster_optimal or others for the modules identification.
		#' @param module_name_prefix default "M"; the prefix of module names; module names are made of the module_name_prefix and numbers;
		#'   numbers are assigned according to the sorting result of node numbers in modules with decreasing trend.
		#' @return res_network with modules, stored in object.
		#' @examples
		#' \donttest{
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
		#' @return None.
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
		#' @return res_network_attr stored in object.
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
		#' Calculate node properties. This function will be deprecated in the next release! Please use get_node_table function!
		#'
		#' @return see the Return part in function get_node_table.
		cal_node_type = function(){
			warning('Please use get_node_table function instead of this! This function will be deprecated in the next release !')
			self$get_node_table(node_roles = TRUE)
		},
		#' @description
		#' Get the node property table. The properties may include the node names, modules allocation, degree, betweenness, abundance, 
		#'   taxonomy, within-module connectivity and among-module connectivity <doi:10.1016/j.geoderma.2022.115866>.
		#'
		#' Authors: Chi Liu, Umer Zeeshan Ijaz
		#'
		#' @param node_roles default TRUE; whether calculate node roles, i.e. Module hubs, Network hubs, Connectors and Peripherals <doi:10.1016/j.geoderma.2022.115866>.
		#' @return res_node_table in object; Abundance expressed as a percentage; z represents within-module connectivity;
		#'   p represents among-module connectivity.		
		#' @examples
		#' \donttest{
		#' t1$get_node_table(node_roles = TRUE)
		#' }
		get_node_table = function(node_roles = TRUE){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			use_abund <- self$use_abund
			node_table <- data.frame(name = V(network)$name) %>% `rownames<-`(.[, 1])
			node_table$degree <- igraph::degree(network)[rownames(node_table)]
			node_table$betweenness <- betweenness(network)[rownames(node_table)]
			# Add abundance info
			sum_abund <- apply(use_abund, 2, function(x) sum(x) * 100/sum(use_abund))
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
					self$use_tax[replace_table[rownames(.), 2], 1:which(colnames(self$use_tax) %in% self$taxa_level), drop = FALSE])
			}else{
				node_table %<>% cbind.data.frame(., self$use_tax[rownames(.), ])
			}
			
			self$res_node_table <- node_table
			message('Result is stored in object$res_node_table ...')
		},
		#' @description
		#' Get the edge property table, including connected nodes, label and weight.
		#'
		#' @return res_edge_table in object.
		#' @examples
		#' \donttest{
		#' t1$get_edge_table()
		#' }
		get_edge_table = function(){
			private$check_igraph()
			private$check_network()
			network <- self$res_network
			# another way:
			# res_edge_table <- as_data_frame(network, what = "edges")
			edges <- t(sapply(1:ecount(network), function(x) ends(network, x)))
			edge_label <- E(network)$label
			if(!is.null(E(network)$weight)){
				edge_weight <- E(network)$weight
			}else{
				edge_weight <- rep(NA, times = length(edge_label))
			}
			res_edge_table <- data.frame(edges, edge_label, edge_weight)
			colnames(res_edge_table) <- c("node1", "node2", "label", "weight")
			self$res_edge_table <- res_edge_table
			message('Result is stored in object$res_edge_table ...')
		},
		#' @description
		#' Get the adjacency matrix from the network graph.
		#'
		#' @param ... parameters passed to as_adjacency_matrix function of igraph package.
		#' @return res_adjacency_matrix in object.
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
		#' Plot the network based on a series of methods from other packages, such as igraph, ggraph and networkD3. 
		#' The networkD3 package provides dynamic network. It is especially useful for a glimpse of the whole network structure and finding 
		#' the interested nodes and edges in a large network. In contrast, the igraph and ggraph methods are suitable for relatively small network.
		#'
		#' @param method default "igraph"; The available options:
		#'   \describe{
		#'     \item{\strong{'igraph'}}{call plot.igraph function in igraph package for a static network; see plot.igraph for the parameters}
		#'     \item{\strong{'ggraph'}}{call ggraph function in ggraph package for a static network}
		#'     \item{\strong{'networkD3'}}{use forceNetwork function in networkD3 package for a dynamic network; see forceNetwork function for the parameters}
		#'   }
		#' @param node_label default "name"; node label shown in the plot for method = "ggraph" or method = "networkD3"; 
		#'   Please see the column names of object$res_node_table, which is the returned table of function object$get_node_table;
		#'   User can select other column names in res_node_table.
		#' @param node_color default NULL; node color assignment for method = "ggraph" or method = "networkD3"; 
		#'   Select a column name of object$res_node_table, such as "module".
		#' @param ggraph_layout default "fr"; for method = "ggraph"; see layout parameter of create_layout function in ggraph package.
		#' @param ggraph_node_size default 2; for method = "ggraph"; the node size.
		#' @param ggraph_text_color default NULL; for method = "ggraph"; a column name of object$res_node_table;
		#'   User can select other column names or change the content of object$res_node_table.
		#' @param ggraph_text_size default 3; for method = "ggraph"; the node label text size.
		#' @param networkD3_node_legend default TRUE; used for method = "networkD3"; logical value to enable node colour legends;
		#'   Please see the legend parameter in networkD3::forceNetwork function.
		#' @param networkD3_zoom default TRUE; used for method = "networkD3"; logical value to enable (TRUE) or disable (FALSE) zooming;
		#'   Please see the zoom parameter in networkD3::forceNetwork function.
		#' @param ... parameters passed to plot.igraph function when method = "igraph" or forceNetwork function when method = "networkD3".
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
				message("Run get_node_table function to get or update the node property table ...")
				self$get_node_table()
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
				g <- g + geom_node_point(aes_string(col = node_color), size = ggraph_node_size, alpha = 0.5) +
					geom_node_text(aes_string(col = ggraph_text_color, label = node_label), size = ggraph_text_size, repel = TRUE) +
					scale_edge_width(range = c(0.5, 2)) +
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
		#' @return res_eigen and res_eigen_expla in object.
		#' @examples
		#' \donttest{
		#' t1$cal_eigen()
		#' }
		cal_eigen = function(){
			private$check_igraph()
			private$check_network()
			use_abund <- self$use_abund
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
		#' @param ... paremeters pass to geom_point.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_taxa_roles(roles_color_background = FALSE)
		#' }
		plot_taxa_roles = function(
			use_type = c(1, 2)[1],
			roles_color_background = FALSE,
			roles_color_values = NULL,
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
		#' @return a new network
		#' @examples
		#' \donttest{
		#' t1$subset_network(node = t1$res_node_table %>% .[.$module == "M1", ] %>% 
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
				stop("Please provide the retained nodes name using node parameter!")
			}
			# whether remove the single node without edges
			if(rm_single == T){
				nodes_raw <- V(sub_network)$name
				edges <- t(sapply(1:ecount(sub_network), function(x) ends(sub_network, x)))
				delete_nodes <- nodes_raw %>% .[! . %in% as.character(c(edges[,1], edges[,2]))]
				if(length(delete_nodes) > 0){
					sub_network %<>% delete_vertices(delete_nodes)
				}
			}
			sub_network
		},
		#' @description
		#' Fit degrees to a power law distribution. First, perform a bootstrapping hypothesis test to determine whether degrees follow a power law distribution.
		#' If the distribution follows power law, then fit degrees to power law distribution and return the parameters.
		#'
		#' @param ... paremeters pass to fit_power_law function in igraph package.
		#' @return res_powerlaw_p and res_powerlaw_fit; see bootstrap_p function in poweRlaw package for the bootstrapping p value details;
		#'   see fit_power_law function in igraph package for the power law fit return details.
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
			bootstrap_res <- poweRlaw::bootstrap_p(resdispl)
			self$res_powerlaw_p <- bootstrap_res
			message('Bootstrap result is stored in object$res_powerlaw_p ...')
			message('Bootstrap p value: ', bootstrap_res$p)
			if(bootstrap_res$p < 0.05){
				message("The p value < 0.05; Degrees do not follow power law distribution ...")
			}else{
				message("Degrees follow power law distribution ...")
				res_powerlaw_fit <- fit_power_law(degree_dis + 1, xmin = est_xmin$xmin, ...)
				message("The estimated alpha: ", res_powerlaw_fit$alpha)
				self$res_powerlaw_fit <- res_powerlaw_fit
				message('Powerlaw fitting result is stored in object$res_powerlaw_fit ...')
			}
		},
		#' @description
		#' Transform classifed features to community-like microtable object for further analysis, such as module-taxa table.
		#'
		#' @param use_col default "module"; which column to use as the 'community'; must be one of the name of res_node_table from function get_node_table.
		#' @return a new \code{\link{microtable}} class.
		#' @examples
		#' \donttest{
		#' t2 <- t1$trans_comm(use_col = "module")
		#' }
		trans_comm = function(use_col = "module"){
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
			abund_table <- self$use_abund
			tax_table <- self$use_tax
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
			cat("trans_network class:\n")
			if(!is.null(self$res_network)){
				cat("network object: \n")
				print.igraph(self$res_network)
				cat("\n")
			}else{
				cat("network object: NULL\n")
			}
			if(!is.null(self$res_network_attr)){
				cat("res_network_attr object: finished\n")
			}else{
				cat("res_network_attr object: NULL\n")
			}
			invisible(self)
		}
		),
	private = list(
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
		rmt = function(cormat, lcor = 0.4, hcor = 0.8){
			s <- seq(0, 3, 0.1)
			pois <- exp(-s)
			ps <- NULL  
			for(i in seq(lcor, hcor, 0.01)){
				cormat1 <- abs(cormat)
				cormat1[cormat1 < i] <- 0  
				eigen_res <- sort(eigen(cormat1)$value)
				ssp <- smooth.spline(eigen_res, control.spar = list(low = 0, high = 3)) 
				nnsd1 <- density(private$nnsd(ssp$y))
				nnsdpois <- density(private$nnsd(pois))
				chival1 <- sum((nnsd1$y - nnsdpois$y)^2/nnsdpois$y/512)
				ps <- rbind(ps, chival1)
			}
			ps <- cbind(ps, c(seq(lcor, hcor, 0.01)))
			tc <- ps[ps[,1] == min(ps[,1]), 2]
			tc
		},
		nnsd = function(x){
			abs(diff(x))
		},
		saveAsGEXF = function(network, filepath = "network.gexf"){
			require("rgexf")
			nodes <- data.frame(cbind(V(network), V(network)$name))
			edges <- get.edges(network, 1:ecount(network))
			vAttrNames <- setdiff(list.vertex.attributes(network), "name")
			nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&", get.vertex.attribute(network, attr))))
			eAttrNames <- setdiff(list.edge.attributes(network), "weight")
			edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&", get.edge.attribute(network, attr))))
			# combine all graph attributes into a meta-data
			graphAtt <- sapply(list.graph.attributes(network), function(attr) sub("&", "&",get.graph.attribute(network, attr)))
			output_gexf <- write.gexf(nodes, edges,
				edgesLabel = as.data.frame(E(network)$label),
				edgesWeight = E(network)$weight,
				nodesAtt = nodesAtt,
				edgesAtt = edgesAtt,
				meta=c(list(creator="trans_network class", description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
			cat(output_gexf$graph, file = filepath)
		},
		network_attribute = function(x){
			res <- data.frame(
				Vertex = round(vcount(x), 0), 
				Edge = round(ecount(x), 0), 
				Average_degree = sum(igraph::degree(x))/length(igraph::degree(x)), 
				Average_path_length = average.path.length(x), 
				Network_diameter = round(diameter(x, directed = FALSE), 0), 
				Clustering_coefficient = transitivity(x), 
				Density = graph.density(x), 
				Heterogeneity = sd(igraph::degree(x))/mean(igraph::degree(x)), 
				Centralization = centr_degree(x)$centralization
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
			amd <- private$among_module_connectivity(comm_graph)
			pc <- private$participation_coeffiecient(amd, td)
			zp <- data.frame(z, pc)
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
			sg1 <- decompose.graph(comm_graph, mode="strong")
			res <- data.frame()
			for(mod in unique(modvs$mod)){
				mod_nodes <- subset(modvs$taxon, modvs$mod == mod)
				neighverts <- unique(unlist(sapply(sg1, FUN = function(s){
					if(any(V(s)$name %in% mod_nodes)){
						V(s)$name
						}else{
						NULL
						}
					})))
				g3 <- induced.subgraph(graph = comm_graph, vids = neighverts)
				mod_degree <- igraph::degree(g3)
				for(i in mod_nodes){
					ki <- mod_degree[which(names(mod_degree) == i)]
					tmp <- data.frame(module = mod, taxa = names(ki), mod_links = ki)
					res <- rbind(res,tmp)
				}
			}
			res
		},
		#compute within-module degree z-score which
		#measures how well-connected a node is to other nodes in the module.
		zscore = function(mod.degree){
			ksi_bar <- aggregate(mod_links ~ module, data = mod.degree, FUN = mean)
			ksi_sigma <- aggregate(mod_links ~ module, data = mod.degree, FUN = sd)
			z <- NULL
			for(i in 1:dim(mod.degree)[1]){
				mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]
				mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]
				z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig
			}
			z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)
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
		participation_coeffiecient = function(mod.degree, total.degree){
			p <- NULL
			for(i in total.degree$taxa){
				ki <- subset(total.degree$total_links, total.degree$taxa == i)
				taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa == i)
				p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)
			}
			p <- as.data.frame(p)
			p
		},
		assign_module_roles = function(zp){
			zp <- na.omit(zp)
			if(nrow(zp) == 0){
				stop("No valid result obtained for this network!")
			}
			zp$taxa_roles <- rep(0, dim(zp)[1])
			outdf <- NULL
			for(i in 1:dim(zp)[1]){
				df <- zp[i, ]
				if(df$z <= 2.5){
					if(df$p <= 0.62){
						df$taxa_roles <- "Peripheral nodes"
					}else{
						df$taxa_roles <- "Connectors"
					}
				}
				else{
					if(df$p <= 0.62){
						df$taxa_roles <- "Module hubs"
					} else {
						df$taxa_roles <- "Network hubs"
					}
				}
				outdf <- rbind(outdf, df)
			}
			outdf
		},
		plot_roles_1 = function(node_roles, roles_color_background, roles_color_values = NULL, module = FALSE, x_lim = c(0, 1), ...){
			if(module == T){
				all_modules <- unique(as.character(node_roles$module))
				node_roles$module <- factor(as.character(node_roles$module), levels = stringr::str_sort(all_modules, numeric = TRUE))
			}
			x1 <- c(0, 0.62, 0, 0.62)
			x2 <- c(0.62, 1, 0.62, 1)
			y1 <- c(-Inf,-Inf, 2.5, 2.5)
			y2 <- c(2.5,2.5, Inf, Inf)
			lab <- c("Peripheral nodes","Connectors" ,"Module hubs","Network hubs")
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

			p <- ggplot(use_data, aes_string(x = use_level, y = "value", color = plot_color, shape = plot_shape, size = plot_size)) + 
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
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
