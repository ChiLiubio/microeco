#' @title
#' Create trans_network object for co-occurrence network analysis.
#'
#' @description
#' This class is a wrapper for a series of network analysis related methods, 
#' including the correlation based <doi:10.1186/1471-2105-13-113>, SpiecEasi <doi:10.1371/journal.pcbi.1004226>,
#' and Probabilistic Graphical Models based <doi:10.1016/j.cels.2019.08.002> network construction approaches, network and node attributes analysis
#' eigengene analysis, network subsetting and other network operations.
#'
#' @export
trans_network <- R6Class(classname = "trans_network",
	public = list(
		#' @description
		#' This function is used to create the trans_network object, store the important intermediate data 
		#'   and calculate correlations if cal_cor parameter is selected.
		#' 
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @param cor_method default "pearson"; "pearson", "spearman" or "kendall"; correlation algorithm, only use for correlation based network.
		#' @param cal_cor default "base"; "base", "WGCNA", "SparCC" or NA; correlation method; NA represent do not calculate correlations, used for non-correlation based network. 
		#' @param taxa_level default "OTU"; taxonomic rank. 
		#' @param filter_thres default 0; the relative abundance threshold. 
		#' @param nThreads default 1; the thread number used for "WGCNA" and SparCC. 
		#' @param SparCC_simu_num default 100; SparCC simulation number for bootstrap. 
		#' @param env_cols default NULL; number or name vector to select the physicochemical data in dataset$sample_table. 
		#' @param add_data default NULL; provide physicochemical table additionally.
		#' @return res_cor_p list.
		#' @examples
		#' \donttest{
		#' data(dataset)
		#' # correlation network
		#' t1 <- trans_network$new(
		#' 		dataset = dataset, 
		#' 		cal_cor = "base", 
		#' 		taxa_level = "OTU", 
		#' 		filter_thres = 0.001)
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
				env_data %<>% .[rownames(dataset1$sample_table), ]
				env_data <- dropallfactors(env_data, unfac2num = TRUE)
				env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
				env_data[] <- lapply(env_data, as.numeric)
				self$env_data <- env_data
			}
			if(taxa_level != "OTU"){
				dataset1 <- dataset1$merge_taxa(taxa = taxa_level)
			}
			# transform each object
			sampleinfo <- dataset1$sample_table
			# store taxonomic table for the following analysis
			self$use_tax <- dataset1$tax_table
			use_abund <- dataset1$otu_table
			use_abund <- use_abund[apply(use_abund, 1, sum)/sum(use_abund) > filter_thres, ]
			
			use_abund <- as.data.frame(t(use_abund))
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
					bootres <- SpiecEasi::sparccboot(use_abund, ncpus = nThreads, R = SparCC_simu_num)
					cor_result <- SpiecEasi::pval.sparccboot(bootres)
					# reshape the results
					use_names <- colnames(bootres$data)
					com_res <- t(combn(use_names, 2))
					res <- cbind.data.frame(com_res, cor_result$cors, cor_result$pvals, stringsAsFactors = FALSE)
					colnames(res) <- c("t1", "t2", "cor", "p")
					res <- rbind.data.frame(
						res, 
						data.frame(t1 = res$t2, t2 = res$t1, cor = res$cor, p = res$p), 
						data.frame(t1 = use_names, t2 = use_names, cor = rep(1, length(use_names)), p = rep(0, length(use_names)))
						)
					res$cor %<>% as.numeric
					res$p %<>% as.numeric
					res_cor <- reshape2::dcast(res, t1~t2, value.var = "cor") %>% 
						`row.names<-`(.[,1]) %>% 
						.[, -1] %>% 
						.[use_names, use_names] %>% 
						as.matrix
					res_p <- reshape2::dcast(res, t1~t2, value.var = "p") %>% 
						`row.names<-`(.[,1]) %>% 
						.[, -1] %>% 
						.[use_names, use_names] %>% 
						as.matrix
					cor_result <- list(cor = res_cor, p = res_p)
				}
				self$res_cor_p <- cor_result
				message('The correlation result list is stored in object$res_cor_p ...')
			}else{
				self$res_cor_p <- NULL
			}
			
			self$use_abund <- use_abund
			self$use_sampleinfo <- sampleinfo
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Calculate network either based on the correlation method or based on SpiecEasi or based on the Probabilistic Graphical Models (PGM) in julia FlashWeave; 
		#' See Deng et al. (2012) <doi:10.1186/1471-2105-13-113> for correlation based method, 
		#' Kurtz et al. (2015) <doi:doi:10.1371/journal.pcbi.1004226> for SpiecEasi method, 
		#' Tackmann et al. (2019) <doi:10.1016/j.cels.2019.08.002> for PGM based method.
		#'
		#' @param network_method default "COR"; "COR", "SpiecEasi" or "PGM"; COR: correlation based method; PGM: Probabilistic Graphical Models based method.
		#' @param p_thres default .01; the p value threshold.
		#' @param COR_weight default TRUE; whether use correlation coefficient as the weight of edges.
		#' @param COR_p_adjust default "fdr"; p.adjust method, see p.adjust.methods.
		#' @param COR_cut default .6; correlation coefficient threshold.
		#' @param COR_low_threshold default .4; the lowest correlation coefficient threshold, use with COR_optimization = TRUE.
		#' @param COR_optimization default FALSE; whether use random matrix theory to optimize the choice of correlation coefficient, see https://doi.org/10.1186/1471-2105-13-113
		#' @param PGM_meta_data default FALSE; whether use env data for the optimization, If TRUE, will automatically find the env_data in the object.
		#' @param PGM_sensitive default "true"; whether use sensitive type in the PGM model.
		#' @param PGM_heterogeneous default "true"; whether use heterogeneous type in the PGM model.
		#' @param SpiecEasi_method default "mb"; either 'glasso' or 'mb';see spiec.easi in package SpiecEasi and https://github.com/zdk123/SpiecEasi.
		#' @param add_taxa_name default "Phylum"; NULL or a taxonomic rank name; used to add taxonomic rank name to network.
		#' @param usename_rawtaxa_when_taxalevel_notOTU default FALSE; whether replace the name of nodes using the taxonomic information.
		#' @param ... paremeters pass to spiec.easi in package SpiecEasi for network_method = "SpiecEasi".
		#' @return res_network in object.
		#' @examples
		#' \donttest{
		#' t1$cal_network(p_thres = 0.01, COR_cut = 0.6)
		#' }
		cal_network = function(
			network_method = c("COR", "SpiecEasi", "PGM")[1],
			p_thres = 0.01,
			COR_weight = TRUE,
			COR_p_adjust = "fdr",
			COR_cut = 0.6,
			COR_low_threshold = 0.4,
			COR_optimization = FALSE,
			PGM_meta_data = FALSE,
			PGM_sensitive = "true",
			PGM_heterogeneous = "true",
			SpiecEasi_method = "mb",
			add_taxa_name = "Phylum",
			usename_rawtaxa_when_taxalevel_notOTU = FALSE,
			...
			){
			if(!require(igraph)){
				stop("igraph package not installed")
			}
			sampleinfo <- self$use_sampleinfo
			taxa_level <- self$taxa_level
			taxa_table <- self$use_tax
			if(!grepl("COR|PGM|SpiecEasi", network_method, ignore.case = TRUE)){
				stop("Unknown network method in network_method parameter!")
			}
			
			message("---------------- ", Sys.time()," : Start ----------------")
			if(grepl("COR", network_method, ignore.case = TRUE)){
				cortable <- self$res_cor_p$cor
				adp <- apply(self$res_cor_p$p, 2, p.adjust, method = COR_p_adjust)
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
				am <- as.matrix(cortable)
				am[abs(cortable) >= tc1] <- 1
				am[adp >= p_thres] <- 0
				am[am != 1] <- 0
				network <- graph.adjacency(am, mode = "undirected")
				edges <- t(sapply(1:ecount(network), function(x) ends(network, x)))
				E(network)$label <- unlist(lapply(seq_len(nrow(edges)), function(x) ifelse(cortable[edges[x, 1], edges[x, 2]] > 0, "+", "-")))
				if(COR_weight == T){
					E(network)$weight <- unlist(lapply(seq_len(nrow(edges)), function(x) abs(cortable[edges[x, 1], edges[x, 2]])))
				}else{
					E(network)$weight <- rep.int(1, ecount(network))
				}
			}
			if(grepl("SpiecEasi", network_method, ignore.case = TRUE)){
				if(!require(SpiecEasi)){
					stop("SpiecEasi package not installed")
				}
				use_abund <- self$use_abund
				use_abund <- as.matrix(use_abund)
				# calculate SpiecEasi network, reference https://github.com/zdk123/SpiecEasi
				network <- spiec.easi(use_abund, method = SpiecEasi_method, ...)
				network <- adj2igraph(getRefit(network))
				V(network)$name <- colnames(use_abund)
				E(network)$label <- unlist(lapply(E(network)$weight, function(x) ifelse(x > 0, "+", "-")))
			}
			if(grepl("PGM", network_method, ignore.case = TRUE)){
				use_abund <- self$use_abund
				# make sure working directory can not be changed by the function when quit.
				oldwd <- getwd()
				on.exit(setwd(oldwd))
				#use_abund <- cbind.data.frame(SampleID = rownames(use_abund), use_abund)
				tem_dir <- tempdir()
				setwd(tem_dir)
				write.table(use_abund, "taxa_table_PGM.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
				L1 <- "using FlashWeave\n"
				L2 <- 'data_path = "taxa_table_PGM.tsv"\n'
				if(PGM_meta_data == T){
					meta_data <- self$env_data
					write.table(meta_data, "meta_table_PGM.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
					L3 <- 'meta_data_path = "meta_table_PGM.tsv"\n'
				}else{
					L3 <- "\n"
				}
				if(PGM_meta_data == T){
					L4 <- paste0("netw_results = learn_network(data_path, meta_data_path, alpha=", p_thres, ", sensitive=", PGM_sensitive, 
						", heterogeneous=", PGM_heterogeneous, ")\n")
				}else{
					L4 <- paste0("netw_results = learn_network(data_path, alpha=", p_thres, ", sensitive=", PGM_sensitive, ", heterogeneous=", 
						PGM_heterogeneous, ")\n")
				}
				L5 <- 'save_network("network_PGM.gml", netw_results)'
				L <- paste0(L1, L2, L3, L4, L5)
				openfile <- file("calculate_network.jl", "wb")
				write(L, file = openfile)
				close(openfile)
				system("julia calculate_network.jl")
				setwd('..')
				message("The temporary files are in ", tem_dir)
				network <- read_graph(paste0(tem_dir, "/network_PGM.gml"), format = "gml")
				network <- set_vertex_attr(network, "name", value = V(network)$label)
				E(network)$label <- unlist(lapply(E(network)$weight, function(x) ifelse(x > 0, "+", "-")))
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
		#' Add network modules to the network.
		#'
		#' @param method default "cluster_fast_greedy"; the method used to find the optimal community structure of a graph;
		#' 	 the following are available functions (options) from igraph package: "cluster_fast_greedy", "cluster_optimal",
		#' 	 "cluster_edge_betweenness", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen",
		#' 	 "cluster_louvain", "cluster_spinglass", "cluster_walktrap".
		#' @param module_name_prefix default "M"; the prefix of module names; module names are made of the module_name_prefix and numbers;
		#'   numbers are assigned according to the sorting result of node numbers in modules with decreasing trend.
		#' @return a network with modules, stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_module()
		#' }
		cal_module = function(method = "cluster_fast_greedy", module_name_prefix = "M"){
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
		#' Save network as gexf style, which can be opened by Gephi <https://gephi.org/>.
		#'
		#' @param filepath default "network.gexf"; file path.
		#' @return None.
		save_network = function(filepath = "network.gexf"){
			if(!require(rgexf)){
				stop("Please install rgexf package")
			}
			private$saveAsGEXF(network = self$res_network, filepath = filepath)
		},
		#' @description
		#' Calculate network properties.
		#'
		#' @return res_network_attr in object.
		#' @examples
		#' \donttest{
		#' t1$cal_network_attr()
		#' }
		cal_network_attr = function(){
			self$res_network_attr <- private$network_attribute(self$res_network)
			message('Result is stored in object$res_network_attr ...')
		},
		#' @description
		#' Calculate node properties.
		#'
		#' @return res_node_type in object.
		#' @examples
		#' \donttest{
		#' t1$cal_node_type()
		#' }
		cal_node_type = function(){
			network <- self$res_network
			node_type <- private$module_roles(network)
			use_abund <- self$use_abund
			# create a replace_table to match the taxa name and marker name when taxa_level is not "OTU"
			if(self$taxa_level != "OTU"){
				replace_table <- data.frame(V(network)$name, V(network)$taxa, stringsAsFactors = FALSE) %>% `row.names<-`(.[,1])
				node_type <- cbind.data.frame(node_type, 
					self$use_tax[replace_table[rownames(node_type), 2], 1:which(colnames(self$use_tax) %in% self$taxa_level), drop = FALSE])
			}else{
				node_type <- cbind.data.frame(node_type, self$use_tax[rownames(node_type), ])
			}
			node_type$degree <- degree(network)[rownames(node_type)]
			node_type$betweenness <- betweenness(network)[rownames(node_type)]
			# Add abundance info
			sum_abund <- apply(use_abund, 2, function(x) sum(x) * 100/sum(use_abund))
			# Same with the above operation to make the names corresponded
			if(self$taxa_level != "OTU"){
				node_type$Abundance <- sum_abund[replace_table[rownames(node_type), 2]]
			}else{
				node_type$Abundance <- sum_abund[rownames(node_type)]
			}
			
			self$res_node_type <- node_type
			message('Result is stored in object$res_node_type ...')
		},
		#' @description
		#' Calculate eigengenes of modules, i.e. the first principal component based on PCA analysis, and the percentage of variance.
		#'
		#' @return res_eigen and res_eigen_expla in object.
		#' @examples
		#' \donttest{
		#' t1$cal_eigen()
		#' }
		cal_eigen = function(){
			use_abund <- self$use_abund
			res_node_type <- self$res_node_type
			# calculate eigengene for each module
			res_eigen <- list()
			res_eigen_expla <- c()
			for(i in unique(as.character(res_node_type$module))){
				tax_names <- rownames(res_node_type[as.character(res_node_type$module) == i, ])
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
		#' Plot the classification and importance of nodes, see object$res_node_type for the variable names used in the parameters.
		#'
		#' @param use_type default 1; 1 or 2; 1 represents taxa roles area plot; 2 represents the layered plot with taxa as x axis.
		#' @param roles_colors default NULL; for use_type 1; colors for each group.
		#' @param plot_module default FALSE; for use_type 1; whether plot the modules information.
		#' @param use_level default "Phylum"; for use_type 2; used taxonomic level in x axis.
		#' @param show_value default c("z", "p"); for use_type 2; used variable in y axis.
		#' @param show_number default 1:10; for use_type 2; showed number in x axis, sorting according to the nodes number.
		#' @param plot_color default "Phylum"; for use_type 2; used variable for color.
		#' @param plot_shape default "taxa_roles"; for use_type 2; used variable for shape.
		#' @param plot_size default "Abundance"; for use_type 2; used for point size; a fixed number (e.g. 5) is also available.
		#' @param color_values default RColorBrewer::brewer.pal(12, "Paired"); for use_type 2; color vector
		#' @param shape_values default c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14); for use_type 2; shape vector, see ggplot2 tutorial for the shape meaning.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' t1$plot_taxa_roles()
		#' }
		plot_taxa_roles = function(
			use_type = c(1, 2)[1],
			roles_colors = NULL,
			plot_module = FALSE,
			use_level = "Phylum",
			show_value = c("z", "p"),
			show_number = 1:10,
			plot_color = "Phylum",
			plot_shape = "taxa_roles",
			plot_size = "Abundance",
			color_values = RColorBrewer::brewer.pal(12, "Paired"),
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14)
			){
			if(use_type == 1){
				res <- private$plot_roles_1(node_roles = self$res_node_type, roles_colors = roles_colors, module = plot_module)
			}
			if(use_type == 2){
				res <- private$plot_roles_2(node_roles = self$res_node_type, 
					plot_color = plot_color,
					plot_shape = plot_shape,
					use_level = use_level, 
					show_value = show_value, 
					show_number = show_number,
					plot_size = plot_size,
					color_values = color_values, 
					shape_values = shape_values
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
		#' t1$subset_network(node = t1$res_node_type %>% .[.$module == "M1", ] %>% 
		#'   rownames, rm_single = TRUE)
		#' # return a sub network that contains all nodes of module M1
		#' }
		subset_network = function(node = NULL, edge = NULL, rm_single = TRUE){
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
		#' Get the table of edges, including connected nodes, labels and weight.
		#'
		#' @return data.frame
		#' @examples
		#' \donttest{
		#' t1$get_edge_table()
		#' }
		get_edge_table = function(){
			network <- self$res_network
			edges <- t(sapply(1:ecount(network), function(x) ends(network, x)))
			edge_label <- E(network)$label
			if(!is.null(E(network)$weight)){
				edge_weight <- E(network)$weight
			}else{
				edge_weight <- rep(NA, times = length(edge_label))
			}
			res_edge_table <- data.frame(edges, edge_label, edge_weight)
			colnames(res_edge_table) <- c("node1", "node2", "label", "weight")
			res_edge_table
		},
		#' @description
		#' Perform a bootstrapping hypothesis test to determine whether degrees follows a power law distribution;
		#' a significant p represents the distribution does not follow power law.
		#'
		#' @param ... paremeters pass to bootstrap_p in poweRlaw package.
		#' @return two lists stored in object; see estimate_xmin and bootstrap_p in poweRlaw package for the details.
		#' @examples
		#' \donttest{
		#' t1$cal_powerlaw_p()
		#' }
		cal_powerlaw_p = function(...){
			network <- self$res_network
			degree_dis <- degree(network)
			if(!require(poweRlaw)){
				stop("Please first install poweRlaw package from CRAN !")
			}
			resdispl <- poweRlaw::displ$new(degree_dis + 1)
			est_xmin <- poweRlaw::estimate_xmin(resdispl)
			resdispl$setXmin(est_xmin)
			res <- poweRlaw::bootstrap_p(resdispl, ...)
			self$res_powerlaw_min <- est_xmin
			message('Estimated lower bound result is stored in object$res_powerlaw_min ...')
			self$res_powerlaw_p <- res
			message('Bootstrap p value is stored in object$res_powerlaw_p ...')			
		},
		#' @description
		#' Fit degrees to a power law distribution.
		#'
		#' @param xmin default NULL; See xmin in fit_power_law function; suggest using the result res_powerlaw_min from cal_powerlaw_p function.
		#' @param ... paremeters pass to fit_power_law function in igraph package.
		#' @return list stored in object; see fit_power_law function for the details explanation.
		#' @examples
		#' \donttest{
		#' t1$cal_powerlaw_fit()
		#' }
		cal_powerlaw_fit = function(xmin = NULL, ...){
			network <- self$res_network
			degree_dis <- degree(network, mode="in")
			self$res_powerlaw_fit <- fit_power_law(degree_dis + 1, xmin = xmin, ...)
			message('Powerlaw fitting result is stored in object$res_powerlaw_fit ...')
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
			if(!is.null(self$res_node_type)){
				cat("res_node_type object: finished\n")
				cat(paste0("	colnames: ", paste0(colnames(self$res_node_type), collapse = ", "), "\n"))
			}else{
				cat("res_node_type object: NULL\n")
			}
			invisible(self)
		}
		),
	private = list(
		cal_corr = function(inputtable, cor_method) {
			N <- ncol(inputtable)
			use_names <- colnames(inputtable)
			com_res <- combn(use_names, 2)
			res <- lapply(seq_len(ncol(com_res)), function(x){
				cor_res <- suppressWarnings(cor.test(inputtable[, com_res[1, x]], inputtable[, com_res[2, x]], method = cor_method));
				c(t1 = com_res[1, x], t2 = com_res[2, x], cor = unname(cor_res$estimate), p = unname(cor_res$p.value))
			})
			res <- do.call(rbind, res) %>% 
				as.data.frame(stringsAsFactors = FALSE)
			res <- rbind.data.frame(
				res, 
				data.frame(t1 = res$t2, t2 = res$t1, cor = res$cor, p = res$p), 
				data.frame(t1 = use_names, t2 = use_names, cor = rep(1, N), p = rep(0, N))
				)
			res$cor %<>% as.numeric
			res$p %<>% as.numeric
			res_cor <- reshape2::dcast(res, t1~t2, value.var = "cor") %>% `row.names<-`(.[,1]) %>% .[, -1] %>% .[use_names, use_names] %>% as.matrix
			res_p <- reshape2::dcast(res, t1~t2, value.var = "p") %>% `row.names<-`(.[,1]) %>% .[, -1] %>% .[use_names, use_names] %>% as.matrix
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
				print(i*100)
			}
			ps <- cbind(ps, c(seq(lcor, hcor, 0.01)))
			tc <- ps[ps[,1] == min(ps[,1]), 2]
			tc
		},
		nnsd = function(sp){
			nns <- NULL
			for(j in 2:length(sp)){
				nn <- abs(sp[j] - sp[j-1])
				nns <- c(nns, nn)
			}
			nns
		},
		saveAsGEXF = function(network, filepath = "network.gexf"){
			require(rgexf)
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
				Average_degree = sum(degree(x))/length(degree(x)), 
				Average_path_length = average.path.length(x), 
				Network_diameter = round(diameter(x, directed = FALSE), 0), 
				Clustering_coefficient = transitivity(x), 
				Density = graph.density(x), 
				Heterogeneity = sd(degree(x))/mean(degree(x)), 
				Centralization = centr_degree(x)$centralization
				)
			res <- base::as.data.frame(t(res))
			colnames(res) <- NULL
			res
		},
		# modified based on microbiomeSeq (http://www.github.com/umerijaz/microbiomeSeq) 
		module_roles = function(comm_graph){
			require(igraph)
			td <- degree(comm_graph) %>% data.frame(taxa = names(.), total_links = ., stringsAsFactors = FALSE)
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
		#compute within-module degree for each of the features
		within_module_degree = function(comm_graph){
			mods <- get.vertex.attribute(comm_graph, "module")
			if(is.null(mods)){
				stop("No modules found! Please first calculate network modules using function cal_module!")
			}
			modvs <- data.frame("taxon"= V(comm_graph)$name, "mod"=mods, stringsAsFactors = FALSE)
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
				mod_degree <- degree(g3)
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
			ksi_bar <- aggregate(mod_links ~ module, data=mod.degree, FUN = mean)
			ksi_sigma <- aggregate(mod_links ~ module, data=mod.degree, FUN = sd)
			z <- NULL
			for(i in 1:dim(mod.degree)[1]){
				mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]
				mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]
				z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig
			}
			z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)
			z
		},
		#calculate the degree (links) of each node to nodes in other modules.
		among_module_connectivity = function(comm_graph){
			modvs <- data.frame(taxon= V(comm_graph)$name, mod=get.vertex.attribute(comm_graph, "module"), stringsAsFactors = FALSE)
			edges <- t(sapply(1:ecount(comm_graph), function(x) ends(comm_graph, x)))
			res <- lapply(modvs$taxon, function(x){
						sapply(unique(c(edges[edges[,1] == x, 2], edges[edges[,2] == x, 1])), function(y){
							c(taxa = x, module = modvs[modvs$taxon==y, "mod"])
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
		plot_roles_1 = function(node_roles, roles_colors = NULL, module = FALSE){
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
			if(is.null(roles_colors)){roles_colors <- rev(c("grey80", RColorBrewer::brewer.pal(3, "Dark2")))}

			p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=lab, alpha = .4))
			p <- p + guides(fill=guide_legend(title="Roles"), alpha = FALSE)
			p <- p + scale_fill_manual(values = roles_colors)
			if(module == T){
				p <- p + geom_point(data=node_roles, aes(x=p, y=z, shape= module)) + theme_bw() + guides(shape=guide_legend(title="Module"))
			}else{
				p <- p + geom_point(data=node_roles, aes(x=p, y=z)) + theme_bw()
			}
			p <- p + theme(strip.background = element_rect(fill = "white")) + 
				xlab("Among-module connectivity") + 
				ylab("Within-module connectivity")
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
			shape_values = c(16, 17, 7, 8, 15, 18, 11, 10, 12, 13, 9, 3, 4, 0, 1, 2, 14)
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
				geom_point(position = "jitter") +
				scale_color_manual(values = color_values) +
				scale_shape_manual(values = shape_values) +				
				facet_grid(variable ~ ., drop = TRUE, scale = "free", space = "fixed") +
				theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 10))
			
			if(plot_color == use_level){
				p <- p + guides(color = FALSE)
			}
			if(plot_shape == use_level){
				p <- p + guides(shape = FALSE)
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
