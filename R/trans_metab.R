#' @title
#' Create \code{trans_metab} object for metabolite analysis.
#'
#' @description
#' This class is a wrapper for a series of metabolite analysis, including origin inference and pathway enrichment.
#'
#' @export
trans_metab <- R6Class(classname = "trans_metab",
	public = list(
		#' @description
		#' Create the \code{trans_metab} object.
		#' 
		#' @param metab default NULL; metabolite data. A \code{\link{microtable}} object or data.frame object.
		#'    If the input is a data.frame object, the function can judge whether it is abundance table and preprocess the data.		
		#' @param microb default NULL; A \code{\link{microtable}} object.
		#' @return \code{data_metab} and \code{data_microb} stored in the object.
		#' @examples
		#' data(soil_metab)
		#' data(soil_microb)
		#' t1 <- trans_metab$new(metab = soil_metab, microb = soil_microb)
		initialize = function(
			metab = NULL,
			microb = NULL
			){
			if(! is.null(metab)){
				if(inherits(metab, "data.frame")){
					if(all(sapply(metab, function(x){inherits(x, "numeric")}))){
						metab_otu_table <- metab
						metab_tax_table <- data.frame(metab_name = rownames(metab))
						rownames(metab_tax_table) <- rownames(metab)
						metab <- microtable$new(otu_table = metab_otu_table, tax_table = metab_tax_table)
					}else{
						metab_tax_table <- metab
						metab <- list()
						metab$tax_table <- metab_tax_table
					}
				}else{
					if(inherits(metab, "microtable")){
						if(is.null(metab$tax_table)){
							metab_tax_table <- data.frame(metab_name = rownames(metab$otu_table))
							rownames(metab_tax_table) <- rownames(metab$otu_table)
							metab$tax_table <- metab_tax_table
						}
					}else{
						stop("Unknown data type of input metab !")
					}
				}
				if(nrow(metab$tax_table) > 3){
					show_metab <- c(rownames(metab$tax_table)[1:3], "...")
				}else{
					show_metab <- rownames(metab$tax_table)
				}
				message("Input metabolites have ", nrow(metab$tax_table), ", including ", paste0(show_metab, collapse = ", "))
			}
			self$data_metab <- metab
			self$data_microb <- microb
		},
		#' @description
		#' Matching compound names against the database from TidyMass2 (DOI: 10.1038/s41467-026-68464-7) to facilitate the acquisition of standardized nomenclature.
		#' An alternative way is to use MetaboAnalyst website (https://www.metaboanalyst.ca/faces/upload/ConvertView.xhtml).
		#' 
		#' @param database_path default "./metorigindb_split_202602"; directory path of the downloaded database. 
		#'	  Please download the pre-collated metorigindb database (RData format) from zenodo (https://zenodo.org/records/18618912) and extract the compressed archive.
		#' @param method default "jw"; method of approximate string matching. Default "jw" is Jaro-Winkler distance.
		#'	  See the \code{method} parameter of \code{amatch} function in stringdist package.
		#' @param maxDist default 0.3; See the \code{maxDist} parameter of \code{amatch} function in stringdist package.
		#' @param ... parameters passed to \code{amatch} function of stringdist package.
		#' @return \code{res_match} table stored in the object.
		#' @examples
		#' \dontrun{
		#' t1$cal_match()
		#' }
		cal_match = function(
			database_path = "./metorigindb_split_202602",
			method = "jw",
			maxDist = 0.3,
			...
			){
			data_metab <- self$data_metab
			
			private$check_input_path(database_path)
			metorigindb_number_path <- file.path(database_path, "metorigindb_number", "metorigindb_number.RData")
			private$check_file_exist(metorigindb_number_path, "metorigindb_number.RData", file.path(database_path, "metorigindb_number"))
			load(metorigindb_number_path, envir = environment())

			user_names <- rownames(data_metab$tax_table)
			standard_names <- metorigindb_number$Compound_name

			user_clean <- private$clean_name(user_names)
			std_clean <- private$clean_name(standard_names)

			# Jaro-Winkler distance, sensitive to the prefix
			matched_index <- stringdist::amatch(
				user_clean,
				table = std_clean,
				method = method,
				maxDist = maxDist,
				...
			)

			res <- data.frame(
				original = user_names,
				cleaned = user_clean,
				matched_standard = ifelse(
					is.na(matched_index), 
					NA, 
					standard_names[matched_index]
				),
				distance = ifelse(
					is.na(matched_index), 
					NA, 
					stringdist::stringdist(user_clean, std_clean[matched_index], method = method)
				)
			)

			self$res_match <- res
			message("Match table is stored in object$res_match ...")
			invisible(self)
		},
		#' @description
		#' Metabolite origin inference based on the preprocessed database from TidyMass2 (DOI: 10.1038/s41467-026-68464-7)
		#' 
		#' @param database_path default "./metorigindb_split_202602"; directory path of the downloaded database. 
		#'	  Please download the pre-collated metorigindb database (RData format) from zenodo (https://zenodo.org/records/18618912) and extract the compressed archive.
		#' @param match_col default "names"; How to match to the data of metorigindb. Default "names" means using the input names of metabolites.
		#'    If the tax_table has other columns like "HMDB_ID" or "KEGG_ID", the user can provide more items, like c("names", "HMDB_ID").
		#'    The program will automatically identify the corresponding columns from the input data based on keywords such as HMDB or KEGG, and then match them with the database.
		#' @param match_names_distance default 0; distance threshold used if the \code{res_match_table} is found in the object, 
		#'    which is calculated from the \code{cal_match} function. Available for the "names" option in \code{match_col} parameter.
		#' @param bac_level default "Genus"; which bacteria level is used to parse the taxa in the \code{data_microb} of object.
		#'    The function can automatically match the taxa those found in the input data.
		#' @return \code{res_origin_rawtable} table and \code{res_origin_list} list stored in the object.
		#'    \code{res_origin_rawtable} is the origin table extracted from the metorigindb database based on the match of name or other ID.
		#'    In \code{res_origin_list}, name is the metabolite; each element is the taxa that may produce the metabolite.
		#' @examples
		#' \dontrun{
		#' t1$cal_origin()
		#' t1$cal_origin(match_col = c("names", "HMDB_ID", "KEGG_ID"))
		#' }
		cal_origin = function(
			database_path = "./metorigindb_split_202602", 
			match_col = "names",
			match_names_distance = 0,
			bac_level = "Genus"
			){
			data_metab <- self$data_metab
			data_microb <- self$data_microb
			
			private$check_input_path(database_path)
			metorigindb_number_path <- file.path(database_path, "metorigindb_number", "metorigindb_number.RData")
			private$check_file_exist(metorigindb_number_path, "metorigindb_number.RData", file.path(database_path, "metorigindb_number"))
			load(metorigindb_number_path, envir = environment())
			
			`%in2%` <- function(x, vec) {!is.na(x) & (x %in% vec)}
			
			length_db <- length(metorigindb_number$Compound_name)
			if("names" %in% match_col){
				res_match_table <- self$res_match
				if(! is.null(res_match_table)){
					message("res_match_table is found in the object. Use the names of matched_standard column to search the database ...")
					select_match_table <- res_match_table[res_match_table$distance == match_names_distance, ]
					name_ID_select <- metorigindb_number$Compound_name %in2% select_match_table$matched_standard
				}else{
					name_ID_select <- metorigindb_number$Compound_name %in2% rownames(data_metab$tax_table)
				}
			}else{
				name_ID_select <- logical(length_db)
			}
			if(any(grepl("HMDB", match_col, ignore.case = TRUE))){
				which_col <- match_col[which(grepl("HMDB", match_col, ignore.case = TRUE))]
				HMDB_ID_select <- metorigindb_number$HMDB_ID_all %>% lapply(., function(x){strsplit(x, "{}", fixed = TRUE)[[1]]}) %>% 
					lapply(., function(x){!all(is.na(x)) & any(x %in% data_metab$tax_table[, which_col])}) %>% 
					unlist
			}else{
				HMDB_ID_select <- logical(length_db)
			}
			if(any(grepl("KEGG", match_col, ignore.case = TRUE))){
				which_col <- match_col[which(grepl("KEGG", match_col, ignore.case = TRUE))]
				KEGG_ID_select <- metorigindb_number$KEGG_ID_all %>% lapply(., function(x){strsplit(x, "{}", fixed = TRUE)[[1]]}) %>% 
					lapply(., function(x){!all(is.na(x)) & any(x %in% data_metab$tax_table[, which_col])}) %>% 
					unlist
			}else{
				KEGG_ID_select <- logical(length_db)
			}

			extract_number_table <- metorigindb_number[name_ID_select | HMDB_ID_select | KEGG_ID_select, ]

			source_path <- file.path(database_path, "metorigindb_source")
			private$check_directory_exist(source_path, "metorigindb_source", database_path)

			source_path_files <- list.files(source_path)
			res <- list()
			ID_pass <- c()

			for(i in extract_number_table$ID){
				if(i %in% ID_pass){
					next
				}
				ID_prefix <- gsub("^N", "", i) %>% as.numeric %>% {./5000} %>% floor %>% {paste0("N", {.*5000 + 1}, "-")}
				ID_path <- source_path_files[grepl(ID_prefix, source_path_files)] %>% {paste0(source_path, "/", .)}
				load(ID_path)
				ID_pass <- c(ID_pass, extract_number_table$ID %>% .[. %in% each_source$ID])
				res[[ID_prefix]] <- each_source[each_source$ID %in% extract_number_table$ID, ]
			}

			source_table <- do.call(rbind, res) %>% as.data.frame

			self$res_origin_rawtable <- source_table
			message("Raw origin table is stored in object$res_origin_rawtable ...")
			
			bac_level_lower <- tolower(bac_level)
			bac_level_prefix <- bac_level_lower %>% substr(., 1, 1)
			
			# split the taxa in source table
			res_origin <- list()
			for(j in seq_len(nrow(source_table))){
				compound_name <- source_table[j, "Compound_name"]
				res_origin[[compound_name]] <- source_table[j, paste0("bacteria_", bac_level_lower)]
			}
			res_origin %<>% lapply(., function(x){strsplit(x, "{}", fixed = TRUE)[[1]]}) %>% 
				lapply(., function(x){x[x != "Unknown"]}) %>%
				lapply(., function(x){gsub("\\(.*\\)", "", x)}) %>%
				lapply(., function(x){trimws(x)})
			
			if(! is.null(data_microb)){
				if(! bac_level %in% colnames(data_microb$tax_table)){
					stop("bac_level: ", bac_level, " is not found in the column names of tax_table in data_microb !")
				}
				all_taxa <- data_microb$tax_table[, bac_level] %>% unique %>% .[. != paste0(bac_level_prefix, "__")] %>% gsub(".__", "", .)
				res_origin %<>% lapply(., function(x){x[x %in% all_taxa]})
			}
			
			self$res_origin_list <- res_origin
			if(! is.null(data_microb)){
				message("Filtered origin taxa at ", bac_level, " level for each metabolite is stored in object$res_origin_list ...")
			}else{
				message("Origin taxa at ", bac_level, " level for each metabolite is stored in object$res_origin_list ... ...")
			}
			invisible(self)
		},
		#' @description
		#' Metabolite-bacteria network based on the \code{res_origin_list} data from the \code{cal_origin} function
		#' 
		#' @return \code{igraph} format network.
		#' @examples
		#' \dontrun{
		#' t1$cal_origin_network()
		#' }
		cal_origin_network = function(){
			res_origin_list <- self$res_origin_list
			if(is.null(res_origin_list)){
					stop("Please first run the cal_origin function !")
			}
			if(!requireNamespace("igraph", quietly = TRUE)){
					stop("Package 'igraph' is required but not installed !")
			}

			# from -> to
			edges_list <- lapply(names(res_origin_list), function(target) {
					sources <- res_origin_list[[target]]
					if (length(sources) == 0) return(NULL)
					data.frame(from = sources, to = target, stringsAsFactors = FALSE)
			})
			# Filter NULL entries before rbind to avoid errors
			edges_list <- edges_list[!sapply(edges_list, is.null)]
			if(length(edges_list) == 0){
					stop("No edges found in res_origin_list. All metabolites have empty origin taxa !")
			}
			edges <- do.call(rbind, edges_list)

			# directed network
			net <- igraph::graph_from_edgelist(as.matrix(edges[, c("from", "to")]), directed = TRUE)

			# add attributes - use the metabolite names as reference for type assignment
			metabolite_names <- names(res_origin_list)
			all_nodes <- igraph::V(net)$name
			igraph::V(net)$type <- ifelse(all_nodes %in% metabolite_names, "metabolite", "bacteria")

			net
		},
		#' @description
		#' Map metabolites to pathways based on a database (KEGG, MetaCyc, Reactome, or custom).
		#' 
		#' @param db_type default "kegg"; pathway database type, one of "kegg", "metacyc", "reactome", "custom".
		#' @param db_path default NULL; pathway database file path (RData format). 
		#'    For KEGG: if NULL, will attempt to fetch from KEGG REST API.
		#'    For MetaCyc/Reactome/Custom: must provide local database file path.
		#' @param id_type default NULL; column name in tax_table containing metabolite IDs.
		#'    If NULL, will be automatically set based on db_type: "KEGG_ID" for kegg, "METACYC_ID" for metacyc, 
		#'    "REACTOME_ID" for reactome, "CUSTOM_ID" for custom.
		#' @param save_db default FALSE; whether save the fetched database (KEGG) to the local folder.
		#' @return \code{res_pathway_map} (data.frame) stored in the object. 
		#'    Contains columns: pathway_id, pathway_name, metab_id, metab_name.
		#' @examples
		#' \dontrun{
		#' # KEGG pathway mapping (online)
		#' t1$cal_pathway()
		#' 
		#' # KEGG pathway mapping (local database)
		#' t1$cal_pathway(db_type = "kegg", db_path = "kegg_pathway_db.RData")
		#' }
		cal_pathway = function(
			db_type = c("kegg", "metacyc", "reactome", "custom")[1],
			db_path = NULL,
			id_type = NULL,
			save_db = FALSE
			){
			data_metab <- self$data_metab
			if(is.null(data_metab)){
					stop("data_metab is not found in the object !")
			}
			if(is.null(data_metab$tax_table)){
					stop("tax_table is not found in data_metab !")
			}

			db_type <- match.arg(db_type, c("kegg", "metacyc", "reactome", "custom"))

			if(is.null(id_type)){
					id_type <- switch(db_type,
							kegg = "KEGG_ID",
							metacyc = "METACYC_ID",
							reactome = "REACTOME_ID",
							custom = "CUSTOM_ID"
					)
			}
			if(! id_type %in% colnames(data_metab$tax_table)){
					stop(id_type, " column is not found in data_metab$tax_table !")
			}

			metab_ids <- data_metab$tax_table[, id_type]
			metab_names <- rownames(data_metab$tax_table)

			if(! is.null(db_path)){
					if(! file.exists(db_path)){
							stop("Database file not found: ", db_path)
					}
					message("Loading pathway database from local file ...")
					db_data <- private$load_pathway_db(db_path, db_type)
			}else{
					if(db_type == "kegg"){
							message("Fetching KEGG pathway information from KEGG REST API ...")
							db_data <- private$fetch_kegg_pathways()
							if(save_db){
								save(db_data, file = "KEGG_pathway_db.RData")
								message("Save KEGG pathway data to the folder: KEGG_pathway_db.RData ...")
							}
					}else{
							stop(db_type, " database requires local file. Please provide db_path parameter !")
					}
			}

			pathways <- db_data$pathways
			metab_to_pathway <- db_data$metab_to_pathway

			metab_valid_ids <- metab_ids[! is.na(metab_ids) & metab_ids != ""]
			metab_valid_names <- metab_names[! is.na(metab_ids) & metab_ids != ""]

			# Use vectorized approach for efficiency instead of row-by-row data.frame creation
			res_pathway_id <- character()
			res_pathway_name <- character()
			res_metab_id <- character()
			res_metab_name <- character()

			for(i in seq_along(metab_valid_ids)){
					metab_id <- metab_valid_ids[i]
					metab_name <- metab_valid_names[i]
					
					pathway_ids <- metab_to_pathway[[metab_id]]
					if(! is.null(pathway_ids) && length(pathway_ids) > 0){
							# Only keep pathway_ids that exist in the pathways dictionary
							valid_pids <- pathway_ids[pathway_ids %in% names(pathways)]
							if(length(valid_pids) > 0){
									res_pathway_id <- c(res_pathway_id, valid_pids)
									res_pathway_name <- c(res_pathway_name, unlist(pathways[valid_pids]))
									res_metab_id <- c(res_metab_id, rep(metab_id, length(valid_pids)))
									res_metab_name <- c(res_metab_name, rep(metab_name, length(valid_pids)))
							}
					}
			}

			if(length(res_pathway_id) == 0){
					stop("No pathway mapping found for the given metabolites. Please check ", id_type, " column !")
			}

			res_pathway_map <- data.frame(
					pathway_id = res_pathway_id,
					pathway_name = res_pathway_name,
					metab_id = res_metab_id,
					metab_name = res_metab_name,
					stringsAsFactors = FALSE
			)

			self$res_pathway_map <- res_pathway_map
			self$param$pathway <- list(db_type = db_type, id_type = id_type, db_path = db_path)
			message("Pathway mapping result is stored in object$res_pathway_map ...")
			message("Total ", nrow(res_pathway_map), " metabolite-pathway mappings found for ", 
					length(unique(res_metab_name)), " metabolites across ", 
					length(unique(res_pathway_id)), " pathways ...")
			invisible(self)
		},
		#' @description
		#' Perform pathway enrichment analysis for target metabolites using Fisher's exact test.
		#' 
		#' @param target_metabs character vector; target metabolite names (row names in tax_table), 
		#'    e.g., differential metabolites. If NULL, will use all metabolites in the pathway mapping.
		#' @param background_metabs character vector; background metabolite names (row names in tax_table). 
		#'    If NULL, will use all metabolites in the pathway mapping result.
		#' @param p_adjust_method default "BH"; p-value adjustment method. See \code{p.adjust} function.
		#' @param p_cutoff default 0.05; significance threshold for adjusted p-value.
		#' @param min_metab_count default 3; minimum number of background metabolites in a pathway for enrichment test.
		#' @return \code{res_pathway_enrich} data.frame stored in the object.
		#'    Contains: pathway_id, pathway_name, metab_count, background_count, 
		#'    enrichment_factor, p_value, p_adjust, significant.
		#' @examples
		#' \dontrun{
		#' t1$cal_pathway()
		#' target <- rownames(t1$data_metab$otu_table)[1:10]
		#' t1$cal_pathway_enrich(target_metabs = target)
		#' }
		cal_pathway_enrich = function(
			target_metabs = NULL,
			background_metabs = NULL,
			p_adjust_method = "BH",
			p_cutoff = 0.05,
			min_metab_count = 3
			){
			res_pathway_map <- self$res_pathway_map
			if(is.null(res_pathway_map)){
					stop("Please first run the cal_pathway function !")
			}

			all_metabs <- unique(res_pathway_map$metab_name)

			if(is.null(background_metabs)){
				background_metabs <- all_metabs
			}else{
				if(! all(background_metabs %in% all_metabs)){
					missing_metabs <- setdiff(background_metabs, all_metabs)
					message("Note: ", length(missing_metabs), " background metabolites not found in pathway mapping, will be ignored ...")
					background_metabs <- intersect(background_metabs, all_metabs)
				}
			}

			if(is.null(target_metabs)){
				target_metabs <- background_metabs
			}else{
				if(! all(target_metabs %in% background_metabs)){
					missing_targets <- setdiff(target_metabs, background_metabs)
					if(length(missing_targets) > 0){
							message("Note: ", length(missing_targets), " target metabolites are not in background, will use intersection ...")
					}
					target_metabs <- intersect(target_metabs, background_metabs)
				}
			}

			if(length(target_metabs) == 0){
				stop("No valid target metabolites !")
			}
			if(length(background_metabs) == 0){
				stop("No valid background metabolites !")
			}

			# Precompute pathway-metabolite lookup for efficiency
			N <- length(background_metabs)
			n <- length(target_metabs)

			# Build target set for fast lookup
			target_set <- unique(target_metabs)
			background_set <- unique(background_metabs)

			pathways <- unique(res_pathway_map$pathway_id)
			res_enrich <- list()

			for(pathway in pathways){
				pathway_metabs <- unique(res_pathway_map[res_pathway_map$pathway_id == pathway, "metab_name"])
				
				metab_in_pathway_bg <- intersect(pathway_metabs, background_set)
				M <- length(metab_in_pathway_bg)
				
				if(M < min_metab_count){
						next
				}

				metab_in_pathway_target <- intersect(pathway_metabs, target_set)
				k <- length(metab_in_pathway_target)

				if(k == 0){
						next
				}

				# Validate contingency table has no negative values
				cell_d <- N - M - n + k
				if(cell_d < 0){
					message("Note: Contingency table has negative cell for pathway ", pathway, 
							". This may indicate inconsistent target/background sets. Skipping ...")
					next
				}

				contingency <- matrix(c(k, M - k, n - k, cell_d), nrow = 2)
				p_value <- fisher.test(contingency, alternative = "greater")$p.value

				# Calculate enrichment factor with safety check
				if(n > 0 && M > 0){
					enrichment_factor <- (k / n) / (M / N)
				}else{
					enrichment_factor <- NA
				}

				pathway_name <- res_pathway_map[res_pathway_map$pathway_id == pathway, "pathway_name"][1]

				res_enrich[[length(res_enrich) + 1]] <- data.frame(
					pathway_id = pathway,
					pathway_name = pathway_name,
					metab_count = k,
					background_count = M,
					enrichment_factor = round(enrichment_factor, 3),
					p_value = p_value,
					stringsAsFactors = FALSE
				)
			}

			if(length(res_enrich) == 0){
				message("No enriched pathways found with minimum metabolite count = ", min_metab_count, " ...")
				self$res_pathway_enrich <- data.frame()
				self$param$pathway_enrich <- list(
						p_adjust_method = p_adjust_method,
						p_cutoff = p_cutoff,
						min_metab_count = min_metab_count
				)
				return(invisible(self))
			}

			res_enrich_df <- do.call(rbind, res_enrich)
			res_enrich_df$p_adjust <- p.adjust(res_enrich_df$p_value, method = p_adjust_method)
			res_enrich_df$significant <- res_enrich_df$p_adjust <= p_cutoff

			res_enrich_df <- res_enrich_df[order(res_enrich_df$p_adjust), ]
			rownames(res_enrich_df) <- NULL

			self$res_pathway_enrich <- res_enrich_df
			self$param$pathway_enrich <- list(
				p_adjust_method = p_adjust_method,
				p_cutoff = p_cutoff,
				min_metab_count = min_metab_count,
				n_target = n,
				N_background = N
			)

			sig_count <- sum(res_enrich_df$significant)
			message("Pathway enrichment analysis completed ...")
			message("Found ", sig_count, " significant enriched pathways (p_adjust <= ", p_cutoff, ") out of ", 
				nrow(res_enrich_df), " tested pathways ...")
			message("Result is stored in object$res_pathway_enrich ...")
			invisible(self)
		},
		#' @description
		#' Visualize pathway enrichment analysis results.
		#' 
		#' @param plot_type default "bubble"; plot type, "bubble" or "bar".
		#' @param top_n default 20; number of top pathways to display.
		#' @param color_by default "p_adjust"; column name for color mapping.
		#' @param size_by default "metab_count"; column name for size mapping (only for bubble plot).
		#' @param order_by default "p_adjust"; column name for ordering pathways.
		#' @param color_gradient_low default "#132B43"; color for low values (e.g., high p-value = not significant).
		#' @param color_gradient_high default "#56B1F7"; color for high values (e.g., low p-value = significant).
		#' @return ggplot2 object.
		#' @examples
		#' \dontrun{
		#' t1$cal_pathway()
		#' target <- rownames(t1$data_metab$otu_table)[1:10]
		#' t1$cal_pathway_enrich(target_metabs = target)
		#' t1$plot_pathway_enrich()
		#' }
		plot_pathway_enrich = function(
			plot_type = c("bubble", "bar")[1],
			top_n = 20,
			color_by = "p_adjust",
			size_by = "metab_count",
			order_by = "p_adjust",
			color_gradient_low = "#132B43",
			color_gradient_high = "#56B1F7"
			){
			res_enrich <- self$res_pathway_enrich
			if(is.null(res_enrich) || nrow(res_enrich) == 0){
				stop("No enrichment results found. Please first run the cal_pathway_enrich function !")
			}

			plot_type <- match.arg(plot_type, c("bubble", "bar"))

			if(! color_by %in% colnames(res_enrich)){
				stop(color_by, " is not found in res_pathway_enrich !")
			}
			if(! order_by %in% colnames(res_enrich)){
				stop(order_by, " is not found in res_pathway_enrich !")
			}

			plot_data <- res_enrich[order(res_enrich[[order_by]]), ]
			if(nrow(plot_data) > top_n){
				plot_data <- plot_data[1:top_n, ]
			}

			plot_data$pathway_label <- paste0(plot_data$pathway_name, " (", plot_data$pathway_id, ")")
			# Order: smallest p_adjust at top of plot
			plot_data$pathway_label <- factor(plot_data$pathway_label, levels = rev(plot_data$pathway_label))

			if(plot_type == "bubble"){
				if(! size_by %in% colnames(plot_data)){
						stop(size_by, " is not found in res_pathway_enrich !")
				}

				g <- ggplot(plot_data, aes(x = pathway_label, y = -log10(.data[[color_by]]), 
								size = .data[[size_by]], color = .data[[color_by]])) +
						geom_point() +
						scale_color_gradient(low = color_gradient_low, high = color_gradient_high) +
						scale_size_continuous(range = c(2, 8)) +
						coord_flip() +
						labs(x = NULL, y = paste0("-log10(", color_by, ")"), size = size_by, color = color_by) +
						theme_bw() +
						theme(axis.text.x = element_text(size = 12, colour = "black"),
								  axis.text.y = element_text(size = 10),
								  panel.grid = element_blank(),
								  panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
								  legend.position = "right")
			}else{
				g <- ggplot(plot_data, aes(x = pathway_label, y = -log10(.data[[color_by]]), 
								fill = .data[[color_by]])) +
						geom_bar(stat = "identity") +
						scale_fill_gradient(low = color_gradient_low, high = color_gradient_high) +
						coord_flip() +
						labs(x = NULL, y = paste0("-log10(", color_by, ")"), fill = color_by) +
						theme_bw() +
						theme(axis.text.x = element_text(size = 12, colour = "black"),
								  axis.text.y = element_text(size = 10),
								  panel.grid = element_blank(),
								  panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
								  legend.position = "right")
			}

			g
		},
		#' @description
		#' Build metabolite-pathway network based on pathway mapping and enrichment results.
		#' 
		#' @param significant_only default TRUE; whether to use only significant enriched pathways.
		#' @param p_cutoff default 0.05; p-value cutoff for significant pathways (used only when significant_only = TRUE).
		#' @param use_enrichment default TRUE; whether to require enrichment results. If FALSE, build network from all pathway mappings.
		#' @return \code{igraph} object representing metabolite-pathway bipartite network.
		#' @examples
		#' \dontrun{
		#' t1$cal_pathway()
		#' target <- rownames(t1$data_metab$otu_table)[1:50]
		#' t1$cal_pathway_enrich(target_metabs = target)
		#' net <- t1$cal_pathway_network()
		#' }
		cal_pathway_network = function(
			significant_only = TRUE,
			p_cutoff = 0.05,
			use_enrichment = TRUE
			){
			res_pathway_map <- self$res_pathway_map

			if(is.null(res_pathway_map)){
				stop("Please first run the cal_pathway function !")
			}

			if(use_enrichment || significant_only){
				res_enrich <- self$res_pathway_enrich
				if(is.null(res_enrich) || nrow(res_enrich) == 0){
					stop("No enrichment results found. Please first run the cal_pathway_enrich function, or set use_enrichment = FALSE !")
				}
			}

			if(significant_only && use_enrichment){
				if("significant" %in% colnames(res_enrich)){
					sig_pathways <- res_enrich[res_enrich$significant == TRUE, "pathway_id"]
				}else{
					sig_pathways <- res_enrich[res_enrich$p_adjust <= p_cutoff, "pathway_id"]
				}
				if(length(sig_pathways) == 0){
					stop("No significant pathways found with p_cutoff = ", p_cutoff, " !")
				}
				res_pathway_map <- res_pathway_map[res_pathway_map$pathway_id %in% sig_pathways, ]
			}

			if(nrow(res_pathway_map) == 0){
				stop("No metabolite-pathway mappings found !")
			}

			if(!requireNamespace("igraph", quietly = TRUE)){
				stop("Package 'igraph' is required but not installed !")
			}

			edges <- data.frame(
					from = res_pathway_map$metab_name,
					to = paste0(res_pathway_map$pathway_id, ": ", res_pathway_map$pathway_name),
					stringsAsFactors = FALSE
			)

			binet <- igraph::graph_from_data_frame(edges, directed = FALSE)

			metabolites <- unique(res_pathway_map$metab_name)
			igraph::V(binet)$type <- ifelse(igraph::V(binet)$name %in% metabolites, "metabolite", "pathway")

			message("Network constructed: ", igraph::ecount(binet), " edges, ", igraph::vcount(binet), " nodes ...")
			message("Metabolites: ", sum(igraph::V(binet)$type == "metabolite"), 
					", Pathways: ", sum(igraph::V(binet)$type == "pathway"), " ...")

			binet
		}
	),
	private = list(
		# Load an RData file and return its content safely.
		# file_path path to the RData file.
		# expected_var optional; name hint for the expected variable (used in error messages).
		load_rdata = function(file_path, expected_var = NULL){
			env <- new.env()
			load(file_path, envir = env)
			var_names <- ls(envir = env)
			if(length(var_names) == 0){
					stop("No objects found in the RData file: ", file_path)
			}
			if(length(var_names) > 1){
					message("Note: Multiple objects found in RData file. Using the first one: ", var_names[1], " ...")
			}
			get(var_names[1], envir = env)
		},
		load_pathway_db = function(db_path, db_type){
			db_data <- private$load_rdata(db_path)
			
			required_fields <- c("pathways", "metab_to_pathway")
			if(!all(required_fields %in% names(db_data))){
					stop("Invalid database format. Required fields: ", 
						 paste(required_fields, collapse = ", "))
			}
			
			message("Loaded ", db_type, " database with ", 
					length(db_data$pathways), " pathways ...")
			
			db_data
		},
		fetch_kegg_pathways = function(){
			message("Attempting to fetch KEGG compound and pathway data via REST API ...")
			message("Note: This may take several minutes for the full KEGG database ...")
			
			kegg_pathways <- list()
			metab_to_pathway <- list()
			
			tryCatch({
					# Fetch pathways
					pathways_url <- "https://rest.kegg.jp/list/pathway"
					pathways_raw <- readLines(pathways_url, warn = FALSE)
					
					pathway_names <- list()
					for(line in pathways_raw){
							parts <- strsplit(line, "\t")[[1]]
							if(length(parts) >= 2){
									pid <- gsub("path:map", "map", parts[1])
									pname <- parts[2]
									pathway_names[[pid]] <- pname
							}
					}
					message("Successfully retrieved ", length(pathway_names), " KEGG pathways ...")
					
					# Fetch all compound-pathway relationships via KEGG link API
					message("Fetching compound-pathway mapping relationships from KEGG ...")
					link_url <- "https://rest.kegg.jp/link/pathway/cpd"
					link_raw <- readLines(link_url, warn = FALSE)
					
					for(line in link_raw){
							parts <- strsplit(line, "\t")[[1]]
							if(length(parts) >= 2){
									if(grepl("path:map", parts[2])){
											cid <- gsub("cpd:", "", parts[1])
											pid <- gsub("path:map", "map", parts[2])
											
											metab_to_pathway[[cid]] <- c(metab_to_pathway[[cid]], pid)
									}
							}
					}
					
					message("Successfully mapped ", length(metab_to_pathway), " compounds to pathways ...")
					
			}, error = function(e){
					stop("Failed to fetch KEGG data: ", conditionMessage(e), 
						 "\nPlease download KEGG database and provide local file path via db_path parameter !")
			})
			
			list(pathways = pathway_names, metab_to_pathway = metab_to_pathway)
		},
		clean_name = function(input) {
			output <- input %>%
				tolower %>%
				gsub("l-|r-", "", .) %>%
				gsub("\\(r\\)|\\(s\\)", "", .) %>%
				gsub("\\d+$", "", .) %>%
				gsub("\\s+", " ", .) %>%
				gsub("degr.*$", " ", .) %>%
				trimws

			output
		},
		check_input_path = function(directory_path){
			if(! dir.exists(directory_path)){
				stop("Provided path is not valid (", directory_path, ")! Please check the input path!")
			}
		},
		check_directory_exist = function(directory_path, directory_name, back_path){
			if(! dir.exists(directory_path)){
				stop(directory_name, " is not found in provided path (", back_path, ")! Please check the downloaded and decompressed folder!")
			}
		},
		check_file_exist = function(file_path, file_name, back_path){
			if(! file.exists(file_path)){
				stop(file_name, " is not found in provided path (", back_path, ")! Please check the downloaded and decompressed folder!")
			}
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
