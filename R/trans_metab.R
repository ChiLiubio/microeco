#' @title
#' Create \code{trans_metab} object for metabolite analysis.
#'
#' @description
#' This class is a wrapper for a series of metabolite analysis, including origin inference.
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
		#' \donttest{
		#' data(soil_metab)
		#' data(soil_microb)
		#' t1 <- trans_env$new(metab = soil_metab, microb = soil_microb)
		#' }
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
							rownames(metab_tax_table) <- rownames(metab)
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
		cal_match = function(
			database_path = "./metorigindb_split_202602",
			method = "jw",
			maxDist = 0.3,
			...
			){
			data_metab <- self$data_metab
			
			metorigindb_number_path <- file.path(database_path, "metorigindb_number", "metorigindb_number.RData")
			if(! file.exists(metorigindb_number_path)){
				stop("metorigindb_number directory and metorigindb_number.RData is not found in provided path (", database_path, ")! Please provide a correct database_path!")
			}
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
		},
		#' @description
		#' Metabolite origin inference based on the preprocessed database from TidyMass2 (DOI: 10.1038/s41467-026-68464-7)
		#' 
		#' @param database_path default "./metorigindb_split_202602"; directory path of the downloaded database. 
		#'	  Please download the pre-collated metorigindb database (RData format) from zenodo (https://zenodo.org/records/18618912) and extract the compressed archive.
		#' @param match_col default "names"; How to match to the data of metorigindb. Default "names" means using the input names of metabolites.
		#'    If the table has other columns like "HMDB_ID" or "KEGG_ID", the user can provide more items, like c("names", "HMDB_ID").
		#' @param match_names_distance default 0; distance threshold used if the \code{res_match_table} is found in the object, 
		#'    which is calculated from the \code{cal_match} function. Available for the "names" option in \code{match_col} parameter.
		#' @param bac_level default "Genus"; which bacteria level is used to parse the taxa in the \code{data_microb} of object.
		#'    The function can automatically match the taxa those found in the input data.		
		cal_origin = function(
			database_path = "./metorigindb_split_202602", 
			match_col = "names",
			match_names_distance = 0,
			bac_level = "Genus"
			){
			data_metab <- self$data_metab
			data_microb <- self$data_microb
			
			metorigindb_number_path <- file.path(database_path, "metorigindb_number", "metorigindb_number.RData")
			if(! file.exists(metorigindb_number_path)){
				stop("metorigindb_number directory and metorigindb_number.RData is not found in provided path (", database_path, ")! Please provide a correct database_path!")
			}
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
			if(! dir.exists(source_path)){
				stop("metorigindb_source directory is not found in provided path (", database_path, ")! Please provide a correct database_path!")
			}

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
		}
	),
	private = list(
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
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
