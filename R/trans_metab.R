#' @title
#' Create \code{trans_metab} object for metabolite analysis.
#'
#' @description
#' This class is a wrapper for a series of metabolite analysis.
#'
#' @export
trans_metab <- R6Class(classname = "trans_metab",
	public = list(
		#' @description
		#' Create the \code{trans_metab} object.
		#' 
		#' @param metab default NULL; metabolite data. A data.frame object or the \code{\link{microtable}} object.
		#' @param microb default NULL; A \code{\link{microtable}} object.
		#' @return data inside the object.
		initialize = function(
			metab = NULL,
			microb = NULL
			){
			self$data_metab <- metab
			self$data_microb <- microb
		},
		#' @description
		#' Metabolite origin inference based on the preprocessed database from TidyMass2 (DOI: 10.1038/s41467-026-68464-7)
		#' 
		#' @param database_path default "./metorigindb_split_202602"; directory path of the downloaded database. 
		#'	  Please first download the pre-collated metorigindb database (RData format) from zenodo (https://zenodo.org/records/18618912) and extract the compressed archive.
		#' @param match_col default "names"; How to match to the data of metorigindb. Default "names" means using the input names of metabolites.
		#'    If the table has other columns like "HMDB_ID" or "KEGG_ID", the user can provide more items, like c("names", "HMDB_ID").
		#' @param bac_level default "Genus"; which bacteria level is used to parse the taxa.
		cal_origin = function(
			database_path = "./metorigindb_split_202602", 
			match_col = "names",
			bac_level = "Genus"
			){
			data_metab <- self$data_metab
			data_microb <- self$data_microb
			
			load(file.path(database_path, "metorigindb_number", "metorigindb_number.RData"))
			
			`%in2%` <- function(x, vec) {!is.na(x) & (x %in% vec)}
			
			length_db <- length(metorigindb_number$Compound_name)
			if("names" %in% match_col){
				name_ID_select <- metorigindb_number$Compound_name %in2% rownames(data_metab$tax_table)
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
				all_taxa <- data_microb$tax_table[, bac_level] %>% unique %>% .[. != paste0(bac_level_prefix, "__")] %>% gsub(".__", "", .)
				res_origin %<>% lapply(., function(x){x[x %in% all_taxa]})
			}
			
			self$res_origin_rawtable <- source_table
			message("Raw origin table is stored in object$res_origin_rawtable ...")
			self$res_origin_list <- res_origin
			if(! is.null(data_microb)){
				message("Filtered origin taxa at ", bac_level, " level for each metabolite is stored in object$res_origin_list ...")
			}else{
				message("Origin taxa at ", bac_level, " level for each metabolite is stored in object$res_origin_list ... ...")
			}
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
