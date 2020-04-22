#' Create trans_func object for functional analysis.
#'
#' This class is a wrapper for a series of functional analysis for communities and species.
#' The functions in this class include \code{\link{cal_spe_func}}, \code{\link{cal_spe_func_perc}}, \code{\link{show_spe_func}},
#' \code{\link{plot_spe_func_perc}}, \code{\link{cal_tax4fun_func}}, \code{\link{cal_biogeo}}.
#'
#' @param dataset the object of \code{\link{microtable}} Class.
#' @return trans_func object with tax_table, otu_table, sample_table.
#' @examples 
#' t1 <- trans_diff$new(dataset = dataset)
#' @export
trans_func <- R6Class(classname = "trans_func",
	public = list(
		initialize = function(dataset = NULL
			){
			self$tax_table <- dataset$tax_table
			self$otu_table <- dataset$otu_table
			self$sample_table <- dataset$sample_table
			},
		cal_biogeo = function(keep_tem = TRUE) {
			otu_file <- self$otu_table
			tax_file <- self$tax_table
			otu_file <- data.frame("#OTU ID" = rownames(otu_file), otu_file, check.names = FALSE, stringsAsFactors = FALSE)
			tax_file <- apply(tax_file, 1, function(x){paste0(x, collapse = ";")})
			otu_file <- data.frame(otu_file, taxonomy = tax_file, check.names = FALSE, stringsAsFactors = FALSE)
			# save to local place
			message("writing the otu_table_for_FAPROTAX.txt for FAPROTAX prediction...")
			write.table(otu_file, "otu_table_for_FAPROTAX.txt", sep = "\t", row.names = FALSE, quote = FALSE)
			code_path <- system.file("extdata", "FAPROTAX_1.2.1", package="microeco")
			use_command <- paste0("python ", code_path, "/collapse_table.py -i otu_table_for_FAPROTAX", ".txt -o ", "FAPROTAX_prediction.tsv -g ", code_path, "/FAPROTAX.txt -d taxonomy --omit_columns 0 --column_names_are_in last_comment_line -f")
			message("run python to predict...")
			system(use_command)
			message("save prediction result FAPROTAX_prediction.tsv ...")
			self$res_biogeo <- read.delim("FAPROTAX_prediction.tsv", check.names = FALSE, row.names = 1)
			if(keep_tem == F){
				message("remove intermediate file otu_table_for_FAPROTAX.txt ...")
				unlink("otu_table_for_FAPROTAX.txt", recursive = FALSE, force = TRUE)
			}
		},
		cal_spe_func = function(){
			data(spe_func)
			# collapse taxonomy
			tax1 <- apply(self$tax_table, 1, function(x){paste0(x, collapse = ";")}) %>% gsub(".__", "", .) %>% gsub(";{1, }$", "", .)
			# reduce computational cost
			tax2 <- unique(tax1)
			# first create result matrix of tax2, then tax1-OTU
			res <- matrix(nrow = length(tax2), ncol = length(spe_func$func_tax))
			rownames(res) <- tax2
			colnames(res) <- names(spe_func$func_tax)
			# match taxa
			for(i in seq_len(ncol(res))){
				res[, i] <- grepl(spe_func$func_tax[i], tax2) %>% as.numeric
			}
			res %<>% as.data.frame
			# identify the inclusion part among groups
			for(i in names(spe_func$func_groups)){
				if(length(spe_func$func_groups[[i]]) == 0){
					next
				}else{
					for(j in seq_along(spe_func$func_groups[[i]])){
						res[, i] <- res[, i] + res[, spe_func$func_groups[[i]][j]]
					}
				}
			}
			# only use 1 represent exist
			res[res > 1] <- 1

			otu_func_table <- res[tax1, ]
			rownames(otu_func_table) <- names(tax1)
			otu_func_table$anaerobic_chemoheterotrophy <- 0
			otu_func_table$anaerobic_chemoheterotrophy[otu_func_table$chemoheterotrophy != 0 & otu_func_table$aerobic_chemoheterotrophy == 0] <- 1			
			self$otu_func_table <- otu_func_table
		},
		cal_spe_func_perc = function(use_community = TRUE, node_type_table = NULL){
			otu_func_table <- self$otu_func_table
			if(use_community){
				x1 <- self$otu_table
				res_spe_func_perc <- sapply(colnames(x1), function(input_samplecolumn){
					z1 <- x1[, input_samplecolumn]
					names(z1) <- rownames(x1)
					# remove species whose abundance is 0
					z1 <- z1[z1 != 0]
					z2 <- unlist(lapply(colnames(otu_func_table), function(x){
						otu_func_table[names(z1), x, drop = TRUE] %>% {sum(. != 0) * 100/length(.)} %>% {round(., 2)}
					}))
					z2
				})
			}else{
				if(is.null(node_type_table)){
					stop("No node_type_table provided! parameter: node_type_table !")
				}else{
					x1 <- node_type_table
					res_spe_func_perc <- sapply(sort(as.character(unique(x1$module))), function(input_samplecolumn){
						z1 <- rownames(x1[x1$module == input_samplecolumn, ])
						z2 <- unlist(lapply(colnames(otu_func_table), function(x){
							otu_func_table[z1, x, drop = TRUE] %>% {sum(. != 0) * 100/length(.)} %>% {round(., 2)}
						}))
						z2
					})
				}
			}
			res_spe_func_perc %<>% t %>% as.data.frame %>% `colnames<-`(colnames(otu_func_table)) %>% .[, apply(., 2, sum) != 0]
			self$res_spe_func_perc <- res_spe_func_perc
		},
		show_spe_func = function(use_func = NULL){
			if(!is.null(use_func)){
				spe_func$func_annotation[[use_func]]
			}
		},
		plot_spe_func_perc = function(filter_func = NULL, group_list = NULL, group_list_default = FALSE){
			plot_data <- self$res_spe_func_perc
			if(!is.null(filter_func)){
				plot_data <- plot_data[, colnames(plot_data) %in% filter_func]
			}
			plot_data <- reshape2::melt(cbind.data.frame(sampname = rownames(plot_data), plot_data), id.vars = "sampname") %>% dropallfactors
			if(group_list_default == T){
				group_list <- private$default_func_group
			}
			# add group and factor
			if(!is.null(group_list)){
				group_table <- reshape2::melt(group_list) %>% dropallfactors
				colnames(group_table) <- c("func", "group")
				filter_func <- group_table$func
				plot_data <- plot_data[plot_data$variable %in% group_table$func, ]
				plot_data <- dplyr::left_join(plot_data, group_table, by = c("variable" = "func"))
				plot_data$group %<>% factor(., levels = names(group_list))
			}
			# make func order better
			if(!is.null(filter_func)){
				plot_data$variable %<>% gsub("_", " ", .) %>% factor(., levels = gsub("_", " ", filter_func))			
			}
			g1 <- ggplot(aes(x=sampname, y=variable, fill=value), data=plot_data) + 
				theme_bw() + 
				geom_tile() + 
				scale_fill_gradient2(low="#00008B", high="#9E0142") +
				scale_y_discrete(position = "right") + 
				labs(y=NULL, x=NULL, fill="Percentage (%)") +
				theme(axis.text.x = element_text(angle = 35, colour = "black", vjust = 1, hjust = 1, size = 14), axis.text.y = element_text(size = 10)) +
				theme(panel.grid = element_blank(), panel.border = element_blank()) +
				theme(panel.spacing = unit(.1, "lines")) + theme(plot.margin=unit(c(1, 0, 0, 1), "cm"))
				
			if(!is.null(group_list)){
				g1 <- g1 + facet_grid(group ~ ., drop=TRUE, scale="free",space="free", switch = "y") +
				theme(strip.background = element_rect(fill = "grey95", colour = "white"), strip.text.y = element_text(angle=180), strip.text=element_text(size=14))
			}
			g1
		},
		cal_tax4fun_func = function(keep_tem = FALSE, folderReferenceData = NULL){
			
			otu_file <- self$otu_table
			tax_file <- self$tax_table
			otu_file <- data.frame("#OTU ID" = rownames(otu_file), otu_file, check.names = FALSE, stringsAsFactors = FALSE)
			tax_file <- apply(tax_file, 1, function(x){paste0(x, collapse = ";")})

			otu_file <- data.frame(otu_file, taxonomy = tax_file, check.names = FALSE, stringsAsFactors = FALSE)
			otu_file$taxonomy %<>% gsub(".__", "", .) %>% paste0(., ";") %>% gsub(";{1, }$", ";", .)

			pathfilename <- "otu_table_filter_tax4fun.txt"
			output <- file(pathfilename, open = "wb")
			# must write this line, otherwise sample line will disappear
			cat("# Constructed from biom file\n", file = output)
			write.table(otu_file, file = output, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
			close(output)
			require(Tax4Fun)
			x1 <- importQIIMEData(pathfilename)
			self$tax4fun_KO <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = TRUE)
			self$tax4fun_path <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = FALSE)
		
			if(keep_tem == F){
				unlink(pathfilename, recursive = FALSE, force = TRUE)
			}
		},
		print = function(){
			cat("trans_func class:\n")
			if(!is.null(self$env_data)){
				cat(paste("otu_func_table have", nrow(self$otu_func_table), "rows and", ncol(self$otu_func_table), "columns\n"))
				cat("\n")
			}
		}
	),
	private = list(
		default_func_group = list(
			"Energy source" = c("aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy", "photoautotrophy", "photoheterotrophy"),
			"C-cycle" = c("chitinolysis", "cellulolysis", "fermentation",  "methanogenesis", "methanotrophy", "methylotrophy"),
			"N-cycle" = c("nitrogen_fixation", "aerobic_ammonia_oxidation","nitrification","aerobic_nitrite_oxidation","nitrate_reduction","nitrate_respiration","nitrite_respiration"),
			"S-cycle" = c("sulfate_respiration", "sulfur_respiration", "sulfite_respiration","dark_sulfide_oxidation", "respiration_of_sulfur_compounds"),
			"Others" = c("dark_hydrogen_oxidation", "iron_respiration", "manganese_oxidation", "fumarate_respiration")
		)
	),
	lock_class = FALSE,
	lock_objects = FALSE
)



#' Confirm species traits.
#'
#' Confirm traits of each OTU by matching the taxonomic assignments to the FAPROTAX database.
#'
#' @return otu_func_table in object.
#' @examples
#' t1$cal_spe_func()
cal_spe_func <- function(){
	dataset$cal_spe_func()
}



#' Calculating the percentages of species with specific trait in communities or modules.
#'
#' The percentages of the OTUs with specific trait can reflect the potential of the corresponding function in the community or the module in the network.
#'
#' @param use_community default TRUE; whether calculate community; if FALSE, use module.
#' @param node_type_table default NULL; If use_community FALSE; provide the node_type_table with the module information, such as the result of \code{\link{cal_node_type}}.
#' @return res_spe_func_perc in object.
#' @examples
#' t1$cal_spe_func_perc(use_community = TRUE)
cal_spe_func_perc <- function(use_community = TRUE, node_type_table = NULL){
	dataset$cal_spe_func_perc()
}

#' Show the basic information for one function.
#'
#'
#' @param use_func default NULL; the function name.
#' @return None.
#' @examples
#' t1$show_spe_func(use_func = "methanotrophy")
show_spe_func <- function(use_community = TRUE, node_type_table = NULL){
	dataset$show_spe_func()
}

#' Plot the percentages of species with specific trait in communities or modules.
#'
#'
#' @param filter_func default NULL; a vector of function names.
#' @param group_list default NULL; a list with group names and the functions in the groups.
#' @param group_list_default default FALSE; whether use the default group list.
#' @return ggplot2.
#' @examples
#' t1$plot_spe_func_perc(group_list_default = TRUE)
plot_spe_func_perc <- function(filter_func = NULL, group_list = NULL, group_list_default = FALSE){
	dataset$plot_spe_func_perc()
}

#' Predict functional potential using tax4fun.
#'
#'
#' @param keep_tem default FALSE; whether keep the intermediate file, that is, the otu table in local place.
#' @param folderReferenceData default NULL; the folder, see http://tax4fun.gobics.de/ and \code{\link{Tax4Fun}} function in Tax4Fun package.
#' @return tax4fun_KO and tax4fun_path in object.
#' @examples
#' t1$cal_tax4fun_func(folderReferenceData = "./SILVA123")
cal_tax4fun_func <- function(keep_tem = FALSE, folderReferenceData = NULL){
	dataset$cal_tax4fun_func()
}

#' Predict functional potential using FAPROTAX.
#'
#'
#' @param code_path default "./FAPROTAX_1.2.1"; the code folder, download from http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Download.
#' @param keep_tem default FALSE; whether keep the intermediate file, that is, the otu_table_for_FAPROTAX.txt in local place.
#' @return res_biogeo in object.
#' @examples
#' t1$cal_biogeo(code_path = "./FAPROTAX_1.2.1")
cal_biogeo <- function(code_path = "./FAPROTAX_1.2.1", keep_tem = TRUE){
	dataset$cal_biogeo()
}

