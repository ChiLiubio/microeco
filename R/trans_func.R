#' @title
#' Create trans_func object for functional analysis.
#'
#' @description
#' This class is a wrapper for a series of functional analysis on species and communities, including the prokaryotes function identification based on Louca et al. (2016) <doi:10.1126/science.aaf4507> or fungi function identification based on Nguyen et al. (2016) <10.1016/j.funeco.2015.06.006>, functional redundancy calculation and metabolic pathway abundance prediction Aßhauer et al. (2015) <10.1093/bioinformatics/btv287>.
#'
#' @export
trans_func <- R6Class(classname = "trans_func",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @return for_what : "prok" or "fungi" or NA, "prok" represent prokaryotes. "fungi" represent fungi. NA represent not identified according to the Kingdom information, 
		#' at this time, if you want to use the functions to identify species traits, you need provide "prok" or "fungi" manually, e.g. dataset$for_what <- "prok".
		#' @examples
		#' data(dataset)
		#' t1 <- trans_func$new(dataset = dataset)
		initialize = function(dataset = NULL
			){
			self$tax_table <- dataset$tax_table
			self$otu_table <- dataset$otu_table
			self$sample_table <- dataset$sample_table
			# confirm the taxonomy automatically, for prokaryotes or fungi
			for_what <- NA
			if(!is.null(self$tax_table$Kingdom)){
				all_Kingdom <- unique(as.character(self$tax_table$Kingdom))
				if(any(grepl("Bacteria|Archaea", all_Kingdom, ignore.case = TRUE))){
					for_what <- "prok"
				}else{
					if(any(grepl("Fungi", all_Kingdom, ignore.case = TRUE))){
						for_what <- "fungi"
					}else{
						message("No Bacteria, Archaea or Fungi found in the Kingdom of tax_table, please set the for_what object use prok or fungi manually!")
					}
				}
			}else{
				message("No Kingdom found in the tax_table, please set the for_what object use prok or fungi manually!")
			}
			self$for_what <- for_what
		},
		#' @description
		#' Confirm traits of each OTU by matching the taxonomic assignments to the functional database;
		#' Prokaryotes: based on the FAPROTAX database, please also cite the original FAPROTAX paper: 
		#' Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272. <doi:10.1126/science.aaf4507>;
		#' Fungi, based on the FUNGuild database or FungalTraits database, please also cite:
		#' FUNGuild: Nguyen, N. H., Song, Z., Bates, S. T., Branco, S., Tedersoo, L., Menke, J., … Kennedy, P. G. (2016). 
		#' FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology, 20(1), 241–248. <doi:10.1016/j.funeco.2015.06.006>
		#' FungalTraits: Põlme, S., Abarenkov, K., Henrik Nilsson, R. et al. FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. 
		#' Fungal Diversity 105, 1–16 (2020). <doi:10.1007/s13225-020-00466-2>
		#'
		#' @param fungi_database default "FUNGuild"; select a fungi trait database for the trait identification, "FUNGuild" or "FungalTraits"; see the description in this function.
		#' @return res_spe_func in object.
		#' @examples
		#' \donttest{
		#' t1$cal_spe_func()
		#' }
		cal_spe_func = function(fungi_database = c("FUNGuild", "FungalTraits")[1]){
			for_what <- self$for_what
			if(is.na(for_what)){
				stop("No for_what object found, please set the for_what object use prok or fungi manually!")
			}
			if(! for_what %in% c("prok", "fungi")){
				stop("Wrong for_what found, please set the for_what object use prok or fungi !")
			}
			if(for_what == "prok"){
				# Copyright (c) 2020, Stilianos Louca
				# All rights reserved.
				# prok_func is a database developed based on the FAPROTAX database (http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Home)
				data("prok_func", envir=environment())
				message("The prokaryotic database is developed based on the FAPROTAX database Version: 1.2.4. Please also cite the original FAPROTAX paper: Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272.\n")
				# collapse taxonomy
				tax1 <- apply(self$tax_table, 1, function(x){paste0(x, collapse = ";")}) %>% gsub(".__", "", .) %>% gsub(";{1, }$", "", .)
				# reduce computational cost
				tax2 <- unique(tax1)
				# first create result matrix of tax2, then tax1-OTU
				res <- matrix(nrow = length(tax2), ncol = length(prok_func$func_tax))
				rownames(res) <- tax2
				colnames(res) <- names(prok_func$func_tax)
				# match taxa
				for(i in seq_len(ncol(res))){
					res[, i] <- grepl(prok_func$func_tax[i], tax2) %>% as.numeric
				}
				res %<>% as.data.frame
				# identify the inclusion part among groups
				for(i in names(prok_func$func_groups)){
					if(length(prok_func$func_groups[[i]]) == 0){
						next
					}else{
						for(j in seq_along(prok_func$func_groups[[i]])){
							res[, i] <- res[, i] + res[, prok_func$func_groups[[i]][j]]
						}
					}
				}
				# only use 1 represent exist
				res[res > 1] <- 1

				otu_func_table <- res[tax1, ]
				rownames(otu_func_table) <- names(tax1)
				otu_func_table$anaerobic_chemoheterotrophy <- 0
				otu_func_table$anaerobic_chemoheterotrophy[otu_func_table$chemoheterotrophy != 0 & otu_func_table$aerobic_chemoheterotrophy == 0] <- 1
			}
			if(for_what == "fungi"){
				if(grepl("FUNGuild", fungi_database, ignore.case = TRUE)){
					data("fungi_func_FUNGuild", envir=environment())
					tax1 <- self$tax_table
					tax1[] <- lapply(tax1, function(x) gsub(".*__", "", x))
					tax1 %<>% dropallfactors
					# add taxon column to store the target taxon
					tax1$taxon <- ""
					# operate the matching for each level that stored in the database
					for(i in c("Phylum", "Order", "Family", "Genus", "Species")){
						use_database <- switch(i, 
							Phylum = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel == "3", ], 
							Order = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel == "7", ], 
							Family = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel == "9", ],
							Genus = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel == "13", ],
							Species = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel %in% c("20", "21", "22", "24"), ])
						# search each OTU even though it has been matched
						for(j in rownames(tax1)){
							if(tax1[j, i] %in% use_database[, "taxon"]){
								tax1[j, "taxon"] <- tax1[j, i]
							}
						}
					}
					# merge two tables
					res_table <- dplyr::left_join(tax1, fungi_func_FUNGuild, by = c("taxon" = "taxon"))
					rownames(res_table) <- rownames(tax1)
					res_table <- res_table[, which(colnames(res_table) %in% "taxon"):ncol(res_table)]
					res_table[is.na(res_table)] <- ""
					# store the raw table similar with the FUNGuild results from python version
					self$res_spe_func_raw_funguild <- res_table
					# generate a data frame store the binary data
					otu_func_table <- res_table[, c("taxon"), drop = FALSE]
					# generate trophicMode binary information
					trophicMode <- c("Pathotroph", "Saprotroph", "Symbiotroph")
					for(i in trophicMode){
						otu_func_table[, i] <- grepl(i, res_table[, "trophicMode"]) %>% as.numeric
					}
					# generate Guild binary information
					Guild <- c("Bryophyte Parasite", "Dung Saprotroph", "Ectomycorrhizal", "Fungal Parasite", "Leaf Saprotroph", "Plant Parasite",
							"Wood Saprotroph", "Animal Pathogen", "Endophyte", "Plant Pathogen", "Lichen Parasite", "Litter Saprotroph", "Soil Saprotroph",
							"Plant Saprotroph", "Epiphyte", "Lichenized", "Arbuscular Mycorrhizal", "Endomycorrhizal", "Ericoid Mycorrhizal", "Orchid Mycorrhizal", 
							"Root Associated Biotroph", "Clavicipitaceous Endophyte", "Animal Endosymbiont", "Wood Saprotrop")
					for(i in Guild){
						otu_func_table[, i] <- grepl(i, res_table[, "guild"]) %>% as.numeric
					}
				}else{
					tax1 <- self$tax_table
					if(!any(grepl("Genus", colnames(tax1), fixed = TRUE))){
						stop("No Genus column found in tax_table of dataset!")
					}
					data("fungi_func_FungalTraits", envir=environment())
					message("This fungi database is developed based on the FungalTraits. Please also cite the original paper: FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. Fungal Diversity 105, 1–16 (2020). <doi:10.1007/s13225-020-00466-2>.\n")
					# remove the redundant data
					fungi_func_FungalTraits %<>% .[- c(3109, 3288, 4741, 7092), ]
					fungi_func_FungalTraits$Decay_type %<>% gsub("white-rot", "white_rot", .)
					fungi_func_FungalTraits$Fruitbody_type %<>% gsub("\\s", "_", .)
					# extract duplicated genus
					duplicated_genus <- fungi_func_FungalTraits$GENUS %>% .[duplicated(.)]
					FungalTraits_undup <- fungi_func_FungalTraits[!fungi_func_FungalTraits$GENUS %in% duplicated_genus, ]
					FungalTraits_undup <- FungalTraits_undup[, -c(1:5)]

					tax1[] <- lapply(tax1, function(x) gsub(".*__", "", x))
					res_table <- dplyr::left_join(tax1, FungalTraits_undup, by = c("Genus" = "GENUS"))
					rownames(res_table) <- rownames(tax1)
					
					# check the duplicated genus
					tax_all_genus <- tax1$Genus %>% .[. != ""] %>% unique
					# if duplicated, use multiple match
					if(any(tax_all_genus %in% duplicated_genus)){
						tax_names <- tax1 %>% .[.$Genus %in% duplicated_genus, ] %>% rownames
						FungalTraits_dup <- fungi_func_FungalTraits[fungi_func_FungalTraits$GENUS %in% duplicated_genus, ]
						FungalTraits_dup$GENUS <- paste0(FungalTraits_dup$Class, "|", FungalTraits_dup$GENUS)
						FungalTraits_dup <- FungalTraits_dup[, -c(1:5)]
						tax_dup <- tax1[tax_names, ]
						tax_dup$Genus <- paste0(tax_dup$Class, "|", tax_dup$Genus)
						res_table_dup <- dplyr::left_join(tax_dup, FungalTraits_dup, by = c("Genus" = "GENUS"))
						rownames(res_table_dup) <- tax_names
						res_table_dup$Genus %<>% gsub(".*\\|", "", .)
						res_table <- res_table[!rownames(res_table) %in% tax_names, ]
						res_table <- rbind.data.frame(res_table, res_table_dup)
						res_table %<>% .[rownames(tax1), ]
					}
					
					# store the raw table for personalized use
					self$res_spe_func_raw_FungalTraits <- res_table					
					# then parse the result for calculation
					filter_data <- res_table
					colnames(filter_data) %<>% gsub("_template", "", .)
					keep_cols <- c("Genus", "primary_lifestyle", "Secondary_lifestyle", "Endophytic_interaction_capability", "Plant_pathogenic_capacity", "Decay_substrate", 
						"Decay_type", "Aquatic_habitat", "Animal_biotrophic_capacity", "Growth_form", "Fruitbody_type", "Hymenium_type", 
						"Ectomycorrhiza_exploration_type", "primary_photobiont")
					filter_data <- filter_data[, keep_cols]
					filter_data[is.na(filter_data)] <- ""

					# generate a data frame store the binary data
					otu_func_table <- filter_data[, c("Genus"), drop = FALSE]
					# generate binary information
					for(i in names(private$default_fungi_func_group$FungalTraits)){
						for(j in private$default_fungi_func_group$FungalTraits[[i]]){
							j_name <- paste0(i, "|", j)
							otu_func_table[, j_name] <- grepl(j, filter_data[, i]) %>% as.numeric
						}
					}
				}
				otu_func_table %<>% .[, -1]
				self$fungi_database <- fungi_database
			}
			self$res_spe_func <- otu_func_table
		},
		#' @description
		#' Calculating the percentages of species with specific trait in communities or modules.
		#' The percentages of the OTUs with specific trait can reflect the potential of the corresponding function in the community or the module in the network.
		#'
		#' @param use_community default TRUE; whether calculate community; if FALSE, use module.
		#' @param node_type_table default NULL; If use_community FALSE; provide the node_type_table with the module information, such as the result of cal_node_type.
		#' @param abundance_weighted default FALSE; whether use abundance. If FALSE, calculate the functional population percentage. If TRUE, calculate the functional individual percentage.
		#' @return res_spe_func_perc in object.
		#' @examples
		#' \donttest{
		#' t1$cal_spe_func_perc(use_community = TRUE)
		#' }
		cal_spe_func_perc = function(use_community = TRUE, node_type_table = NULL, abundance_weighted = FALSE){
			res_spe_func <- self$res_spe_func
			otu_table <- self$otu_table
			if(use_community){
				res_spe_func_perc <- sapply(colnames(otu_table), function(input_samplecolumn){
					sample_otu <- otu_table[, input_samplecolumn]
					names(sample_otu) <- rownames(otu_table)
					# remove species whose abundance is 0
					sample_otu <- sample_otu[sample_otu != 0]
					res_table <- unlist(lapply(colnames(res_spe_func), function(x){
						if(abundance_weighted){
							(res_spe_func[names(sample_otu), x, drop = TRUE] * sample_otu) %>% 
								{sum(.) * 100/sum(sample_otu)} %>% 
								{round(., 2)}
						}else{
							res_spe_func[names(sample_otu), x, drop = TRUE] %>% 
								{sum(. != 0) * 100/length(.)} %>% 
								{round(., 2)}
						}
					}))
					res_table
				})
			}else{
				if(is.null(node_type_table)){
					stop("No node_type_table provided! parameter: node_type_table !")
				}else{
					if(abundance_weighted){
						otu_total_abund <- apply(otu_table, 1, sum)
					}
					# calculate for each module
					res_spe_func_perc <- sapply(sort(unique(as.character(node_type_table$module))), function(input_module){
							otu_names <- rownames(node_type_table[node_type_table$module == input_module, ])
							res_table <- unlist(lapply(colnames(res_spe_func), function(x){
								if(abundance_weighted){
									(res_spe_func[otu_names, x, drop = TRUE] * otu_total_abund[otu_names]) %>% 
										{sum(.) * 100/sum(otu_total_abund[otu_names])} %>% 
										{round(., 2)}
								}else{
									res_spe_func[otu_names, x, drop = TRUE] %>% 
										{sum(. != 0) * 100/length(.)} %>% 
										{round(., 2)}
								}
							}))
						res_table
					})
				}
			}
			res_spe_func_perc %<>% t %>% as.data.frame %>% `colnames<-`(colnames(res_spe_func)) %>% .[, apply(., 2, sum) != 0]
			self$res_spe_func_perc <- res_spe_func_perc
		},
		#' @description
		#' Show the basic information for a specific function of prokaryotes.
		#'
		#'
		#' @param use_func default NULL; the function name.
		#' @return None.
		#' @examples
		#' \donttest{
		#' t1$show_prok_func(use_func = "methanotrophy")
		#' }
		show_prok_func = function(use_func = NULL){
			data("prok_func", envir=environment())
			if(!is.null(use_func)){
				prok_func$func_annotation[[use_func]]
			}
		},
		#' @description
		#' Plot the percentages of species with specific trait in communities or modules.
		#'
		#' @param filter_func default NULL; a vector of function names used to show in the plot.
		#' @param use_group_list default TRUE; If TRUE, use default group list; If use personalized group list, 
		#'    first set trans_func$func_group_list object with a list of group names and functions.
		#' @param add_facet default TRUE; whether use group names as the facets in the plot, see trans_func$func_group_list object.
		#' @param select_samples default NULL; character vector, select partial samples to show
		#' @return ggplot2.
		#' @examples
		#' \donttest{
		#' t1$plot_spe_func_perc(use_group_list = TRUE)
		#' }
		plot_spe_func_perc = function(filter_func = NULL, use_group_list = TRUE, add_facet = TRUE, select_samples = NULL){
			plot_data <- self$res_spe_func_perc
			if(!is.null(filter_func)){
				plot_data <- plot_data[, colnames(plot_data) %in% filter_func]
			}
			plot_data <- reshape2::melt(cbind.data.frame(sampname = rownames(plot_data), plot_data), id.vars = "sampname") %>% dropallfactors
			# add group and factor
			if(use_group_list){
				if(self$for_what == "fungi"){
					if(grepl("FUNGuild", self$fungi_database, ignore.case = TRUE)){
						func_group_list_use <- self$func_group_list[["FUNGuild"]]
					}else{
						func_group_list_use_raw <- self$func_group_list[["FungalTraits"]]
						func_group_list_use <- lapply(names(func_group_list_use_raw), function(x){
							paste0(x, "|", func_group_list_use_raw[[x]])
						})
						names(func_group_list_use) <- names(func_group_list_use_raw)
					}
				}else{
					func_group_list_use <- self$func_group_list
				}
				
				group_table <- reshape2::melt(func_group_list_use) %>% dropallfactors
				colnames(group_table) <- c("func", "group")
				filter_func <- group_table$func
				plot_data <- plot_data[plot_data$variable %in% group_table$func, ]
				plot_data <- dplyr::left_join(plot_data, group_table, by = c("variable" = "func"))
				plot_data$group %<>% factor(., levels = names(func_group_list_use))
			}
			# make func order better
			if(!is.null(filter_func)){
				plot_data$variable %<>% gsub("_", " ", .) %>% factor(., levels = gsub("_", " ", filter_func))			
			}
			if(!is.null(select_samples)){
				plot_data %<>% .[.$sampname %in% select_samples, , drop = FALSE]
				plot_data$sampname %<>% factor(., levels = select_samples)
			}
			if(self$for_what == "fungi"){
				if(grepl("FungalTraits", self$fungi_database, ignore.case = TRUE)){
					plot_data$variable %<>% gsub(".*\\|", "", .)
				}
			}
			g1 <- ggplot(aes(x=sampname, y=variable, fill=value), data=plot_data) + 
				theme_bw() + 
				geom_tile() + 
				scale_fill_gradient2(low="#00008B", high="#9E0142") +
				scale_y_discrete(position = "right") + 
				labs(y=NULL, x=NULL, fill="Percentage (%)") +
				theme(axis.text.x = element_text(angle = 35, colour = "black", vjust = 1, hjust = 1, size = 14), axis.text.y = element_text(size = 10)) +
				theme(panel.grid = element_blank(), panel.border = element_blank()) +
				theme(panel.spacing = unit(.1, "lines")) + 
				theme(plot.margin=unit(c(1, 0, 0, 1), "cm"))
				
			if(use_group_list & add_facet){
				g1 <- g1 + facet_grid(group ~ ., drop=TRUE, scale="free",space="free", switch = "y") +
					theme(strip.background = element_rect(fill = "grey95", colour = "white"), strip.text.y.left = element_text(angle=360), strip.text=element_text(size=14))
			}
			g1
		},
		#' @description
		#' Predict functional potential of communities using FAPROTAX.
		#' please also cite the original FAPROTAX paper: 
		#' Louca, S., Parfrey, L. W., & Doebeli, M. (2016). Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272. <doi:10.1126/science.aaf4507>;
		#'
		#' @param keep_tem default FALSE; whether keep the intermediate file, that is, the otu_table_for_FAPROTAX.txt in local place.
		#' @param Ref_folder default "./FAPROTAX_1.2.1"; see http://www.loucalab.com/archive/FAPROTAX
		#' @return res_FAPROTAX in object.
		cal_FAPROTAX = function(keep_tem = TRUE, Ref_folder = "./FAPROTAX_1.2.1") {
			message("This is FAPROTAX database 1.2.1 with python 2.7. The newer versions may exist in http://www.loucalab.com/archive/FAPROTAX ")
			otu_file <- self$otu_table
			tax_file <- self$tax_table
			otu_file <- data.frame("#OTU ID" = rownames(otu_file), otu_file, check.names = FALSE, stringsAsFactors = FALSE)
			tax_file <- apply(tax_file, 1, function(x){paste0(x, collapse = ";")})
			otu_file <- data.frame(otu_file, taxonomy = tax_file, check.names = FALSE, stringsAsFactors = FALSE)
			# save to local place
			pathfilename <- "otu_table_for_FAPROTAX"
			pathfilename <- tempfile(pathfilename, fileext = ".txt")
			message("writing the otu_table_for_FAPROTAX to temporary file ", pathfilename, " ...")
			write.table(otu_file, pathfilename, sep = "\t", row.names = FALSE, quote = FALSE)
			code_path <- Ref_folder
			use_command <- paste0("python ", code_path, "/collapse_table.py -i ", pathfilename, " -o ", "FAPROTAX_prediction.tsv -g ", 
				code_path, "/FAPROTAX.txt -d taxonomy --omit_columns 0 --column_names_are_in last_comment_line -f")
			message("run python to predict...")
			system(use_command)
			message("save prediction result in FAPROTAX_prediction.tsv ...")
			self$res_FAPROTAX <- read.delim("FAPROTAX_prediction.tsv", check.names = FALSE, row.names = 1)
			if(keep_tem == F){
				message("remove intermediate file otu_table_for_FAPROTAX.txt ...")
				unlink(pathfilename, recursive = FALSE, force = TRUE)
			}
		},
		#' @description
		#' Predict functional potential of communities using tax4fun.
		#' please also cite: Aßhauer, K. P., Wemheuer, B., Daniel, R., & Meinicke, P. (2015). 
		#' Tax4Fun: Predicting functional profiles from metagenomic 16S rRNA data. Bioinformatics, 31(17), 2882–2884. <doi:10.1093/bioinformatics/btv287>
		#'
		#' @param keep_tem default FALSE; whether keep the intermediate file, that is, the otu table in local place.
		#' @param folderReferenceData default NULL; the folder, see http://tax4fun.gobics.de/ and Tax4Fun function in Tax4Fun package.
		#' @return tax4fun_KO and tax4fun_path in object.
		cal_tax4fun = function(keep_tem = FALSE, folderReferenceData = NULL){
			if(is.null(folderReferenceData)){
				stop("No folderReferenceData provided! Please see the help document!")
			}
			if(!require(Tax4Fun)){
				stop("Tax4Fun package not installed! see http://tax4fun.gobics.de/ ")
			}
			otu_file <- self$otu_table
			tax_file <- self$tax_table
			otu_file <- data.frame("#OTU ID" = rownames(otu_file), otu_file, check.names = FALSE, stringsAsFactors = FALSE)
			tax_file <- apply(tax_file, 1, function(x){paste0(x, collapse = ";")})

			otu_file <- data.frame(otu_file, taxonomy = tax_file, check.names = FALSE, stringsAsFactors = FALSE)
			otu_file$taxonomy %<>% gsub(".__", "", .) %>% paste0(., ";") %>% gsub(";{1, }$", ";", .)

			pathfilename <- "otu_table_filter_tax4fun"
			pathfilename <- tempfile(pathfilename, fileext = ".txt")
			output <- file(pathfilename, open = "wb")
			# must write this line, otherwise sample line will disappear
			cat("# Constructed from biom file\n", file = output)
			write.table(otu_file, file = output, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
			close(output)

			x1 <- importQIIMEData(pathfilename)
			self$tax4fun_KO <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = TRUE)
			self$tax4fun_path <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = FALSE)
		
			if(keep_tem == F){
				unlink(pathfilename, recursive = FALSE, force = TRUE)
			}
		},
		#' @description
		#' Print the trans_func object.
		print = function(){
			cat("trans_func class:\n")
			if(!is.null(self$env_data)){
				cat(paste("otu_func_table have", nrow(self$otu_func_table), "rows and", ncol(self$otu_func_table), "columns\n"))
				cat("\n")
			}
		}
	),
	active = list(
		#### Active ----
		
		#' @field func_group_list store and show the function group list
		func_group_list = function(func_list) {
			if (missing(func_list)){
				switch(self$for_what, prok = private$default_prok_func_group, 
					fungi = private$default_fungi_func_group, NA)
			}else{
				if(self$for_what == "prok"){
					private$default_prok_func_group <- func_list
				}else{
					if(self$for_what == "fungi"){
						private$default_fungi_func_group <- func_list
					}
				}
			}
		}
	),
	private = list(
		default_prok_func_group = list(
			"Energy source" = c("aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy", "photoautotrophy", "photoheterotrophy"),
			"C-cycle" = c("chitinolysis", "cellulolysis", "fermentation",  "methanogenesis", "methanotrophy", "methylotrophy"),
			"N-cycle" = c("nitrogen_fixation", "aerobic_ammonia_oxidation","nitrification","aerobic_nitrite_oxidation", 
				"nitrate_reduction","nitrate_respiration","nitrite_respiration"),
			"S-cycle" = c("sulfate_respiration", "sulfur_respiration", "sulfite_respiration","dark_sulfide_oxidation", "respiration_of_sulfur_compounds"),
			"Others" = c("dark_hydrogen_oxidation", "iron_respiration", "manganese_oxidation", "fumarate_respiration")
		),
		default_fungi_func_group = list(
			FUNGuild = list(
				"Trophic Mode" = c("Pathotroph", "Saprotroph", "Symbiotroph"),
				"Guild" = c("Bryophyte Parasite", "Dung Saprotroph", "Ectomycorrhizal", "Fungal Parasite", "Leaf Saprotroph", "Plant Parasite",
					"Wood Saprotroph", "Animal Pathogen", "Endophyte", "Plant Pathogen", "Lichen Parasite", "Litter Saprotroph", "Soil Saprotroph",
					"Plant Saprotroph", "Epiphyte", "Lichenized", "Arbuscular Mycorrhizal", "Endomycorrhizal", "Ericoid Mycorrhizal", "Orchid Mycorrhizal", 
					"Root Associated Biotroph", "Clavicipitaceous Endophyte", "Animal Endosymbiont", "Wood Saprotrop")
			),
			FungalTraits = list(
				"primary_lifestyle" = c("algal_ectosymbiont", "algal_parasite", "algivorous/protistivorous", "animal-associated", "animal_endosymbiont", "animal_parasite", 
					"arbuscular_mycorrhizal", "arthropod-associated", "bacterivorous", "dung_saprotroph", "ectomycorrhizal", "epiphyte", "foliar_endophyte", 
					"lichen_parasite", "lichenized", "litter_saprotroph", "moss_symbiont", "mycoparasite", "nectar/tap_saprotroph", "plant_pathogen", "pollen_saprotroph", 
					"protistan_parasite", "root_endophyte", "soil_saprotroph", "sooty_mold", "unspecified_pathotroph", "unspecified_saprotroph", 
					"unspecified_symbiotroph", "wood_saprotroph"),
				"Secondary_lifestyle" = c("algal_decomposer", "algal_parasite", "algal_symbiont", "animal-associated", "animal_decomposer", "animal_parasite", 
					"arbuscular_mycorrhizal", "arthropod-associated", "arthropod_parasite", "bryophilous", "coral-associated", "dung_saprotroph", "epiphyte", 
					"ericoid_mycorrhizal", "fatty_acid_producer", "fish_parasite", "foliar_endophyte", "fungal_decomposer", "insect-associated", 
					"invertebrate-associated", "invertebrate_parasite", "lichen_parasite", "litter_saprotroph", "liverwort-associated", "moss_parasite", 
					"mycoparasite", "myxomycete_decomposer", "nectar/tap_saprotroph", "nematophagous", "plant_pathogen", "pollen_saprotroph", "protistan_parasite", 
					"resin_saprotroph", "rock-inhabiting", "root-associated", "root_endophyte", "root_endophyte_dark_septate", "soil_saprotroph", "sooty_mold", 
					"termite_symbiont", "unsepcified_saprotroph", "unspecified_saprotroph", "unspecified_symbiotroph", "vertebrate-associated", "wood_saprotroph"),
				"Endophytic_interaction_capability" = c("class1_clavicipitaceous_endophyte", "foliar_endophyte", 
					"moss-associated", "no_endophytic_capacity", "root-associated", "root_endophyte", "root_endophyte_dark_septate"),
				"Plant_pathogenic_capacity" = c("algal_parasite", "leaf/fruit/seed_pathogen", "moss-associated", "moss_parasite", "other_plant_pathogen", "root-associated", 
					"root_pathogen", "wood_pathogen"),
				"Decay_substrate" = c("algal_material", "animal_material", "pollen", "dung", "protist_material", "soil", "sugar-rich_substrates", "fungal_material", "leaf/fruit/seed", 
					"burnt_material", "hydrocarbon-rich_substrates_(fuel)", "roots", "wood", "resins"),
				"Decay_type" = c("blue-staining", "brown_rot", "chitinolytic", "keratinolytic", "mold", "other", "soft_rot", "white_rot"),
				"Aquatic_habitat" = c("aquatic", "freshwater", "marine", "non-aquatic", "partly_aquatic", "partly_freshwater_(partly_non-aquatic)", "partly_marine_(partly_non-aquatic)"),
				"Animal_biotrophic_capacity" = c("animal-associated", "opportunistic_human_parasite", "animal_endosymbiont", "animal_parasite", 
					 "arthropod-associated", "arthropod_ectosymbiont", "arthropod_parasite", "coral-associated", "fish_parasite", "invertebrate-associated", 
					 "invertebrate_parasite", "nematophagous", "opportunistic_animal_parasite", "opportunistic_animal_parasite_G._destructans", 
					 "opportunistic_human_parasite_Pythium_insidiosum", "termite_symbiont", "vertebrate-associated", "vertebrate_parasite"),
				"Growth_form" = c("amoeboid-biflagellate", "biflagellate", "biflagellate-rhizomycelial", "dimorphic_yeast", "filamentous_mycelium", 
					"labyrinthulid", "other", "thallus_photosynthetic", "unicellular-aflagellate_(non-yeast)", "yeast", "zoosporic-plasmodium_(aphelid)", 
					"zoosporic-rhizomycelial_(chytrid-like)", "zoosporic_only"),
				"Fruitbody_type" = c("agaricoid", "apothecium_(hymenium_on_surface)", "chasmothecium_(tiny_initially_closed_ascome_rupture_opening)_", 
					"clathroid", "clavarioid", "cleistothecium", "cleistothecium_(closed,_spherical)", "corticioid", "cyphelloid", "gasteroid", "gasteroid-hypogeous", 
					"hysterothecium_(opening_slith-like)", "mazaedium_(pushpin-like)", "other", "perithecium_(hymenium_hidden,_narrow_opening)", 
					"phalloid", "polyporoid", "rust", "smut", "thyrothecium_(tiny_flattened_disc,_aperture_in_the_middle)", "tremelloid", "zoosporangium", "zygosporangium"),
				"Hymenium_type" = c("closed", "gel", "gills", "hydnoid", "none", "other", "poroid", "smooth"),
				"Ectomycorrhiza_exploration_type" = c("contact", "long-distance", "mat", "medium-distance_fringe", "medium-distance_smooth", 
					"short-distance_coarse", "short-distance_delicate", "unknown"),
				"primary_photobiont" = c("chlorococcoid", "cyanobacterial", "trentepohlioid")
			)
		)
	),
	lock_class = FALSE,
	lock_objects = FALSE
)



