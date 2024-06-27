#' @title
#' Create \code{trans_func} object for functional prediction.
#'
#' @description
#' This class is a wrapper for a series of functional prediction analysis on species and communities, including the prokaryotic trait prediction based on 
#' Louca et al. (2016) <doi:10.1126/science.aaf4507> and Lim et al. (2020) <10.1038/s41597-020-0516-5>, 
#' or fungal trait prediction based on Nguyen et al. (2016) <10.1016/j.funeco.2015.06.006> and Polme et al. (2020) <doi:10.1007/s13225-020-00466-2>; 
#' functional redundancy calculation and metabolic pathway abundance prediction Abhauer et al. (2015) <10.1093/bioinformatics/btv287>.
#'
#' @export
trans_func <- R6Class(classname = "trans_func",
	public = list(
		#' @description
		#' Create the \code{trans_func} object. This function can identify the data type for Prokaryotes or Fungi automatically.
		#' 
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @return \code{for_what}: "prok" or "fungi" or NA, "prok" represent prokaryotes. "fungi" represent fungi. NA stand for unknown according to the Kingdom information. 
		#' In this case, if the user still want to use the function to identify species traits, please provide "prok" or "fungi" manually, 
		#' e.g. \code{t1$for_what <- "prok"}.
		#' @examples
		#' data(dataset)
		#' t1 <- trans_func$new(dataset = dataset)
		initialize = function(dataset = NULL
			){
			check_microtable(dataset)
			self$tax_table <- dataset$tax_table
			self$otu_table <- dataset$otu_table
			self$sample_table <- dataset$sample_table
			self$rep_fasta <- dataset$rep_fasta
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
						message("No Bacteria, Archaea or Fungi found in the Kingdom of tax_table in the input dataset! ",
							"Please assign for_what object using prok or fungi manually, such as object$for_what <- 'fungi' !")
					}
				}
			}else{
				message("No Kingdom column found in the tax_table! Please assign for_what object using prok or fungi manually, such as object$for_what <- 'fungi'")
			}
			self$for_what <- for_what
		},
		#' @description
		#' Identify traits of each feature by matching taxonomic assignments to functional database.
		#'
		#' @param prok_database default "FAPROTAX"; \code{"FAPROTAX"} or \code{"NJC19"}; select a prokaryotic trait database:
		#'   \describe{
		#'     \item{\strong{'FAPROTAX'}}{FAPROTAX; Reference: Louca et al. (2016). Decoupling function and taxonomy in the global ocean microbiome. 
		#'     	  Science, 353(6305), 1272. <doi:10.1126/science.aaf4507>}
		#'     \item{\strong{'NJC19'}}{NJC19: Lim et al. (2020). Large-scale metabolic interaction network of the mouse and human gut microbiota. 
		#'     	  Scientific Data, 7(1). <10.1038/s41597-020-0516-5>. 
		#'     	  Note that the matching in this database is performed at the species level, 
		#'     	  hence utilizing it demands a higher level of precision in regards to the assignments of species in the taxonomic information table.}
		#'   }
		#' @param fungi_database default "FUNGuild"; \code{"FUNGuild"} or \code{"FungalTraits"}; select a fungal trait database:
		#'   \describe{
		#'     \item{\strong{'FUNGuild'}}{Nguyen et al. (2016) FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild.
		#'     	  Fungal Ecology, 20(1), 241-248, <doi:10.1016/j.funeco.2015.06.006>}
		#'     \item{\strong{'FungalTraits'}}{version: FungalTraits_1.2_ver_16Dec_2020V.1.2; Polme et al. 
		#'     	  FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles.  
		#'     	  Fungal Diversity 105, 1-16 (2020). <doi:10.1007/s13225-020-00466-2>}
		#'   }
		#' @param FUNGuild_confidence default c("Highly Probable", "Probable", "Possible"). 
		#'    Selected 'confidenceRanking' when \code{fungi_database = "FUNGuild"}.
		#' @return \code{res_spe_func} stored in object.
		#' @examples
		#' \donttest{
		#' t1$cal_spe_func(prok_database = "FAPROTAX")
		#' t1$cal_spe_func(fungi_database = "FungalTraits")
		#' }
		cal_spe_func = function(
			prok_database = c("FAPROTAX", "NJC19")[1], 
			fungi_database = c("FUNGuild", "FungalTraits")[1],
			FUNGuild_confidence = c("Highly Probable", "Probable", "Possible")
			){
			for_what <- self$for_what
			if(is.na(for_what) | is.null(for_what)){
				stop("No for_what object found, please assign the for_what object use prok or fungi manually!")
			}
			if(! for_what %in% c("prok", "fungi")){
				stop("Wrong for_what, please make sure for_what object is one of prok and fungi !")
			}
			tax1 <- self$tax_table
			if(for_what == "prok"){
				if(grepl("FAPROTAX", prok_database, ignore.case = TRUE)){
					# Copyright (c) 2023, Stilianos Louca. All rights reserved.
					# developed based on the FAPROTAX database (http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Home)
					data("prok_func_FAPROTAX", envir = environment())
					message("FAPROTAX v", prok_func_FAPROTAX$ver, ". Please also cite the original FAPROTAX paper: Louca et al. (2016).")
					message("Decoupling function and taxonomy in the global ocean microbiome. Science, 353(6305), 1272.\n")
					# collapse taxonomy
					tax1 <- apply(tax1, 1, function(x){paste0(x, collapse = ";")}) %>% 
						gsub(".__", "", .) %>% 
						gsub(";{1, }$", "", .)
					# reduce computational cost
					tax2 <- unique(tax1)
					# first create result matrix of tax2, then tax1-OTU
					res <- matrix(nrow = length(tax2), ncol = length(prok_func_FAPROTAX$func_tax))
					rownames(res) <- tax2
					colnames(res) <- names(prok_func_FAPROTAX$func_tax)
					# match taxa
					for(i in seq_len(ncol(res))){
						res[, i] <- grepl(prok_func_FAPROTAX$func_tax[i], tax2) %>% as.numeric
					}
					res %<>% as.data.frame
					# identify the inclusion part among groups
					for(i in names(prok_func_FAPROTAX$func_add_groups)){
						if(length(prok_func_FAPROTAX$func_add_groups[[i]]) == 0){
							next
						}else{
							for(j in seq_along(prok_func_FAPROTAX$func_add_groups[[i]])){
								res[, i] <- res[, i] + res[, prok_func_FAPROTAX$func_add_groups[[i]][j]]
							}
						}
					}
					# identify the subtraction part among groups
					for(i in names(prok_func_FAPROTAX$func_subtract_groups)){
						if(length(prok_func_FAPROTAX$func_subtract_groups[[i]]) == 0){
							next
						}else{
							for(j in seq_along(prok_func_FAPROTAX$func_subtract_groups[[i]])){
								res[, i] <- res[, i] - res[, prok_func_FAPROTAX$func_subtract_groups[[i]][j]]
							}
						}
					}
					# only use 1 to represent existence
					res[res > 1] <- 1
					res[res < 0] <- 0
					otu_func_table <- res[tax1, ]
					rownames(otu_func_table) <- names(tax1)
					otu_func_table$anaerobic_chemoheterotrophy <- 0
					otu_func_table$anaerobic_chemoheterotrophy[otu_func_table$chemoheterotrophy != 0 & otu_func_table$aerobic_chemoheterotrophy == 0] <- 1
				}else{
					tax1[] <- lapply(tax1, function(x) gsub(".*__", "", x))
					if(!any(grepl("Species", colnames(tax1)))){
						stop("No Species column in the tax_table !")
					}
					if(all(tax1$Species == "")){
						stop("No species name exist in the Species column!")
					}
					if(!any(grepl("\\s", tax1$Species))){
						stop("No blank space found in the Species column! Please check whether species names are right!")
					}
					
					data("prok_func_NJC19_list", envir=environment())
					# list to frame
					frame1 <- data.frame()
					for(i in names(prok_func_NJC19_list)){
						x1 <- prok_func_NJC19_list[[i]]
						for(j in names(x1)){
							x2 <- data.frame(i, x1[[j]], j, stringsAsFactors = FALSE)
							colnames(x2) <- c("Species", "trait", "Type")
							frame1 <- rbind(frame1, x2)
						}
					}
					frame1$use_trait <- paste(frame1[, 2], frame1[, 3], sep = " -- ")
					frame1$value <- 1
					frame1 <- frame1[, -c(2:3)]
					frame2 <- suppressMessages(reshape2::dcast(frame1, Species ~ use_trait, value.var="value"))

					otu_func_table <- dplyr::left_join(tax1, frame2, by = c("Species" = "Species"))
					rownames(otu_func_table) <- rownames(tax1)
					otu_func_table <- otu_func_table[, -c(1:which(colnames(otu_func_table) == "Species"))]
					otu_func_table[is.na(otu_func_table)] <- 0
				}
				self$database <- prok_database
			}
			if(for_what == "fungi"){
				if(grepl("FUNGuild", fungi_database, ignore.case = TRUE)){
					data("fungi_func_FUNGuild", envir=environment())
					message("Please also cite the original paper: Nguyen et al. (2016).")
					message("FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology, 20(1), 241-248. \n")
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
							Species = fungi_func_FUNGuild[fungi_func_FUNGuild$taxonomicLevel %in% c("20", "21", "22", "24"), ]
						)
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
					message('Mapped raw FUNGuild result is stored in object$res_spe_func_raw_funguild ...')
					
					if(! any(FUNGuild_confidence %in% c("Highly Probable", "Probable", "Possible"))){
						stop("FUNGuild_confidence parameter must comes from 'Highly Probable', 'Probable' and 'Possible' !")
					}
					res_table[! res_table$confidenceRanking %in% FUNGuild_confidence, ] <- ""
					if(all(res_table$confidenceRanking == "")){
						if(all(FUNGuild_confidence %in% c("Highly Probable", "Probable", "Possible"))){
							message("No available results are predicted ...")
						}else{
							stop("No available result! Please change the FUNGuild_confidence parameter!")
						}
					}
					# generate a data frame store the binary data
					otu_func_table <- res_table[, c("taxon"), drop = FALSE]
					# generate trophicMode binary information
					trophicMode <- c("Pathotroph", "Saprotroph", "Symbiotroph")
					for(i in trophicMode){
						otu_func_table[, i] <- grepl(i, res_table[, "trophicMode"]) %>% as.numeric
					}
					# generate Guild binary information
					Guild <- private$default_fungi_func_group$FUNGuild[["Guild"]]
					for(i in Guild){
						otu_func_table[, i] <- grepl(i, res_table[, "guild"]) %>% as.numeric
					}
				}else{
					if(!any(grepl("Genus", colnames(tax1), fixed = TRUE))){
						stop("No Genus column found in tax_table of dataset!")
					}
					data("fungi_func_FungalTraits", envir=environment())
					message("Please also cite: ")
					message("FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. Fungal Diversity 105, 1-16 (2020).\n")
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
					message('Mapped raw FungalTraits result is stored in object$res_spe_func_raw_FungalTraits ...')
					# then parse the result for calculation
					filter_data <- res_table
					colnames(filter_data) %<>% gsub("_template", "", .)
					keep_cols <- c(
						"Genus", 
						"primary_lifestyle", 
						"Secondary_lifestyle", 
						"Endophytic_interaction_capability", 
						"Plant_pathogenic_capacity", 
						"Decay_substrate", 
						"Decay_type", 
						"Aquatic_habitat", 
						"Animal_biotrophic_capacity", 
						"Growth_form", 
						"Fruitbody_type", 
						"Hymenium_type", 
						"Ectomycorrhiza_exploration_type", 
						"primary_photobiont")
					filter_data <- filter_data[, keep_cols]
					filter_data[is.na(filter_data)] <- ""

					# generate data frame store the binary data
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
				self$database <- fungi_database
			}
			if(any(is.na(otu_func_table))){
				warning("NA found in the final table! Convert NA to 0 ...")
				otu_func_table[is.na(otu_func_table)] <- 0
			}
			self$res_spe_func <- otu_func_table
			message('The functional binary table is stored in object$res_spe_func ...')
		},
		#' @description
		#' Calculating the percentages of species with specific trait in communities.
		#' The percentages of the taxa with specific trait can reflect corresponding functional potential in the community.
		#' So this method is one representation of functional redundancy (FR) without the consideration of phylogenetic distance among taxa.
		#' The FR is defined:
		#'      \deqn{FR_{kj}^{unweighted} = \frac{N_{j}}{N_{k}}}
		#'      \deqn{FR_{kj}^{weighted} = \frac{\sum_{i=1}^{N_{j}} A_{i}}{\sum_{i=1}^{N_{k}} A_{i}}}
		#' where \eqn{FR_{kj}} denotes the FR for sample k and function j. \eqn{N_{k}} is the species number in sample k.
		#' \eqn{N_{j}} is the number of species with function j in sample k.
		#' \eqn{A_{i}} is the abundance (counts) of species i in sample k.
		#' 
		#' @param abundance_weighted default FALSE; whether use abundance of taxa. If FALSE, calculate the functional population percentage. 
		#' 	  If TRUE, calculate the functional individual percentage.
		#' @param perc default TRUE; whether to use percentages in the result. If TRUE, value is bounded between 0 and 100.
		#' 	  If FALSE, the result is relative proportion (`abundance_weighted = FALSE`) or relative abundance (`abundance_weighted = TRUE`) bounded between 0 and 1.
		#' @param dec default 2; remained decimal places.
		#' @return \code{res_spe_func_perc} stored in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_spe_func_perc(abundance_weighted = TRUE)
		#' }
		cal_spe_func_perc = function(abundance_weighted = FALSE, perc = TRUE, dec = 2){
			if(is.null(self$res_spe_func)){
				stop("Please first run cal_spe_func function !")
			}
			bound_value <- ifelse(perc, 100, 1)
			res_spe_func <- self$res_spe_func
			otu_table <- self$otu_table
			res_spe_func_perc <- sapply(colnames(otu_table), function(input_samplecolumn){
				sample_otu <- otu_table[, input_samplecolumn]
				names(sample_otu) <- rownames(otu_table)
				# remove species whose abundance is 0
				sample_otu <- sample_otu[sample_otu != 0]
				res_table <- unlist(lapply(colnames(res_spe_func), function(x){
					if(abundance_weighted){
						(res_spe_func[names(sample_otu), x, drop = TRUE] * sample_otu) %>% 
							{sum(.) * bound_value/sum(sample_otu)} %>% 
							{round(., dec)}
					}else{
						res_spe_func[names(sample_otu), x, drop = TRUE] %>% 
							{sum(. != 0) * bound_value/length(.)} %>% 
							{round(., dec)}
					}
				}))
				res_table
			})
			res_spe_func_perc %<>% t %>% 
				as.data.frame %>% 
				`colnames<-`(colnames(res_spe_func)) %>% 
				.[, apply(., 2, sum) != 0]
			self$res_spe_func_perc <- res_spe_func_perc
			message('The result table is stored in object$res_spe_func_perc ...')
		},
		#' @description
		#' Show the annotation information for a function of prokaryotes from FAPROTAX database.
		#'
		#'
		#' @param use_func default NULL; the function name.
		#' @return None.
		#' @examples
		#' \donttest{
		#' t1$show_prok_func(use_func = "methanotrophy")
		#' }
		show_prok_func = function(use_func = NULL){
			data("prok_func_FAPROTAX", envir=environment())
			if(!is.null(use_func)){
				prok_func_FAPROTAX$func_annotation[[use_func]]
			}
		},
		#' @description
		#' Transform the \code{res_spe_func_perc} table to the long table format for the following visualization.
		#' Also add the group information if the database has hierarchical groups.
		#' 
		#' @return \code{res_spe_func_perc_trans} stored in the object.
		#' @examples
		#' \donttest{
		#' t1$trans_spe_func_perc()
		#' }
		trans_spe_func_perc = function(){
			if(is.null(self$res_spe_func_perc)){
				stop("Please first run cal_spe_func and cal_spe_func_perc function !")
			}
			res_spe_func_perc <- self$res_spe_func_perc

			trans_perc <- reshape2::melt(cbind.data.frame(sampname = rownames(res_spe_func_perc), res_spe_func_perc), id.vars = "sampname") %>% dropallfactors
			# add group and factor
			if(self$for_what == "fungi"){
				if(grepl("FUNGuild", self$database, ignore.case = TRUE)){
					func_group_list_use <- self$func_group_list[["FUNGuild"]]
				}else{
					func_group_list_use_raw <- self$func_group_list[["FungalTraits"]]
					func_group_list_use <- lapply(names(func_group_list_use_raw), function(x){
						paste0(x, "|", func_group_list_use_raw[[x]])
					})
					names(func_group_list_use) <- names(func_group_list_use_raw)
				}
			}else{
				if(grepl("FAPROTAX", self$database, ignore.case = TRUE)){
					func_group_list_use <- self$func_group_list[["FAPROTAX"]]
				}else{
					func_group_list_use <- NULL
				}
			}
			if(!is.null(func_group_list_use)){
				group_table <- reshape2::melt(func_group_list_use) %>% dropallfactors
				colnames(group_table) <- c("func", "group")
				trans_perc <- dplyr::left_join(trans_perc, group_table, by = c("variable" = "func"))
				if(grepl("FAPROTAX", self$database, ignore.case = TRUE)){
					trans_perc %<>% .[.$variable != "chemoheterotrophy", ]
					if(any(is.na(trans_perc$group))){
						trans_perc$group[is.na(trans_perc$group)] <- "Other"
						trans_perc$group %<>% factor(., levels = c(names(func_group_list_use), "Other"))
					}else{
						trans_perc$group %<>% factor(., levels = names(func_group_list_use))
					}
				}else{
					trans_perc$group %<>% factor(., levels = names(func_group_list_use))
				}
			}
			if(self$for_what == "fungi"){
				if(grepl("FungalTraits", self$database, ignore.case = TRUE)){
					trans_perc$variable %<>% gsub(".*\\|", "", .)
				}
			}
			self$res_spe_func_perc_trans <- trans_perc
			message('Transformed long format table is stored in object$res_spe_func_perc_trans ...')
		},
		#' @description
		#' Plot the percentages of species with specific trait in communities.
		#'
		#' @param add_facet default TRUE; whether use group names as the facets in the plot, see \code{trans_func$func_group_list} object.
		#' @param order_x default NULL; character vector; to sort the x axis text; can be also used to select partial samples to show.
		#' @param color_gradient_low default "#00008B"; the color used as the low end in the color gradient.
		#' @param color_gradient_high default "#9E0142"; the color used as the high end in the color gradient.
		#' @return ggplot2.
		#' @examples
		#' \donttest{
		#' t1$plot_spe_func_perc()
		#' }
		plot_spe_func_perc = function(
			add_facet = TRUE, 
			order_x = NULL,
			color_gradient_low = "#00008B",
			color_gradient_high = "#9E0142"
			){
			if(is.null(self$res_spe_func_perc_trans)){
				message("The res_spe_func_perc_trans object is not found. Run the trans_spe_func_perc function to get it ...")
				self$trans_spe_func_perc()
			}
			plot_data <- self$res_spe_func_perc_trans
			
			if(!is.null(order_x)){
				plot_data %<>% .[.$sampname %in% order_x, , drop = FALSE]
				plot_data$sampname %<>% factor(., levels = order_x)
			}
			
			g1 <- ggplot(aes(x = sampname, y = variable, fill = value), data = plot_data) + 
				theme_bw() + 
				geom_tile() + 
				scale_fill_gradient2(low = color_gradient_low, high = color_gradient_high) +
				scale_y_discrete(position = "right") + 
				labs(y = NULL, x = NULL, fill = "Percentage (%)") +
				theme(axis.text.x = element_text(angle = 35, colour = "black", vjust = 1, hjust = 1, size = 14), axis.text.y = element_text(size = 10)) +
				theme(panel.grid = element_blank(), panel.border = element_blank()) +
				theme(panel.spacing = unit(.1, "lines")) + 
				theme(plot.margin=unit(c(1, 0, 0, 1), "cm"))
				
			if(add_facet){
				if("group" %in% colnames(plot_data)){
					g1 <- g1 + facet_grid(group ~ ., drop=TRUE, scale="free",space="free", switch = "y") +
						theme(strip.background = element_rect(fill = "grey95", colour = "white"), strip.text.y.left = element_text(angle = 360), strip.text = element_text(size = 14))
				}
			}
			g1
		},
		#' @description
		#' Predict functional potential of communities using \code{tax4fun} package.
		#' please cite: Tax4Fun: Predicting functional profiles from metagenomic 16S rRNA data. Bioinformatics, 31(17), 2882-2884, <doi:10.1093/bioinformatics/btv287>.
		#' Note that this function requires a standard prefix in taxonomic table with double underlines (e.g. 'g__') .
		#'
		#' @param keep_tem default FALSE; whether keep the intermediate file, that is, the feature table in local place.
		#' @param folderReferenceData default NULL; the folder, see \href{http://tax4fun.gobics.de/}{http://tax4fun.gobics.de/}  and Tax4Fun function in Tax4Fun package.
		#' @return \code{tax4fun_KO} and \code{tax4fun_path} in object.
		cal_tax4fun = function(keep_tem = FALSE, folderReferenceData = NULL){
			if(is.null(folderReferenceData)){
				stop("No folderReferenceData provided! Please see the help document!")
			}
			if(!require("Tax4Fun")){
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
			suppressWarnings(write.table(otu_file, file = output, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE))
			close(output)

			x1 <- importQIIMEData(pathfilename)
			self$tax4fun_KO <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = TRUE)
			message('The KO abundance result is stored in object$tax4fun_KO ...')
			self$tax4fun_path <- Tax4Fun(x1, folderReferenceData = folderReferenceData, fctProfiling = FALSE)
			message('The pathway abundance result is stored in object$tax4fun_path ...')

			if(keep_tem == F){
				unlink(pathfilename, recursive = FALSE, force = TRUE)
			}
		},
		#' @description
		#' Predict functional potential of communities with Tax4Fun2 method. 
		#'   The function was adapted from the raw Tax4Fun2 package to make it compatible with the microtable object.
		#'   Pleas cite: 
		#'   Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences. Environmental Microbiome 15, 11 (2020).
		#' 	 <doi:10.1186/s40793-020-00358-7>
		#'
		#' @param blast_tool_path default NULL; the folder path, e.g., ncbi-blast-2.5.0+/bin ; blast tools folder downloaded from 
		#'   "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+"  ; e.g., ncbi-blast-2.5.0+-x64-win64.tar.gz  for windows system; 
		#'   if blast_tool_path is NULL, search the tools in the environmental path variable.
		#' @param path_to_reference_data default "Tax4Fun2_ReferenceData_v2"; the path that points to files used in the prediction; 
		#'   The directory must contain the Ref99NR or Ref100NR folder; 
		#'   download Ref99NR.zip from "https://cloudstor.aarnet.edu.au/plus/s/DkoZIyZpMNbrzSw/download" or 
		#'   Ref100NR.zip from "https://cloudstor.aarnet.edu.au/plus/s/jIByczak9ZAFUB4/download".
		#' @param path_to_temp_folder default NULL; The temporary folder to store the logfile, intermediate file and result files; if NULL, 
		#' 	 use the default temporary in the computer system.
		#' @param database_mode default 'Ref99NR'; "Ref99NR" or "Ref100NR"; Ref99NR: 99\% clustering reference database; Ref100NR: no clustering.
		#' @param normalize_by_copy_number default TRUE; whether normalize the result by the 16S rRNA copy number in the genomes. 
		#' @param min_identity_to_reference default 97; the sequences identity threshold used for finding the nearest species.
		#' @param use_uproc default TRUE; whether use UProC to functionally anotate the genomes in the reference data.
		#' @param num_threads default 1; the threads used in the blastn.
		#' @param normalize_pathways default FALSE; Different to Tax4Fun, when converting from KEGG functions to KEGG pathways, 
		#' 	 Tax4Fun2 does not equally split KO gene abundances between pathways a functions is affiliated to. The full predicted abundance is affiliated to each pathway. 
		#'   Use TRUE to split the abundances (default is FALSE).
		#' @return \code{res_tax4fun2_KO} and \code{res_tax4fun2_pathway} in object.
		#' @examples
		#' \dontrun{
		#' t1$cal_tax4fun2(blast_tool_path = "ncbi-blast-2.5.0+/bin", 
		#'     path_to_reference_data = "Tax4Fun2_ReferenceData_v2")
		#' }
		cal_tax4fun2 = function(
			blast_tool_path = NULL,
			path_to_reference_data = "Tax4Fun2_ReferenceData_v2",
			path_to_temp_folder = NULL,
			database_mode = 'Ref99NR',
			normalize_by_copy_number = T,
			min_identity_to_reference = 97, 
			use_uproc = T,
			num_threads = 1,
			normalize_pathways = F
			){
			if(is.null(path_to_temp_folder)){
				path_to_temp_folder <- tempdir()
				message("The intermediate files are saved in the temporary directory --- ", path_to_temp_folder)
			}
			if(!dir.exists(path_to_temp_folder)){
				stop("Temporay folder--", path_to_temp_folder, " is not found! Please first create it!")
			}
			if(!dir.exists(path_to_reference_data)){
				stop("Tax4Fun2 ReferenceData folder--", path_to_reference_data, " is not existed!")
			}
			if(is.null(self$rep_fasta)){
				stop("The rep_fasta is missing in your microtable object! It is necessary for Tax4Fun2! Use help(microtable) to see the rep_fasta description!")
			}
			# first check whether blastn tool is available
			if(is.null(blast_tool_path)){
				blast_bin <- "blastn"
			}else{
				blast_bin <- file.path(blast_tool_path, "blastn")
				if (tolower(Sys.info()[["sysname"]]) != "windows"){
					system(paste("chmod +x", blast_bin))
				}else{
					blast_bin <- file.path(blast_tool_path, "blastn.exe")
				}
			}
			res <- system(command = paste(blast_bin, "-help"), intern = T)
			if(length(res) == 0){
				blast_bin <- "blastn"
				res <- system(command = paste(blast_bin, "-help"), intern = T)
				if(length(res) == 0){
					stop(blast_bin, " not found! Please check the file path!")
				}
			}
			# Choose which refernence data set is used
			if(! database_mode %in% c("Ref99NR", "Ref100NR")){
				stop('Please provide valid database_mode! Must be Ref99NR or Ref100NR!')
			}else{
				path_to_ref_db <- file.path(path_to_reference_data, paste0(database_mode, "/", database_mode,'.fasta'))
				path_to_ref_dir <- file.path(path_to_reference_data, database_mode)
			}
			
			self$res_tax4fun2_database_mode <- database_mode
			message('database_mode is stored in object$res_tax4fun2_database_mode.')
			self$res_tax4fun2_path_to_ref_dir <- path_to_ref_dir

			if(!file.exists(path_to_ref_db)){
				stop("Reference database--", path_to_ref_db, " not found!")
			}
			# check whether the blastdb is available!
			check_files <- paste0(database_mode, '.fasta.', c("ndb", "nhr", "nin", "not", "nsq", "ntf", "nto"))
			check_files_path <- file.path(path_to_reference_data, paste0(database_mode, "/", check_files))
			if(!all(file.exists(check_files_path))){
				# use makeblastdb
				if(is.null(blast_tool_path)){
					makeblastdb_bin <- "makeblastdb"
				}else{
					makeblastdb_bin <- file.path(blast_tool_path, "makeblastdb")
					if (tolower(Sys.info()[["sysname"]]) != "windows"){
						system(paste("chmod +x", makeblastdb_bin))
					}else{
						makeblastdb_bin <- file.path(blast_tool_path, "makeblastdb.exe")
					}
				}
				res <- system(command = paste(makeblastdb_bin, "-help"), intern = T)
				if(length(res) == 0){
					makeblastdb_bin <- "makeblastdb"
					res <- system(command = paste(makeblastdb_bin, "-help"), intern = T)
					if(length(res) == 0){
						stop(makeblastdb_bin, " not found! Please check the file path!")
					}
				}
				# use refernence db
				message('Generate blast reference database.')
				cmd <- paste(makeblastdb_bin, '-dbtype nucl -in', path_to_ref_db)
				if (tolower(Sys.info()[["sysname"]]) == "windows"){
					system(cmd, show.output.on.console = F)
				}else{
					system(cmd, ignore.stdout = T, ignore.stderr = T)
				}
				message('Reference database finished.')
			}
			# write the fasta file
			rep_fasta_path <- file.path(path_to_temp_folder, "rep_fasta.tmp.fasta")
			if(inherits(self$rep_fasta, "list")){
				seqinr::write.fasta(self$rep_fasta, names = names(self$rep_fasta), file.out = rep_fasta_path)
			}else{
				if(inherits(self$rep_fasta, "DNAStringSet")){
					Biostrings::writeXStringSet(x = self$rep_fasta, filepath = rep_fasta_path)
				}else{
					stop("Unknown fasta format! Must be either list (from read.fasta of seqinr package) or DNAStringSet (from readDNAStringSet of Biostrings package)!")
				}
			}
			
			res_blast_path <- file.path(path_to_temp_folder, 'ref_blast.txt')
			message('Blast start.')
			cmd <- paste(blast_bin, '-db', path_to_ref_db, '-query', rep_fasta_path, '-evalue 1e-20 -max_target_seqs 100 -outfmt 6 -out', res_blast_path, '-num_threads', num_threads)
			if (tolower(Sys.info()[["sysname"]]) == "windows"){
				system(cmd, show.output.on.console = F)
			}else{
				system(cmd, ignore.stdout = T, ignore.stderr = T)
			}
			message('Blast finished.')

			# Write to log file
			path_to_log_file = file.path(path_to_temp_folder, 'logfile.txt')
			write(x = "RefBlast", file = path_to_log_file, append = F)
			write(x = date(), file = path_to_log_file, append = T)
			write(x = database_mode, file = path_to_log_file, append = T)
			write(x = rep_fasta_path, file = path_to_log_file, append = T)

			# filter redundant data
			private$blastTableReducer(path_to_blast_file = res_blast_path)
			message('Cleanup finished.')

			# parse file and make the prediction
			write(x = "Functional prediction", file = path_to_log_file, append = T)
			if(min_identity_to_reference < 90){
				warning("Minimum identity of less than 90% will likly results in inaccurate predictions!")
			}
			message(paste0("Using minimum identity cutoff of ", min_identity_to_reference, "% to nearest neighbor."))
			ref_blast_result <- read.delim(res_blast_path, h = F)
			ref_blast_result_reduced <- ref_blast_result[which(ref_blast_result$V3 >= min_identity_to_reference), 1:2]

			# Reading and filtering the otu table
			otu_table <- self$otu_table %>% tibble::rownames_to_column()
			# for the taxa mapping
			raw_otu_table_reduced <- merge(x = ref_blast_result_reduced, y = otu_table, by.x = "V1", by.y = "rowname")
			otu_table_reduced <- raw_otu_table_reduced[, -1]
			# for the calculation
			otu_table_reduced_aggregated <- aggregate(x = otu_table_reduced[, -1, drop = FALSE], by = list(otu_table_reduced[,1]), sum)
			self$res_tax4fun2_otu_table_reduced_aggregated <- otu_table_reduced_aggregated

			# Write unknown fraction to log file
			if((ncol(otu_table) - 1) == 1){
				unknown_fraction1 = as.data.frame(round(1 - sum(ifelse(otu_table_reduced[,-1]>0,1,0)) / sum(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
				write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
				write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
				unknown_fraction2 = as.data.frame(round(1 - sum(otu_table_reduced_aggregated[,-1]) / sum(otu_table[,-1]), digits = 5))
				write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
				write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
			} else {
				unknown_fraction1 = as.data.frame(round(1 - colSums(ifelse(otu_table_reduced[,-1]>0,1,0)) / colSums(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
				write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
				write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
				unknown_fraction2 = as.data.frame(round(1 - colSums(otu_table_reduced_aggregated[,-1]) / colSums(otu_table[,-1]), digits = 5))
				write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
				write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
			}

			# normalize or not normalize is the question
			n = 1
			if(use_uproc) n = 3
			if(normalize_by_copy_number) n = n + 1

			# Generate reference profile
			message('Generating reference profile.')
			reference_profile = NULL
			for(reference_id in otu_table_reduced_aggregated$Group.1){
				reference_file_path <- file.path(path_to_ref_dir, paste0(reference_id, '.tbl.gz'))
				reference_file <- read.delim(file = reference_file_path)
				reference_profile <- rbind(reference_profile, as.numeric(reference_file[, n]))
			}
			self$res_tax4fun2_reference_profile <- reference_profile
			
			# all the required KEGG files are stored in Tax4Fun2_KEGG Rdata
			data("Tax4Fun2_KEGG", envir=environment())
			ko_list <- Tax4Fun2_KEGG$ko_list

			map_reference_profile <- data.frame(id = otu_table_reduced_aggregated$Group.1, reference_profile)
			res_tax4fun2_reference_profile <- dplyr::left_join(raw_otu_table_reduced[, 1:2], map_reference_profile, by = c("V2" = "id"))
			colnames(res_tax4fun2_reference_profile) <- c("Taxa", "id", ko_list$ko)
			write.table(res_tax4fun2_reference_profile, file.path(path_to_temp_folder, 'res_tax4fun2_reference_profile.tsv'), sep = "\t")
			message("Reference profile file is saved in ", file.path(path_to_temp_folder, "res_tax4fun2_reference_profile.tsv"), ".")
			
			# Calculate functional profiles sample-wise
			message('Generating functional profile for samples.')
			functional_prediction = NULL
			for(sample_num in 2:ncol(otu_table_reduced_aggregated)){
				functional_prediction_sample <- reference_profile * as.numeric(otu_table_reduced_aggregated[, sample_num])
				functional_prediction_sample <- colMeans(functional_prediction_sample)
				functional_prediction_sample <- functional_prediction_sample / sum(functional_prediction_sample)
				if(is.na(sum(functional_prediction_sample))){
					functional_prediction_sample[1:nrow(ko_list)] <- 0
				}
				functional_prediction <- cbind(functional_prediction, functional_prediction_sample)
			}

			colnames(functional_prediction) <- names(otu_table)[2:ncol(otu_table_reduced_aggregated)]
			functional_prediction_final <- data.frame(KO = ko_list$ko, functional_prediction, description = ko_list$description)
			if(ncol(functional_prediction) >= 2) keep <- which(rowSums(functional_prediction) > 0)
			if(ncol(functional_prediction) == 1) keep <- which(functional_prediction > 0)
			if (length(keep) == 0){
				stop("No functional prediction possible!\nEither no nearest neighbor found or your table is empty!")
			}
			functional_prediction_final <- functional_prediction_final[keep,]
			write.table(x = functional_prediction_final, file = file.path(path_to_temp_folder, 'functional_prediction.txt'), append = F, quote = F, 
				sep = "\t", row.names = F, col.names = T)
			self$res_tax4fun2_KO <- functional_prediction_final
			message('Result KO abundance is stored in object$res_tax4fun2_KO.')

			# Converting the KO profile to a profile of KEGG pathways
			message('Converting functions to pathways.')
			ko2ptw <- Tax4Fun2_KEGG$ko2ptw
			if(normalize_pathways){
				functional_prediction_norm <- functional_prediction / ko_list$pathway_count
			}
			pathway_prediction <- aggregate(x = functional_prediction[ko2ptw$nrow,], by = list(ko2ptw$ptw), sum)
			if(ncol(pathway_prediction) >= 3){
				col_sums <- colSums(pathway_prediction[,-1])
				col_sums[col_sums == 0] <- 1
				pathway_prediction[,-1] <- t(t(pathway_prediction[,-1]) / col_sums)
				keep <- which(rowSums(pathway_prediction[,-1]) > 0)
			} else {
				pathway_prediction[,-1] <- t(t(pathway_prediction[,-1]) / sum(pathway_prediction[,-1]))
				keep <- which(pathway_prediction[,2] > 0)
			}
			if(sum(pathway_prediction[, -1]) == 0) stop("Conversion to pathway failed!")
			pathway_prediction %<>% .[keep, ]
			names(pathway_prediction) <- names(otu_table)
			rownames(pathway_prediction) <- pathway_prediction[, 1]
			pathway_prediction <- pathway_prediction[, -1, drop = FALSE]
			
			self$res_tax4fun2_pathway <- pathway_prediction
			message('Pathway abundance table is stored in object$res_tax4fun2_pathway.')			
			ptw_desc <- Tax4Fun2_KEGG$ptw_desc
			pathway_prediction_final <- data.frame(pathway_prediction, ptw_desc[rownames(pathway_prediction), ])
			pathway_prediction_final <- data.frame(pathway = rownames(pathway_prediction_final), pathway_prediction_final)
			write.table(x = pathway_prediction_final, file = file.path(path_to_temp_folder, 'pathway_prediction.txt'), 
				append = F, quote = F, sep = "\t", row.names = F, col.names = T)
		},
		#' @description
		#' Calculate (multi-) functional redundancy index (FRI) of prokaryotic community with Tax4Fun2 method.
		#' This function is used to calculating aFRI and rFRI use the intermediate files generated by the function cal_tax4fun2().
		#' please also cite: 
		#' Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences. 
		#' 	 Environmental Microbiome 15, 11 (2020). <doi:10.1186/s40793-020-00358-7>
		#'
		#' @return res_tax4fun2_aFRI and res_tax4fun2_rFRI in object.
		#' @examples
		#' \dontrun{
		#' t1$cal_tax4fun2_FRI()
		#' }
		cal_tax4fun2_FRI = function(){
			if(is.null(self$res_tax4fun2_reference_profile) | is.null(self$res_tax4fun2_database_mode) | is.null(self$res_tax4fun2_path_to_ref_dir)){
				stop("Please first run cal_tax4fun2 !")
			}else{
				reference_profile <- self$res_tax4fun2_reference_profile
				database_mode <- self$res_tax4fun2_database_mode
				path_to_ref_dir <- self$res_tax4fun2_path_to_ref_dir
			}

			path_to_ref_tree <- file.path(path_to_ref_dir, paste0(database_mode, ".tre"))
			otu_table_reduced_aggregated <- self$res_tax4fun2_otu_table_reduced_aggregated
			# Reading reference tree
			reference_tree <- ape::read.tree(path_to_ref_tree)
			distance_matrix <- cophenetic(reference_tree)

			rownumber <- NULL
			for(reference_id in otu_table_reduced_aggregated$Group.1){
				rownumber <- c(rownumber, which(row.names(distance_matrix) == reference_id))
			}
			distance_matrix_reduced <- distance_matrix[rownumber, rownumber]
			distance_matrix_mean <- mean(as.dist(distance_matrix))
			distance_matrix_reduced_mean <- mean(as.dist(distance_matrix_reduced))
			rm(distance_matrix)
			
			data("Tax4Fun2_KEGG", envir=environment())
			ko_list <- Tax4Fun2_KEGG$ko_list

			# Calculate functional redundancy sample-wise
			abs_functional_redundancy_tab <- NULL
			rel_functional_redundancy_tab <- NULL
			for(sample_use in 2:ncol(otu_table_reduced_aggregated)){
				print(paste('Calculate functional redundancy for sample', names(otu_table_reduced_aggregated[sample_use])))
				functional_prediction_sample <- reference_profile * as.numeric(otu_table_reduced_aggregated[, sample_use])
				functional_prediction_sample_mod <- ifelse(functional_prediction_sample >= 1, 1, 0)

				abs_functional_redundancy_sample <- NULL
				rel_functional_redundancy_sample <- NULL
				for(i in 1:nrow(ko_list)){
				  ko_count <- functional_prediction_sample_mod[, i]
				  aFRI <- (mean(as.dist(distance_matrix_reduced * ko_count)) * (sum(ko_count) / length(ko_count))) / distance_matrix_mean
				  rFRI <- (mean(as.dist(distance_matrix_reduced * ko_count)) * (sum(ko_count) / length(ko_count))) / distance_matrix_reduced_mean
				  abs_functional_redundancy_sample <- c(abs_functional_redundancy_sample, aFRI)
				  rel_functional_redundancy_sample <- c(rel_functional_redundancy_sample, rFRI)
				}
				abs_functional_redundancy_tab <- cbind(abs_functional_redundancy_tab, abs_functional_redundancy_sample)
				rel_functional_redundancy_tab <- cbind(rel_functional_redundancy_tab, rel_functional_redundancy_sample)
			}

			abs_functional_redundancy_tab %<>% as.data.frame
			rel_functional_redundancy_tab %<>% as.data.frame

			colnames(abs_functional_redundancy_tab) <- colnames(rel_functional_redundancy_tab) <- colnames(otu_table_reduced_aggregated)[2:ncol(otu_table_reduced_aggregated)]

			abs_functional_redundancy_final <- data.frame(KO = ko_list$ko, abs_functional_redundancy_tab, description = ko_list$description)
			rel_functional_redundancy_final <- data.frame(KO = ko_list$ko, rel_functional_redundancy_tab, description = ko_list$description)

			self$res_tax4fun2_aFRI <- abs_functional_redundancy_final
			message('Absolute functional redundancy is stored in object$res_tax4fun2_aFRI')
			self$res_tax4fun2_rFRI <- rel_functional_redundancy_final
			message('Relative functional redundancy is stored in object$res_tax4fun2_rFRI')
		},
		#' @description
		#' Print the trans_func object.
		print = function(){
			cat("trans_func object:\n")
			cat(paste("Functional analysis for", self$for_what, ".\n"))
			if(!is.null(self$sample_table)){
				cat("sample_table is available.\n")
			}
			if(!is.null(self$otu_table)){
				cat("otu_table is available.\n")
			}
			if(!is.null(self$tax_table)){
				cat("tax_table is available.\n")
			}
			if(!is.null(self$rep_fasta)){
				cat("rep_fasta is available.\n")
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
			FAPROTAX = list(
				"Energy source" = c("aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy", "photoautotrophy", "photoheterotrophy"),
				"C-cycle" = c("cellulolysis", "xylanolysis", "chitinolysis", "ligninolysis", "fermentation",  "methanogenesis", "methanotrophy", "methylotrophy", 
					"hydrocarbon_degradation", "oil_bioremediation"),
				"N-cycle" = c("nitrogen_fixation", "nitrification", "aerobic_ammonia_oxidation", "aerobic_nitrite_oxidation", 
					"nitrate_reduction", "nitrate_respiration","nitrite_respiration", "nitrogen_respiration",
					"denitrification", "nitrite_denitrification", "nitrate_denitrification", "nitrous_oxide_denitrification", "ureolysis"),
				"S-cycle" = c("sulfate_respiration", "sulfur_respiration", "sulfite_respiration","dark_sulfide_oxidation", "respiration_of_sulfur_compounds",
					"thiosulfate_respiration", "dark_oxidation_of_sulfur_compounds")
			)
		),
		default_fungi_func_group = list(
			FUNGuild = list(
				"Trophic Mode" = c("Pathotroph", "Saprotroph", "Symbiotroph"),
				"Guild" = c("Bryophyte Parasite", "Dung Saprotroph", "Ectomycorrhizal", "Fungal Parasite", "Leaf Saprotroph", "Plant Parasite",
					"Wood Saprotroph", "Animal Pathogen", "Endophyte", "Plant Pathogen", "Lichen Parasite", "Litter Saprotroph", "Soil Saprotroph",
					"Plant Saprotroph", "Epiphyte", "Lichenized", "Arbuscular Mycorrhizal", "Endomycorrhizal", "Ericoid Mycorrhizal", "Orchid Mycorrhizal", 
					"Root Associated Biotroph", "Clavicipitaceous Endophyte", "Animal Endosymbiont")
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
				"Decay_substrate" = c("algal_material", "animal_material", "pollen", "dung", "protist_material", "soil", 
					"sugar-rich_substrates", "fungal_material", "leaf/fruit/seed", 
					"burnt_material", "hydrocarbon-rich_substrates_(fuel)", "roots", "wood", "resins"),
				"Decay_type" = c("blue-staining", "brown_rot", "chitinolytic", "keratinolytic", "mold", "other", "soft_rot", "white_rot"),
				"Aquatic_habitat" = c("aquatic", "freshwater", "marine", "non-aquatic", "partly_aquatic", "partly_freshwater_(partly_non-aquatic)", 
					"partly_marine_(partly_non-aquatic)"),
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
		),
		blastTableReducer = function(path_to_blast_file = '') {
			if(file.size(path_to_blast_file) == 0) stop('Blast result file empty!')
			id1 = ""
			file_in = file(description = path_to_blast_file, open = "r")
			file_out = file(description = paste0(path_to_blast_file, ".tmp"), open = "w")
			while (TRUE) 
			{
			  line = readLines(con = file_in, n = 1)
			  if (length(line) == 0) break
			  id2 = strsplit(x = line, split = "\t", fixed = T)[[1]][1]
			  if(id1 != id2)
			  {
				id1 = id2
				write(x = line, file = file_out, append = T)
			  }
			}
			close(file_in)
			close(file_out)
			# Remove old file and rename tmp file
			file.remove(path_to_blast_file)
			file.rename(from = paste0(path_to_blast_file, ".tmp"), to = path_to_blast_file)
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)
