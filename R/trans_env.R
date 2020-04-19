
#' Create trans_env object for the analysis of the effects of physicochemical indexes on communities.
#'
#' This class is a wrapper for a series of operations associated with physicochemical measurements.
#' The functions in this class include \code{\link{cal_rda}}, \code{\link{trans_rda}}, \code{\link{plot_rda}}, \code{\link{cal_mantel}},
#' \code{\link{cal_cor}}, \code{\link{plot_corr}}
#'
#' @param dataset the object of \code{\link{microtable}} Class.
#' @param env_cols default NULL; a vector; If the physicochemical data is in sample_table, use this parameter to indicate the columns.
#' @param add_data default NULL; provide the physicochemical data frame individually.
#' @param complete_na default FALSE; Whether fill the NA in the physicochemical data.
#' @return env_data and dataset in trans_env object.
#' @examples
#' t1 <- trans_env$new(dataset = dataset, add_data = env_data)
#' @export
trans_env <- R6Class(classname = "trans_env",
	public = list(
		initialize = function(dataset = NULL, env_cols = NULL, add_data = NULL, complete_na = FALSE
			){
			if(is.null(add_data)){
				env_data <- dataset$sample_table[, env_cols, drop = FALSE]
			}else{
				env_data <- add_data[rownames(add_data) %in% rownames(dataset$sample_table), ]
			}
			if(!is.null(dataset)){
				dataset1 <- clone(dataset)
				dataset1$sample_table %<>% base::subset(rownames(.) %in% rownames(env_data))
				dataset1$tidy_dataset(main_data = FALSE)
				env_data %<>% .[rownames(dataset1$sample_table), ]
			}
			if(complete_na == T){
				env_data[env_data == ""] <- NA
				env_data <- dropallfactors(env_data, unfac2num = TRUE)
				env_data[] <- lapply(env_data, function(x){if(is.character(x)) as.factor(x) else x})
				env_data %<>% mice::mice(print = FALSE) %>% mice::complete(., 1)
			}
			self$env_data <- env_data
			self$dataset <- dataset1
			},
		cal_rda = function(use_dbrda = TRUE, add_matrix = NULL, use_measure = NULL, feature_sel = FALSE, taxa_level = NULL, taxa_filter_thres = NULL){
			env_data <- self$env_data
			if(use_dbrda == T){
				if(is.null(self$dataset$beta_diversity) & is.null(add_matrix)){
					stop("No distance matrix provided; please use set add_matrix parameter")
				}
				if(!is.null(self$dataset$beta_diversity)){
					if(!is.null(use_measure)){
						use_matrix <- self$dataset$beta_diversity[[use_measure]]
					}else{
						use_matrix <- self$dataset$beta_diversity[[1]]
					}
				}else{
					use_matrix <- add_matrix
				}
				use_data <- use_matrix[rownames(env_data), rownames(env_data)] %>% as.dist
			}else{
				if(is.null(self$dataset)){
					stop("No abundance dataset provided; please set dataset parameter in creating Class")
				}
				if(is.null(taxa_level)){
					cat("No taxa_level provided, Genus used automatically!")
					taxa_level <- "Genus"
				}
				newdat <- self$dataset$merge_taxa(taxa_level)
				use_abund <- newdat$otu_table
				if(!is.null(taxa_filter_thres)){
					use_abund <- use_abund[apply(use_abund, 1, sum)/sum(use_abund) > taxa_filter_thres, ]
				}
				use_data <- as.data.frame(t(use_abund))
			}
			if(feature_sel == T){
				if(use_dbrda == T){
					mod0 <- dbrda(use_data ~ 1, env_data)
					mod1 <- dbrda(use_data ~ ., env_data)
				}else{
					mod0 <- rda(use_data ~ 1, env_data)
					mod1 <- rda(use_data ~ ., env_data)					
				}
				forward_res <- ordiR2step(mod0, scope = formula(mod1), direction="forward", perm.max = 999)
				res_sign <- gsub("+ ", "", rownames(data.frame(forward_res$anova)), fixed = TRUE)
				if(length(res_sign) == 0){
					stop("Non variables obtained after selection according to model. Check method and data!")
				}
				res_sign <- res_sign[1:(length(res_sign) - 1)]
				env_data <- env_data[,c(res_sign),drop=FALSE]
			}
			self$use_dbrda <- use_dbrda
			self$taxa_level <- taxa_level
			if(use_dbrda == T){
				self$res_rda <- dbrda(use_data ~ ., env_data)
			}else{
				self$res_rda <- rda(use_data ~ ., env_data)
			}
		},
		trans_rda = function(show_taxa = 10, adjust_arrow_length = FALSE, min_perc_env = 1, max_perc_env = 100, min_perc_tax = 1, max_perc_tax = 100){
			res_rda <- self$res_rda
			scrs <- scores(res_rda ,choices = c(1, 2), display = c("sp", "wa", "cn"))
			scrs$biplot <- scores(res_rda, choices=c(1, 2), "bp", scaling="sites")
			df_sites <- cbind.data.frame(scrs$sites, self$dataset$sample_table[rownames(scrs$sites), ])
			colnames(df_sites)[1:2] <- c("x","y")
			multiplier <- vegan:::ordiArrowMul(scrs$biplot)
			df_arrows<- scrs$biplot * multiplier
			colnames(df_arrows)<-c("x","y")
			df_arrows=as.data.frame(df_arrows)
			eigval <- res_rda$CCA$eig/sum(res_rda$CCA$eig)
			eigval <- round(100 * eigval, 1)
			eigval[1] <- paste0("RDA1", " [", eigval[1], "%]")
			eigval[2] <- paste0("RDA2", " [", eigval[2], "%]")

			if(self$use_dbrda == F){
				scrs$biplot_spe <- scores(res_rda, choices=c(1, 2), "sp", scaling="species")
				df_species <- scrs$species
				colnames(df_species)[1:2] <- c("x","y")
				multiplier_spe <- vegan:::ordiArrowMul(scrs$biplot_spe)
				df_arrows_spe <- scrs$biplot_spe * multiplier_spe
				colnames(df_arrows_spe)<-c("x","y")
				df_arrows_spe = dropallfactors(cbind.data.frame(df_arrows_spe, self$dataset$tax_table[rownames(df_arrows_spe), self$taxa_level, drop = FALSE]))
				df_arrows_spe %<>% .[!grepl("__$|__uncultured|sp$", .[, 3]), ]
				df_arrows_spe <- df_arrows_spe %>% {.[,1]^2 + .[,2]^2} %>% `names<-`(rownames(df_arrows_spe)) %>% 
					sort(., decreasing = TRUE) %>% .[1:10] %>% names %>% df_arrows_spe[., ]

			}else{
				df_species = NULL
				df_arrows_spe = NULL
			}
			if(adjust_arrow_length == T){
				df_arrows[,1:2] <- private$stand_fun(df_arrows[,1:2], min_perc = min_perc_env, max_perc = max_perc_env)
				if(self$use_dbrda == F){
					df_arrows_spe[,1:2] <- private$stand_fun(df_arrows_spe[,1:2], min_perc = min_perc_tax, max_perc = max_perc_tax)
				}
			}
			
			self$res_rda_trans = list(df_sites = df_sites, df_arrows = df_arrows, eigval = eigval, df_species = df_species, df_arrows_spe = df_arrows_spe)
		},
		plot_rda = function(plot_color = NULL, plot_shape = NULL, color_values = RColorBrewer::brewer.pal(8, "Dark2"), taxa_text_color = "firebrick1", 
			taxa_text_type = "italic"){
			p <- ggplot()
			p <- p + theme_bw()
			p <- p + theme(panel.grid=element_blank())
			p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey80")
			p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "grey80")
			p <- p + geom_point(data=self$res_rda_trans$df_sites,aes_string("x", "y", colour = plot_color, shape = plot_shape),size= 3.5)
			# plot arrows
			p <- p + geom_segment(data=self$res_rda_trans$df_arrows, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(0.2, "cm")), color = "grey30")
			p <- p + ggrepel::geom_text_repel(data=as.data.frame(self$res_rda_trans$df_arrows*1), 
				aes(x, y, label = gsub("`", "", rownames(self$res_rda_trans$df_arrows))),size=3.7, color = "black", segment.color = "white")
			if(!is.null(plot_color)){
				p <- p + scale_color_manual(values = color_values)
			}
			p <- p + xlab(self$res_rda_trans$eigval[1]) + ylab(self$res_rda_trans$eigval[2])
			
			if(self$use_dbrda == F){
				df_arrows_spe1 <- self$res_rda_trans$df_arrows_spe
				p <- p + geom_segment(data=df_arrows_spe1, aes(x = 0, y = 0, xend = x, yend = y), 
					arrow = arrow(length = unit(0.2, "cm")), color = "firebrick1", alpha = .6)
				df_arrows_spe1[, self$taxa_level] %<>% gsub(".*__", "", .) %>% gsub("Candidatus ", "", .) 
				if(taxa_text_type == "italic"){
					df_arrows_spe1[, self$taxa_level] %<>%  paste0("italic('", .,"')")
				}
				p <- p + ggrepel::geom_text_repel(data=df_arrows_spe1, aes_string("x", "y", label = self$taxa_level), size=3, color = taxa_text_color, segment.alpha = .01, parse = TRUE)
			}
			p
		},
		cal_mantel = function(select_env_data = NULL, partial_mantel = FALSE, add_matrix = NULL, use_measure = NULL, method = "pearson", ...){
			if(is.null(self$dataset$beta_diversity) & is.null(add_matrix)){
				stop("No distance matrix provided; please use set add_matrix parameter or use provide dataset in creating Class")
			}
			if(is.null(add_matrix)){
				if(!is.null(use_measure)){
					use_matrix <- self$dataset$beta_diversity[[use_measure]]
				}else{
					use_matrix <- self$dataset$beta_diversity[[1]]
				}
			}else{
				use_matrix <- add_matrix
			}
			env_data <- self$env_data
			# if no selection, automatically select the numeric columns
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric))]
			}else{
				env_data <- env_data[, select_env_data]
			}
			veg.dist <- as.dist(use_matrix[rownames(env_data), rownames(env_data)])
			variable_name <- c()
			corr_res <- c()
			p_res <- c()

			for(i in 1:ncol(env_data)){
				env.dist <- vegdist(scale(env_data[, i, drop=FALSE]), "euclid")
				if(partial_mantel == T){
					zdis <- vegdist(scale(env_data[, -i, drop=FALSE]), "euclid")
					man1 <- mantel.partial(veg.dist, env.dist, zdis)
				}else{
					man1 <- mantel(veg.dist, env.dist)
				}
				variable_name <- c(variable_name, colnames(env_data)[i])
				corr_res <- c(corr_res, man1$statistic)
				p_res <- c(p_res, man1$signif)
			}
			cor_method <- rep(method, length(p_res))
			significance <- cut(p_res, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
			res_mantel <- data.frame(variable_name, cor_method, corr_res, p_res, significance)
			if(partial_mantel == T){
				self$res_mantel_partial <- res_mantel
			}else{
				self$res_mantel <- res_mantel
			}
		},
		cal_cor = function(use_data = c("Genus", "all", "other")[1], select_env_data = NULL, cor_method = c("pearson", "spearman", "kendall")[1],
			p_adjust_type = c("Type", "Taxa", "Env")[3], p_adjust_method = "fdr", add_abund_table = NULL, by_group = NULL,
			other_taxa = NULL, group_use = NULL, group_select = NULL,
			taxa_name_full = TRUE
			){
			env_data <- self$env_data
			if(is.null(select_env_data)){
				env_data <- env_data[, unlist(lapply(env_data, is.numeric))]
			}
			if(!is.null(add_abund_table)){
				abund_table <- add_abund_table
			}else{
				if(use_data %in% taxonomic_ranks){
					abund_table <- self$dataset$taxa_abund[[use_data]]
				}
				if(grepl("all|other", use_data, ignore.case = TRUE)){
					abund_table <- do.call(rbind, unname(self$dataset$taxa_abund))
					if(use_data == "other"){
						if(is.null(other_taxa)){
							stop("You select other, but no other_taxa provided!")
						}
						abund_table <- abund_table[other_taxa, ]
					}
				}
				abund_table %<>% .[!grepl("__$|__uncultured$", rownames(.)), ]
				abund_table <- as.data.frame(t(abund_table))
			}
			# filter samples by one group
			if(!is.null(group_use)){
				if(is.null(group_select)){
					stop("You select group_use parameter, but no group_select parameter provided!")
				}
				sel_sample_names <- self$dataset$sample_table %>% .[.[, group_use] %in% group_select, ] %>% rownames
				abund_table <- abund_table[sel_sample_names, ]
			}
			env_data %<>% .[rownames(.) %in% rownames(abund_table), ]
			abund_table <- abund_table[rownames(env_data), ]
			if(is.null(by_group)){
				groups <- rep("All", nrow(env_data))
			}else{
				groups <- self$dataset$sample_table[, by_group] %>% as.character
				message("Calculate the corr by the groups in ", by_group, " of sample_table, respectively")
			}
			comb_names <- expand.grid(unique(groups), colnames(abund_table), colnames(env_data)) %>% t %>% as.data.frame(stringsAsFactors = FALSE)
			df1 <- sapply(comb_names, function(x){
				suppressWarnings(cor.test(abund_table[groups == x[1], x[2]], env_data[groups == x[1], x[3]], method = cor_method)) %>%
				{c(x, Correlation = unname(.$estimate), Pvalue = unname(.$p.value))}
			})
			df1 %<>% t %>% as.data.frame(stringsAsFactors = FALSE)
			colnames(df1) <- c("Type", "Taxa", "Env", "Correlation","Pvalue")
			df1$Pvalue %<>% as.numeric
			df1$Correlation %<>% as.numeric
			df1$AdjPvalue <- rep(0, nrow(df1))
			choose_col <- which(c("Type", "Taxa", "Env") %in% p_adjust_type)
			comb_names2 <- comb_names[choose_col, ] %>% t %>% as.data.frame %>% unique %>% t %>% as.data.frame(stringsAsFactors = FALSE)
			invisible(lapply(comb_names2, function(x){
				row_sel <- unlist(lapply(as.data.frame(t(df1[, choose_col, drop = FALSE])), function(y) all(y %in% x)));
				df1$AdjPvalue[row_sel] <<- p.adjust(df1[row_sel, "Pvalue"], method = p_adjust_method)
			}))
			df1$Significance <- cut(df1$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
			df1<-df1[complete.cases(df1), ]
			df1$Env <- factor(df1$Env, levels = unique(as.character(df1$Env)))
			if(taxa_name_full == F){
				df1$Taxa %<>% gsub(".*__(.*?)$", "\\1", .)
			}
			self$res_cor <- df1
			self$cor_method <- cor_method
		},
		plot_corr = function(pheatmap = FALSE, ylab_type_italic = FALSE, 
			keep_full_name = FALSE, keep_prefix = TRUE,
			color_vector = c("#00008B", "#102D9B", "#215AAC", "#3288BD", "#66C2A5",  "#E6F598", "#FFFFBF", "#FED690", "#FDAE61", "#F46D43", "#D53E4F"),
			plot_x_size = 9, mylabels_x = NULL, font_family = NULL
			){
			use_data <- self$res_cor
			if(keep_full_name == F){
				use_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(keep_prefix == F){
				use_data$Taxa %<>% gsub(".*__", "", .)
			}
			clu_data_1 <- reshape2::dcast(use_data, Taxa~Env, value.var = "Correlation") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
			sig_data <- reshape2::dcast(use_data, Taxa~Env, value.var = "Significance") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
			if(pheatmap == T){
				if(!require(pheatmap)){
					stop("pheatmap package not installed")
				}
				if(ylab_type_italic == T){
					eval(parse(text = paste0("mylabels_y <- c(", paste0("expression(italic(", paste0('"', rownames(clu_data_1),'"'), "))", collapse = ","),")", collapse = "")))
				}else{
					mylabels_y <- rownames(clu_data_1)
				}
				color_vector_use <- colorRampPalette(color_vector)(100)
				p <- pheatmap(clu_data_1, clustering_distance_row = "correlation", clustering_distance_cols= "correlation",border_color = NA, scale = "none",
						 fontsize=plot_x_size,cluster_cols=TRUE,cluster_rows=TRUE, labels_row = mylabels_y, labels_col = mylabels_x, 
						 display_numbers = sig_data, number_color = "black",  fontsize_number = plot_x_size*1.2,
						 color = color_vector_use)
				p$gtable
			}else{
				lim_y <- hclust(dist(clu_data_1)) %>% {.$labels[.$order]}
				lim_x <- hclust(dist(t(clu_data_1))) %>% {.$labels[.$order]}

				p <- ggplot(aes(x=Env, y=Taxa, fill=Correlation), data = use_data) +
					theme_bw() + 
					geom_tile() + 
					scale_fill_gradientn(colours = color_vector) +
					scale_y_discrete(limits=lim_y, position = "left") + 
					scale_x_discrete(limits=lim_x) +
					geom_text(aes(label=Significance), color="black", size=4) + 
					labs(y = NULL, x = "Measure", fill = self$cor_method) +
					theme(axis.text.x = element_text(angle = 40, colour = "black", vjust = 1, hjust = 1, size = 10)) +
					theme(strip.background = element_rect(fill = "grey85", colour = "white")) +
					theme(strip.text=element_text(size=11), panel.border = element_blank(), panel.grid = element_blank())
				if(ylab_type_italic == T){
					p <- p + theme(axis.text.y = element_text(face = 'italic'))
				}
				if(!is.null(font_family)){
					p <- p + theme(text = element_text(family = font_family))
				}
				p
			}
		},
		print = function(){
			cat("trans_env class:\n")
			if(!is.null(self$env_data)){
				cat(paste0("Env table have ", ncol(self$env_data), " variables: ", paste0(colnames(self$env_data), collapse = ",")))
				cat("\n")
			}
		}
	),
	private = list(
		stand_fun = function(x, min_perc = 1, max_perc = 10) {
			# x must be a two column data.frame or matrix
			t1 <- x[,1]^2 + x[,2]^2
			a <- min_perc * max(t1)
			b <- max_perc * max(t1)
			Ymax <- max(t1)
			Ymin <- min(t1)
			k <- (b-a)/(Ymax-Ymin) 
			norY <- a + k*(t1-Ymin)
			per <- abs(x[,1]/x[,2])
			newx <- (((norY*per^2) / (per^2 + 1)) ^ (1/2)) * sapply(x[, 1], function(y) ifelse(y > 0, 1, -1))
			newy <- (newx/per) * sapply(x[, 2], function(y) ifelse(y > 0, 1, -1))
			res <- data.frame(newx, newy)
			colnames(res) <- colnames(x)
			res
		}
	),
	lock_class = FALSE,
	lock_objects = FALSE
)



#' calculate RDA.
#'
#' @param use_dbrda default TRUE; whether use db-RDA, if FALSE, use RDA.
#' @param add_matrix default NULL; additional distance matrix provided, if you donot want to use the beta diversity matrix in the dataset.
#' @param use_measure default NULL; name of beta diversity matrix. If necessary and not provided, use the first beta diversity matrix.
#' @param feature_sel default FALSE; whether perform the feature selection.
#' @param taxa_level default NULL; If use RDA, provide the taxonomic rank.
#' @param taxa_filter_thres default NULL; If want to filter taxa, provide the relative abundance threshold.
#' @return res_rda in object.
#' @examples
#' t1$cal_rda(use_dbrda = TRUE, use_measure = "bray")
cal_rda <- function(use_dbrda = TRUE, add_matrix = NULL, use_measure = NULL, feature_sel = FALSE, taxa_level = NULL, taxa_filter_thres = NULL){
	dataset$cal_rda()
}




#' transform RDA result for the following plotting.
#'
#' @param show_taxa default 10; taxa number shown in the plot.
#' @param adjust_arrow_length default FALSE; whether adjust the arrow length to be clear
#' @param min_perc_env default 1; minimum value for env arrow, relatively.
#' @param max_perc_env default 100; maximum value for env arrow, relatively.
#' @param min_perc_tax default 1; minimum value for tax arrow, relatively.
#' @param max_perc_tax default 100; maximum value for tax arrow, relatively.
#' @return res_rda_trans in object.
#' @examples
#' t1$trans_rda(adjust_arrow_length = TRUE, max_perc_env = 10)
trans_rda <- function(show_taxa = 10, adjust_arrow_length = FALSE, min_perc_env = 1, max_perc_env = 100, min_perc_tax = 1, max_perc_tax = 100){
	dataset$trans_rda()
}


#' plot RDA result.
#'
#' @param plot_color default NULL; group used for color.
#' @param plot_shape default NULL; group used for shape.
#' @param color_values color pallete.
#' @param taxa_text_color default "firebrick1"; taxa text colors.
#' @param taxa_text_type default "italic"; taxa text style; better to use "italic" for Genus, use "normal" for others.
#' @return ggplot object.
#' @examples
#' t1$plot_rda(plot_color = "Group")
plot_rda <- function(plot_color = NULL, plot_shape = NULL, color_values = RColorBrewer::brewer.pal(8, "Dark2"),
				taxa_text_color = "firebrick1", taxa_text_type = "italic"){
	dataset$plot_rda()
}


#' Mantel test between beta diversity matrix and physicochemical data.
#'
#' @param select_env_data default NULL; numeric or character vector to select columns in env_data; if not provided, automatically select the columns with numeric attributes.
#' @param partial_mantel default FALSE; whether use partial mantel test.
#' @param add_matrix default NULL; additional distance matrix provided, if you donot want to use the beta diversity matrix in the dataset.
#' @param use_measure default NULL; name of beta diversity matrix. If necessary and not provided, use the first beta diversity matrix.
#' @param method default "pearson"; one of c("pearson", "spearman", "kendall"); correlation method.
#' @return res_mantel in object.
#' @examples
#' t1$cal_mantel(use_measure = "bray")
cal_mantel <- function(select_env_data = NULL, partial_mantel = FALSE, add_matrix = NULL, use_measure = NULL, method = "pearson", ...){
	dataset$cal_mantel()
}

#' Correlations between abundance data table and physicochemical data table.
#'
#' Calculating any correlation between two columns from two tables.
#'
#' @param use_data default "Genus"; one of c("Genus", "all", "other"); Genus: genus abundance, all: all taxa abundance, other: provide additional data with other_taxa parameter.
#' @param select_env_data default NULL; numeric or character vector to select columns in env_data; if not provided, automatically select the columns with numeric attributes.
#' @param method default "pearson"; one of c("pearson", "spearman", "kendall"); correlation method.
#' @param p_adjust_method default "fdr"; p.adjust method.
#' @param p_adjust_type default "Env"; one of c("Type", "Taxa", "Env"); p.adjust type; Env: physicochemical data; Taxa: taxa data; Type: group used.
#' @param add_abund_table default NULL; additional data table to be used. Samples must be rows.
#' @param by_group default NULL; numeric or character vector to select one column in env_data; correlations separately for groups.
#' @param other_taxa default NULL; provide additional taxa, see use_data parameter.
#' @param group_use default NULL; numeric or character vector to select one column in env_data for selecting samples; together with group_select.
#' @param group_select default NULL; the group name used; will retain samples within the group.
#' @param taxa_name_full default TRUE; Whether retain the complete taxonomic name of taxa.
#' @return res_cor in object.
#' @examples
#' t2 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", rf_taxa_level = "Genus")
#' t1 <- trans_env$new(dataset = dataset, add_data = env_data[, 4:11])
#' t1$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t2$res_rf$Taxa[1:40])

cal_cor <- function(use_data = c("Genus", "all", "other")[1], select_env_data = NULL, cor_method = c("pearson", "spearman", "kendall")[1],
			p_adjust_method = "fdr", p_adjust_type = c("Env", "Taxa", "Type")[1], add_abund_table = NULL, by_group = NULL,
			other_taxa = NULL, group_use = NULL, group_select = NULL,
			taxa_name_full = TRUE){
	dataset$cal_cor()
}


#' Plot correlations heatmap.
#'
#'
#' @param color_vector color pallete.
#' @param pheatmap default FALSE; whether use heatmap with clustering plot.
#' @param ylab_type_italic default FALSE; whether use italic type for y lab text.
#' @param keep_full_name default FALSE; whether use the complete taxonomic name.
#' @param keep_prefix default TRUE; whether retain the taxonomic prefix.
#' @param plot_x_size default 9; x axis text size.
#' @param mylabels_x default NULL; provide x axis text labels additionally; only available when pheatmap = TRUE.
#' @param font_family default NULL; font family used in ggplot2; only available when pheatmap = FALSE.
#' @return plot.
#' @examples
#' t1$plot_corr()

plot_corr <- function(color_vector = c("#00008B", "#102D9B", "#215AAC", "#3288BD", "#66C2A5",  "#E6F598", "#FFFFBF", "#FED690", "#FDAE61", "#F46D43", "#D53E4F"),
			pheatmap = FALSE, ylab_type_italic = FALSE, 
			keep_full_name = FALSE, keep_prefix = TRUE,
			plot_x_size = 9, mylabels_x = NULL, font_family = NULL){
	dataset$plot_corr()
}












