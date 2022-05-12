
#' Copy an R6 class object completely
#'
#' @param x R6 class object
#' @param deep default TRUE; deep copy
#' @return identical but unrelated R6 object.
#' @examples
#' data("dataset")
#' clone(dataset)
#' @export
clone <- function(x, deep = TRUE){
	y <- x$clone(deep = deep)
	y
}

#' Remove all factors in a data frame
#'
#' @param x data frame
#' @param unfac2num default FALSE; whether try to convert all character to numeric; if FALSE, only try to convert column with factor attribute.
#'   Note that this can only transform the columns that may be transformed to numeric without using factor.
#' @param char2num default FALSE; whether force all the character to be numeric class by using factor as an intermediate.
#' @return data frame without factor
#' @examples
#' data("taxonomy_table_16S")
#' taxonomy_table_16S[, 1] <- as.factor(taxonomy_table_16S[, 1])
#' str(dropallfactors(taxonomy_table_16S))
#' @export
dropallfactors <- function(x, unfac2num = FALSE, char2num = FALSE){
	# check x class
	if(!is.data.frame(x)){
		stop("input data must be data.frame class")
	}
	if(unfac2num == T){
		x[] <- lapply(x, function(x) trycharnum(x))
	}else{
		x[] <- lapply(x, function(x) if(is.factor(x)) trycharnum(x) else x)
	}
	if(char2num == T){
		x[] <- lapply(x, function(x) if(is.character(x)) as.factor(x) else x)
		x[] <- lapply(x, function(x) as.numeric(x))
	}
	x
}

# inner function
trycharnum <- function(x){
	if(suppressWarnings(sum(is.na(as.numeric(as.character(x)))) != sum(is.na(x)))) {
		x <- as.character(x)
	} else {
		x <- as.numeric(as.character(x))
	}
	x
}

#' Clean up the taxonomic table to make taxonomic assignments consistent.
#'
#' @param taxonomy_table a data.frame with taxonomic information.
#' @param column default "all"; "all" or a number; 'all' represents cleaning up all the columns; a number represents cleaning up this column.
#' @param pattern default see the function parameter; the characters (regular expression) to be cleaned up or replaced; cleaned up when parameter replacement = "", 
#'   replaced when parameter replacement has something; Note that the capital and small letters are not distinguished.
#' @param replacement default ""; the characters used to replace the character in pattern parameter.
#' @param ignore.case default TRUE; if FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.
#' @param na_fill default ""; used to replace the NA.
#' @return taxonomic table.
#' @format \code{\link{data.frame}} object.
#' @examples
#' data("taxonomy_table_16S")
#' tidy_taxonomy(taxonomy_table_16S)
#' @export
tidy_taxonomy <- function(taxonomy_table, 
	column = "all",
	pattern = c(".*uncultur.*", ".*unknown.*", ".*unidentif.*", ".*unclassified.*", ".*No blast hit.*", ".*sp\\.$",
		".*metagenome.*", ".*cultivar.*", ".*archaeon$", "__synthetic.*", ".*\\sbacterium$", ".*bacterium\\s.*"),
	replacement = "",
	ignore.case = TRUE,
	na_fill = ""
	){
	if(column == "all"){
		taxonomy_table[] <- lapply(seq_len(ncol(taxonomy_table)), 
			function(x) tidy_taxonomy_column(taxonomy_table, i = x, pattern = pattern, replacement = replacement, ignore.case = ignore.case, na_fill = na_fill))
	}else{
		if(!inherits(column, "numeric")){
			stop("The input column is not numeric class !")
		}
		taxonomy_table[, column] <- tidy_taxonomy_column(taxonomy_table, i = column, pattern = pattern, 
			replacement = replacement, ignore.case = ignore.case, na_fill = na_fill)
	}
	taxonomy_table
}

# inner function
tidy_taxonomy_column <- function(taxonomy_table, i, pattern, replacement, ignore.case, na_fill){
	taxonomy_table[, i] <- gsub(paste0(pattern, collapse = "|"), replacement, taxonomy_table[, i], ignore.case = ignore.case)
	# delete the blank space in beginning and end
	taxonomy_table[, i] <- gsub("^\\s+|\\s+$", "", taxonomy_table[, i])
	# delete any " in the text
	taxonomy_table[, i] <- gsub('"', "", taxonomy_table[, i], fixed = TRUE)
	# some data have single underline, so first double, then single
	taxonomy_table[, i] <- gsub("^.*__", "", taxonomy_table[, i])
	taxonomy_table[, i] <- gsub("^._", "", taxonomy_table[, i])
	# check the missing data
	taxonomy_table[, i][is.na(taxonomy_table[, i])] <- na_fill
	# paste the final result with double underlines
	taxonomy_table[, i] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__", taxonomy_table[, i])
	taxonomy_table[, i]
}

# inner function
summarySE_inter = function(usedata = NULL, measurevar, groupvars = NULL, na.rm = TRUE) {
	length2 <- function(x, na.rm=TRUE) ifelse(na.rm, sum(!is.na(x)), length(x))
	datac <- usedata %>% 
			dplyr::grouped_df(groupvars) %>% 
			dplyr::summarise(N = length2(!!sym(measurevar), na.rm=na.rm), Mean = mean(!!sym(measurevar), na.rm=na.rm), SD = stats::sd(!!sym(measurevar), na.rm=na.rm)) %>%
			as.data.frame
	datac$SE <- datac$SD / sqrt(datac$N)
	datac
}

#' A color palette for 20 elements.
#'
#' This is one palette option for users who have <= 20 elements to plot.
#' The palletes of RColorBrewer package provide at most 12 discrete colors, such as "Set3" and "Paired".
#' This palette is adapted from D3.js library and can be used as a supplement.
#'
#' @export
color_palette_20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5", 
	"#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")


# inner function
# Add correlation or regression statictics to a scatter plot
# developed based on the stat_cor of ggpubr package
stat_corlm <- function(mapping = NULL, data = NULL, 
	type = c("cor", "lm")[1], cor_method = "pearson", label_sep = ";",
	pvalue_trim = 4, cor_coef_trim = 3, lm_fir_trim = 2, lm_sec_trim = 2, lm_squ_trim = 2,
	label.x.npc = "left", label.y.npc = "top", label.x = NULL, label.y = NULL,
	geom = "text", position = "identity", na.rm = FALSE, show.legend = NA,
	inherit.aes = TRUE, ...
	){
	layer(
		stat = StatCorLm, data = data, mapping = mapping, geom = geom,
		position = position, show.legend = show.legend, inherit.aes = inherit.aes,
		params = list(label.x.npc  = label.x.npc , label.y.npc  = label.y.npc, label.x = label.x, label.y = label.y,
			type = type, cor_method = cor_method, label_sep = label_sep,
			pvalue_trim = pvalue_trim, cor_coef_trim = cor_coef_trim, lm_fir_trim = lm_fir_trim, lm_sec_trim = lm_sec_trim, lm_squ_trim = lm_squ_trim,
			parse = TRUE, na.rm = na.rm, ...)
	)
}

StatCorLm <- ggproto("StatCorLm", Stat,
	required_aes = c("x", "y"),
	default_aes = aes(hjust = ..hjust.., vjust = ..vjust..),
	compute_group = function(data, scales, type, cor_method, label_sep,
		pvalue_trim, cor_coef_trim, lm_fir_trim, lm_sec_trim, lm_squ_trim,
		label.x.npc, label.y.npc, label.x, label.y
		){
		if(length(unique(data$x)) < 2){
			return(data.frame())
		}
		# Returns a data frame with estimate, p.value, label, method
		.test <- .corlm_test(
			data$x, data$y, type = type, cor_method = cor_method, label_sep = label_sep,
			pvalue_trim = pvalue_trim,
			cor_coef_trim = cor_coef_trim, 
			lm_fir_trim = lm_fir_trim, 
			lm_sec_trim = lm_sec_trim, 
			lm_squ_trim = lm_squ_trim
		)
		# Returns a data frame with label: x, y, hjust, vjust
		.label.pms <- ggpubr:::.label_params(data = data, scales = scales,
			label.x.npc = label.x.npc, label.y.npc = label.y.npc,
			label.x = label.x, label.y = label.y ) %>%
			dplyr::mutate(hjust = 0)
		cbind(.test, .label.pms)
	}
)

.corlm_test <- function(
	x, y, 
	type, cor_method = "pearson", label_sep = ";",
	pvalue_trim = 4, cor_coef_trim = 3, lm_fir_trim = 2, lm_sec_trim = 2, lm_squ_trim = 2
	){
	label_sep_use <- paste0("*`", label_sep, "`~")
	if(type == "cor"){
		fit <- suppressWarnings(stats::cor.test(x, y, method = cor_method, alternative = "two.sided"))
		pvalue <- fit$p.value
		estimate <- fit$estimate
		cor_var <- list(
			cor_coef = unname(round(estimate, digits = cor_coef_trim)), 
			cor_p = ifelse(pvalue < 0.0001, " < 0.0001", paste0(" = ", round(pvalue, digits = pvalue_trim)))
		)
		cor_coef_exp <- substitute(italic(R)~"="~cor_coef, cor_var)
		cor_p_exp <- substitute(~italic(P)*cor_p, cor_var)
		res <- paste0(c(cor_coef_exp, label_sep_use, cor_p_exp), collapse = "")
	}else{
		fit <- stats::lm(y ~ x)
		pvalue <- stats::anova(fit)$`Pr(>F)`[1]
		inte <- round(unname(stats::coef(fit))[1], digits = lm_sec_trim)
		lm_var <- list(
			lm_a = ifelse(inte < 0, paste0(" - ", abs(inte)), paste0(" + ", as.character(inte))),
			lm_b = round(unname(stats::coef(fit))[2], digits = lm_fir_trim),
			lm_r2 = round(summary(fit)$r.squared, digits = lm_squ_trim),
			lm_p = ifelse(pvalue < 0.0001, " < 0.0001", paste0(" = ", round(pvalue, digits = pvalue_trim)))
			)
		lm_ab_exp <- substitute(italic(y) == lm_b %.% italic(x)*lm_a, lm_var)
		lm_r2_exp <- substitute(~italic(R)^2~"="~lm_r2, lm_var)
		lm_p_exp <- substitute(~italic(P)*lm_p, lm_var)
		
		res <- paste0(c(lm_ab_exp, label_sep_use, lm_r2_exp, label_sep_use, lm_p_exp), collapse = "")
	}
	data.frame(label = as.character(as.expression(res)))
}



