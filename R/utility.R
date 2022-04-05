
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
	taxonomy_table[, i] %<>% gsub(paste0(pattern, collapse = "|"), replacement, ., ignore.case = ignore.case)
	# delete the blank space in beginning and end
	taxonomy_table[, i] %<>% gsub("^\\s+|\\s+$", "", .)
	# delete any " in the text
	taxonomy_table[, i] %<>% gsub('"', "", ., fixed = TRUE)
	# some data have single underline, so first double, then single
	taxonomy_table[, i] %<>% gsub("^.*__", "", .) %>% gsub("^._", "", .)
	# check the missing data
	taxonomy_table[, i][is.na(taxonomy_table[, i])] <- na_fill
	# paste the final result with double underlines
	taxonomy_table[, i] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__", taxonomy_table[, i])
	taxonomy_table[, i]
}

# inner function
summarySE_inter = function(usedata=NULL, measurevar, groupvars=NULL, na.rm=TRUE) {
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






