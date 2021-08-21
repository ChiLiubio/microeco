
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

#' Clear up the taxonomic table to make taxonomic assignments consistent.
#'
#' @param taxonomy_table a data.frame with taxonomic information.
#' @return taxonomic table.
#' @format \code{\link{data.frame}} object.
#' @examples
#' data("taxonomy_table_16S")
#' tidy_taxonomy(taxonomy_table_16S)
#' @export
tidy_taxonomy <- function(taxonomy_table){
	taxonomy_table[] <- lapply(seq_len(ncol(taxonomy_table)), function(x) tidy_taxonomy_column(taxonomy_table, x))
	taxonomy_table
}

# inner function
tidy_taxonomy_column <- function(taxonomy_table, i){
	taxonomy_table[,i] <- gsub(".*No blast hit.*|.*Unknown.*|.*unidentif.*|.*sp\\.$|.*Unclassified.*", "", taxonomy_table[,i], ignore.case = TRUE)
	taxonomy_table[,i] <- gsub(".*metagenome.*|.*uncultur.*|.*cultivar.*|D_6__synthetic.*|.*archaeon$", "", taxonomy_table[,i], ignore.case = TRUE)
	taxonomy_table[,i] <- gsub('"', "", taxonomy_table[,i], fixed = TRUE)
	taxonomy_table[,i] <- gsub("^\\s+|\\s+$|.*\\sbacterium$|.*bacterium\\s.*", "", taxonomy_table[,i])
	taxonomy_table[,i] <- gsub("^.*__", "", taxonomy_table[,i])
	# some data have single underline
	taxonomy_table[,i] <- gsub("^._", "", taxonomy_table[,i])
	# check the missing data
	taxonomy_table[,i][is.na(taxonomy_table[,i])] <- ""
	# paste the final result with double underlines
	taxonomy_table[,i] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__", taxonomy_table[,i])
	taxonomy_table[,i]
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

