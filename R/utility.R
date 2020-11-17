


#' Transform microtable object in microeco package to the phyloseq object in phyloseq package.
#'
#' @param dataset a microtable object.
#' @return phyloseq object.
#' @examples
#' \donttest{
#' data("dataset")
#' meco2phyloseq(dataset)
#' }
#' @export
meco2phyloseq <- function(dataset){
	otu_table_trans <- dataset$otu_table
	sample_table_trans <- dataset$sample_table
	tax_table_trans <- dataset$tax_table
	phylo_tree_trans <- dataset$phylo_tree
	# OTU table for phyloseq
	otu_table_trans <- phyloseq::otu_table(otu_table_trans, taxa_are_rows = TRUE)
	sampledata <- phyloseq::sample_data(sample_table_trans)
	rownames(sampledata) <- rownames(sample_table_trans)
	# Taxonomic table for phyloseq
	if(!is.null(tax_table_trans)){
		tax_table_trans <- phyloseq::tax_table(as.matrix(tax_table_trans))
	}
	physeq <- phyloseq::phyloseq(otu_table_trans, tax_table_trans, sampledata, phylo_tree_trans)
	physeq
}


#' Transform the phyloseq object in phyloseq package to microtable object in microeco package.
#'
#' @param physeq a phyloseq object.
#' @return microtable object.
#' @examples
#' \donttest{
#' library(phyloseq)
#' data("GlobalPatterns")
#' phyloseq2meco(GlobalPatterns)
#' }
#' @export
phyloseq2meco <- function(physeq){
	if(physeq@otu_table@taxa_are_rows){
		otu_table_trans <- as.data.frame(physeq@otu_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
	}else{
		otu_table_trans <- as.data.frame(t(physeq@otu_table@.Data), check.names = FALSE, stringsAsFactors = FALSE)
	}
	sample_table_trans <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE, stringsAsFactors = FALSE)
	tax_table_trans <- as.data.frame(physeq@tax_table@.Data, check.names = FALSE, stringsAsFactors = FALSE)
	tax_table_trans %<>% tidy_taxonomy
	phylo_tree_trans <- physeq@phy_tree
	
	dataset <- microtable$new(sample_table = sample_table_trans, otu_table = otu_table_trans, tax_table = tax_table_trans, phylo_tree = phylo_tree_trans)
	dataset
}



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
#' @return data frame without factor
#' @examples
#' data("taxonomy_table_16S")
#' taxonomy_table_16S[, 1] <- as.factor(taxonomy_table_16S[, 1])
#' dropallfactors(taxonomy_table_16S)
#' @export
dropallfactors <- function(x, unfac2num = FALSE){
	# check x class
	if(!is.data.frame(x)){
		stop("input data must be data.frame class")
	}
	trycharnum <- function(x){
		if(suppressWarnings(sum(is.na(as.numeric(as.character(x)))) != sum(is.na(x)))) {
			x <- as.character(x)
		} else {
			x <- as.numeric(as.character(x))
		}
		x
	}
	if(unfac2num == T){
		x[] <- lapply(x, function(x) trycharnum(x))
	}else{
		x[] <- lapply(x, function(x) if(is.factor(x)) trycharnum(x) else x)
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
	taxonomy_table[,i] <- gsub(".*metagenome.*|.*uncultur.*|.*cultivar.*|D_6__synthetic.*|.*archaeon$", "", taxonomy_table[,i], ignore.case = TRUE)
	taxonomy_table[,i] <- gsub('"', "", taxonomy_table[,i], fixed = TRUE)
	taxonomy_table[,i] <- gsub("^\\s+|\\s+$|.*\\sbacterium$|.*bacterium\\s.*", "", taxonomy_table[,i])
	taxonomy_table[,i] <- gsub("^.*__", paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__"), taxonomy_table[,i])
	taxonomy_table[,i][is.na(taxonomy_table[,i])] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__")
	taxonomy_table[,i][grepl("^$",taxonomy_table[,i])] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__")
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



##################################################################################
##################################################################################

# metastat code from White et al. (2009) <doi:10.1371/journal.pcbi.1000352>.
#************************************************************************
# ************************** SUBROUTINES ********************************
#************************************************************************

#*****************************************************************************************************
#  calc two sample two statistics
#  g is the first column in the matrix representing the second condition
#*****************************************************************************************************
calc_twosample_ts <- function(Pmatrix, g, nrows, ncols)
{
	C1 <- array(0, dim=c(nrows,3));  # statistic profiles
	C2 <- array(0, dim=c(nrows,3)); 
	Ts <- array(0, dim=c(nrows,1));

	if (nrows == 1){
		C1[1,1] = mean(Pmatrix[1:g-1]);
		C1[1,2] = stats::var(Pmatrix[1:g-1]); # variance
		C1[1,3] = C1[1,2]/(g-1);    # std err^2

		C2[1,1] = mean(Pmatrix[g:ncols]);
		C2[1,2] = stats::var(Pmatrix[g:ncols]);  # variance
		C2[1,3] = C2[1,2]/(ncols-g+1); # std err^2
	}else{
		# generate statistic profiles for both groups
		# mean, var, stderr
		for (i in 1:nrows){ # for each taxa
			# find the mean of each group
			C1[i,1] = mean(Pmatrix[i, 1:g-1]);  
			C1[i,2] = stats::var(Pmatrix[i, 1:g-1]); # variance
			C1[i,3] = C1[i,2]/(g-1);    # std err^2

			C2[i,1] = mean(Pmatrix[i, g:ncols]);  
			C2[i,2] = stats::var(Pmatrix[i, g:ncols]);  # variance
			C2[i,3] = C2[i,2]/(ncols-g+1); # std err^2
		}
	}

	# permutation based t-statistics
	for (i in 1:nrows){ # for each taxa
		xbar_diff = C1[i,1] - C2[i,1]; 
		denom = sqrt(C1[i,3] + C2[i,3]);
		Ts[i] = xbar_diff/denom;  # calculate two sample t-statistic 
	}

	return (Ts);

}


#*****************************************************************************************************
#  function to calculate qvalues.
#  takes an unordered set of pvalues corresponding the rows of the matrix
#*****************************************************************************************************
calc_qvalues <- function(pvalues)
{
	nrows = length(pvalues);

	# create lambda vector
	lambdas <- seq(0,0.95,0.01);
	pi0_hat <- array(0, dim=c(length(lambdas)));

	# calculate pi0_hat
	for (l in 1:length(lambdas)){ # for each lambda value
		count = 0;
		for (i in 1:nrows){ # for each p-value in order
			if (pvalues[i] > lambdas[l]){
				count = count + 1; 	
			}
			pi0_hat[l] = count/(nrows*(1-lambdas[l]));
		}
	}

	f <- unclass(stats::smooth.spline(lambdas,pi0_hat,df=3));
	f_spline <- f$y;
	pi0 = f_spline[length(lambdas)];   # this is the essential pi0_hat value

	# order p-values
	ordered_ps <- order(pvalues);
	pvalues <- pvalues;
	qvalues <- array(0, dim=c(nrows));
	ordered_qs <- array(0, dim=c(nrows));

	ordered_qs[nrows] <- min(pvalues[ordered_ps[nrows]]*pi0, 1);
	for(i in (nrows-1):1) {
		p = pvalues[ordered_ps[i]];
		new = p*nrows*pi0/i;

		ordered_qs[i] <- min(new,ordered_qs[i+1],1);
	}

	# re-distribute calculated qvalues to appropriate rows
	for (i in 1:nrows){
		qvalues[ordered_ps[i]] = ordered_qs[i];
	}

	################################
	# plotting pi_hat vs. lambda
	################################
	# plot(lambdas,pi0_hat,xlab=expression(lambda),ylab=expression(hat(pi)[0](lambda)),type="p");
	# lines(f);

	return (qvalues);
}


#*****************************************************************************************************
#  function to calculate permuted pvalues from Storey and Tibshirani(2003)
#  B is the number of permutation cycles
#  g is the first column in the matrix of the second condition 
#*****************************************************************************************************
permuted_pvalues <- function(Imatrix, tstats, B, g, Fmatrix)
{
	# B is the number of permutations were going to use!
	# g is the first column of the second sample
	# matrix stores tstats for each taxa(row) for each permuted trial(column)

	M = nrow(Imatrix);
	ps <- array(0, dim=c(M)); # to store the pvalues
	if (is.null(M) || M == 0){
		return (ps);
	}
	permuted_ttests <- array(0, dim=c(M, B));
	ncols = ncol(Fmatrix);
	# calculate null version of tstats using B permutations.
	for (j in 1:B){
		trial_ts <- permute_and_calc_ts(Imatrix, sample(1:ncol(Imatrix)), g);
		permuted_ttests[,j] <- abs(trial_ts); 
	}

	# calculate each pvalue using the null ts
	if ((g-1) < 8 || (ncols-g+1) < 8){
		# then pool the t's together!
		# count how many high freq taxa there are
		hfc = 0;
		for (i in 1:M){                   # for each taxa
			if (sum(Fmatrix[i,1:(g-1)]) >= (g-1) || sum(Fmatrix[i,g:ncols]) >= (ncols-g+1)){
				hfc = hfc + 1;
			}
		}
		# the array pooling just the frequently observed ts  
		cleanedpermuted_ttests <- array(0, dim=c(hfc,B));
		hfc = 1;
		for (i in 1:M){
			if (sum(Fmatrix[i,1:(g-1)]) >= (g-1) || sum(Fmatrix[i,g:ncols]) >= (ncols-g+1)){
				cleanedpermuted_ttests[hfc,] = permuted_ttests[i,];
				hfc = hfc + 1;
			}
		}
		#now for each taxa
		for (i in 1:M){  
			ps[i] = (1/(B*hfc))*sum(cleanedpermuted_ttests > abs(tstats[i]));
		}
	}else{
		for (i in 1:M){
			ps[i] = (1/(B+1))*(sum(permuted_ttests[i,] > abs(tstats[i]))+1);
		}
	}

	return (ps);
}


#*****************************************************************************************************
# takes a matrix, a permutation vector, and a group division g.
# returns a set of ts based on the permutation.
#*****************************************************************************************************
permute_and_calc_ts <- function(Imatrix, y, g)
{
	nr = nrow(Imatrix);
	nc = ncol(Imatrix);
	# first permute the rows in the matrix
	Pmatrix <- Imatrix[,y[1:length(y)]];
	Ts <- calc_twosample_ts(Pmatrix, g, nr, nc);

	return (Ts);
}



# http://metastats.cbcb.umd.edu/detect_DA_features.r
#*****************************************************************************************************
#*****************************************************************************************************
#  Last modified: 4/14/2009 
#  
#  Author: james robert white, whitej@umd.edu, Center for Bioinformatics and Computational Biology.
#  University of Maryland - College Park, MD 20740
#
#  This software is designed to identify differentially abundant features between two groups
#  Input is a matrix of frequency data. Several thresholding options are available.
#  See documentation for details.
#*****************************************************************************************************
#*****************************************************************************************************

#*****************************************************************************************************
#  detect_differentially_abundant_features:
#  the major function - inputs an R object "jobj" containing a list of feature names and the 
#  corresponding frequency matrix, the argument g is the first column of the second group. 
#  
#  -> set the pflag to be TRUE or FALSE to threshold by p or q values, respectively
#  -> threshold is the significance level to reject hypotheses by.
#  -> B is the number of bootstrapping permutations to use in estimating the null t-stat distribution.
#*****************************************************************************************************
#*****************************************************************************************************


detect_differentially_abundant_features <- function(jobj, g, pflag = NULL, threshold = NULL, B = NULL){
	#**********************************************************************************
	# ************************ INITIALIZE COMMAND-LINE ********************************
	# ************************        PARAMETERS       ********************************
	#**********************************************************************************
	qflag = FALSE;
	if (is.null(B)){
		B = 1000;
	}
	if (is.null(threshold)){
		threshold = 0.05;
	}
	if (is.null(pflag)){
		pflag = TRUE;
		qflag = FALSE;
	}
	if (pflag == TRUE){
		qflag = FALSE;
	}
	if (pflag == FALSE){
		qflag = TRUE;
	}

	#********************************************************************************
	# ************************ INITIALIZE PARAMETERS ********************************
	#********************************************************************************

	#*************************************
	Fmatrix <- jobj$matrix;                   # the feature abundance matrix
	taxa <- jobj$taxa;                        # the taxa/(feature) labels of the TAM
	nrows = nrow(Fmatrix);                   
	ncols = ncol(Fmatrix);
	Pmatrix <- array(0, dim=c(nrows,ncols));  # the relative proportion matrix
	C1 <- array(0, dim=c(nrows,3));           # statistic profiles for class1 and class 2
	C2 <- array(0, dim=c(nrows,3));           # mean[1], variance[2], standard error[3]   
	T_statistics <- array(0, dim=c(nrows,1)); # a place to store the true t-statistics 
	pvalues <- array(0, dim=c(nrows,1));      # place to store pvalues
	qvalues <- array(0, dim=c(nrows,1));      # stores qvalues
	#*************************************

	#*************************************
	#  convert to proportions
	#  generate Pmatrix
	#*************************************
	totals <- array(0, dim=c(ncol(Fmatrix)));
	for (i in 1:ncol(Fmatrix)) {
		# sum the ith column 
		totals[i] = sum(Fmatrix[,i]);
	}

	for (i in 1:ncols) {   # for each subject
		for (j in 1:nrows) { # for each row
			Pmatrix[j,i] = Fmatrix[j,i]/totals[i];
		}
	}

	#********************************************************************************
	# ************************** STATISTICAL TESTING ********************************
	#********************************************************************************

	if (ncols == 2){  # then we have a two sample comparison
		#************************************************************
		#  generate p values using chisquared or fisher's exact test
		#************************************************************
		for (i in 1:nrows){           # for each feature
			f11 = sum(Fmatrix[i,1]);
			f12 = sum(Fmatrix[i,2]);
			f21 = totals[1] - f11;
			f22 = totals[2] - f12;
			C1[i,1] = f11/totals[1];                       # proportion estimate
			C1[i,2] = (C1[i,1]*(1-C1[i,1]))/(totals[1]-1); # sample variance
			C1[i,3] = sqrt(C1[i,2]);                       # sample standard error
			C2[i,1] = f12/totals[2];
			C2[i,2] = (C2[i,1]*(1-C2[i,1]))/(totals[2]-1);
			C2[i,3] = sqrt(C2[i,2]); 

			#  f11  f12
			#  f21  f22  <- contigency table format
			contingencytable <- array(0, dim=c(2,2));
			contingencytable[1,1] = f11;
			contingencytable[1,2] = f12;
			contingencytable[2,1] = f21;
			contingencytable[2,2] = f22;

			if (f11 > 20 && f22 > 20){
				csqt <- stats::chisq.test(contingencytable);
				pvalues[i] = csqt$p.value;
			}else{
				ft <- stats::fisher.test(contingencytable, workspace = 8e6, alternative = "two.sided", conf.int = FALSE);
				pvalues[i] = ft$p.value;
			}

		}

		#*************************************
		#  calculate q values from p values
		#*************************************
		qvalues <- calc_qvalues(pvalues);

	}else{ # we have multiple subjects per population

		#*************************************
		#  generate statistics mean, var, stderr    
		#*************************************
		for (i in 1:nrows){ # for each taxa
			# find the mean of each group
			C1[i,1] = mean(Pmatrix[i, 1:g-1]);  
			C1[i,2] = stats::var(Pmatrix[i, 1:g-1]); # variance
			C1[i,3] = C1[i,2]/(g-1);    # std err^2 (will change to std err at end)

			C2[i,1] = mean(Pmatrix[i, g:ncols]);  
			C2[i,2] = stats::var(Pmatrix[i, g:ncols]);  # variance
			C2[i,3] = C2[i,2]/(ncols-g+1); # std err^2 (will change to std err at end)
		}

		#*************************************
		#  two sample t-statistics
		#*************************************
		for (i in 1:nrows){                   # for each taxa
			xbar_diff = C1[i,1] - C2[i,1]; 
			denom = sqrt(C1[i,3] + C2[i,3]);
			T_statistics[i] = xbar_diff/denom;  # calculate two sample t-statistic
		}

		#*************************************
		# generate initial permuted p-values
		#*************************************
		pvalues <- permuted_pvalues(Pmatrix, T_statistics, B, g, Fmatrix);

		#*************************************
		#  generate p values for sparse data 
		#  using fisher's exact test
		#*************************************
		for (i in 1:nrows){                   # for each taxa
			if (sum(Fmatrix[i,1:(g-1)]) < (g-1) && sum(Fmatrix[i,g:ncols]) < (ncols-g+1)){
				# then this is a candidate for fisher's exact test
				f11 = sum(Fmatrix[i,1:(g-1)]);
				f12 = sum(Fmatrix[i,g:ncols]);
				f21 = sum(totals[1:(g-1)]) - f11;
				f22 = sum(totals[g:ncols]) - f12;
				#  f11  f12
				#  f21  f22  <- contigency table format
				contingencytable <- array(0, dim=c(2,2));
				contingencytable[1,1] = f11;
				contingencytable[1,2] = f12;
				contingencytable[2,1] = f21;
				contingencytable[2,2] = f22;
				ft <- fisher.test(contingencytable, workspace = 8e6, alternative = "two.sided", conf.int = FALSE);
				pvalues[i] = ft$p.value; 
			}  
		}

		#*************************************
		#  calculate q values from p values
		#*************************************
		qvalues <- calc_qvalues(pvalues);

		#*************************************
		#  convert stderr^2 to std error
		#*************************************
		for (i in 1:nrows){
			C1[i,3] = sqrt(C1[i,3]);
			C2[i,3] = sqrt(C2[i,3]);
		}
	}

	#*************************************
	#  threshold sigvalues and print
	#*************************************
	sigvalues <- array(0, dim=c(nrows,1));
	if (pflag == TRUE){  # which are you thresholding by?
		sigvalues <- pvalues;
	}else{
		sigvalues <- qvalues;
	}
	s = sum(sigvalues <= threshold);
	Differential_matrix <- array(0, dim=c(s,9));

	dex = 1;
	for (i in 1:nrows){
		if (sigvalues[i] <= threshold){
			Differential_matrix[dex,1]   = jobj$taxa[i];
			Differential_matrix[dex,2:4] = C1[i,];
			Differential_matrix[dex,5:7] = C2[i,];
			Differential_matrix[dex,8]   = pvalues[i];  
			Differential_matrix[dex,9]   = qvalues[i];
			dex = dex+1;  
		}
	}

	# show(Differential_matrix);

	Total_matrix <- array(0, dim=c(nrows,9));
	for (i in 1:nrows){
		Total_matrix[i,1]   = jobj$taxa[i];
		Total_matrix[i,2:4] = C1[i,];
		Total_matrix[i,5:7] = C2[i,];
		Total_matrix[i,8]   = pvalues[i];
		Total_matrix[i,9]   = qvalues[i];
	}

	colnames(Total_matrix) <- c("taxa", "mean(group1)", "variance(group1)", "std.err(group1)", "mean(group2)", "variance(group2)", "std.err(group2)", "pvalue", "qvalue")
	Total_matrix <- Total_matrix[order(Total_matrix[, "pvalue"], decreasing = FALSE), ]
	Total_matrix
}

#*****************************************************************************************************
#  load up the frequency matrix from a file
#*****************************************************************************************************
# Note sep
load_frequency_matrix <- function(input){
	# load names 
	subjects <- array(0,dim=c(ncol(input)));
	for(i in 1:length(subjects)) {
		subjects[i] <- as.character(colnames(input)[i]);
	}
	# load taxa
	taxa <- array(0,dim=c(nrow(input)));
	for(i in 1:length(taxa)) {
		taxa[i] <- as.character(rownames(input)[i]);
	}

	# load remaining counts
	matrix <- array(0, dim=c(length(taxa),length(subjects)));
	for(i in 1:length(taxa)){
		for(j in 1:length(subjects)){ 
			matrix[i,j] <- as.numeric(input[i,j]);
		}
	}
	jobj <- list(matrix=matrix, taxa=taxa)
	return(jobj)
}



calculate_metastat <- function(inputdata, g, pflag = FALSE, threshold = NULL, B = NULL){
	trans_data <- load_frequency_matrix(input = inputdata)
	res <- detect_differentially_abundant_features(jobj = trans_data, g = g, pflag = pflag, threshold = threshold, B = B)
	res
}

