#' @title
#' Feature abundance normalization/transformation.
#'
#' @description
#' Feature abundance normalization/transformation for a microtable object or data.frame object.
#'
#' @export
trans_norm <- R6Class(classname = "trans_norm",
	public = list(
		#' @description
		#' Get a transposed abundance table if the input is microtable object. In the table, rows are samples, and columns are features. 
		#'   This can make the further operations same with the traditional ecological methods.
		#' @param dataset the \code{\link{microtable}} object or \code{data.frame} object. 
		#' 	 If it is \code{data.frame} object, please make sure that rows are samples, and columns are features.
		#' @return data_table, stored in the object. 
		#' @examples
		#' library(microeco)
		#' data(dataset)
		#' t1 <- trans_norm$new(dataset = dataset)
		initialize = function(dataset = NULL)
			{
			if(inherits(dataset, "microtable")){
				dataset$tidy_dataset()
				abund_table <- dataset$otu_table
				abund_table <- t(abund_table)
			}else{
				abund_table <- dataset
				if(!inherits(abund_table, "matrix")){
					abund_table <- as.matrix(abund_table)
				}
				if(ncol(abund_table) < nrow(abund_table)){
					message("Please make sure that the rows in input data.frame object are samples ...")
				}
			}
			message("Sample number: ", nrow(abund_table), "; Feature number: ", ncol(abund_table), " ...")
			message("Prepared abundance table is stored in object$data_table ...")
			self$data_table <- abund_table
			self$dataset <- dataset
		},
		#' @description
		#' Normalization/transformation methods.
		#' @param method default "rarefy"; See the following available options. \cr 
		#' \cr 
		#' Methods for normalization:
		#' \itemize{
		#'   \item \code{"rarefy"}: classic rarefaction based on R sample function.
		#'   \item \code{"SRS"}: scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <doi:10.7717/peerj.9593>.
		#'   \item \code{"clr"}: Centered log-ratio normalization <ISBN:978-0-412-28060-3> <doi: 10.3389/fmicb.2017.02224>. 
		#' 	   	 It is defined:  \deqn{clr_{ki} = \log\frac{x_{ki}}{g(x_i)}}
		#' 	   	 where \eqn{x_{ki}} is the abundance of \eqn{k}th feature in sample \eqn{i}, \eqn{g(x_i)} is the geometric mean of abundances for sample \eqn{i}.
		#' 	   	 A pseudocount need to be added to deal with the zero. For more information, please see the 'clr' method in \code{decostand} function of vegan package.
		#'   \item \code{"rclr"}: Robust centered log-ratio normalization <doi:10.1128/msystems.00016-19>.
		#' 	   	 It is defined:  \deqn{rclr_{ki} = \log\frac{x_{ki}}{g(x_i > 0)}}
		#' 	   	 where \eqn{x_{ki}} is the abundance of \eqn{k}th feature in sample \eqn{i}, \eqn{g(x_i > 0)} is the geometric mean of abundances (> 0) for sample \eqn{i}.
		#' 	   	 In rclr, zero values are kept as zeroes, and not taken into account.
		#'   \item \code{"GMPR"}: Geometric mean of pairwise ratios <doi: 10.7717/peerj.4600>. 
		#' 	   	 For a given sample \eqn{i}, the size factor \eqn{s_i} is defined:
		#' 	   	     \deqn{s_i =  \biggl( {\displaystyle\prod_{j=1}^{n} Median_{k|c_{ki}c_{kj} \ne 0} \lbrace \dfrac{c_{ki}}{c_{kj}} \rbrace} \biggr) ^{1/n}}
		#' 	   	 where \eqn{k} denotes all the features, and \eqn{n} denotes all the samples. 
		#' 	   	 For sample \eqn{i}, \eqn{GMPR = \frac{x_{i}}{s_i}}, where \eqn{x_i} is the feature abundances of sample \eqn{i}.
		#'   \item \code{"CSS"}: Cumulative sum scaling normalization based on the \code{metagenomeSeq} package <doi:10.1038/nmeth.2658>.
		#' 	   	 For a given sample \eqn{j}, the scaling factor \eqn{s_{j}^{l}} is defined:
		#' 	   	     \deqn{s_{j}^{l} = {\displaystyle\sum_{i|c_{ij} \leqslant q_{j}^{l}} c_{ij}}}
		#' 	   	 where \eqn{q_{j}^{l}} is the \eqn{l}th quantile of sample \eqn{j}, that is, in sample \eqn{j} there are \eqn{l} features with counts smaller than \eqn{q_{j}^{l}}.
		#' 	   	 \eqn{c_{ij}} denotes the count (abundance) of feature i in sample \eqn{j}.
		#' 	   	 For \eqn{l} = 0.95\eqn{m} (feature number), \eqn{q_{j}^{l}} corresponds to the 95th percentile of the count distribution for sample \eqn{j}.
		#' 	   	 Normalized counts \eqn{\tilde{c_{ij}} = (\frac{c_{ij}}{s_{j}^{l}})(N)}, where \eqn{N} is an appropriately chosen normalization constant.
		#'   \item \code{"TSS"}: Total sum scaling. Abundance is divided by the sequencing depth.
		#' 	   	 For a given sample \eqn{j}, normalized counts is defined:
		#' 	   	     \deqn{\tilde{c_{ij}} = \frac{c_{ij}}{\sum_{i=1}^{N_{j}} c_{ij}}}
		#' 	   	 where \eqn{c_{ij}} is the counts of feature \eqn{i} in sample \eqn{j}, and \eqn{N_{j}} is the feature number of sample \eqn{j}.
		#'   \item \code{"eBay"}: Empirical Bayes approach to normalization <10.1186/s12859-020-03552-z>. 
		#' 	   	 The implemented method is not tree-related. In the output, the sum of each sample is 1.
		#'   \item \code{"TMM"}: Trimmed mean of M-values method based on the \code{normLibSizes} function of \code{edgeR} package <doi: 10.1186/gb-2010-11-3-r25>.
		#'   \item \code{"DESeq2"}: Median ratio of gene counts relative to geometric mean per gene based on the DESeq function of \code{DESeq2} package <doi: 10.1186/s13059-014-0550-8>.
		#' 	   	 This option can invoke the \code{trans_diff} class and extract the normalized data from the original result.
		#' 	   	 Note that either \code{group} or \code{formula} should be provided.
		#' 	   	 The scaling factor is defined:
		#' 	   	     \deqn{s_{j} = Median_{i} \frac{c_{ij}}{\bigl( {\prod_{j=1}^{n} c_{ij}} \bigr) ^{1/n}}}
		#' 	   	 where \eqn{c_{ij}} is the counts of feature \eqn{i} in sample \eqn{j}, and \eqn{n} is the total sample number.
		#'   \item \code{"Wrench"}: Group-wise and sample-wise compositional bias factor <doi: 10.1186/s12864-018-5160-5>.
		#' 	   	 Note that condition parameter is necesary to be passed to \code{condition} parameter in \code{wrench} function of Wrench package.
		#' 	   	 As the input data must be microtable object, so the input condition parameter can be a column name of \code{sample_table}.
		#' 	   	 The scaling factor is defined:
		#' 	   	     \deqn{s_{j} = \frac{1}{p} \sum_{ij} W_{ij} \frac{X_{ij}}{\overline{X_{i}}}}
		#' 	   	 where \eqn{X_{ij}} represents the relative abundance (proportion) for feature \eqn{i} in sample \eqn{j},
		#' 	   	 \eqn{\overline{X_{i}}} is the average proportion of feature \eqn{i} across the dataset,
		#' 	   	 \eqn{W_{ij}} represents a weight specific to each technique, and \eqn{p} is the feature number in sample.
		#'   \item \code{"RLE"}: Relative log expression. 
		#' }
		#' Methods based on \code{decostand} function of vegan package:
		#' \itemize{
		#'   \item \code{"total"}: divide by margin total (default MARGIN = 1, i.e. rows - samples).
		#'   \item \code{"max"}: divide by margin maximum (default MARGIN = 2, i.e. columns - features).
		#'   \item \code{"normalize"}:  make margin sum of squares equal to one (default MARGIN = 1).
		#'   \item \code{"range"}: standardize values into range 0...1 (default MARGIN = 2). If all values are constant, they will be transformed to 0.
		#'   \item \code{"standardize"}: scale x to zero mean and unit variance (default MARGIN = 2).
		#'   \item \code{"pa"}: scale x to presence/absence scale (0/1).
		#'   \item \code{"log"}: logarithmic transformation.
		#' }
		#' Other methods for transformation:
		#' \itemize{
		#'   \item \code{"AST"}: Arc sine square root transformation.
		#' }
		#' @param sample.size default NULL; libray size for rarefaction when method = "rarefy" or "SRS". If not provided, use the minimum number across all samples. 
		#'    For "SRS" method, this parameter is passed to \code{Cmin} parameter of \code{SRS} function of SRS package.
		#' @param rngseed default 123; random seed. Available when method = "rarefy" or "SRS".
		#' @param replace default TRUE; see \code{\link{sample}} for the random sampling; Available when \code{method = "rarefy"}.
		#' @param pseudocount default 1; add pseudocount for those features with 0 abundance when \code{method = "clr"}.
		#' @param intersect.no default 10; the intersecting taxa number between paired sample for \code{method = "GMPR"}.
		#' @param ct.min default 1; the minimum number of counts required to calculate ratios for \code{method = "GMPR"}.
		#' @param condition default NULL; Only available when \code{method = "Wrench"}. 
		#'    This parameter is passed to the \code{condition} parameter of \code{wrench} function in Wrench package
		#'    It must be a column name of \code{sample_table} or a vector with same length of samples.
		#' @param MARGIN default NULL; 1 = samples, and 2 = features of abundance table; only available when method comes from \code{decostand} function of vegan package.
		#'    If MARGIN is NULL, use the default value in decostand function.
		#' @param logbase default 2; The logarithm base.
		#' @param ... parameters pass to \code{vegan::decostand}, or \code{metagenomeSeq::cumNorm} when method = "CSS", 
		#'    or \code{edgeR::normLibSizes} when method = "TMM" or "RLE", 
		#'    or \code{trans_diff} class when method = "DESeq2",
		#'    or \code{wrench} function of Wrench package when method = "Wrench".
		#' 
		#' @return new microtable object or data.frame object.
		#' @examples
		#' newdataset <- t1$norm(method = "clr")
		#' newdataset <- t1$norm(method = "log")
		norm = function(method = "rarefy", sample.size = NULL, rngseed = 123, replace = TRUE, pseudocount = 1, intersect.no = 10, 
			ct.min = 1, condition = NULL, MARGIN = NULL, logbase = 2, ...)
			{
			abund_table <- self$data_table
			if(is.null(method)){
				stop("Please select a method!")
			}
			method <- tolower(method)
			method <- match.arg(method, c("rarefy", "srs", "gmpr", "clr", "rclr", "css", "tss", "ebay", "tmm", "deseq2", "rle", "wrench", "ast", 
				"total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log"))
			
			if(method %in% c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")){
				if(is.null(MARGIN)){
					MARGIN <- switch(method, total = 1, max = 2, frequency = 2, normalize = 1, range = 2, rank = 1, standardize = 2, chi.square = 1, NULL)
				}
				res_table <- vegan::decostand(x = abund_table, method = method, MARGIN = MARGIN, logbase = logbase, ...)
			}
			if(method %in% c("rarefy", "srs")){
				newotu <- as.data.frame(t(abund_table))
				suppressMessages(tmpobj <- microtable$new(newotu))
				set.seed(rngseed)

				if(is.null(sample.size)){
					sample.size <- min(tmpobj$sample_sums())
					message("Use the minimum number across samples: ", sample.size)
				}
				if(length(sample.size) > 1){
					stop("Input sample.size had more than one value!")
				}
				if(sample.size <= 0){
					stop("sample.size less than or equal to zero. Need positive sample size to work!")
				}
				if (max(tmpobj$sample_sums()) < sample.size){
					stop("sample.size is larger than the maximum of sample sums, pleasure check input sample.size!")
				}
				if (min(tmpobj$sample_sums()) < sample.size) {
					rmsamples <- tmpobj$sample_names()[tmpobj$sample_sums() < sample.size]
					message(length(rmsamples), " samples are removed because of fewer reads than input sample.size ...")
					feature_num_raw <- length(tmpobj$taxa_names())
					tmpobj$sample_table <- base::subset(tmpobj$sample_table, ! tmpobj$sample_names() %in% rmsamples)
					tmpobj$tidy_dataset()
					feature_num_filter <- length(tmpobj$taxa_names())
					if(feature_num_filter < feature_num_raw){
						message((feature_num_raw - feature_num_filter), " features with 0 abundance are removed after filtering samples ...")
					}
				}
				newotu <- tmpobj$otu_table
				if(method == "rarefy"){
					newotu <- as.data.frame(apply(newotu, 2, private$rarefaction_subsample, sample.size = sample.size, replace = replace))
				}else{
					newotu <- SRS::SRS(newotu, Cmin = sample.size, set_seed = TRUE, seed = rngseed)
				}
				rownames(newotu) <- rownames(tmpobj$otu_table)
				tmpobj$otu_table <- newotu
				rmtaxa <- apply(newotu, 1, sum) %>% .[. == 0]
				if(length(rmtaxa) > 0){
					message(length(rmtaxa), " features are removed because they are no longer present in any sample after random subsampling ...")
					tmpobj$tidy_dataset()
				}
				res_table <- t(tmpobj$otu_table)
			}
			if(method == "deseq2"){
				if(!inherits(self$dataset, "microtable")){
					stop("For DESeq2 method, the input dataset must be microtable object when creating the object!")
				}
				use_dataset <- clone(self$dataset)
				use_dataset$tax_table <- data.frame(OTU = rownames(use_dataset$otu_table))
				rownames(use_dataset$tax_table) <- use_dataset$tax_table$OTU
				res_obj <- suppressMessages(trans_diff$new(dataset = use_dataset, method = "DESeq2", taxa_level = "OTU", ...))
				res_raw <- res_obj$res_diff_raw
				res_table <- t(counts(res_raw, normalized = TRUE))
			}
			if(method == "gmpr"){
				transposed_table <- t(abund_table)
				size_factor <- private$GMPR(transposed_table, intersect.no = intersect.no, ct.min = ct.min)
				message("Sample size factor is stored in object$size_factor ...")
				self$size_factor <- size_factor
				res_table <- t(transposed_table) / size_factor
			}
			if(method %in% c("clr", "rclr")){
				if(is.null(MARGIN)){
					MARGIN <- 1
				}
				if(method == "clr"){
					if(any(abund_table == 0)){
						abund_table <- abund_table + pseudocount
					}
				}
				res_table <- vegan::decostand(x = abund_table, method = method, MARGIN = MARGIN, logbase = logbase, ...)
			}
			if(method == "css"){
				obj <- metagenomeSeq::newMRexperiment(t(abund_table))
				obj_1 <- metagenomeSeq::cumNorm(obj, ...)
				res_table <- t(metagenomeSeq::MRcounts(obj_1, norm = TRUE))
			}
			if(method == "tss"){
				res_table <- apply(abund_table, 1, function(x){x/sum(x)}) %>% t
			}
			if(method == "ebay"){
				res_table <- private$ebay(abund_table)
			}
			if(method == "wrench"){
				if(!require("Wrench")){
					stop('Please first install Wrench package: BiocManager::install("Wrench")')
				}
				if(is.null(condition)){
					stop("Please provide the condition parameter, which is necessary for Wrench method!")
				}
				if(length(condition) == 1){
					if(!inherits(self$dataset, "microtable")){
						stop("Please provide a microtable object when input condition parameter is a column name!")
					}
					use_dataset <- clone(self$dataset)
					check_table_variable(use_dataset$sample_table, condition, "condition", "sample_table")
					condition <- use_dataset$sample_table[, condition]
				}else{
					if(length(condition) != nrow(abund_table)){
						stop("Provided condition must have a length same with samples!")
					}
				}				
				transposed_table <- t(abund_table)
				res_raw <- wrench(mat = transposed_table, condition = condition, ...)
				self$res_wrench_raw <- res_raw
				size_factor <- res_raw$nf
				res_table <- t(transposed_table) / size_factor
			}
			if(method %in% c("tmm", "rle")){
				abund_table %<>% t
				upmethod <- toupper(method)
				libsize <- edgeR::normLibSizes(abund_table, method = upmethod, ...)
				effec_libsize <- colSums(abund_table) * libsize
				ref_libsize <- mean(effec_libsize)
				res_table <- sweep(abund_table, MARGIN = 2, effec_libsize, "/") * ref_libsize
				res_table %<>% t
			}
			if(method == "ast"){
				res_table <- private$AST(abund_table)
			}
			if(inherits(self$dataset, "microtable")){
				res_dataset <- clone(self$dataset)
				res_dataset$otu_table <- as.data.frame(t(res_table))
				res_dataset$tidy_dataset()
				res_dataset
			}else{
				res_table
			}
		}
	),
	private = list(
		rarefaction_subsample = function(x, sample.size, replace=FALSE){
			# Adapted from the rarefy_even_depth in phyloseq package.
			rarvec <- numeric(length(x))
			if(sum(x) <= 0){
				# Protect against, and quickly return an empty vector, 
				return(rarvec)
			}
			if(replace){
				suppressWarnings(subsample <- sample(1:length(x), sample.size, replace = TRUE, prob = x))
			} else {
				# resample without replacement
				obsvec <- apply(data.frame(OTUi = 1:length(x), times = x), 1, function(x){
					rep_len(x["OTUi"], x["times"])
				})
				obsvec <- unlist(obsvec, use.names=FALSE)
				suppressWarnings(subsample <- sample(obsvec, sample.size, replace = FALSE))
			}
			sstab <- table(subsample)
			# Assign the tabulated random subsample values to the species vector
			rarvec[as(names(sstab), "integer")] <- sstab
			return(rarvec)
		},
		AST = function(x){
			sign(x) * asin(sqrt(abs(x)))
		},
		# modified from https://github.com/jchen1981/GMPR
		GMPR = function(comm, intersect.no, ct.min) {
			comm[comm < ct.min] <- 0	
			comm.no <- numeric(ncol(comm))
			output <- sapply(1:ncol(comm),  function(i) {
						if (i %% 100 == 0) {
							cat(i, '\n')
						}
						x <- comm[, i]
						# Compute the pairwise ratio
						pr <- x / comm
						# Handling of the NA, NaN, Inf
						pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
						# Counting the number of non-NA, NaN, Inf
						incl.no <- colSums(!is.na(pr))		
						# Calculate the median of PR
						pr.median <- matrixStats::colMedians(pr, na.rm=TRUE)
						# Record the number of samples used for calculating the GMPR
						comm.no[i] <- sum(incl.no >= intersect.no)
						# Geometric mean of PR median
						if (comm.no[i] > 1) {
							exp(mean(log(pr.median[incl.no >= intersect.no])))
						} else {
							NA
						}
					}
			)
			if (sum(is.na(output))) {
				warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(output)], collapse='\n'), 
						'\ndo not share at least ', intersect.no, ' common taxa with other samples! ',
						'For these samples, their size factors are set to be NA! \n', 
						'You may consider removing these samples since they are potentially outliers or negative controls!\n',
						'You may also consider decreasing the minimum number of intersecting taxa (intersect.no parameter) and rerun the procedure!\n'))
			}
			names(output) <- names(comm.no) <- colnames(comm)
			attr(output, 'NSS') <- comm.no
			output
		},
		# modified from https://github.com/liudoubletian/eBay
		ebay = function(comm) {
			sample_n <- nrow(comm)

			B_e <- MGLM::MGLMreg(comm~1, dist="DM")@coefficients

			gr <- matrix(rep(1, sample_n))
			alpha_e <- exp(gr%*%B_e)

			exp_norm <- comm

			for (n in 1:sample_n) {
				exp_norm[n, ] <- unlist(comm[n, ] + alpha_e[n, ]) / (sum(comm[n, ]) + sum(alpha_e[n, ]))
			}
			exp_norm
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
