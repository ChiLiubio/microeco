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
		#' @param method default NULL; See the following details and available options. \cr 
		#' \cr 
		#' Methods for normalization:
		#' \itemize{
		#'   \item \code{GMPR}: Geometric mean of pairwise ratios <doi: 10.7717/peerj.4600>. 
		#'   \item \code{clr}: Centered log-ratio normalization <doi: 10.3389/fmicb.2017.02224>. 
		#'   \item \code{rclr}: Robust centered log-ratio normalization <doi: doi:10.1128/msystems.00016-19>.
		#'   \item \code{CCS}: Cumulative sum scaling normalization based on the \code{metagenomeSeq} package <doi:10.1038/nmeth.2658>.
		#'   \item \code{TSS}: Total sum scaling, divided by the sequencing depth.
		#'   \item \code{TMM}: Trimmed mean of M-values method based on the \code{normLibSizes} function of \code{edgeR} package.
		#'   \item \code{SRS}: scaling with ranked subsampling method based on the SRS package provided by Lukas Beule and Petr Karlovsky (2020) <DOI:10.7717/peerj.9593>.
		#' }
		#' Methods based on \code{\link{decostand}} function:
		#' \itemize{
		#'   \item \code{total}: divide by margin total (default MARGIN = 1, i.e. rows - samples).
		#'   \item \code{max}: divide by margin maximum (default MARGIN = 2, i.e. columns - features).
		#'   \item \code{normalize}:  make margin sum of squares equal to one (default MARGIN = 1).
		#'   \item \code{range}: standardize values into range 0...1 (default MARGIN = 2). If all values are constant, they will be transformed to 0.
		#'   \item \code{standardize}: scale x to zero mean and unit variance (default MARGIN = 2).
		#'   \item \code{pa}: scale x to presence/absence scale (0/1).
		#'   \item \code{log}: logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of the logarithm; zeros are left as zeros. Higher bases give less weight to quantities and more to presences, and logbase = Inf gives the presence/absence scaling. Please note this is not log(x+1). Anderson et al. (2006) suggested this for their (strongly) modified Gower distance (implemented as method = "altGower" in vegdist), but the standardization can be used independently of distance indices.
		#' }
		#' Other methods for transformation:
		#' \itemize{
		#'   \item \code{AST}: Arc sine square root transformation.
		#' }
		#' @param MARGIN default NULL; 1 = samples, and 2 = features of abundance table; only available when method comes from \code{\link{decostand}} function.
		#'    If MARGIN is NULL, use the default value in decostand function.
		#' @param logbase default exp(1); The logarithm base.
		#' @param Cmin default NULL; see Cmin parameter in \code{SRS::SRS} function; Only available when \code{method = "SRS"}.
		#'    If not provided, use the minimum number across all the samples.
		#' @param pseudocount default 1; add pseudocount for those features with 0 abundance when \code{method = "clr"}.
		#' @param intersect.no default 10; the intersecting taxa number between paired sample for \code{method = "GMPR"}.
		#' @param ct.min default 1; the minimum number of counts required to calculate ratios for \code{method = "GMPR"}.
		#' @param ... parameters pass to \code{\link{decostand}} or \code{metagenomeSeq::cumNorm} when method = "CCS" or 
		#'    \code{edgeR::normLibSizes} when method = "TMM" or GMPR function (intersect.no and ct.min parameters) when method = "GMPR" (https://github.com/jchen1981/GMPR).
		#' 
		#' @return new microtable object or data.frame object.
		#' @examples
		#' newdataset <- t1$norm(method = "log")
		#' newdataset <- t1$norm(method = "clr")
		norm = function(method = NULL, MARGIN = NULL, logbase = 2, Cmin = NULL, pseudocount = 1, intersect.no = 10, ct.min = 1, ...)
			{
			abund_table <- self$data_table
			if(is.null(method)){
				stop("Please select a method!")
			}
			method <- tolower(method)
			method <- match.arg(method, c("gmpr", "clr", "rclr", "ccs", "tss", "tmm", "srs", "ast", 
				"total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log"))
			
			if(method %in% c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")){
				if(is.null(MARGIN)){
					MARGIN <- switch(method, total = 1, max = 2, frequency = 2, normalize = 1, range = 2, rank = 1, standardize = 2, chi.square = 1, NULL)
				}
				res_table <- vegan::decostand(x = abund_table, method = method, MARGIN = MARGIN, logbase = logbase, ...)
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
			if(method == "ccs"){
				obj <- metagenomeSeq::newMRexperiment(t(abund_table))
				obj_1 <- metagenomeSeq::cumNorm(obj, ...)
				res_table <- t(metagenomeSeq::MRcounts(obj_1, norm = TRUE))
			}
			if(method == "tss"){
				res_table <- apply(abund_table, 1, function(x){x/sum(x)}) %>% t
			}
			if(method == "srs"){
				newotu <- as.data.frame(t(abund_table))
				if(is.null(Cmin)){
					Cmin <- min(colSums(newotu))
				}
				res_table <- SRS::SRS(newotu, Cmin = Cmin, set_seed = TRUE, seed = 123)
				res_table <- t(res_table)
				colnames(res_table) <- colnames(abund_table)
			}
			if(method == "tmm"){
				libsize <- edgeR::normLibSizes(abund_table, method = "TMM", ...)
				effec_libsize <- colSums(abund_table) * libsize
				ref_libsize <- mean(effec_libsize)
				res_table <- sweep(abund_table, MARGIN = 2, effec_libsize, "/") * ref_libsize
			}
			if(method == "ast"){
				res_table <- private$AST(abund_table)
			}
			if(inherits(self$dataset, "microtable")){
				res_dataset <- clone(self$dataset)
				res_dataset$otu_table <- as.data.frame(t(res_table))
				res_dataset
			}else{
				res_table
			}
		}
	),
	private = list(
		AST = function(x){
			sign(x) * asin(sqrt(abs(x)))
		},
		# modified based on https://github.com/jchen1981/GMPR
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
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
