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
		#' Normalization/transformation methods include CLR, CCS, TSS, TMM, AST and those based on the \code{\link{decostand}} function in vegan package.
		#' @param method default NULL; See the following details and available options. \cr 
		#' \cr 
		#' Methods for normalization:
		#' \itemize{
		#'   \item \code{CLR}: Centered log-ratio normalization. 
		#'   \item \code{CCS}: Cumulative sum scaling normalization based on the \code{metagenomeSeq} package.
		#'   \item \code{TSS}: Total sum scaling, dividing counts by the sequencing depth.
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
		#' @param logbase default exp(1); The logarithm base used in method = "log" or "CLR".
		#' @param Cmin default NULL; see Cmin parameter in \code{SRS::SRS} function; Only available when \code{method = "SRS"}.
		#'    If not provided, use the minimum number across all the samples.
		#' @param pseudocount default 1; add pseudocount for those features with 0 abundance when \code{method = "CLR"}.
		#' @param ... parameters pass to \code{\link{decostand}} or \code{metagenomeSeq::cumNorm} when method = "CCS" or 
		#'    \code{edgeR::normLibSizes} when method = "TMM".
		#' 
		#' @return new microtable object or data.frame object.
		#' @examples
		#' newdataset <- t1$norm(method = "log")
		#' newdataset <- t1$norm(method = "CLR")
		norm = function(method = NULL, MARGIN = NULL, logbase = exp(1), Cmin = NULL, pseudocount = 1, ...)
			{
			abund_table <- self$data_table
			method <- match.arg(method, c("CLR", "CCS", "TSS", "TMM", "SRS", "AST", 
				"total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log"))
			
			if(method %in% c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")){
				if(is.null(MARGIN)){
					MARGIN <- switch(method, total = 1, max = 2, frequency = 2, normalize = 1, range = 2, rank = 1, standardize = 2, chi.square = 1, NULL)
				}
				res_table <- vegan::decostand(x = abund_table, method = method, MARGIN = MARGIN, logbase = logbase, ...)
			}
			if(method == "CLR"){
				if(any(abund_table == 0)){
					abund_table <- abund_table + pseudocount
				}
				if(is.null(MARGIN)){
					MARGIN <- 1
				}
				res_table <- vegan::decostand(x = abund_table, method = "clr", MARGIN = MARGIN, logbase = logbase)
			}
			if(method == "CCS"){
				obj <- metagenomeSeq::newMRexperiment(t(abund_table))
				obj_1 <- metagenomeSeq::cumNorm(obj, ...)
				res_table <- t(metagenomeSeq::MRcounts(obj_1, norm = TRUE))
			}
			if(method == "TSS"){
				res_table <- apply(abund_table, 1, function(x){x/sum(x)}) %>% t
			}
			if(method == "SRS"){
				newotu <- as.data.frame(t(abund_table))
				if(is.null(Cmin)){
					Cmin <- min(colSums(newotu))
				}
				res_table <- SRS::SRS(newotu, Cmin = Cmin, set_seed = TRUE, seed = 123)
				res_table <- t(res_table)
				colnames(res_table) <- colnames(abund_table)
			}
			if(method == "TMM"){
				libsize <- edgeR::normLibSizes(abund_table, method = "TMM", ...)
				effec_libsize <- colSums(abund_table) * libsize
				ref_libsize <- mean(effec_libsize)
				res_table <- sweep(abund_table, MARGIN = 2, effec_libsize, "/") * ref_libsize
			}
			if(method == "AST"){
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
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
