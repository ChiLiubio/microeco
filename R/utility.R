
#' Copy an R6 class object completely
#'
#' @param x R6 class object
#' @param deep default TRUE; deep copy
#' @return identical but unrelated R6 object.
#' @examples
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
#' dropallfactors(x)
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
#' tidy_taxonomy(taxonomy_table)
#' @export
tidy_taxonomy <- function(taxonomy_table){
	tidy_taxonomy_raw <- function(taxonomy_table, i){
		taxonomy_table[,i] %<>% 
			gsub(".*No blast hit.*", "", ., ignore.case = TRUE) %>%
			gsub(".*Unknown.*", "", ., ignore.case = TRUE) %>%
			gsub(".*unidentif.*", "", ., ignore.case = TRUE) %>%
			gsub(".*sp\\.$", "", ., ignore.case = TRUE) %>%
			gsub(".*Unclassified.*", "", ., ignore.case = TRUE) %>%
			gsub(".*metagenome.*", "", ., ignore.case = TRUE) %>%
			gsub(".*uncultur.*|.*cultivar.*", "", ., ignore.case = TRUE) %>%
			gsub("D_6__synthetic.*", "", ., ignore.case = TRUE) %>%
			gsub(".*archaeon$", "", ., ignore.case = TRUE) %>%
			gsub("^\\s+", "", .) %>%
			gsub("\\s+$", "", .) %>%
			gsub('"', "", ., fixed = TRUE) %>%
			gsub(".*\\sbacterium$|.*bacterium\\s.*", "", .)
		taxonomy_table[,i] <- gsub("^.*__", paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__"), taxonomy_table[,i])
		taxonomy_table[,i][is.na(taxonomy_table[,i])] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__")
		taxonomy_table[,i][grepl("^$",taxonomy_table[,i])] <- paste0(tolower(substr(colnames(taxonomy_table)[i], 1, 1)), "__")
		taxonomy_table[,i]
	}
	taxonomy_table[] <- lapply(seq_len(ncol(taxonomy_table)), function(x) tidy_taxonomy_raw(taxonomy_table, x))
	taxonomy_table
}



##################################################################################


# inner function
summarySE_inter = function(usedata=NULL, measurevar, groupvars=NULL, na.rm=TRUE) {
	length2 <- function(x, na.rm=TRUE) ifelse(na.rm, sum(!is.na(x)), length(x))
	datac <- usedata %>% 
			dplyr::grouped_df(groupvars) %>% 
			dplyr::summarise(N = length2(!!sym(measurevar), na.rm=na.rm), Mean = mean(!!sym(measurevar), na.rm=na.rm), SD = sd(!!sym(measurevar), na.rm=na.rm)) %>%
			as.data.frame
	datac$SE <- datac$SD / sqrt(datac$N)
	datac
}

##################################################################################
##################################################################################
##################################################################################
##################################################################################

#' sparcc wrapper
#'
#' An implementation of SparCC R code from https://github.com/zdk123/SpiecEasi
#' @param data Community count data matrix
#' @param iter Number of iterations in the outer loop
#' @param inner_iter Number of iterations in the inner loop
#' @param th absolute value of correlations below this threshold are considered zero by the inner SparCC loop.
#' @seealso \code{\link{sparccboot}}
#' @export
sparcc <- function(data, iter=20, inner_iter=10, th=.1) {
#  without all the 'frills'
    sparccs <- lapply(1:iter, function(i)
                      sparccinner(t(apply(data, 1, norm_diric)),
                                  iter=inner_iter, th=th))
    # collect
    cors <- array(unlist(lapply(sparccs, function(x) x$Cor)),
                 c(ncol(data),ncol(data),iter))
    corMed <- apply(cors, 1:2, median)
    covs <- array(unlist(lapply(sparccs, function(x) x$Cov)),
                 c(ncol(data),ncol(data),iter))
    covMed <- apply(covs, 1:2, median)
    covMed <- cor2cov(corMed, sqrt(diag(covMed)))
    list(Cov=covMed, Cor=corMed)
}

#' Bootstrap SparCC
#'
#' Get bootstrapped estimates of SparCC correlation coefficients. To get empirical p-values, pass this output to \code{pval.sparccboot}. 
#'
#' @param data Community count data
#' @param sparcc.params named list of parameters to pass to \code{sparcc}
#' @param statisticboot function which takes data and bootstrap sample indices and results the upper triangle of the bootstapped correlation matrix
#' @param statisticperm function which takes data and permutated sample indices and results the upper triangle of the null correlation matrix
#' @param R number of bootstraps
#' @param ncpus number of cores to use for parallelization
#' @param ... additional arguments that are passed to \code{boot::boot}
#' @export
sparccboot <- function(data, sparcc.params=list(),
                        statisticboot=function(data, indices) tril(do.call("sparcc",
                      c(list(data[indices,,drop=FALSE]), sparcc.params))$Cor),
                statisticperm=function(data, indices) tril(do.call("sparcc",  c(list(apply(data[indices,], 2, sample)), sparcc.params))$Cor),
                      R, ncpus=1, ...) {

    if (!requireNamespace('boot', quietly=TRUE))
      stop('\'boot\' package is not installed')

    res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=ncpus, ...)
    null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus)
    class(res) <- 'list'
    structure(c(res, list(null_av=null_av)), class='sparccboot')
}

#' SparCC p-vals
#'
#' Get empirical p-values from bootstrap SparCC output, modified based on https://github.com/zdk123/SpiecEasi
#'
#' @param x output from \code{sparccboot}
#' @param sided type of p-value to compute. Only two sided (sided="both") is implemented.
#' @export
pval.sparccboot <- function(x, sided='both') {
# calculate 1 or 2 way pseudo p-val from boot object
# Args: a boot object
    if (sided != "both") stop("only two-sided currently supported")
    nparams  <- ncol(x$t)
    tmeans   <- colMeans(x$null_av$t)
#    check to see whether correlations are unstable -- confirm
#    that sample correlations are in 95% confidence interval of
#    bootstrapped samples
    niters   <- nrow(x$t)
    ind95    <- max(1,round(.025*niters)):round(.975*niters)
    boot_ord <- apply(x$t, 2, sort)
    boot_ord95 <- boot_ord[ind95,]
    outofrange <- unlist(lapply(1:length(x$t0), function(i) {
            aitvar <- x$t0[i]
            range  <- range(boot_ord95[,i])
            range[1] > aitvar || range[2] < aitvar
        }))
    # calc whether center of mass is above or below the mean
    bs_above <- unlist(lapply(1:nparams, function(i)
                    length(which(x$t[, i] > tmeans[i]))))
    is_above <- bs_above > x$R/2
    cors <- x$t0
#    signedAV[is_above] <- -signedAV[is_above]
    pvals    <- ifelse(is_above, 2*(1-bs_above/x$R), 2*bs_above/x$R)
    pvals[pvals > 1]  <- 1
    pvals[outofrange] <- NaN
	# reshape
	use_names <- colnames(x$data)
	com_res <- t(combn(use_names, 2))
	res <- cbind.data.frame(com_res, cors, pvals, stringsAsFactors = FALSE)
	colnames(res) <- c("t1", "t2", "cor", "p")
	res <- rbind.data.frame(res, data.frame(t1 = res$t2, t2 = res$t1, cor = res$cor, p = res$p), 
		data.frame(t1 = use_names, t2 = use_names, cor = rep(1, length(use_names)), p = rep(0, length(use_names))))
	res$cor %<>% as.numeric
	res$p %<>% as.numeric
	res_cor <- reshape2::dcast(res, t1~t2, value.var = "cor") %>% `row.names<-`(.[,1]) %>% .[, -1] %>% .[use_names, use_names] %>% as.matrix
	res_p <- reshape2::dcast(res, t1~t2, value.var = "p") %>% `row.names<-`(.[,1]) %>% .[, -1] %>% .[use_names, use_names] %>% as.matrix
	res <- list(cor = res_cor, p = res_p)
	res
}



#' @keywords internal
sparccinner <- function(data.f, T=NULL, iter=10, th=0.1) {
    if (is.null(T))   T  <- av(data.f)
    res.bv <- basis_var(T)
    Vbase  <- res.bv$Vbase
    M      <- res.bv$M
    cbase  <- C_from_V(T, Vbase)
    Cov    <- cbase$Cov
    Cor    <- cbase$Cor

    ## do iterations here
    excluded <- NULL
    for (i in 1:iter) {
        res.excl <- exclude_pairs(Cor, M, th, excluded)
        M <- res.excl$M
        excluded <- res.excl$excluded
        if (res.excl$break_flag) break
        res.bv <- basis_var(T, M=M, excluded=excluded)
        Vbase  <- res.bv$Vbase
        M      <- res.bv$M
        K <- M
        diag(K) <- 1
        cbase  <- C_from_V(T, Vbase)
        Cov    <- cbase$Cov
        Cor    <- cbase$Cor
    }
    list(Cov=Cov, Cor=Cor, i=i, M=M, excluded=excluded)
}

#' @keywords internal
exclude_pairs <- function(Cor, M, th=0.1, excluded=NULL) {
# exclude pairs with high correlations
    break_flag <- FALSE
    C_temp <- abs(Cor - diag(diag(Cor)) )  # abs value / remove diagonal
    if (!is.null(excluded)) C_temp[excluded] <- 0 # set previously excluded correlations to 0
    exclude <- which(abs(C_temp - max(C_temp)) < .Machine$double.eps*100)[1:2]
    if (max(C_temp) > th)  {
        i <- na.exclude(arrayInd(exclude, c(nrow(M), ncol(M)))[,1])
        M[i,i] <- M[i,i] - 1
        excluded_new <- c(excluded, exclude)
    } else {
        excluded_new <- excluded
        break_flag   <- TRUE
    }
    list(M=M, excluded=excluded_new, break_flag=break_flag)
}

#' @keywords internal
basis_cov <- function(data.f) {
# data.f -> relative abundance data
# OTUs in columns, samples in rows (yes, I know this is transpose of normal)
    # first compute aitchison variation
    T <- av(data.f)
    res.bv <- basis_var(T)
    Vbase  <- res.bv$Vbase
    M      <- res.bv$M
    cbase  <- C_from_V(T, Vbase)
    Cov    <- cbase$Cov
    Cor    <- cbase$Cor
    list(Cov=Cov, M=M)
}

#' @keywords internal
basis_var <- function(T, CovMat = matrix(0, nrow(T), ncol(T)),
                      M = matrix(1, nrow(T), ncol(T)) + (diag(ncol(T))*(ncol(T)-2)),
                      excluded = NULL, Vmin=1e-4) {

    if (!is.null(excluded)) {
        T[excluded] <- 0
     #   CovMat[excluded] <- 0
    }
    Ti     <- matrix(rowSums(T))
    CovVec <- matrix(rowSums(CovMat - diag(diag(CovMat)))) # row sum of off diagonals
    M.I <- tryCatch(solve(M), error=function(e) MASS::ginv(M))
    Vbase <- M.I %*% (Ti + 2*CovVec)
    Vbase[Vbase < Vmin] <- Vmin
    list(Vbase=Vbase, M=M)
}

#' @keywords internal
C_from_V <- function(T, Vbase) {
    J      <- matrix(1, nrow(T), ncol(T))
    Vdiag  <- diag(c(Vbase))
    CovMat <- .5*((J %*% Vdiag) + (Vdiag %*% J) - T)
    CovMat <- (CovMat + t(CovMat))/2  # enforce symmetry
    # check that correlations are within -1,1
    CorMat <- cov2cor(CovMat)
    CorMat[abs(CorMat) > 1] <- sign(CorMat[abs(CorMat) > 1])
    CovMat <- cor2cov(CorMat, sqrt(as.vector(Vbase)))
    list(Cov=CovMat, Cor=CorMat)
}

#' @keywords internal
av <- function(data) {
    cov.clr <- cov(clr(data))
    J <- matrix(1, ncol(data), ncol(data))
    (J %*% diag(diag(cov.clr))) + (diag(diag(cov.clr)) %*% J) - (2*cov.clr)
}

#' @keywords internal
norm_diric   <- function(x, rep=1) {
    dmat <- VGAM::rdiric(rep, x+1)
    norm_to_total(colMeans(dmat))
}

#' Centered log-ratio functions
#' @export
clr <- function(x, ...) {
  UseMethod('clr', x)
}

#' @method clr default
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
  apply(x.f, mar, clr, ...)
}

#' @method clr data.frame
#' @export
clr.data.frame <- function(x.f, mar=2, ...) {
  clr(as.matrix(x.f), mar, ...)
}

#' Total sum normalization of a (presumably positive) data vector
norm_to_total <- function(x) x/sum(x)

#' @keywords internal
tril <- function(x) x[lower.tri(x)]

#' Convert a symmetric correlation matrix to a covariance matrix
#' given the standard deviation
#'
#' @param cor a symmetric correlation matrix
#' @param sds standard deviations of the resulting covariance.
#' @return Covariance matrix of sample dimension as cor
cor2cov <- function(cor, sds) {
  if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
  cor * sds * rep(sds, each=nrow(cor))
}




##################################################################################
##################################################################################
##################################################################################
##################################################################################

# add significance label
#' @keywords internal
stat_compare_signif <- function(mapping = NULL, data = NULL,
                     method = "wilcox.test", paired = FALSE, ref.group = NULL,
                     comparisons = NULL, hide.ns = TRUE, label.sep = ", ",
                     label = NULL, label.x.npc = "left", label.y.npc = "top",
                     label.x = NULL, label.y = NULL, increase_signif = 0.15, tip_length = 0.008,
                     geom = "text", position = "identity",  na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, map_signif_level = TRUE, ...) {

    method.info <- .method_info(method)
    method <- method.info$method

    test.args <- list(paired = paired)
    if(method == "wilcox.test")
      test.args$exact <- FALSE

    pms <- list(...)
    size <- ifelse(is.null(pms$size), 0.3, pms$size)
    color <- ifelse(is.null(pms$color), "black", pms$color)

    step_increase <- ifelse(is.null(label.y), increase_signif, 0)
    geom_signif(comparisons = comparisons, y_position = label.y,
                          test = method, test.args = test.args, tip_length = tip_length,
                          step_increase = step_increase, size = size, color = color,
                          map_signif_level = map_signif_level)

}

#' @keywords internal
geom_signif <- function(mapping = NULL, data = NULL, stat = "signif",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, comparisons=NULL, test="wilcox.test", test.args=NULL,
                        annotations=NULL, map_signif_level=TRUE, y_position=NULL,xmin=NULL, xmax=NULL,
                        margin_top=0.1, step_increase=0, tip_length=0.008,
                        size=0.5, textsize = 3.88, family="", vjust = 0,
                        ...) {
  params <- list(na.rm = na.rm, ...)
  if (identical(stat, "signif")) {
    params <- c(params, list(comparisons=comparisons, test=test, test.args=test.args,
                   annotations=annotations, map_signif_level=map_signif_level,
                   y_position=y_position,xmin=xmin, xmax=xmax,
                   margin_top=margin_top, step_increase=step_increase,
                   tip_length=tip_length))
  }
  ggplot2::layer(
    stat = StatSignif, geom = GeomSignif, mapping = mapping,  data = data,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = params
  )
}

#' @keywords internal
GeomSignif <- ggplot2::ggproto("GeomSignif", ggplot2::Geom,
                           required_aes = c("x", "xend", "y", "yend", "annotation"),
                           default_aes = ggplot2::aes(shape = 19, colour = "black", textsize = 3.88, angle = 0, hjust = 0.5,
                                             vjust = 0, alpha = NA, family = "", fontface = 1, lineheight = 1.2, linetype=1, size=0.5),
                           draw_key = ggplot2::draw_key_point,

                           draw_group = function(data, panel_params, coord) {
                             coords <- coord$transform(data, panel_params)
                             grid::gList(
                               grid::textGrob(
                                 label=as.character(coords$annotation),
                                 x=mean(c(coords$x[1], tail(coords$xend, n=1))),
                                 y=max(c(coords$y, coords$yend)) - 0.045,
                                 default.units = "native",
                                 hjust = coords$hjust, vjust = coords$vjust,
                                 rot = coords$angle,
                                 gp = grid::gpar(
                                   col = scales::alpha(coords$colour, coords$alpha),
                                   fontsize = coords$textsize * ggplot2::.pt,
                                   fontfamily = coords$family,
                                   fontface = coords$fontface,
                                   lineheight = coords$lineheight
                                 )
                               ),
                               grid::segmentsGrob(
                                 coords$x, coords$y-0.05,
                                 default.units = "native",
                                 coords$xend, coords$yend-0.05,
                                 gp = grid::gpar(
                                   col = scales::alpha(coords$colour, coords$alpha),
                                   lty = coords$linetype,
                                   lwd = coords$size * ggplot2::.pt
                                 )
                               )
                             )
                           }
)





#' @keywords internal
StatSignif <- ggplot2::ggproto("StatSignif", ggplot2::Stat,
                  required_aes = c("x", "y", "group"),
                  setup_params = function(data, params) {
                    # if(any(data$group == -1)|| any(data$group != data$x)){
                    if(any(data$group == -1)){
                      stop("Can only handle data with groups that are plotted on the x-axis")
                    }

                    if (is.character(params$test)) params$test <- match.fun(params$test)
                    params$complete_data <- data
                    if(is.null(params$xmin) != is.null(params$xmax) || length(params$xmin) != length(params$xmax))
                      stop("If xmin or xmax is set, the other one also needs to be set and they need to contain the same number of values")
                    if(!is.null(params$xmin) && !is.null(params$comparisons))
                      stop("Set either the xmin, xmax values or the comparisons")
                    if(!is.null(params$xmin) && is.null(params$y_position))
                      stop("If xmin, xmax are defined also define y_position")
                    if(! is.null(params$y_position) && length(params$y_position) == 1)
                      params$y_position <- rep(params$y_position, max(length(params$comparisons), length(params$xmin)))
                    if(length(params$margin_top) == 1) params$margin_top <- rep(params$margin_top, max(length(params$comparisons),length(params$xmin)))
                    if(length(params$step_increase) == 1) params$step_increase <- rep(params$step_increase, max(length(params$comparisons),length(params$xmin)))
                    if(length(params$tip_length) == 1) params$tip_length <- rep(params$tip_length, max(length(params$comparisons),length(params$xmin)) * 2)
                    if(length(params$tip_length) == length(params$comparisons)) params$tip_length <- rep(params$tip_length, each=2)
                    if(length(params$tip_length) == length(params$xmin)) params$tip_length <- rep(params$tip_length, each=2)
                    if(! is.null(params$annotations) && length(params$annotations) == 1)
                      params$annotations <- rep(params$annotations, max(length(params$comparisons),length(params$xmin)))
                    if(! is.null(params$annotations) && length(params$annotations) != max(length(params$comparisons),length(params$xmin)))
                      stop(paste0("annotations contains a different number of elements (", length(params$annotations),
                                  ") than comparisons or xmin (", max(length(params$comparisons),length(params$xmin)), ")."))

                    if(all(params$map_signif_level == TRUE)){
                      params$map_signif_level <- c("***"=0.001, "**"=0.01, "*"=0.05)
                    }else if(is.numeric(params$map_signif_level)){
                      if(is.null(names(params$map_signif_level)) ){
                        if(length(params$map_signif_level) <= 3){
                          names(params$map_signif_level) <- tail(c("***", "**", "*"), n=length(params$map_signif_level))
                        }else{
                          stop('Cannot handle un-named map for significance values, please provide in the following format: c("***"=0.001, "**"=0.01, "*"=0.05)')
                        }
                      }
                    }
                    return(params)
                  },
                  compute_group = function(data, scales, comparisons, test, test.args, complete_data,
                                           annotations, map_signif_level, y_position, xmax, xmin,
                                           margin_top, step_increase, tip_length) {
                    if(! is.null(comparisons)){
                      i <- 0
                      result <- lapply(comparisons, function(comp){
                        i <<- i + 1
                        # All entries in group should be the same
                        if(scales$x$map(comp[1]) == data$group[1]){
                          test_result <- if(is.null(annotations)){
                            group_1 <- complete_data$y[complete_data$x == scales$x$map(comp[1]) & complete_data$PANEL == data$PANEL[1]]
                            group_2 <- complete_data$y[complete_data$x == scales$x$map(comp[2]) & complete_data$PANEL == data$PANEL[1]]
                            p_value <- do.call(test, c(list(group_1, group_2), test.args))$p.value
                            if(is.numeric(map_signif_level)){
                              temp_value <- names(which.min(map_signif_level[which(map_signif_level > p_value)]))
                              if(is.null(temp_value)){
                                "ns"
                              }else{
                                temp_value
                              }
                            }else{
                              if(p_value < 2.2e-16){
                                "< 2.2e-16"
                              }else{
                                as.character(signif(p_value, digits=2))
                              }

                            }
                          }else{
                            annotations[i]
                          }
                          y_scale_range <- (scales$y$range$range[2] - scales$y$range$range[1])
                          y_pos <- if(is.null(y_position)){
                            scales$y$range$range[2] + y_scale_range * margin_top[i] + y_scale_range * step_increase[i] * (i-1)
                          }else{
                            y_pos <- y_position[i]
                          }
                          data.frame(x=c(min(comp[1],comp[2]),min(comp[1],comp[2]),max(comp[1],comp[2])),
                                     xend=c(min(comp[1],comp[2]),max(comp[1],comp[2]),max(comp[1],comp[2])),
                                     y=c(y_pos - y_scale_range*tip_length[(i-1)*2+1], y_pos, y_pos),
                                     yend=c(y_pos, y_pos, y_pos-y_scale_range*tip_length[(i-1)*2+2]),
                                     annotation=test_result, group=paste(c(comp, i), collapse = "-"))
                        }
                      })
                      do.call(rbind, result)
                    }else{
                      if(data$x[1] == min(complete_data$x) & data$group[1] == min(complete_data$group)){
                        y_scale_range <- (scales$y$range$range[2] - scales$y$range$range[1])
                        data.frame(x=c(xmin, xmin, xmax),
                                   xend=c(xmin, xmax, xmax),
                                   y=c(y_position - y_scale_range*tip_length[seq_len(length(tip_length))%% 2 == 1], y_position, y_position),
                                   yend=c(y_position, y_position, y_position-y_scale_range*tip_length[seq_len(length(tip_length))%% 2 == 0]),
                                   annotation=rep(annotations, times=3), group=rep(seq_along(xmin), times=3))
                      }
                    }
                  }
)




# Check and get test method info
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# return a list(method, name)
.method_info <- function(method){

  if(is.null(method))
    method = "wilcox.test"

  allowed.methods <- list(
    t = "t.test", t.test = "t.test", student = "t.test",
    wiloxon = "wilcox.test", wilcox = "wilcox.test", wilcox.test = "wilcox.test",
    anova = "anova", aov = "anova",
    kruskal = "kruskal.test", kruskal.test = "kruskal.test")

  method.names <- list(
    t.test = "T-test", wilcox.test = "Wilcoxon",
    anova = "Anova", kruskal.test = "Kruskal-Wallis")

  if(!(method %in% names(allowed.methods)))
    stop("Non-supported method specified. Allowed methods are one of: ",
         .collapse(allowed.methods, sep =", "))
  method <- allowed.methods[[method]]
  method.name <- method.names[[method]]

  list(method = method, name = method.name)
}



# Collapse one or two vectors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
.collapse <- function(x, y = NULL, sep = "."){
  if(missing(y))
    paste(x, collapse = sep)
  else if(is.null(x) & is.null(y))
    return(NULL)
  else if(is.null(x))
    return (as.character(y))
  else if(is.null(y))
    return(as.character(x))
  else
    paste0(x, sep, y)
}




##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

# metastat code from Statistical methods for detecting differentially abundant features in clinical metagenomic samples
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
		C1[1,2] = var(Pmatrix[1:g-1]); # variance
		C1[1,3] = C1[1,2]/(g-1);    # std err^2

		C2[1,1] = mean(Pmatrix[g:ncols]);
		C2[1,2] = var(Pmatrix[g:ncols]);  # variance
		C2[1,3] = C2[1,2]/(ncols-g+1); # std err^2
	}else{
		# generate statistic profiles for both groups
		# mean, var, stderr
		for (i in 1:nrows){ # for each taxa
			# find the mean of each group
			C1[i,1] = mean(Pmatrix[i, 1:g-1]);  
			C1[i,2] = var(Pmatrix[i, 1:g-1]); # variance
			C1[i,3] = C1[i,2]/(g-1);    # std err^2

			C2[i,1] = mean(Pmatrix[i, g:ncols]);  
			C2[i,2] = var(Pmatrix[i, g:ncols]);  # variance
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

	f <- unclass(smooth.spline(lambdas,pi0_hat,df=3));
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
				csqt <- chisq.test(contingencytable);
				pvalues[i] = csqt$p.value;
			}else{
				ft <- fisher.test(contingencytable, workspace = 8e6, alternative = "two.sided", conf.int = FALSE);
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
			C1[i,2] = var(Pmatrix[i, 1:g-1]); # variance
			C1[i,3] = C1[i,2]/(g-1);    # std err^2 (will change to std err at end)

			C2[i,1] = mean(Pmatrix[i, g:ncols]);  
			C2[i,2] = var(Pmatrix[i, g:ncols]);  # variance
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


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


