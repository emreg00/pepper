
#' Get z score
#'
#' Calculates z-score for each gene in each sample
#' using the mean and sd over provided control samples.
#' @param expr Expression matrix (probes vs samples).
#' @param samples.background Names of the background (control) samples.
#' @param method Method to calculate the z score, defaults to mean and sd, use median for med and mad.
#' @param samples.to.exclude Names of the samples to be excluded from background (for CV).
#' @return Data frame containing z-scores
#' @keywords internal
#' @export
get.z.score<-function(expr, samples.background, method="mean", samples.to.exclude=NULL) {
    if(!is.null(samples.to.exclude)) { # For CV 
	samples.background = setdiff(samples.background, samples.to.exclude)
    }
    if(method == "mean") {
	#n = length(samples.background)
	vals.mean = apply(expr[,samples.background], 1, mean)
	vals.sd = apply(expr[,samples.background], 1, sd) #* sqrt((n-1)/n) for uncorrected sd 
    } else if(method == "median") {
	vals.mean = apply(expr[,samples.background], 1, median)
	vals.sd = apply(expr[,samples.background], 1, mad) 
    } else {
	stop(sprintf("Warning: Unknown method ", method))
    }
    z = (expr - vals.mean) / vals.sd
    return(z)
}


