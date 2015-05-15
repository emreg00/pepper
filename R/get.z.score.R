
#' Get z score
#'
#' Calculates z-score for each gene in each sample
#' using the mean and sd over provided control samples.
#' @param expr Expression matrix (probes vs samples).
#' @param samples.background Names of the background (control) samples.
#' @param method Method to calculate the z score, defaults to mean and sd, use median for med and mad.
#' @return Data frame containing z-scores
#' @keywords internal
#' @export
get.z.score<-function(expr, samples.background, method="mean") {
    if(method == "mean") {
	#n = length(samples.background)
	vals.mean = apply(expr[,samples.background], 1, mean)
	vals.sd = apply(expr[,samples.background], 1, sd) #* sqrt((n-1)/n) for uncorrected sd 
    } else if(method == "median") {
	vals.mean = apply(expr[,samples.background], 1, median)
	vals.sd = apply(expr[,samples.background], 1, mad) 
    } else {
	print(c("Warning: Unknown method", method))
    }
    samples = setdiff(colnames(expr), samples.background)
    z = (expr[, samples] - vals.mean) / vals.sd
    for(sample in samples.background) {
	samples = setdiff(samples.background, sample)
	if(method == "mean") {
	    vals.mean = apply(expr[,samples], 1, mean)
	    vals.sd = apply(expr[,samples], 1, sd) 
	} else if(method == "median") {
	    vals.mean = apply(expr[,samples], 1, median)
	    vals.sd = apply(expr[,samples], 1, mad) 
	}
	z[, sample] = (expr[,sample] - vals.mean) / vals.sd
    }
    return(z)
}



