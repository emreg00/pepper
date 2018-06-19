
#' Get peeps from z-score matrix
#'
#' Returns peeps for each sample in a given z-score matrix.
#' @param z Matrix containing z-scores (genes vs samples), tab separated.
#' @param cutoff Threshold for deciding peeps, either a z-score or adjusted p-value (if convert.to.pvalues=T).
#' @param convert.to.pvalues Flag to convert z-scores to p-values. 
#'	  If TRUE, the z-scores are converted to P-values which are then corrected for multiple hypothesis testing.
#' @param adjust.method P-value correction method (see p.adjust).
#' @return Data frame containing sample name and geneid of genes in the peeps
#' @export
#' @examples
#' peeps <- get.peeps.from.z.matrix(z, cutoff=0.05, convert.to.pvalues=T) 
get.peeps.from.z.matrix<-function(z, cutoff, convert.to.pvalues, adjust.method="BH") {
    if (!requireNamespace("plyr", quietly=TRUE)) {
	stop("This function requires plyr package to be installed!")
    }
    # Get peep geneids for each sample
    e = abs(z) 
    if(convert.to.pvalues) {
	# The result of below will be matrix (not data frame, loosing rownames)
	e = apply(e, 2, function(x) { convert.z.score(x) })
	e = apply(e, 2, function(x) { p.adjust(x, method=adjust.method) })
    }
    if(convert.to.pvalues) {
	indices = apply(e, 2, function(x) { which(x <= cutoff)})
    } else {
	indices = apply(e, 2, function(x) { which(x >= cutoff)})
    }
    #print(indices) 
    geneids = lapply(indices, function(x) { rownames(z)[x] })
    d = ldply(geneids, cbind)
    # below gives error in batch run, although okey in interactive run
    colnames(d) = c("sample", "geneid")
    return(d)
}

