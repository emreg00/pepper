
#' Convert z score
#'
#' Convert z-scores to P-values.
#' @param z z-scores.
#' @param one.sided if conversion should take the sign of the z-scores (two.sided if NULL).
#' @return P-values.
#' @keywords internal
#' @export
convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
	#pval = pnorm(-abs(z), mean=mean(z), sd=sd(z))
	pval = pnorm(-abs(z))
	pval = 2 * pval
	pval[pval > 1] = 1
    } else if(one.sided=="-") {
	pval = pnorm(z)
    } else {
	pval = pnorm(-z)
    }
    return(pval)
}
