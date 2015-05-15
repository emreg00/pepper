
#' Get z matrix
#'
#' Returns z-score matrix for a given GEO data set.
#' The z-scores are calculated for each gene in each sample
#' using the mean and sd over provided control samples.
#' @param expr Expression matrix (genes vs samples), tab separated.
#' @param sample.mapping Sample - condition mapping.
#' @param states.control Label of control (background) samples. 
#'        If NULL sample.mapping is assumed to include the following types: "case" "control".
#' @param states.case Label of case (disease) samples.
#' @param method Method to calculate the z score, defaults to mean and sd, use median for med and mad.
#' @param out.file Output file for writing z score matrix.
#' @return Data frame containing z-scores
#' @export
#' @examples
#' gds.data = fetch.expression.data("GDS4966", do.log2=F, probe.conversion="Gene ID")
#' expr = gds.data$expr
#' sample.mapping = gds.data$sample.mapping
#' z = get.z.matrix(expr, sample.mapping, 
#' 	states.control = c("healthy donor"), 
#' 	states.case = c("tuberculosis", "latent tuberculosis infection"))
get.z.matrix<-function(expr, sample.mapping, states.control=NULL, states.case=NULL, method="mean", out.file=NULL) {
    if(!is.null(out.file)) {
	if (file.exists(out.file)) {
	    z = read.table(out.file, sep="\t", header=T, check.names=F) 
	    return(z)
	}
    } 
    if(!is.null(states.control) & !is.null(states.case)){
	sample.mapping = convert.sample.mapping.to.case.control(sample.mapping, states.case, states.control)
    }
    # Get case/control samples
    samples = as.vector(sample.mapping[sample.mapping$type == "case", "sample"])
    samples.background = as.vector(sample.mapping[sample.mapping$type == "control", "sample"])
    expr = expr[,c(samples, samples.background)]
    z = get.z.score(expr, samples.background, method)
    # Ignoring NAs (genes with missing expression) to avoid future problems in counting
    z = na.omit(z)
    if(!is.null(out.file)) {
	write.table(z, file=out.file, row.names=T, quote=F, sep="\t") 
    }
    return(z)
}


