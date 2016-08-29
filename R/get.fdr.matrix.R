
#' Get fdr matrix 
#'
#' Calculates FDRs from z-scores for each gene in each sample.
#' @param z Data frame containing z-scores (probes vs samples).
#' @param adjust.method P-value correction method (see p.adjust).
#' @param out.file Output file for writing z score matrix.
#' @return Data frame containing FDR values.
#' @keywords internal
#' @export
get.fdr.matrix<-function(z, adjust.method, out.file=NULL) {
    if(!is.null(out.file)) {
	if (file.exists(out.file)) {
	    e = read.table(out.file, sep="\t", header=T, check.names=F) 
	    return(e)
	}
    } 
    e = abs(z) 
    # The result of below will be matrix (not data frame, loosing rownames)
    e = apply(e, 2, function(x) { convert.z.score(x) })
    e = apply(e, 2, function(x) { p.adjust(x, method=adjust.method) })
    colnames(e) = colnames(z)
    rownames(e) = rownames(z)
    if(!is.null(out.file)) {
	write.table(e, file=out.file, row.names=T, quote=F, sep="\t") 
    }
    return(e)
}

