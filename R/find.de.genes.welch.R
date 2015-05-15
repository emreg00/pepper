
#' Find differentially expressed genes using Welch 
#' (t-test w/ unequal variance) test
#'
#' @param expr expression matrix.
#' @param sample.mapping Sample - condition mapping.
#' @param states Conditions to be considered as case.
#' @param out.file File to write output. If NULL (default) not used.
#' @param state.background Condition to be considered as control.
#' @param adjust.method Multiple hypothesis testing correction method. 
#'    Defaults to BH.
#' @param cutoff Adjust p-value cutoff. Defaults to 0.2.
#' @return Data frame with results
find.de.genes.welch<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    if(is.null(state.background)) {
	state.background = "case"
    }
    expr = expr[,as.vector(sample.mapping$sample)]
    samples.background = which(sample.mapping$type == state.background)
    samples = setdiff(1:ncol(expr), samples.background)
    p.values = apply(expr, 1, function(x) { b = t.test(x[samples], x[samples.background]); return(b$p.value) })
    fc = apply(expr, 1, function(x) { log(mean(x[samples]) / mean(x[samples.background])) })

    d = data.frame(GeneID=names(p.values), logFC=as.vector(fc), P.Value=as.vector(p.values), adj.P.Val=p.adjust(p.values, method=adjust.method))
    d = d[order(d$adj.P.Val),]

    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
    }
    return(d) 
}

