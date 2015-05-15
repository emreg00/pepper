
#' Find differentially expressed genes using SAM
#'
#' @param expr expression matrix.
#' @param sample.mapping Sample - condition mapping.
#' @param states Conditions to be considered as case.
#' @param out.file File to write output. If NULL (default) not used.
#' @param state.background Condition to be considered as control.
#' @param adjust.method Multiple hypothesis testing correction method. 
#'    Defaults to BH.
#' @param cutoff Adjust p-value cutoff. Defaults to 0.2
#' @return Data frame with results.
find.de.genes.sam<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    if (!requireNamespace("samr", quietly=TRUE)) {
	stop("Differential expression using SAM requires samr package to be installed")
    }
    if(is.null(state.background)) {
	state.background = "case"
    }
    expr = expr[,as.vector(sample.mapping$sample)]
    x = as.matrix(expr)
    y = ifelse(sample.mapping$type == state.background, 1, 2)

    samfit = samr::SAM(x, y, resp.type="Two class unpaired", fdr.output=cutoff)
    #samfit = invisible(samr::SAM(x, y, resp.type="Two class unpaired", fdr.output=cutoff))
    d = rbind(samfit$siggenes.table$genes.up, samfit$siggenes.table$genes.lo)
    d = data.frame(GeneID=d[,2], logFC=as.double(d[,6]), adj.P.Val=as.double(d[,7])/100)
    d = d[order(d$adj.P.Val),]

    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
    }
    return(d) 
}

