
#' Convert probe level expression to gene level
#'
#' @param expr Expression matrix.
#' @param gene.mapping Probe-gene mapping.
#' @param selection.function Function to handle multiple probes 
#'	corresponding to the same gene. Defaults to NULL (choses
#'	probe with max expression).
#' @return a data frame containin gene level expression data
#' @keywords internal
#' @export 
#' @examples
#' #expr.gene = convert.probe.to.gene.expression(expr, gene.mapping)
convert.probe.to.gene.expression<-function(expr, gene.mapping, selection.function=NULL) {
    if(nrow(gene.mapping) != nrow(expr)) {
	print("Warning: dimension inconsisitency in gene annotation!")
	gene.mapping = as.data.frame(gene.mapping[gene.mapping$Probe %in% rownames(expr),]) #gene.mapping[rownames(gene.mapping) %in% rownames(expr),])
	expr = expr[rownames(expr) %in% gene.mapping$Probe,] #expr[rownames(gene.mapping),]
    }
    gene.mapping = factor(gene.mapping[,"Gene"])
    if(is.null(selection.function)) {
	selection.function<-function(asample){ # max
	   return(tapply(abs(asample), gene.mapping, max)) #mean)) 
	}
	#variances = apply(expr, 1, var)
	#get.max<-function(e) {
	#    idx = which.max(e[,2])
	#    return(e[idx,1])
	#}
	#selection.function<-function(asample) { # max.variance
	#    return(by(data.frame(asample, variances), gene.mapping, get.max))
	#}
    }
    expr.gene<-apply(expr, 2, selection.function)
    expr.gene = expr.gene[rownames(expr.gene) != "", ] # Filter the one from no gene probes
    return(expr.gene)
}

