
#' Convert probe level expression to gene level
#'
#' @param expr Expression matrix.
#' @param gene.mapping Probe-gene mapping.
#' @param selection.method How to handle multiple probes 
#'	corresponding to the same gene. Defaults to iqr (probe with
#'	highest IQR). Other options var (highest variance), med (median
#'	of the probes)
#' @return a data frame containin gene level expression data
#' @keywords internal
#' @export 
#' @examples
#' #expr.gene = convert.probe.to.gene.expression(expr, gene.mapping)
convert.probe.to.gene.expression<-function(expr, gene.mapping, selection.method="iqr") {
    if(nrow(gene.mapping) != nrow(expr)) {
	print("Warning: dimension inconsisitency in gene annotation!")
	gene.mapping = as.data.frame(gene.mapping[gene.mapping$Probe %in% rownames(expr),]) #gene.mapping[rownames(gene.mapping) %in% rownames(expr),])
	expr = expr[rownames(expr) %in% gene.mapping$Probe,] #expr[rownames(gene.mapping),]
    }
    gene.mapping = factor(gene.mapping[,"Gene"])
    if(selection.method == "med") {
	selection.function<-function(asample) { 
	    # Max abs value of the probes
	    #return(tapply(asample, gene.mapping, function(x) { x[which.max(abs(x))] })) 
	    # Median of probes 
	    return(tapply(asample, gene.mapping, median)) 
	}
    } else {
	if(selection.method == "iqr") {
	    # Probe with highest IQR
	    values = apply(expr, 1, IQR)
	} else if(selection.method == "var") {
	    # Probe with max variance
	    values = apply(expr, 1, var)
	} else {
	    stop("Uknown selection method")
	}
	get.max<-function(e) {
	    idx = which.max(e[,2])
	    return(e[idx,1])
	}
	selection.function<-function(asample) { 
	    # Max variance or inter quartile region
	    return(by(data.frame(asample, values), gene.mapping, get.max))
	}
    }
    expr.gene<-apply(expr, 2, selection.function)
    expr.gene = expr.gene[rownames(expr.gene) != "", ] # Filter the one from no gene probes
    return(expr.gene)
}

