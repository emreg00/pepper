
#' Find differentially expressed genes using LIMMA
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
find.de.genes.limma<-function(expr, sample.mapping, states, out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.2) {
    if (!requireNamespace("limma", quietly=TRUE)) {
	stop("Differential expression using LIMMA requires limma package to be installed")
    }
    #if(ncol(expr) != nrow(sample.mapping)) {
    #	print(c("Warning: inconsistent dimenstions!", ncol(expr), nrow(sample.mapping)))
    #}
    expr = expr[,as.vector(sample.mapping$sample)]
    design = model.matrix(~ 0 + sample.mapping$type)
    colnames(design) = gsub(" ", "_", states) #colnames(design))
    fit = limma::lmFit(expr, design)
    #contrast = unlist(sapply(states, function (x) { if(ref != x) { paste(ref, x, sep = "-") } }))
    contrast = c()
    if(!is.null(state.background)) {
	# W.r.t. background
	for(i in 1:length(states)) { 
	    if(states[i] != state.background)
		contrast = c(contrast, paste(gsub(" ", "_", states[i]), gsub(" ", "_", state.background), sep = "-"))
	}
    } else {
	# All vs all DE test
	for(i in 1:length(states)) { 
	    for(j in 1:length(states)) { 
		if(i<j) { 
		    contrast = c(contrast, paste(gsub(" ", "_", states[i]), gsub(" ", "_", states[j]), sep = "-"))
		}
	    }
	}
    }
    cont.matrix = limma::makeContrasts(contrasts=contrast, levels=design)
    fit2 = limma::contrasts.fit(fit, cont.matrix)
    fit2 = limma::eBayes(fit2)
    #topTable(fit2, coef=1, adjust="BH") #adjust="fdr", sort.by="B",
    d = limma::topTable(fit2, coef=1, adjust=adjust.method, sort.by="B", number=100000)
    d$GeneID = rownames(d)
    n = ncol(d)
    d = d[,c(n, 1:(n-1))]
    if(!is.null(out.file)) {
	write.table(d, file=out.file, row.names=F, quote=F, sep="\t")
	#results = decideTests(fit2, adjust.method=adjust.method, p.value=cutoff) # none
	#png(paste(out.file, ".png", sep=""))
	#vennDiagram(results)
	#dev.off()
    }
    return(d) #(d[d$adj.P.Val<=cutoff,])
}


