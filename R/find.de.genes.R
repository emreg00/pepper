
#' Find differentially expressed genes
#'
#' @param expr expression matrix.
#' @param sample.mapping Sample - condition mapping.
#' @param states Conditions to be considered as case.
#' @param method Differential expression analysis method: limma (Default) | sam | welch.
#' @param out.file File to write output. If NULL (default) not used.
#' @param state.background Condition to be considered as control.
#' @param adjust.method Multiple hypothesis testing correction method. 
#'    Defaults to BH.
#' @param cutoff Adjust p-value cutoff. Defaults to 0.05
#' @param functional.enrichment GO or KEGG based functional enrichment analysis
#' @return data frame with results
#' @export
#' @examples
#' gds.data = fetch.expression.data("GDS4966", do.log2=F, probe.conversion="Gene ID")
#' expr = gds.data$expr
#' sample.mapping = gds.data$sample.mapping
#' sample.mapping = convert.sample.mapping.to.case.control(sample.mapping,
#'  		    states.control = c("healthy donor"), 
#' 		    states.case = c("tuberculosis", "latent tuberculosis infection"))
#' d = find.de.genes(expr, sample.mapping, c("case", "control"), method="limma")
find.de.genes<-function(expr, sample.mapping, states, method="limma", out.file=NULL, state.background=NULL, adjust.method='BH', cutoff=0.05, functional.enrichment=NULL) {
    d = NULL
    if(method == "all") {
	# SAM
	d = find.de.genes(expr, sample.mapping, states, method="sam", out.file=paste(out.file, ".sam", sep=""), state.background=state.background, adjust.method=adjust.method, cutoff=cutoff) 
	e.sam = d[d$adj.P.Val<=cutoff,]
	
	# Welch test
	d = find.de.genes(expr, sample.mapping, states, method="welch", out.file=paste(out.file, ".welch", sep=""), state.background=state.background, adjust.method=adjust.method, cutoff=cutoff) 
	e.welch = d[d$adj.P.Val<=cutoff,]

	# LIMMA
	d = find.de.genes(expr, sample.mapping, states, method="limma", out.file=paste(out.file, ".limma", sep=""), state.background=state.background, adjust.method=adjust.method, cutoff=cutoff) 
	e.limma = d[d$adj.P.Val<=cutoff,]

	d = union(union(e.sam, e.welch), e.limma)
	#print(c("limma vs sam", length(intersect(e.limma$GeneID, e.sam$GeneID))))
	#print(c("limma vs t", length(intersect(e.limma$GeneID, e.welch$GeneID))))
	#print(c("sam vs t", length(intersect(e.sam$GeneID, e.welch$GeneID))))
    } else if(method == "limma") {
	d = find.de.genes.limma(expr, sample.mapping, states, out.file, state.background, adjust.method, cutoff)
    } else if(method == "sam") {
	d = find.de.genes.sam(expr, sample.mapping, states, out.file, state.background, adjust.method, cutoff)
    } else if(method == "welch") {
	d = find.de.genes.welch(expr, sample.mapping, states, out.file, state.background, adjust.method, cutoff)
    }
    if(!is.null(functional.enrichment)) {
	if (!requireNamespace("limma", quietly=TRUE)) {
	    stop("Functional enrichment analysis requires limma package to be installed")
	}
	genes = rownames(d)
	genes = unlist(lapply(strsplit(genes, " /// "), function(x) {x[1]}))
	if(length(genes) > 0) {
	    if(!is.null(out.file)) {
		file.name = paste0(out.file, ".", functional.enrichment)
	    }
	    if(functional.enrichment == "go") {
		limma::go = goana(genes)
		a = limma::topGO(go, number=500)
	    } else if(functional.enrichment == "kegg") {
		kegg = limma::kegga(genes)
		a = limma::topKEGG(kegg, number=Inf)
	    } else {
		stop("Unrecognized functional enrichment type!")
	    }
	    if(is.null(out.file)) {
		print(a)
	    } else {
		write.table(a, file=file.name, row.names=T, quote=F, sep="\t")
	    }
	} else {
	    print("No DE genes at given cutoff!")
	}
    }
    return(d)
}


