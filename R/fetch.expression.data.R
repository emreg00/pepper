
#' Fetch expression data from GEO / or given folder
#'
#' @param geo.id GEO id.
#' @param do.log2 Apply log2 transformation to the expression values. Defaults to TRUE.
#' @param output.dir Directory where all files will be written.
#' @param probe.conversion Convert probe level expression to gene level using provided 
#' annotation label (uses platform specific annotations downloaded from GEO). 
#' Defaults to NULL (no conversion). In case of multiple probes, probe with absolute max
#' value is chosen.
#' @return A list containing 3 data frames: expression matrix, sample mapping, gene mapping.
#' @export
#' @examples
#' gds.data = fetch.expression.data("GDS4966", do.log2=F, probe.conversion="Gene ID")
#' expr = gds.data$expr
#' sample.mapping = gds.data$sample.mapping
fetch.expression.data<-function(geo.id, sample.mapping.column="characteristics_ch1", do.log2=NULL, probe.conversion=NULL, output.dir=paste(geo.id, "/", sep="")) {
    if (!file.exists(output.dir)){
	dir.create(file.path(output.dir))
    } 
    out.file = paste(output.dir, "expr.dat", sep="")
    out.file.2 = paste(output.dir, "mapping.dat", sep="")
    out.file.3 = paste(output.dir, "annotation.dat", sep="")
    # Load expression file if exists
    if(!all(file.exists(out.file), file.exists(out.file.2), file.exists(out.file.3))) {
	# Get data set
	data.set = get.data.set(geo.id, output.dir)
	# Using Biobase + GEOquery methods
	if(substr(geo.id,1,3) == "GDS") {
	    eset = GEOquery::GDS2eSet(data.set) #, do.log2=do.log2)
	    # Get sample - phenotype mapping
	    sample.mapping = Biobase::pData(eset)[,c(1:2)]
	    colnames(sample.mapping) = c("sample", "type")
	} else { # GSE
	    if (length(data.set) > 1) idx = grep(geo.id, attr(data.set, "names")) else idx = 1 
	    eset = data.set[[idx]]
	    sample.mapping = data.frame(sample=rownames(Biobase::pData(eset)), type=Biobase::pData(eset)[, sample.mapping.column])
	}
	write.table(sample.mapping, out.file.2, sep="\t", quote=F, row.names=F)
	expr = Biobase::exprs(eset)
	if(is.null(do.log2)) {
	    # Check whether data is already log transformed (from GEO2R)
	    qx = as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	    LogC = (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	    if (LogC) { 
		do.log2 = T
		expr[which(expr <= 0)] = NA
		print("Expression will be log2 transformed")
	    } else {
		do.log2 = F
	    }
	}
	if(do.log2) {
	    expr = log2(expr)
	}
	# This was below before causing to ignore genes with NA probes (i.e. results for joerg)
	expr = na.omit(expr)
	# Get probe - gene mapping
	if(!is.null(probe.conversion)) {
	    if(is.null(nrow(Biobase::fData(eset)))) {
		if(substr(geo.id,1,3) == "GDS") {
		    geo.id.platform = GEOquery::Meta(data.set)$platform
		} else {
		    geo.id.platform = Biobase::annotation(eset)
		}
		gene.mapping = get.platform.annotation(geo.id.platform, probe.conversion, output.dir)
	    } else {
		gene.mapping = data.frame(Probe = Biobase::fData(eset)[,"ID"], Gene = Biobase::fData(eset)[,probe.conversion]) 
	    }
	    expr = convert.probe.to.gene.expression(expr, gene.mapping) 
	} else {
	    gene.mapping = data.frame(Probe=rownames(expr), Gene=rownames(expr))
	}
	write.table(gene.mapping, out.file.3, sep="\t", quote=F, row.names=F)
	#expr = na.omit(expr) 
	write.table(expr, out.file, sep="\t", quote=F)
    }
    expr = read.table(out.file, sep="\t", header=T, check.names=F) 
    sample.mapping = read.table(out.file.2, sep="\t", header=T, quote='"') 
    gene.mapping = read.table(out.file.3, sep="\t", header=T, quote='"', check.names=F) #, row.names=1) 
    return(list(expr=expr, sample.mapping=sample.mapping, gene.mapping=gene.mapping))
}

