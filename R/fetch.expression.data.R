
#' Fetch expression data from GEO / or given folder
#'
#' @param geo.id GEO id.
#' @param do.log2 Apply log2 transformation to the expression values. Defaults to TRUE.
#' @param output.dir Directory where all files will be written.
#' @param probe.conversion Convert probe level expression to gene level using provided 
#' annotation label (uses platform specific annotations downloaded from GEO). 
#' Defaults to NULL (no conversion). In case of multiple probes, probe with absolute max
#' value is chosen.
#' @param conversion.mapping Mapping of platform specific ids to user provided ids
#' @param conversion.mapping.function Function to process the mapped name such that it matches
#' with the ids provided in conversion.map
#' @param geo.id.sub GEO id for the sub-set (e.g., specific to the platform).
#' @return A list containing 3 data frames: expression matrix, sample mapping, gene mapping.
#' @export
#' @examples
#' gds.data = fetch.expression.data("GDS4966", do.log2=F, probe.conversion="Gene ID")
#' expr = gds.data$expr
#' sample.mapping = gds.data$sample.mapping
fetch.expression.data<-function(geo.id, sample.mapping.column="characteristics_ch1", do.log2=NULL, probe.conversion=NULL, conversion.mapping=NULL, conversion.mapping.function=NULL, output.dir=paste(geo.id, "/", sep=""), geo.id.sub=NULL, reprocess=NULL) {
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
	    idx = which(colnames(Biobase::pData(eset)) == sample.mapping.column)
	    if(length(idx) == 0) { # sample.mapping.column == "characteristics_ch1") {
		idx = 2 
	    } 
	    sample.mapping = Biobase::pData(eset)
	    write.table(sample.mapping, paste0(out.file.2, ".full"), sep="\t", quote=F, row.names=F)
	    sample.mapping = Biobase::pData(eset)[,c(1, idx)] #[,c(1:2)]
	    colnames(sample.mapping) = c("sample", "type")
	} else { # GSE
	    if (length(data.set) > 1) {
		if(is.null(geo.id.sub)) {
		    idx = grep(geo.id, attr(data.set, "names")) 
		} else {
		    idx = grep(geo.id.sub, attr(data.set, "names")) 
		}
	    } else {
		idx = 1 
	    }
	    eset = data.set[[idx]]
	    sample.mapping = Biobase::pData(eset)
	    sample.mapping$sample = rownames(Biobase::pData(eset))
	    write.table(sample.mapping, paste0(out.file.2, ".full"), sep="\t", quote=F, row.names=F)
	    sample.mapping = data.frame(sample=rownames(Biobase::pData(eset)), type=Biobase::pData(eset)[, sample.mapping.column])
	}
	write.table(sample.mapping, out.file.2, sep="\t", quote=F, row.names=F)
	# Get expression data
	if(!is.null(reprocess)) {
	    output.dir.raw = paste0(output.dir, "raw/")
	    if (!file.exists(output.dir.raw)){
		file.paths = GEOquery::getGEOSuppFiles(geo.id, baseDir = paste0(output.dir, "../")) 
		untar(paste0(output.dir, geo.id, "_RAW.tar"), exdir=output.dir.raw)
	    }
	    if(substr(reprocess, 1, nchar("affy")) == "affy") {
		if (!requireNamespace("affy", quietly=TRUE)) {
		    stop("Reprocessing of Affymetrix arrays requires affy package to be installed")
		}
		arguments = unlist(strsplit(reprocess, "\\|"))
		if(length(arguments) == 1) {
		    d = affy::ReadAffy(celfile.path = output.dir.raw) 
		    Biobase::sampleNames(d) <- sub("(_|\\.).*CEL\\.gz","",  Biobase::sampleNames(d))
		    eset.raw = affy::rma(d) 
		    expr = Biobase::exprs(eset.raw)
		} else {
		    if (!requireNamespace("makecdfenv", quietly=TRUE)) {
		    	stop("Reprocessing of Affymetrix arrays with custom annotation requires makecdfenv package to be installed")
		    }
		    # CDF environment is not recognized in rma (potentially due to namespace issues)
		    # Load manually created expr instead
		    if(file.exists(paste0(out.file, ".probe"))) {
			expr = read.table(paste0(out.file, ".probe"), sep="\t", header=T, check.names=F) 
		    } else {
			cdf = makecdfenv::make.cdf.env(arguments[2], cdf.path = output.dir.raw, compress=T) 
			environment(cdf) = asNamespace('affy')
			d = affy::ReadAffy(celfile.path = output.dir.raw, cdfname="cdf")
			Biobase::sampleNames(d) <- sub("(_|\\.).*CEL\\.gz","",  Biobase::sampleNames(d))
			eset.raw = affy::rma(d) 
			expr = Biobase::exprs(eset.raw)
			#write.table(expr, paste0(out.file, ".probe"), sep="\t", quote=F)
		    }
		}
	    } else if(substr(reprocess, 1, nchar("illumina")) == "illumina") {
		#if (!requireNamespace("beadarray", quietly=TRUE)) {
		#    stop("Reprocessing of Illumina arrays requires beadarray package to be installed")
		#}
		if (!requireNamespace("limma", quietly=TRUE) | !requireNamespace("preprocessCore", quietly=TRUE)) {
		    stop("Reprocessing of Illumina arrays requires limma and preprocessCore packages to be installed")
		}
		arguments = unlist(strsplit(reprocess, "\\|"))
		d = read.csv(paste0(output.dir, "filelist.txt"), sep="\t")
		data.set = limma::read.ilmn(d[d$Type=="TXT","Name"], path=output.dir.raw, probeid = arguments[2], expr = arguments[3]) # "ProbeSet_name" "Beadstudio_intensity"
		expr = apply(data.set$E, 2, as.numeric)
		expr = preprocessCore::normalize.quantiles(expr)
		rownames(expr) = rownames(data.set$E)
		colnames(expr) = sub("(_|\\.).*txt\\.gz","", d[d$Type=="TXT","Name"])
		# ctrlfiles = "AsuragenMAQC -controls.txt", annotation = "TargetID", other.columns = c("Detection  Pval") 
		#d = beadarray::readBeadSummaryData(paste0(output.dir, geo.id, "_RAW.tar"), skip=0, ProbeID = "ProbeSet_name", columns=list(exprs = "Beadstudio_intensity")) # AVG_Signal
		#d = readIllumina(output.dir.raw, useImages=FALSE, illuminaAnnotation = "Humanv3")
		#expr = Biobase::exprs(d)
		#eset = beadarray::normaliseIllumina(expr, method="quantile", transform="log2")
	    } else {
		stop("Unrecognized reprocessing type")
	    }
	} else {
	    expr = Biobase::exprs(eset)
	}
	if(is.null(do.log2)) {
	    # Check whether data is already log transformed (from GEO2R)
	    qx = as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	    LogC = (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	    if (LogC) { 
		do.log2 = T
	    } else {
		do.log2 = F
	    }
	}
	if(do.log2) {
	    print("Expression will be log2 transformed")
	    # To avoid leaving out probes with negative values
	    #expr[which(expr <= 0)] = NA
	    if(min(expr, na.rm=T) < 0) {
		expr = expr - min(expr, na.rm=T)
	    }
	    expr = log2(expr+1)
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
	    #print(head(gene.mapping, 20))
	    if(!is.null(conversion.mapping.function)) {
		gene.mapping$Gene = unlist(lapply(as.character(gene.mapping$Gene), conversion.mapping.function))
	    }
	    #print(head(gene.mapping, 20)) 
	    if(!is.null(conversion.mapping)) {
		# Likely to convert accession numbers to geneids - get rid of version (trailing dots and digits)
		gene.mapping = data.frame(Probe = gene.mapping$Probe, Gene = as.character(conversion.mapping[sub("\\.[0-9]+", "", gene.mapping$Gene)]))
		gene.mapping$Gene = factor(gene.mapping$Gene, levels=c(levels(gene.mapping$Gene), ""))
		gene.mapping[gene.mapping$Gene %in% c("NULL", "---", "-", " "), "Gene"] = ""
		#gene.mapping[gene.mapping$Gene == "NULL", "Gene"] = NA
		#gene.mapping = na.omit(gene.mapping) 
	    }
	    #print(head(gene.mapping, 20)) 
	    if(length(levels(factor(gene.mapping$Gene))) == 1) {
		stop("Problem with gene mapping!")
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

