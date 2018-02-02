
#' Get expression data set 
#'
#' @param geo.id GEO id.
#' @param output.dir Output directory to write / look for files.
#' @return data set 
#' @export
get.data.set<-function(geo.id, output.dir, is.annotation=F) {
    data.set = NULL
    #a = list.files(path = output.dir, pattern = geo.id)
    #if(length(a) > 0) {
    #	file.name = paste0(output.dir, a[1])
    #}
    #data.set = GEOquery::getGEO(filename=file.name) # does not get annotation info
    if(is.annotation) {
	for(suffix in c("annot", "soft", "annot.gz", "soft.gz")) {
	    if(is.null(data.set)) {
		file.name = file.path(output.dir, paste0(geo.id, ".", suffix))
		if(file.exists(file.name)) {
		    data.set = GEOquery::getGEO(filename=file.name)
		} 
	    } else {
		break()
	    }
	}
    }
    if(is.null(data.set)) {
	# Get GEO file
	if(substr(geo.id,1,3) == "GSE") {
	    data.set = GEOquery::getGEO(geo.id, destdir=output.dir, GSEMatrix=T)
	} else if(substr(geo.id,1,3) == "GDS") {
	    data.set = GEOquery::getGEO(geo.id, destdir=output.dir, GSEMatrix=T, AnnotGPL=T)
	} else { # Assumes ArrayExpress file
	    if (!requireNamespace("ArrayExpress", quietly=TRUE)) {
		stop("Retreiving ArrayExpress data sets requires ArrayExpress package to be installed!")
	    }
	    data.set = ArrayExpress::getAE(geo.id, path=output.dir, type="full", local=T)
	}
    }
    return(data.set)
}


