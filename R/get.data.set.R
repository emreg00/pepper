
#' Get expression data set 
#'
#' @param geo.id GEO id.
#' @param output.dir Output directory to write / look for files.
#' @return data set 
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
		file.name = paste0(output.dir, geo.id, ".", suffix) 
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
	} else {
	    data.set = GEOquery::getGEO(geo.id, destdir=output.dir, GSEMatrix=T, AnnotGPL=T)
	}
    }
    return(data.set)
}


