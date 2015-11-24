
#' Get expression data set 
#'
#' @param geo.id GEO id.
#' @param output.dir Output directory to write / look for files.
#' @return data set 
get.data.set<-function(geo.id, output.dir) {
    data.set = NULL
    file.name = paste(output.dir, geo.id, ".annot", sep="")
    file.name2 = paste(output.dir, geo.id, ".soft", sep="")
    if(file.exists(file.name)) {
	data.set = GEOquery::getGEO(filename=file.name)
    } else if(file.exists(file.name2)) {
	data.set = GEOquery::getGEO(filename=file.name2)
    } else {
	file.name = paste(output.dir, geo.id, ".annot.gz", sep="")
	file.name2 = paste(output.dir, geo.id, ".soft.gz", sep="")
	if(file.exists(file.name)) {
	    data.set = GEOquery::getGEO(filename=file.name)
	} else if(file.exists(file.name2)) {
	    data.set = GEOquery::getGEO(filename=file.name2)
	} else {
	    # Get GEO file
	    if(substr(geo.id,1,3) == "GSE") {
		data.set = GEOquery::getGEO(geo.id, destdir=output.dir, GSEMatrix=T)
	    } else {
		data.set = GEOquery::getGEO(geo.id, destdir=output.dir, AnnotGPL=T)
	    }
	}
    }
    return(data.set)
}


