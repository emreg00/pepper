
#' Get expression data set 
#'
#' @param gds.id GEO data set id.
#' @param output.dir Output directory to write / look for files.
#' @return data set 
get.data.set<-function(gds.id, output.dir) {
    data.set = NULL
    file.name = paste(output.dir, gds.id, ".annot.gz", sep="")
    file.name2 = paste(output.dir, gds.id, ".soft.gz", sep="")
    if(file.exists(file.name)) {
	data.set = GEOquery::getGEO(filename=file.name)
    } else if(file.exists(file.name2)) {
	data.set = GEOquery::getGEO(filename=file.name2)
    } else {
	file.name = paste(output.dir, gds.id, ".annot", sep="")
	file.name2 = paste(output.dir, gds.id, ".soft", sep="")
	if(file.exists(file.name)) {
	    data.set = GEOquery::getGEO(filename=file.name)
	} else if(file.exists(file.name2)) {
	    data.set = GEOquery::getGEO(filename=file.name2)
	} else {
	    # Get GDS file
	    data.set = GEOquery::getGEO(gds.id, destdir=output.dir, AnnotGPL=T)
	}
    }
    return(data.set)
}


