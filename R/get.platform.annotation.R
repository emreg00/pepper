
#' Get array platform annotation for probes
#'
#' @param geo.id GEO id.
#' @param probe.conversion Id of the column used for probe gene mapping.
#' @param output.dir Directory to store files.
#' @return Data frame containing gene mapping.
#' @keywords internal
get.platform.annotation<-function(geo.id, probe.conversion, output.dir) {
    data.set = get.data.set(geo.id, output.dir, is.annotation=T)
    d = GEOquery::Table(data.set)
    #print(colnames(d)) 
    #print(probe.conversion)
    #print(length(d[,probe.conversion]))
    gene.mapping = data.frame(Probe = d[,"ID"], Gene = d[,probe.conversion]) #row.names = d[,"ID"],  "Gene.Symbol",
    return(gene.mapping)
    # Alternative using bioconductor annotation packages
    #source("http://www.bioconductor.org/biocLite.R")
    #biocLite("hgu133a")
    #mget("121_at",hgu133aSYMBOL)
    #mget("121_at",hgu133aUNIGENE)
}


