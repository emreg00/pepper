
#' Get array platform annotation for probes
#'
#' @param data.set GEO data set.
#' @param probe.conversion Id of the column used for probe gene mapping.
#' @param output.dir Directory to store files.
#' @return Data frame containing gene mapping.
get.platform.annotation<-function(data.set, probe.conversion, output.dir) {
    gds.id = GEOquery::Meta(data.set)$platform
    data.set = get.data.set(gds.id, output.dir)
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


