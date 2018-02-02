
#' Convert sample mapping to case control from conditions
#'
#' @param sample.mapping Sample-condition mapping for the data set.
#' @param states.case Conditions to be used as case
#' @param states.control Conditions to be used as background
#' @param out.file File name to write the converted mapping
#' @return a data frame containin sample mapping
#' @export
#' @examples
#' #sample.mapping = convert.sample.mapping.to.case.control(sample.mapping,
#' # 		    states.control = c("healthy donor"), 
#' #		    states.case = c("tuberculosis", "latent tuberculosis infection"))
convert.sample.mapping.to.case.control<-function(sample.mapping, states.case, states.control, out.file=NULL) {
    #sample.mapping2 = data.frame(sample=sample.mapping$sample, type=sapply(sample.mapping$type, as.character), stringsAsFactors=F)
    sample.mapping2 = data.frame(sample=sample.mapping$sample, type=NA, stringsAsFactors=F)
    for(state in states.case) {
	sample.mapping2[sample.mapping$type == state, "type"] = "case"
    }
    for(state in states.control) {
	sample.mapping2[sample.mapping$type == state, "type"] = "control"
    }
    sample.mapping2 = na.omit(sample.mapping2)
    sample.mapping2$type = factor(sample.mapping2$type)
    if(!is.null(out.file)) {
	write.table(sample.mapping2, out.file, quote=F, row.names=F, sep="\t")
    }
    return(sample.mapping2)
}

