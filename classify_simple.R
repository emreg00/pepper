library(ROCR)
library(caret)
library(plyr)
library(PEPPER)

# GLOBAL definitions
data.dir = "../data/"
seed = 142341
set.seed(seed) 
adjust.method = 'BH' # "bonferroni"
fdr.cutoff = 0.05
z.cutoff = 2.5
n.pool.cutoff = NULL
use.fdr.cutoff = T 
n.fold = 5 
k.n.n = 10 
n.repetition = 100


main<-function() {
    classify("GSE17755", "RepeatedCVWithinBatch", "peep") # "LOOCVWithinBatch"
}


classify<-function(data.set, method, classifier) { 
	if(classifier == "peep") {
	    # Peep gene pool overlap based scores
	    get.fraction.based.classification(data.set, c(data.set), method) 
	} else if(classifier == "knn") {
	    # k-NN based approach
	    get.kNN.based.classification(data.set, c(data.set), method) 
	}
}


get.fraction.based.classification<-function(data.set, data.sets, method="") {
    output.dir = paste0(data.dir, data.set, "/")
    file.name = paste0(output.dir, "mapping_case.dat")
    mapping = read.csv(file.name, sep="\t")
    file.name = paste0(output.dir, "expr.dat") 
    expr = read.csv(file.name, sep="\t", check.names=F)
    # Get number of background genes
    n.gene = NULL
    if(grepl("WithinBatch$", method)) {
	if(is.null(n.pool.cutoff) | grepl("^RepeatedCV", method)) {
	    n.gene = nrow(expr)
	}
	data.sets = c(data.set)
    }
    if(use.fdr.cutoff) {
	cutoff = fdr.cutoff
    } else {
	cutoff = z.cutoff
    }

    if(method == "LOOCVWithinBatch") {
	n.fold = -1
	n.repetition = 1
    }
    if(n.fold == -1) { # leave-one-out
	n.fold = nrow(mapping)
	print(n.fold) 
    }
    mapping$flag = ifelse(mapping$type == "case", 1, 0)
    d = data.frame(sample = 1:(n.repetition * nrow(mapping)), score=NA, flag=NA, fold=NA, rep=NA)
    i = 1
    #values = vector("list", n.repetition)
    values.pool.cutoff = c()
    for(repetition in 1:n.repetition) { # Repeat n.repetition times
	cv = createFolds(factor(mapping$flag), k = n.fold, list = TRUE, returnTrain=F) 
	#scores = vector("list", n.fold) 
	#labels = vector("list", n.fold)
	for(k in 1:length(cv)) { 
	    indices = cv[[k]]
	    samples.to.exclude = mapping[indices, "sample"]
	    #print(sprintf("%s %d", samples.to.exclude, mapping[indices, "flag"])) 
	    e = get.peeps(expr, mapping, samples.to.exclude, cutoff=cutoff, convert.to.pvalues=use.fdr.cutoff)
	    vals = get.pool(e, n.cutoff=n.pool.cutoff, n.gene=n.gene)
	    pool = vals$pool
	    n.cutoff = vals$n.cutoff
	    values.pool.cutoff = c(values.pool.cutoff, n.cutoff)
	    #scores.inner = c()
	    #labels.inner = c()
	    for(sample in samples.to.exclude) {
		peep = e[e$sample == sample, "geneid"]
		score = length(intersect(peep, pool)) / length(pool)
		#print(c(length(intersect(peep, pool)), length(peep), length(pool))) 
		if(is.na(score)) {
		    score = 0
		}
		label = mapping[mapping$sample == sample, "flag"]
		d[i, "sample"] = sample
		d[i, "score"] = score
		d[i, "flag"] = label
		d[i, "fold"] = k
		d[i, "rep"] = repetition
		#print(c(sample, d[i,]))
		i = i + 1
		#scores.inner = c(scores.inner, score)
		#labels.inner = c(labels.inner, label)
	    }
	    #scores[[k]] = scores.inner
	    #labels[[k]] = labels.inner
	}
	#print(scores)
	#print(labels) 
	#pred = prediction(scores, labels)
	#perf = performance(pred, "auc")
	#a = 100*perf@y.values[[1]]
	#values[[repetition]] = a
    }
    #d = na.omit(d)
    #print(rbind(head(d), tail(d)))
    if(n.fold != nrow(mapping)) {
	# AUC per cv per repetition 
	a = ddply(d, .(rep, fold), function(x) { pred=prediction(x$score, x$flag); perf=performance(pred, "auc"); 100*perf@y.values[[1]] })
	values = ddply(a, .(rep), function(x) { mean(x$V1) } )
	#print(values)
	print(sprintf("AUC: %.1f (std: %.1f)", mean(values$V1), sd(values$V1)))
    } else {
	# AUC per repetition over folds jointly 
	values = ddply(d, .(rep), function(x) { pred=prediction(x$score, x$flag); perf=performance(pred, "auc"); 100*perf@y.values[[1]] })
	#print(values)
	print(sprintf("AUC: %.1f", mean(values$V1)))
    }
    print(table(values.pool.cutoff)) 
    
    if(use.fdr.cutoff) { 
	out.file = paste0(output.dir, "classifier_fraction", method, "_fdr", cutoff, "_n", n.pool.cutoff, ".dat")
    } else {
	out.file = paste0(output.dir, "classifier_fraction", method, "_z", cutoff, "_n", n.pool.cutoff, ".dat")
    }
    write.table(d, out.file, sep="\t", quote=F, row.names=F)
    return(mean(values$V1))
}


get.kNN.based.classification<-function(data.set, data.sets, method="") {
    output.dir = paste0(data.dir, data.set, "/")
    file.name = paste0(output.dir, "mapping_case.dat")
    mapping = read.csv(file.name, sep="\t")
    file.name = paste0(output.dir, "expr.dat") 
    expr = read.csv(file.name, sep="\t", check.names=F)
    if(grepl("WithinBatch$", method)) {
	data.sets = c(data.set)
    } else {
	stop("Not implemented!")
    }

    if(method == "LOOCVWithinBatch") {
	n.fold = -1
	n.repetition = 1
    }
    if(n.fold == -1) { # leave-one-out
	n.fold = nrow(mapping)
	print(n.fold) 
    }
    mapping$flag = ifelse(mapping$type == "case", 1, 0)
    d = data.frame(sample = 1:(n.repetition * nrow(mapping)), score=NA, flag=NA, fold=NA, rep=NA)
    i = 1
    for(repetition in 1:n.repetition) { # Repeat n.repetition times
	cv = createFolds(factor(mapping$flag), k = n.fold, list = TRUE, returnTrain=F) 
	for(k in 1:length(cv)) { 
	    indices = cv[[k]]
	    samples.to.exclude = mapping[indices, "sample"]
	    #print(sprintf("%s %d", samples.to.exclude, mapping[indices, "flag"])) 
	    samples = mapping[-indices, "sample"]
	    for(sample in samples.to.exclude) {
		values = c()
		values.label = c() #ifelse(mapping[-indices, "flag"]==1, 1, 0)
		for(sample2 in samples) { 
		    if(sample != sample2) {
			values = c(values, cor(expr[,sample], expr[, sample2]))
			values.label = c(values.label, mapping[mapping$sample == sample2, "flag"]) #ifelse(mapping[mapping$sample == sample2, "type"] == "case", 1, 0)) 
		    }
		}
		idx = rev(order(values))[1:k.n.n]
		score = sum(values[idx] * values.label[idx])
		d[d$sample == sample, "score"] = score
		label = mapping[mapping$sample == sample, "flag"]
		d[i, "sample"] = sample
		d[i, "score"] = score
		d[i, "flag"] = label
		d[i, "fold"] = k
		d[i, "rep"] = repetition
		#print(c(sample, d[i,]))
		i = i + 1
	    }
	}
    }
    #d = na.omit(d)
    if(n.fold != nrow(mapping)) {
	# AUC per cv per repetition 
	a = ddply(d, .(rep, fold), function(x) { pred=prediction(x$score, x$flag); perf=performance(pred, "auc"); 100*perf@y.values[[1]] })
	values = ddply(a, .(rep), function(x) { mean(x$V1) } )
	#print(values)
	print(sprintf("AUC: %.1f (std: %.1f)", mean(values$V1), sd(values$V1)))
    } else {
	# AUC per repetition over folds jointly 
	values = ddply(d, .(rep), function(x) { pred=prediction(x$score, x$flag); perf=performance(pred, "auc"); 100*perf@y.values[[1]] })
	#print(values)
	print(sprintf("AUC: %.1f", mean(values$V1)))
    }

    out.file = paste0(output.dir, "classifier_knn", method, "_k", k.n.n, ".dat")
    write.table(d, out.file, sep="\t", quote=F, row.names=F)
}


get.peeps<-function(expr, sample.mapping, samples.to.exclude, cutoff, convert.to.pvalues) {
    # Get case/control samples
    samples = as.vector(sample.mapping[sample.mapping$type == "case", "sample"])
    samples.background = as.vector(sample.mapping[sample.mapping$type == "control", "sample"])
    #expr = expr[,c(samples, samples.background)]
    z = get.z.score(expr, samples.background, method="mean", samples.to.exclude=samples.to.exclude)
    # Ignoring NAs (i.e. genes with 0 sd) 
    #z = na.omit(z)
    e = get.peeps.from.z.matrix(z, cutoff, convert.to.pvalues) 
    e$flag = ifelse(e$sample %in% samples, 1, 0)
    e$exclude = ifelse(e$sample %in% samples.to.exclude, 1, 0)
    return(e) # sample geneid flag exclude
}


get.pool<-function(e, n.cutoff=NULL, n.gene=NULL) {
    e = e[e$flag==1 & e$exclude==0,]
    # Get n.cutoff (number of samples) based on binomial model
    if(is.null(n.cutoff)) {
	if(is.null(n.gene)) {
	    stop("For calculating pool n.cutoff, n.gene in the data set should be provided!")
	}
	n.cutoff = get.pool.cutoff(e, n.gene)
	#print(sprintf("n.cutoff: %d", n.cutoff)) 
    }
    f = count(e, "geneid")
    pool = as.vector(f[f$freq >= n.cutoff,"geneid"])
    return(list(pool=pool, n.cutoff=n.cutoff))
}


get.pool.cutoff<-function(e, n.gene) {
    g = count(e, "sample")
    n.sample = length(levels(factor(g$sample)))
    prob = mean(as.vector(g$freq)) / n.gene
    #print(as.vector(g$freq))
    #print(sprintf("n.sample: %d n.peep: %.1f n.gene: %d prob: %e", n.sample, mean(g$freq), n.gene, prob)) 
    a = 1-cumsum(dbinom(0:n.sample, n.sample, prob))
    n.cutoff = min(which(a < 1/n.gene))
    #print(a) 
    #n.cutoff = n.cutoff + 1 
    return(n.cutoff)
}


main()

