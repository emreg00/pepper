# PePPeR - PErsonalized Perturbation ProfilER

PePPeR (Personalized Perturbation ProfileR) is an R package providing methods to fetch expression data sets from the GEO database, identify per-individual or group-wise differentially expressed (DE) genes and construct individual perturbation profiles. 

## Requirements

Depends:
-    R (>= 3.0.2)

Imports:
-    Biobase
-    GEOquery

Suggests:
-    limma
-    affy
-    samr
-    preprocessCore
-    makecdfenv


## Installation

### Directly from GitHub
You can use `install_github()` method in [devtools](https://github.com/hadley/devtools) package

```R
install.packages("devtools")
library(devtools)
install_github("emreg00/pepper", dependencies = T)
```


### From an archieve file 
Create a tarball containing the files in the repository

```bash
tar cvzf pepper.tgz --exclude .git pepper/
```

Install it using R

```R
R CMD INSTALL pepper.tgz
```

## Usage

```R
> library(PEPPER)
> 
```

See the documentation on the following functions for their use (?function.name):
- `fetch.expression.data`: for fetching information from NCBI GEO data base 
- `find.de.genes`: for getting group-wise DE genes
- `get.z.matrix` and `get.z.score` and `get.peeps.from.z.matrix`: for getting per-individual DE genes

See `classify.simple.R` function for an example on how to use PeePs for sample classification.


## Case study: Getting DE genes in a reprocessed GEO Parkinson data set (GSE7621)

- Fetch the Parkinson data set available at [GEO](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7621), reprocess it using affy package and map probes to gene ids

```R
# Data set specific parameters
geo.id <- "GSE7621" # GEO id of the data set
probe.conversion <- "ENTREZ_GENE_ID" # column name for gene id mapping
conversion.map <- NULL # probe to gene mapping, if NULL uses the mapping in the data set
conversion.mapping.function <- NULL # modify probe names using this function 
sample.mapping.column <- "characteristics_ch1" # column to use for sample mapping
geo.id.sub <- NULL # the platform to use if there are multiple platform annotations
reprocess <- "affy" # reprocessing type for raw data
output.dir <- "./"

# Get the expression and sample mapping info from the reprocessed data set
d <- fetch.expression.data(geo.id, sample.mapping.column = sample.mapping.column, do.log2 = NULL, probe.conversion = probe.conversion, conversion.map = conversion.map, conversion.mapping.function = conversion.mapping.function, output.dir = output.dir, geo.id.sub = geo.id.sub, reprocess = reprocess)
expr <- d$expr
sample.mapping <- d$sample.mapping
```

- Get group-wise DE
```R
states.case <- c("Parkinson's Disease")
states.control <- c("Old Control") 
out.file <- "case_mapping.dat"
sample.mapping <- convert.sample.mapping.to.case.control(sample.mapping, states.case, states.control, out.file = out.file)
adjust.method <- 'BH'
fdr.cutoff <- 0.05
out.file <- "de.dat"
de <- find.de.genes(expr, sample.mapping, c("case", "control"), method="limma", out.file, adjust.method=adjust.method, cutoff=fdr.cutoff, functional.enrichment="kegg") 
de <- de[abs(de$logFC)>=1,]
```

- Get per-individual DE
```R
# Get z scores
out.file <- "z.dat"
cutoff <- 2.5
z = get.z.matrix(expr, sample.mapping, method="mean", out.file=out.file)
indices <- apply(abs(z), 2, function(x) { which(x >= cutoff)})
geneids <- lapply(indices, function(x) { rownames(z)[x] })
# Alteratively you can use get.peeps.from.z.matrix
peeps <- get.peeps.from.z.matrix(z, cutoff=2, convert.to.pvalues=F) 
```


## Citation

Menche J et al., Integrating personalized gene expression profiles into predictive disease-associated gene pools. Npj Systems Biology and Applications 2017;3:10 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28649437)

