
# cell #2
memory.limit(size=3500)
library(cluster)
library(ReactomePA)

# cell #3
mrnaNorm <- read.table("/data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
            header = F, fill = T, skip = 2)
mrnaIDs <- read.table("/data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
            header = F, fill = T, nrows = 1)
mrnaIDs <- mrnaIDs[, -1][, -1]

# cell #3
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
sampleType <- as.data.frame(samp)
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaNorm)
# 20531 1213 columns are patients (except the 1st for gene name) rows are expression levels for each gene
dim(mrnaIDs)
# 1 1213   the first column is the gene name, the others are one patient per row
dim(mrnaClass)
# 1 1212 one patients per row   1 = tumor, 0 = normal
table(unlist(sampClass))
#   0    1 
# 112 1100                     112 normals and 1100 tumor
sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum) 

# cell #4
geneNames <- mrnaNorm[1] # extract the gene names from mrnaNorm as its first column
dim(geneNames)
# 20531 genes

# cell #5
mrnaData = t(mrnaNorm[, -1]) # remove first column of mrnaData and transpose it to have genes as columns
rm(samp)
rm(sampClass)
rm(mrnaNorm)
gc()

# cell #6
bssWssFast <- function (X, givenClassArr, numClass=2)
# between squares / within square feature selection
{
	classVec <- matrix(0, numClass, length(givenClassArr))
	for (k in 1:numClass) {
		temp <- rep(0, length(givenClassArr))
		temp[givenClassArr == (k - 1)] <- 1
		classVec[k, ] <- temp
	}
	classMeanArr <- rep(0, numClass)
	ratio <- rep(0, ncol(X))
	for (j in 1:ncol(X)) {
		overallMean <- sum(X[, j]) / length(X[, j])
		for (k in 1:numClass) {
			classMeanArr[k] <- 
				sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
		}
	  classMeanVec <- classMeanArr[givenClassArr + 1]
	  bss <- sum((classMeanVec - overallMean)^2)
	  wss <- sum((X[, j] - classMeanVec)^2)
	  ratio[j] <- bss/wss
	}
      sort(ratio, decreasing = TRUE, index = TRUE)
}

# cell #7
# select features
dim(mrnaData)
# 1212 20531  matrix
dim(mrnaClass)
# 1 1212
dim(mrnaClassNum)
# 1 1212
dim(geneNames)
# 20531 genes
bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)
mrnaDataReduced <- mrnaData[,bss$ix[1:100]]
dim(mrnaDataReduced)
# 1212  100
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

# cell #8
set.seed(20)
kmeans.clusters <- kmeans(mrnaDataReduced, 2, nstart = 20)
table(kmeans.clusters$cluster)

# cell #9
table(kmeans.clusters$cluster, trainClasses)

# cell #10
clusplot(mrnaDataReduced, kmeans.clusters$cluster, color=TRUE, shade=TRUE, 
  	labels=2, lines=0)

# cell #11
pam.clusters <- pam(mrnaDataReduced, 2)
table(pam.clusters$clustering)

# cell #12
table(pam.clusters$clustering, trainClasses)

# cell #13
clusplot(mrnaDataReduced, pam.clusters$clustering, color=TRUE, shade=TRUE, 
  	labels=2, lines=0)

# cell #14
genes <- geneNames[bss$ix[1:100],1]
genes <- as.character(genes)
hugoNames <- lapply(genes, function(t) substr(t, 1, regexpr("\\|", t) - 1))   # extract HUGO names
entrezNames <- lapply(genes, function(t) substr(t, regexpr("\\|", t) + 1, nchar(t)))  # extract Entrez names
paths = enrichPathway(unlist(entrezNames), pvalueCutoff=1)
head(summary(paths))

# cell #15
enrichMap(paths, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
