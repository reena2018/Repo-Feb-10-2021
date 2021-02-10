
# cell #1
mrnaNorm <- read.table("/data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
            header = F, fill = T, skip = 2)
class(mrnaNorm)

# cell #2
mrnaIDs <- read.table("/data/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
            header = F, fill = T, nrows = 1)
class(mrnaIDs)
dim(mrnaIDs)
mrnaIDs <- mrnaIDs[, -1][, -1]  # remove the first solumn "Hybridization REF" which spans two columns because of the space between the two words
dim(mrnaIDs)

# cell #3
dim(mrnaNorm)
dim(mrnaIDs)
mrnaIDs5 = mrnaIDs[,c(1,2,3,4,5)] # extract first 5 patientIDs
mrnaNorm5x5 = mrnaNorm[1:5, 1:5] # first 5 rows and columns 
head(mrnaIDs5, 2) # display first two rows
head(mrnaNorm, 2) # display first two rows
summary(mrnaNorm5x5) # summary statistics

# cell #4
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
# extract the sample type (tumor / normal), excluding the first column since it is a gene ID
sampleType <- as.data.frame(samp)
dim(sampleType)
# 1 1212

# cell #5
unique(sampleType)
# extracts how many unique objects there are 

tab <- table(unlist(sampleType))
tab
# count how many of each type
#   01   06   11 
#  1093   7  112 

# cell #6
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaClass)
table(unlist(sampClass))
#   0    1 
# 112 1100 

# cell #7
sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum) 

# cell #8
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

# cell #9
bss <- bssWssFast(t(mrnaNorm[, -1]), t(mrnaClassNum), 2)
# returns a list with two elements: a list of bss/wss ratios in decreasing order bss[1] and a list of features bss[2]

bss$ix[1:51]
# show list of 50 ranked gene indexes

# cell #10
genes <- mrnaNorm[bss$ix[1:51],1]
genes
# show list of 50 tanked gene names
# [1] FHL1|2273        ADAMTS4|9507     ABAT|18          PAG1|55824       ARHGAP1|392     
# [6] TP53TG3B|729355  SCRN3|79634      TIE1|7075        SNRPA|6626       SCCPDH|51097    
#[11] PCYT1B|9468      DSE|29940        LGI4|163175      LIX1L|128077     LCE3D|84648     
#[16] POTEE|445582     GIMAP2|26157     CABYR|26256      PLGLA|285189     LRRC57|255252   
#[21] C21orf33|8209    CDC14A|8556      CD300LF|146722   CEP57|9702       MAGEC2|51438    
#[26] LGALS3|3958      CAT|847          LOC647979|647979 ABCA8|10351      PPIL3|53938     
#[31] EFCAB1|79645     CCDC45|90799     ABCA5|23461      CAV1|857         TP53I3|9540     
#[36] ADRB1|153        SOD1|6647        LDLRAP1|26119    LYPLA1|10434     KLHL10|317719   
#[41] HLA-DPB1|3115    SGCG|6445        MLL2|8085        FGD4|121512      CD33|945        
#[46] BMS1|9790        MARCKSL1|65108   BTNL8|79908      SS18L2|51188     ABCA7|10347   


# cell #11
# merge class values into mrnaNorm as the last row  - after the genes and select only top 100 genes
mrnaSetClass <- rbind(mrnaNorm[bss$ix[1:100],-1], setNames(mrnaClassNum, names(mrnaNorm[,-1])))
dim(mrnaSetClass)
# 101  1212       101 genes as rows and 1212 patients as columns
# to select by class, transpose the matrix
transpose <- as.data.frame(t(mrnaSetClass))
colnames(transpose)[101] <- "class"
# select normals
normals <- subset(transpose, transpose$class == 0)
dim(normals)
#[1] 112 101
# select tumors
tumors <- subset(transpose, transpose$class == 1)
dim(tumors)
#[1] 1100 101
# combine them together
both <- rbind(normals, tumors[1:112,])
dim(both)
# select a color scheme - red for tumor and blue for normal
color.map <- function(class) { if (class==0) "#FF0000" else "#0000FF" } # red for tumor and blue for normal
# create heatmap
groupColors <- unlist(lapply(both$class, color.map))
heatmap(as.matrix(t(both)), scale = "row", col=topo.colors(100), ColSideColors=groupColors)


