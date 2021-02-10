
library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(ggplot2)

# cell #1
dataset <- read.table("readonly/heart-h.txt", header = TRUE, sep=",", quote="")
class(dataset)

# cell #2
nrow(dataset) # number of rwos
ncol(dataset) # number of columns
dim(dataset) # number of rows and columns
head(dataset, 0) # returns the list of variables
summary(dataset) # descriptive statistics

# cell #3
md.pattern(dataset)  # (mice) displays all the missing values, NA for missing values

# cell #4
mice_plot <- aggr(dataset, col=c("green","red"),
                    numbers=TRUE, sortVars=TRUE,
                    labels=names(dataset), cex.axis=.7, ylab=c("Missing data","Pattern")) # (VIM) display graphically missing values

# cell #5
dataset[,"chol_imputed"] <- with(dataset, impute(chol, mean))  # (Hmisc) impute with mean
summary(dataset)

# cell #6
#dataset_n <- Normalize(chol_imputed ~., data = dataset) # (RWeka) normalizes all numeric except chol_imputed
#dataset_n <- Normalize(~ chol_imputed , data = dataset) # normalizes only chol_imputed
dataset_n <- Normalize(~. , data = dataset) # normalizes all variables
summary(dataset_n)

# cell #7
dataset$chol_bin <- as.numeric(cut2(dataset$chol_imputed, g=3)) # (Hmisc) create 3 bins g = 3 quantile bins equal-depth
summary(dataset)
#head(dataset)

# cell #8
dataset[,"chol_bin"] <- as.numeric(cut(dataset[,"chol_imputed"], 5)) # create bins of same width
summary(dataset)
#head(dataset)

# cell #9
ggplot(dataset,aes(x=age,y=chol, color=num)) + geom_point(size = 4) # scatterplot (age, chol)

# cell #10
ggplot(dataset, aes(chest_pain, fill=factor(num))) + geom_bar() # stacked histogram for categorical variable

# cell #11
ggplot(dataset, aes(x=age, fill=factor(num))) + geom_bar() # stacked histogram for numeric variable


