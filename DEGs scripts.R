#set working directory

setwd("D:/Project 1st/R program/TB/TB2/HIV")
getwd()

#called library
#install.packages("DESeq2")
#BiocManager::install("DESeq2")
library("DESeq2")


#Read csv file
Count_data = read.table(file = "counts_HIV.csv", header = T, sep = ",", row.names = 1)

dim(Count_data)  #dimension of dataframe
dim(Count_data)  #dimension of dataframe



Col_data <- read.table(file = "HIV_col_data.csv", header = T, sep = ",",row.names = 1)
dim(Col_data)  #dimensions of dataframe

all(rownames(Col_data)==colnames(Count_data)) #check the row and col name equal or not
#count the number of NA value in the matrix
class(Count_data)
class(Col_data)
which(is.na(Count_data),arr.ind=TRUE)

# --------- IF there is no NA value then you don't have run this block of code----------

#install.packages("zoo")
library(zoo)
Count_data[]<-t(na.aggregate(t(Count_data))) #for replacing all Na value
Count_data [is.na(Count_data)] = 0

which(is.na(Count_data),arr.ind=TRUE) #recheck 

# ---------------------------------------------------------------- 

dds = DESeqDataSetFromMatrix(countData = round(Count_data),
                             colData = Col_data,
                             design = ~ condition) # we're testing for the different condidtion

dds$condition <- relevel(dds$condition, ref = 'Control') #relevel commend for define the reference parameter
dds
dds <- DESeq(dds) #it will normalize the data in background and calculate the fold change foe each gene
res1 <- results(dds) # get the result 
summary(res1) #summary of the result
res1




###   keep only sig results, padj<0.05 and log2FoldChange >1
resSigUp <- subset(res1, pvalue < 0.05 & log2FoldChange >1)
summary(resSigUp)
write.csv(resSigUp, "Upregulated_HIV.csv")

###keep only sig results, padj<0.05 and log2FoldChange < -1
resSigDown <- subset(res1, pvalue < 0.05 & log2FoldChange < -1)
summary(resSigDown)
write.csv(resSigDown, "downgulated_HIV.csv")

###keep UP and Down in one file with padj<0.05

resSig <- subset(res1, log2FoldChange >1 & pvalue <0.05 | log2FoldChange < -1 & pvalue < 0.05)
summary(resSig)
write.csv(resSig, "DE_HIV.csv")

##Write for volcano plot
##write.table(res1, file="Volcano", row.names=F, sep="\t")

#install library
#BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)

EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                title = "",
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )

