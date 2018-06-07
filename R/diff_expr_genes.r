#Differentially expressed genes with EdgeR

library(edgeR)
library(locfit)
library(splines)

#Working directory
setwd("C:/Users/yhueber.CGIARAD/Desktop/RNAseq_new_genotypes/")
#Filename
x <- read.delim("all_indiv_RNAseq.counts.txt", row.names="Gene")


#xx <- x[,49:54]
xx <- x[, c('C1','C2','C3','GM47','GM46','GM54')]

#xx <- x[,c(1, 4:8, 10)]
#xx <- x[, c('C1','C2','C3','C4','C5','C6')]

#gros michel 'GM47','GM46','GM54','GM56','GM58','GM61'
#poyo 'P73','P87','P89','P91','P93','P94'
#red dacca 'RD89','RD93','RD94','RD79','RD82','RD91'

group <- c(1,1,1,2,2,2) #biological replicates
y <- DGEList(counts=xx, group=group)
colnames(y)
dim(y) #total number of unique tags
y$samples #number of tags ("genes") per sample

#Filtering
keep <- rowSums(cpm(y)>1) >= 3

#keep <- rowSums(cpm(y)>1) >= 3

#y <- y[keep,]
#dim(y) #total number of unique tags
#y$samples$lib.size <- colSums(y$counts) #recompute library sizes

y <- y[keep, , keep.lib.sizes=FALSE]
y$samples #number of tags ("genes") per sample

dim(y) #total number of unique tags

#Normalization
y <- calcNormFactors(y, method="RLE")
y$samples #number of tags ("genes") per sample

#plotMDS(y)

#design
#design <- model.matrix(~group)

#Pairwise comparisons
#1. Estimating dispersions
y <- estimateDisp(y) #classic mode
#y <- estimateDisp(y,design)

#2. Test for DE (differentially expressed)
et <- exactTest(y, pair=c(1,2))
topTags(et)

#total number of genes at 5% FDR
#summary(de <- decideTestsDGE(et, p.value=0.05, adjust.method="none")) 
summary(de <- decideTestsDGE(et, p.value=0.05, adjust.method="BH"))
#summary(de <- decideTestsDGE(et, p.value=0.05, adjust.method="fdr"))
#summary(de <- decideTestsDGE(et, p.value=0.05, adjust.method="BY"))
#summary(de <- decideTestsDGE(et, p.value=0.05, adjust.method="holm"))

#plot smear
#detags <- rownames(y)[as.logical(de)]
#plotSmear(et, de.tags=detags)
#abline(h=c(-1, 1), col="blue")

#Print list of DE genes with pvalue=0.05
top <- topTags(et, n=8375)
write.table(top,file="genesDE_edgR_v3.12.1_RLE_Cachaco_vs_Gros_Michel_control_vs_control_FDR05.txt",col.names=TRUE,sep="\t",append=FALSE)




