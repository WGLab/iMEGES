###########

library(pcaMethods)

TP_score <- read.table("Final_variants_score.txt",header=F)

dat1 <- data.frame(TP_score)
dat1[dat1=='.'] <- NA

dat1[,15][is.na(dat1[,15])] <- 0.01

dat2 =dat1[rowSums(is.na(dat1[,10:14]))!=5,]

dat3 <- apply(as.matrix(dat2[, 10:20]), 2, as.numeric)

dat4 <- dat3[rowSums(is.na(dat3))!=5, ]

dat1_pca <- pca(dat4, nPcs=3, method="bpca")


dat1_imputed <- completeObs(dat1_pca)
dat1_imputed <- dat1_imputed[,1:11]
labels <- rep(0,dim(dat1_imputed)[[1]])

Final_lables1 <- cbind(dat1_imputed,labels)

Final_lables2  <- cbind(dat2[,c(2,4:9)],Final_lables1)

colnames(Final_lables2) <- c("Chr", "Start", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "Eigen", "CADD","dann", "GWAVA_unmatched_score", "FATHMM_noncoding", "gnomAD_genome_ALL","Brain.known.eQTLs", "H3K27Ac", "H3K27me3","H3K4me1","H3K4me3","labels")

Final_lables2[,7][is.na(Final_lables2[,7])] <- "."

write.table(Final_lables2, file="SCHI_TN.txt", col.names=TRUE, row.names=F,  quote = F, sep = "\t", append=F)

