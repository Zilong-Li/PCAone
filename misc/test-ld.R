## creat a fake assoc File
## CHR BP A1 A2 P


bim <- read.table("example/merge.bim")
dim(bim)

ord <- sort(sample(nrow(bim), 0.02*nrow(bim)))
out <- bim[ord,c(1,4,5,6)]
colnames(out) <- c("CHR", "BP", "A1", "A2")
out$P <- runif(nrow(out), 0, 0.1)

sum(out$P < 0.0001)
head(out)

write.table(out, file = "example/clump.assoc.file", quote = F, row.names=F, sep = "\t")
