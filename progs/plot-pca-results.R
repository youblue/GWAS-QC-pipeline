setwd("${OUTPUT_DIR}")

pcs=read.table("data_PCA.eigenvec",h=F,skip=1)

pdf("PCA.plot.pdf")
plot(pcs$V3,pcs$V4, xlab="PC1", ylab="PC2", pch=19)
dev.off()
