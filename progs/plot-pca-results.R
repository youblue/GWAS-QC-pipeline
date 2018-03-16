setwd("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105")

pcs=read.table("data_PCA.eigenvec",h=F,skip=1)

pdf("PCA.plot.pdf")
plot(pcs$V3,pcs$V4, xlab="PC1", ylab="PC2", pch=19)
dev.off()

pdf("PCA2.plot.pdf")
plot(pcs$V3,pcs$V5, xlab="PC1", ylab="PC3", pch=19)
dev.off()

pdf("PCA3.plot.pdf")
plot(pcs$V4,pcs$V5, xlab="PC2", ylab="PC3", pch=19)
dev.off()