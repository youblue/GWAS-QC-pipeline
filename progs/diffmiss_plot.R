#Load SNP differential missingness file and generate distribution

setwd("${OUTPUT_DIR}")

b.frq <- read.table("clean_inds_data_test_missing.missing",header=T)
b.frq$logP = log10(b.frq$P)
pdf("diffmiss_plot.pdf")
plot(ecdf(b.frq$logP), xlim=c(-10,0),ylim=c(0,1),pch=20, main="Distribution of differential missingness P-values", xlab="logP Differential Missingness", ylab="Fraction of SNPs",col="red",axes=T)

dev.off()
