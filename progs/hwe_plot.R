#Load HWE P-value file and generate frequency_distribution

setwd("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105")

b.frq <- read.table("clean_inds_data_hwe.hwe",header=T)
pdf("hwe_plot.pdf")
b.frq$logP = log10(b.frq$P)
plot(ecdf(b.frq$logP), xlim=c(-10,0),ylim=c(0,0.80),pch=20, main="HWE P-value",xlab="logP (HWE)", ylab="Fraction of SNPs",axes=T)

dev.off()