#Load SNP frequency file and generate histogram

setwd("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105")

b.frq <- read.table("clean_inds_data_missing.lmiss",header=T)
pdf("snpmiss_plot.pdf")
plot(ecdf(b.frq$F_MISS),xlim=c(0,0.10),ylim=c(0,1),pch=20, main="SNP Missingness Distribution", xlab="Missingness Frequency", ylab="Fraction of SNPs",col="blue",axes=T)

dev.off()