#Load SNP frequency file and generate histogram

setwd("${OUTPUT_DIR}")

b.frq <- read.table("clean_inds_data_missing.lmiss",header=T)
pdf("snpmiss_plot.pdf")
plot(ecdf(b.frq$F_MISS),xlim=c(0,0.10),ylim=c(0,1),pch=20, main="SNP Missingness Distribution", xlab="Missingness Frequency", ylab="Fraction of SNPs",col="blue",axes=T)

dev.off()
