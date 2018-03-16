#Load SNP frequency file and generate cumulative freequency distribution

setwd("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105")

b.frq <- read.table("clean_inds_data_freq.frq",header=T)
pdf("maf_plot.pdf")
plot(ecdf(b.frq$MAF), xlim=c(0,1),ylim=c(0,1),pch=20, main="MAF cumulative distribution",xlab="Minor allele frequency (MAF)", ylab="Fraction of SNPs",axes=T)

dev.off()