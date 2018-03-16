setwd("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105")

pcs=read.table("data_PCA.eigenvec",h=F,skip=1)
pc1 = pcs$V3

range_6sd <- c(mean(pc1)-6*sd(pc1), mean(pc1)+6*sd(pc1))
failed_id = pcs[which(pc1 < range_6sd[1] | pc1 > range_6sd[2]), ]

write.table(failed_id, "fail_PCA.txt", row.names=F, col.names=F, quote=F)