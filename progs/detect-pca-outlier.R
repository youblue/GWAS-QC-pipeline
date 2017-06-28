setwd("${OUTPUT_DIR}")

pcs=read.table("data_PCA.eigenvec",h=F,skip=1)
pc1 = pcs$V3

range_3sd <- c(mean(pc1)-3*sd(pc1), mean(pc1)+3*sd(pc1))
failed_id = pcs[which(pc1 < range_3sd[1] | pc1 > range_3sd[2]), ]

write.table(failed_id, "fail_PCA.txt", row.names=F, col.names=F, quote=F)
