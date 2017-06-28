#--INSPECT MISSINGNESS PATTERNS--#

#IMPORT PLINK FILES WITH MISSINGNESS INFORMATION
#requires the files data_miss.imiss and data_het.het to be present in the script folder

setwd("${OUTPUT_DIR}")

imiss <- read.table("data_miss.imiss",header=T)
het <- read.table("data_het.het",header=T)

#CALCULATE CALL RATE, LOG10(F_FMISS) and mean heterozygosity
imiss$CALL_RATE <- 1-imiss$F_MISS
imiss$logF_MISS = log10(imiss[,6])
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
het$meanHet <- ifelse(het$meanHet=="NaN", c(0),c(het$meanHet))
imiss.het <- merge(het,imiss,by=c("FID","IID"))

#GENERATE CALL RATE BY HETEROZYGOSITY PLOT
colors  <- densCols(imiss$logF_MISS,het$meanHet)
pdf("pairs.imiss-vs-het.pdf")
plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-3,0),ylim=c(0.17,0.21),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
axis(2,at=seq(0.17,0.21, 0.01),tick=T)
axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
#Heterozygosity thresholds (Horizontal Line)
abline(h=mean(het$meanHet)-(3*sd(het$meanHet)),col="RED",lty=2)
abline(h=mean(het$meanHet)+(3*sd(het$meanHet)),col="RED",lty=2)
#Missing Data Thresholds (Vertical Line)
#abline(v=-1.30103, col="BLUE", lty=2) #THRESHOLD=0.05
#abline(v=-1.522879, col="RED", lty=2) #THRESHOLD=0.03

abline(v=-2, col="RED", lty=2) #THRESHOLD=0.01


dev.off()

output_3sd <- c(mean(het$meanHet)-3*sd(het$meanHet), mean(het$meanHet)+3*sd(het$meanHet))
write.csv(output_3sd, "het_output_3sd.csv", row.names = FALSE)
