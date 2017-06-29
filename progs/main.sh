#!/bin/bash

###################################################################
#                     GWAS QC WORKFLOW                            #
#                Wangshu Zhang 2017-06-25                         #
###################################################################

# Current Work Path
export DIR=${DIR}
cd $DIR

# Set Data folder, Program folder, Output folder
export DATA_DIR=$DIR/data     
export PROGS_DIR=$DIR/progs 
export OUTPUTS_DIR=$DIR/outputs

# Load R, Plink, EIGENSOFT for PCA
module load r/3.3.1
export PATH="$PROGS_DIR/plink-1.07-x86_64:$PATH"
export PATH="$PROGS_DIR/EIG-6.1.4/bin:$PATH"

################## t321-0105 train data only ##################
export TIER=${TIER}                            # Data tier
export OUTPUTS_DIR=$OUTPUTS_DIR/QC\_train\_$TIER # Update Output folder
mkdir -p $OUTPUTS_DIR
export DATA=$DATA_DIR/train/$TIER/train\_$TIER   # Update Data folder

# Read Original data
plink --bfile ${DATA} --noweb

#########################################################
#			QC SYNTAX			#
#			Sample QC			#
#########################################################

# 1.Identify individuals with discordant sex information
plink --bfile ${DATA} --check-sex --out ${OUTPUTS_DIR}/sexstat --noweb

# 2.Select individuals with Status=PROBLEM in the file sexstat.sexcheck
grep "PROBLEM" ${OUTPUTS_DIR}/sexstat.sexcheck > ${OUTPUTS_DIR}/fail_sex_check.txt

# 3. Calculate missingness score each individual
plink --bfile ${DATA} --missing --out ${OUTPUTS_DIR}/data_miss --noweb

# 4.Calculate heterozygosity score for each individual  
plink --bfile ${DATA} --het --out ${OUTPUTS_DIR}/data_het --noweb

# 5. Plot the distribution of missingness and heterozygosity scores using the script het_miss.R
R CMD BATCH ${PROGS_DIR}/miss_het_plot.R

# 6. Based on preselected cutoff identify individuals with high missingness and/or outlier heterozygosity
perl ${PROGS_DIR}/select_miss_het.pl

##. DELETE: Join QC failed individuals
cat ${OUTPUTS_DIR}/fail_sex_check.txt $OUTPUTS_DIR/fail_miss_het.txt | sort -k1 | uniq > ${OUTPUTS_DIR}/fail_data_inds.txt 

# 7. Prune SNPs for IBS matrix so that no pair of SNPs within a given window (50kb) has an r2 greater that a given threshold (0.2)
plink --bfile ${DATA} --indep-pairwise 50 5 0.2 --out ${OUTPUTS_DIR}/data_exclude_highLD --noweb # Actual use: Not exclude anything

# 8. Generate pair-wise IBS for all pairs of individuals in the study based on the reduced marker set
plink --bfile ${DATA} --extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in --genome --out ${OUTPUTS_DIR}/pairwiseIBS --noweb

# 9. Identify all pairs of individuals with an IBD > 0.185
perl ${PROGS_DIR}/run-IBD-QC.pl

# 10. Plot pairwise IBD
Rscript ${PROGS_DIR}/plot-IBD.R ${OUTPUTS_DIR}/pairwiseIBS ${OUTPUTS_DIR} ${DATA}.fam

# 11. Extract HapMap SNPs
plink --bfile ${DATA} --extract ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt --make-bed --out ${OUTPUTS_DIR}/extract-hapmap-snps --noweb

# 12. Merge HapMap - phase 1
#plink --bfile ${OUTPUTS_DIR}/extract-hapmap-snps \
#--bmerge ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam \
#--extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in --make-bed --out ${OUTPUTS_DIR}/merge-hapmap --noweb

# 13. Merge HapMap - flip SNPs
#plink --bfile ${OUTPUTS_DIR}/extract-hapmap-snps \
#--extract ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt \
#--flip ${OUTPUTS_DIR}/merge-hapmap.missnp \
#--make-bed \
#--out ${OUTPUTS_DIR}/extract-hapmap-flip-snps \
#--noweb

# 14. Merge HapMap - phase 2
#plink --bfile ${OUTPUTS_DIR}/extract-hapmap-flip-snps \
#--bmerge ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam \
#--extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in \
#--make-bed \
#--out ${OUTPUTS_DIR}/merge-hapmap \
#--noweb

# 15. Corrupted SNPs, unable to flip, exluce corrupt SNPs (differing alleles between HapMap data and my data)
#plink --bfile ${OUTPUTS_DIR}/extract-hapmap-flip-snps \
#--exclude ${OUTPUTS_DIR}/merge-hapmap.missnp \
#--extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in \
#--make-bed \
#--out ${OUTPUTS_DIR}/extract-hapmap-flip-snps-exclude-corrupt \
#--noweb
        	                
# 16. Merge HapMap - phase 3
#plink --bfile ${OUTPUTS_DIR}/extract-hapmap-flip-snps-exclude-corrupt \
#--bmerge ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim \
#	 ${DATA_DIR}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam \
#--extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in \
#--make-bed \
#--out ${OUTPUTS_DIR}/05-merge-hapmap \
#--noweb

# 17. Use PLINK 1.9 to do PCA, top 5
${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/pairwiseIBS --pca 5 --out ${OUTPUTS_DIR}/data_PCA --noweb

# 18. Plot the first two components
Rscript ${PROGS_DIR}/plot-pca-results.R

# 19. Detect outliers from top PCs (out of 3SD)
Rscript ${PROGS_DIR}/detect-pca-outlier.R

# 20. Remove QC failed individuals from data
cat ${OUTPUTS_DIR}/fail* | awk '{ print $1 }' | sort -k1 | uniq > ${OUTPUTS_DIR}/fail_qc_inds.txt # Join QC failed individuals
plink --bfile ${DATA} --remove ${OUTPUTS_DIR}/fail_qc_inds.txt --make-bed --out ${OUTPUTS_DIR}/clean_inds_data --noweb


#########################################################
#			QC SYNTAX			#
#			SNP QC				#
#########################################################

export DATA_SAMPLEQC=$OUTPUTS_DIR/clean_inds_data # data after Sample QC

# 9. Calculate minor allele frequencies
plink --bfile ${DATA_SAMPLEQC} --freq --out ${OUTPUTS_DIR}/clean_inds_data_freq --noweb

# 10. Plot the distribution of MAF values using the script maf_plot.R to decide a MAF cutoff
R CMD BATCH ${PROGS_DIR}/maf_plot.R

# 11. Calculate SNP missingness
plink --bfile ${DATA_SAMPLEQC} --missing --out ${OUTPUTS_DIR}/clean_inds_data_missing --noweb

# 12. Plot the distribution of missingness values using the script snpmiss_plot.R to decide a missingness cutoff
Rscript ${PROGS_DIR}/snpmiss_plot.R

# 13. Calculate differential missingness
plink --bfile ${DATA_SAMPLEQC} --test-missing --out ${OUTPUTS_DIR}/clean_inds_data_test_missing --noweb

# 14. Plot the distribution of differential missingness P-values using the script snpmiss_plot.R to decide a differential missingness P-value cutoff
Rscript ${PROGS_DIR}/diffmiss_plot.R

# 15. Select SNPs showing extreme differential missingness 
perl ${PROGS_DIR}/select_diffmiss.pl 

# 16. Identify SNPs with extreme HWE deviations
plink --bfile ${DATA_SAMPLEQC} --hardy --out ${OUTPUTS_DIR}/clean_inds_data_hwe --noweb

#17. Select Unaffected only for HWE Plot
head -1 ${OUTPUTS_DIR}/clean_inds_data_hwe.hwe > ${OUTPUTS_DIR}/clean_inds_data_hweu.hwe | grep "UNAFF" ${OUTPUTS_DIR}/clean_inds_data_hwe.hwe >> ${OUTPUTS_DIR}/clean_inds_data_hweu.hwe

# 18. Plot the distribution of HWE P-values(in controls) using the script hwe_plot.R to decide a HWE P-value cutoff
Rscript ${PROGS_DIR}/hwe_plot.R

# 19 Remove SNPs failing QC
plink --bfile ${DATA_SAMPLEQC} --maf 0.05 --geno 0.05 --exclude ${OUTPUTS_DIR}/fail_diffmiss_data.txt --hwe 0.00001 --make-bed --out ${OUTPUTS_DIR}/clean_data --noweb

# 20. Remove X chr SNPs
#plink --bfile ${OUTPUTS_DIR}/clean_data --chr 23 --make-bed --out ${OUTPUTS_DIR}/xsnps --noweb
#plink --noweb --bfile ${OUTPUTS_DIR}/clean_data --exclude ${OUTPUTS_DIR}/xsnps.bim --make-bed --out ${OUTPUTS_DIR}/qced_data

#################################################################
#		  CLEAN UP		    		       	#
#   move log and intermediate files to another directory        #
#################################################################

cd $OUTPUTS_DIR

mv *.log Logfiles/.
mv *.nof logfiles/.
mv clean_inds*.* logfiles/.
mv data_*.* logfiles/.
mv clean_data*.* logfiles/.
mv fail_*.* logfiles/.
mv sexstat*.* logfiles/.
mv xsnps.* logfiles/.
mv *.Rout logfiles/.










