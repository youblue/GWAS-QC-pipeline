#!/bin/bash

###################################################################
#                GNS STUDY using GWAS QC WORKFLOW                 #
#                Wangshu Zhang 2017-06-25                         #
###################################################################

# Current Work Path
export DIR=/dc2/wzhang01/Adarsh
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
export TIER=t321-0105                            # Data tier
export OUTPUTS_DIR=$OUTPUTS_DIR/QC\_train\_$TIER # Update Output folder
mkdir -p $OUTPUTS_DIR
export DATA=$DATA_DIR/train/$TIER/train\_$TIER   # Update Data folder

# Read Original data
#plink --bfile ${DATA} --noweb

#########################################################
#			QC SYNTAX			#
#			Sample QC			#
#########################################################

# 1.Identify individuals with discordant sex information
#plink --bfile ${DATA} --check-sex --out ${OUTPUTS_DIR}/sexstat --noweb

# 2.Select individuals with Status=PROBLEM in the file sexstat.sexcheck
#grep "PROBLEM" ${OUTPUTS_DIR}/sexstat.sexcheck > ${OUTPUTS_DIR}/fail_sex_check.txt

# 3. Calculate missingness score each individual
#plink --bfile ${DATA} --missing --out ${OUTPUTS_DIR}/data_miss --noweb
#plink --bfile ${DATA} --chr X --missing --out ${OUTPUTS_DIR}/data_miss_chrX --noweb # only for Chr X

# 4.Calculate heterozygosity score for each individual  
#plink --bfile ${DATA} --het --out ${OUTPUTS_DIR}/data_het --noweb
#plink --bfile ${DATA} --chr 23 --het --out ${OUTPUTS_DIR}/data_het_chrX --noweb


# 5. Plot the distribution of missingness and heterozygosity scores using the script het_miss.R
#R CMD BATCH ${PROGS_DIR}/miss_het_plot.R
#R CMD BATCH ${PROGS_DIR}/miss_het_plot_chrX.R

# 6. Based on preselected cutoff identify individuals with high missingness and/or outlier heterozygosity
#perl ${PROGS_DIR}/select_miss_het.pl

# 7. Prune SNPs for IBS matrix so that no pair of SNPs within a given window (50kb) has an r2 greater that a given threshold (0.1)
#plink --bfile ${DATA} --maf 0.05 --indep-pairwise 50 5 0.1 --out ${OUTPUTS_DIR}/data_exclude_highLD --noweb

# 8. Generate pair-wise IBS for all pairs of individuals in the study based on the reduced marker set
#plink --bfile ${DATA} --extract ${OUTPUTS_DIR}/data_exclude_highLD.prune.in --genome --out ${OUTPUTS_DIR}/pairwiseIBS --noweb

# 9. Identify all pairs of individuals with an IBD > 0.185
#perl ${PROGS_DIR}/run-IBD-QC.pl

# 10. Plot pairwise IBD
#Rscript ${PROGS_DIR}/plot-IBD.R ${OUTPUTS_DIR}/pairwiseIBS ${OUTPUTS_DIR} ${DATA}.fam

# 11. Use PLINK 1.9 to do PCA, top 5
#${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/pairwiseIBS --pca 5 --out ${OUTPUTS_DIR}/data_PCA --noweb

# 12. Plot the first two components
#Rscript ${PROGS_DIR}/plot-pca-results.R

# 13. Detect outliers from top PCs (out of 3SD)
#Rscript ${PROGS_DIR}/detect-pca-outlier.R

# 14. Remove QC failed individuals from data
#cat ${OUTPUTS_DIR}/fail* | sort -k1 | uniq > ${OUTPUTS_DIR}/fail_qc_inds.txt # Join QC failed individuals
#plink --bfile ${DATA} --remove ${OUTPUTS_DIR}/fail_qc_inds.txt --make-bed --out ${OUTPUTS_DIR}/clean_inds_data --noweb


#########################################################
#			QC SYNTAX			#
#			SNP QC				#
#########################################################

#export DATA_SAMPLEQC=$OUTPUTS_DIR/clean_inds_data # data after Sample QC
export DATA_SAMPLEQC=$OUTPUTS_DIR/clean_inds_data_uniqueSNPs

# Print duplicated SNPs for same pos
#${PROGS_DIR}/plink --bfile ${DATA_SAMPLEQC} --list-duplicate-vars ids-only --out ${OUTPUTS_DIR}/duplicatedSNPs --noweb

# 9. Calculate minor allele frequencies
#plink --bfile ${DATA_SAMPLEQC} --freq --out ${OUTPUTS_DIR}/clean_inds_data_freq --noweb

# 10. Plot the distribution of MAF values using the script maf_plot.R to decide a MAF cutoff
#R CMD BATCH ${PROGS_DIR}/maf_plot.R

# 11. Calculate SNP missingness
#plink --bfile ${DATA_SAMPLEQC} --missing --out ${OUTPUTS_DIR}/clean_inds_data_missing --noweb


# Sort SNPs by missingnesss
#sort -k 5 ${OUTPUTS_DIR}/clean_inds_data_missing.lmiss | awk '{print $2}' > ${OUTPUTS_DIR}/sorted_dupSNP

# extract data again using sorted SNPs
#plink --bfile ${DATA_SAMPLEQC} --extract ${OUTPUTS_DIR}/sorted_dupSNP --make-bed --out ${OUTPUTS_DIR}/clean_inds_data_sorted --noweb

#export DATA_SAMPLEQC=$OUTPUTS_DIR/clean_inds_data_sorted



# Exclude duplicated SNPs
#plink --bfile ${DATA_SAMPLEQC} --exclude ${OUTPUTS_DIR}/duplicatedSNPs.txt --make-bed --out ${OUTPUTS_DIR}/clean_inds_data_uniqueSNPs --noweb
#export DATA_SAMPLEQC=$OUTPUTS_DIR/clean_inds_data_uniqueSNPs


# 12. Plot the distribution of missingness values using the script snpmiss_plot.R to decide a missingness cutoff
#Rscript ${PROGS_DIR}/snpmiss_plot.R

# 13. Calculate differential missingness
#plink --bfile ${DATA_SAMPLEQC} --test-missing --out ${OUTPUTS_DIR}/clean_inds_data_test_missing --noweb

# 14. Plot the distribution of differential missingness P-values using the script snpmiss_plot.R to decide a differential missingness P-value cutoff
#Rscript ${PROGS_DIR}/diffmiss_plot.R

# 15. Select SNPs showing extreme differential missingness 
#perl ${PROGS_DIR}/select_diffmiss.pl 

# 16. Identify SNPs with extreme HWE deviations
#plink --bfile ${DATA_SAMPLEQC} --hardy --out ${OUTPUTS_DIR}/clean_inds_data_hwe --noweb

#17. Select Unaffected only for HWE Plot
#head -1 ${OUTPUTS_DIR}/clean_inds_data_hwe.hwe > ${OUTPUTS_DIR}/clean_inds_data_hweu.hwe | grep "UNAFF" ${OUTPUTS_DIR}/clean_inds_data_hwe.hwe >> ${OUTPUTS_DIR}/clean_inds_data_hweu.hwe

# 18. Plot the distribution of HWE P-values(in controls) using the script hwe_plot.R to decide a HWE P-value cutoff
#Rscript ${PROGS_DIR}/hwe_plot.R

# 19 Remove SNPs failing QC

#plink --bfile ${DATA_SAMPLEQC} --indep-pairwise 500K 5 0.8 --out ${OUTPUTS_DIR}/dataSampleQC_pruned --noweb
#plink --bfile ${DATA_SAMPLEQC} --extract ${OUTPUTS_DIR}/dataSampleQC_pruned.prune.in --maf 0.05 --geno 0.05 --exclude ${OUTPUTS_DIR}/fail_diffmiss_data.txt --hwe 0.000001 --make-bed --out ${OUTPUTS_DIR}/clean_data1 --noweb

#plink --bfile ${DATA_SAMPLEQC} --indep-pairwise 500K 5 0.2 --out ${OUTPUTS_DIR}/dataSampleQC_pruned0.2 --noweb
#plink --bfile ${DATA_SAMPLEQC} --extract ${OUTPUTS_DIR}/dataSampleQC_pruned0.2.prune.in --maf 0.05 --geno 0.05 --exclude ${OUTPUTS_DIR}/fail_diffmiss_data.txt --hwe 0.000001 --make-bed --out ${OUTPUTS_DIR}/clean_data1_r0.2 --noweb

#plink --bfile ${DATA_SAMPLEQC} --indep-pairwise 500K 5 0.5 --out ${OUTPUTS_DIR}/dataSampleQC_pruned0.5 --noweb
#plink --bfile ${DATA_SAMPLEQC} --extract ${OUTPUTS_DIR}/dataSampleQC_pruned0.5.prune.in --maf 0.05 --geno 0.05 --exclude ${OUTPUTS_DIR}/fail_diffmiss_data.txt --hwe 0.000001 --make-bed --out ${OUTPUTS_DIR}/clean_data1_r0.5 --noweb





#### Not used:
#plink --bfile ${DATA_SAMPLEQC} --maf 0.05 --geno 0.05 --exclude ${OUTPUTS_DIR}/fail_diffmiss_data.txt --hwe 0.000001 --make-bed --out ${OUTPUTS_DIR}/filtered_data --noweb
#plink --bfile ${OUTPUTS_DIR}/filtered_data --indep-pairwise 500K 5 0.8 --out ${OUTPUTS_DIR}/dataSampleQC_pruned_final --noweb
#plink --bfile ${OUTPUTS_DIR}/filtered_data --extract ${OUTPUTS_DIR}/dataSampleQC_pruned_final.prune.in --make-bed --out ${OUTPUTS_DIR}/clean_data --noweb


# 20. Remove X chr SNPs
#plink --bfile ${OUTPUTS_DIR}/clean_data --chr 23 --make-bed --out ${OUTPUTS_DIR}/xsnps --noweb
#plink --noweb --bfile ${OUTPUTS_DIR}/clean_data --exclude ${OUTPUTS_DIR}/xsnps.bim --make-bed --out ${OUTPUTS_DIR}/qced_data

#################################################################
#		  CLEAN UP		    		       	#
#   move log and intermediate files to another directory        #
#################################################################

#cd $OUTPUTS_DIR

#mv *.log Logfiles/.
#mv *.nof logfiles/.
#mv clean_inds*.* logfiles/.
#mv data_*.* logfiles/.
#mv clean_data*.* logfiles/.
#mv fail_*.* logfiles/.
#mv sexstat*.* logfiles/.
#mv xsnps.* logfiles/.
#mv *.Rout logfiles/.




############### for PCA
# 8. Generate pair-wise IBS for all pairs of individuals in the study based on the reduced marker set
#plink --bfile ${OUTPUTS_DIR}/clean_data1 --genome --make-bed --out ${OUTPUTS_DIR}/clean_data_pairwiseIBS --noweb
#plink --bfile ${OUTPUTS_DIR}/clean_data1_r0.2 --genome --make-bed --out ${OUTPUTS_DIR}/clean_data_pairwiseIBS_r0.2 --noweb

# 11. Use PLINK 1.9 to do PCA, top 5
#${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/clean_data_pairwiseIBS --pca 5 --out ${OUTPUTS_DIR}/clean_data_PCA --noweb
#${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/clean_data_pairwiseIBS_r0.2 --pca 5 --out ${OUTPUTS_DIR}/clean_data_PCA_r0.2 --noweb



######### merge and pca #########
#plink --bfile ${OUTPUTS_DIR}/clean_data1 --bmerge /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1.bed /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1.bim /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1.fam --make-bed --out ${OUTPUTS_DIR}/merge
#plink --bfile ${OUTPUTS_DIR}/clean_data1_r0.2 --bmerge /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.2.bed /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.2.bim /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.2.fam --make-bed --out ${OUTPUTS_DIR}/merge_r0.2
#plink --bfile ${OUTPUTS_DIR}/clean_data1_r0.5 --bmerge /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.5.bed /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.5.bim /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data1_r0.5.fam --make-bed --out ${OUTPUTS_DIR}/merge_r0.5






#plink --bfile ${OUTPUTS_DIR}/merge --genome --make-bed --out ${OUTPUTS_DIR}/merge_data_pairwiseIBS --noweb
#plink --bfile ${OUTPUTS_DIR}/merge_r0.2 --genome --make-bed --out ${OUTPUTS_DIR}/merge_data_pairwiseIBS_r0.2 --noweb
#${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/merge_data_pairwiseIBS --pca 5 --out ${OUTPUTS_DIR}/merge_data_PCA --noweb
#${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/merge_data_pairwiseIBS_r0.2 --pca 5 --out ${OUTPUTS_DIR}/merge_data_PCA_r0.2 --noweb





#################### only extract known snps ####################
#plink --bfile ${DATA} --extract ${DATA_DIR}/snp_known_list.txt --make-bed --out ${OUTPUTS_DIR}/clean_data_knownSnps --noweb
#plink --bfile ${OUTPUTS_DIR}/clean_data_knownSnps --bmerge /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data_knownSnps.bed /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data_knownSnps.bim /dc2/wzhang01/Adarsh/outputs/QC_train_t321-0106/clean_data_knownSnps.fam --make-bed --out ${OUTPUTS_DIR}/merge_knownSNPs
plink --bfile ${OUTPUTS_DIR}/merge_knownSNPs --genome --make-bed --out ${OUTPUTS_DIR}/merge_data_knownSNPs_pairwiseIBS --noweb
${PROGS_DIR}/plink --bfile ${OUTPUTS_DIR}/merge_data_knownSNPs_pairwiseIBS --pca 5 --out ${OUTPUTS_DIR}/merge_data_knownSNPs_PCA --noweb




























