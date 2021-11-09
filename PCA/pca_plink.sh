#build pca for high coverage uncorrected individuals
#convert vcf to plink format
plink2 --vcf highcov-nofilter-biallelic-AN-MM-pcc.vcf --make-bed --allow-extra-chr 

#run PCs
plink2 --bfile highcov-nofilter-biallelic-AN-MM-pcc --allow-extra-chr --pca 10 --out highcov-nofilter-biallelic-AN-MM-pcc --threads 4


#build pca for all uncorrected individuals

#build pca for all corrected individuals
