##First run vcftools relatedness on high coverage individuals, split into corrected (after pca) species groups
##Input should be files already converted to only autosomal files

##First split full file into high and low coverage files
bcftools view highcov-lowcov-biallelic-AN-MM-pcc.vcf.gz -S highcov-nofilter.txt -O z -o highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-lowcov-biallelic-AN-MM-pcc.vcf.gz -S lowcov-nofilter.txt -O z -o lowcov-nofilter-biallelic-AN-MM-pcc.vcf.gz 


##Split high coverage files into specific populations
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S bengals-corr.txt -O z -o bengal-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S amurs-corr.txt -O z -o amur-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S generic-corr.txt -O z -o generic-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S indochinese-corr.txt -O z -o indochinese-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S southchina-corr.txt -O z -o southchina-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S sumatran-corr.txt -O z -o sumatran-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz
bcftools view highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz -S malayan-corr.txt -O z -o malayan-highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz


##Then run relatedness measures on each vcf file:
for file in *.vcf.gz;do
  vcftools --vcf ${bn}.vcf.gz --relatedness --out ${bn}
  truffle --vcf ${bn}.vcf.gz --cpu 10 --segments --out ${bn}.truffle
done

