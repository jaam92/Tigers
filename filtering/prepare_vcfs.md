First, process high coverage vcf file from initial genotyping

Select only SNPs
```
bcftools view -m2 -M2 -v snps genotypes.highcov.vcf.gz | bgzip -c > genotypes.highcov.biallelic.vcf.gz
```

Index with tabix
```
tabix -p vcf genotypes.highcov.biallelic.vcf.gz
```

Look at summary stats
```
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%DP\n' genotypes.highcov.biallelic.vcf.gz | bgzip -c > AFs.tab.gz
```
Rename putative cat chromosomes 
```
zcat genotypes.highcov.biallelic.vcf.gz | sed -e 's/;HRSCAF=//g' > genotypes.highcov.biallelic2.vcf.gz
zcat genotypes.highcov.biallelic2.vcf.gz | sed -e 's/Scaffold_120HRSCAF207/chrA1/g' | sed -e 's/Scaffold_13HRSCAF59/chrA2/g' | sed -e 's/Scaffold_11HRSCAF30/chrA3/g' | sed -e 's/Scaffold_1HRSCAF1/chrB1/g' | sed -e 's/Scaffold_7HRSCAF14/chrB2/g' | sed -e 's/Scaffold_3HRSCAF9/chrB3/g' | sed -e 's/Scaffold_430HRSCAF546/chrB4/g' | sed -e 's/Scaffold_22HRSCAF95/chrC1/g' | sed -e 's/Scaffold_6HRSCAF13/chrC2/g' | sed -e 's/Scaffold_9HRSCAF19/chrD1/g' | sed -e 's/Scaffold_8HRSCAF18/chrD2/g' | sed -e 's/Scaffold_14HRSCAF68/chrD3/g' | sed -e 's/Scaffold_10HRSCAF21/chrD4/g' | sed -e 's/Scaffold_119HRSCAF206/chrE1/g' | sed -e 's/Scaffold_90HRSCAF172/chrE2/g' | sed -e 's/Scaffold_18HRSCAF89/chrE3/g' | sed -e 's/Scaffold_15HRSCAF71/chrF1/g' | sed -e 's/Scaffold_2HRSCAF5/chrF2/g' > genotypes.highcov.biallelic.renameChroms.vcf.gz
```
Do some filtering to get to get only sites where 90% of the individuals have a call
filter to 90% coverage, 177 individuals so to get the number of GTs for 90 percent coverage -> (177*2)*.9 ~ 318 which is also median AN
```
bcftools filter -i 'AN > 317' genotypes.highcov.biallelic.renameChroms.vcf | bgzip -c > genotypes.highcov.biallelic.renameChroms.vcf.gz
tabix -p vcf genotypes.highcov.biallelic.renameChroms.vcf.gz
```
Pull out putative cat chromosomes 
```
bcftools view -r chrA1,chrA2,chrA3,chrB1,chrB2,chrB3,chrB4,chrC1,chrC2,chrD1,chrD2,chrD3,chrD4,chrE1,chrE2,chrE3,chrF1,chrF2 genotypes.highcov.biallelic.renameChroms.vcf.gz | bgzip -c > genotypes.highcov.biallelic.putativeCatChroms.vcf.gz
tabix -p vcf genotypes.highcov.biallelic.putativeCatChroms.vcf.gz
bcftools filter -i 'AN > 317' genotypes.highcov.biallelic.putativeCatChroms.vcf.gz | bgzip -c > genotypes.highcov.biallelic.putativeCatChroms.ANfilt.vcf.gz
```
Merge with low coverage file
```
bcftools merge genotypes.highcov.biallelic.putativeCatChroms.ANfilt.vcf.gz lowcov-all86-impute.biallelic.putativeCatChroms.vcf.gz -O z -o highcov-lowcov-biallelic-AN-pCC.vcf.gz 
```
Filter out sites that do not have missing calls for at least 90% of individuals (this will clean up the imputed individuals)
```
vcftools --gzvcf highcov-lowcov-biallelic-AN-pCC.vcf.gz --max-missing 0.9 --recode-INFO-all --recode --out highcov-lowcov-biallelic-AN-MM-pcc.vcf
```
Rezip
```
bgzip -c highcov-lowcov-biallelic-AN-MM-pcc.vcf.recode.vcf > highcov-lowcov-biallelic-AN-MM-pcc.vcf.gz
```
Index & proceed to next steps
```
tabix -p vcf highcov-lowcov-biallelic-AN-MM-pcc.vcf.gz
```
