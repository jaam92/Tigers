Here, we are combining the phased reference files so that we can try out local ancestry calling.
Files in 1MB chunks and need to do some renaming soooo
```
for file in $(find . -type f -name "*.vcf.gz");do
  bn = `basename $file .vcf.gz`
  zcat $file | grep -v '^#' > ${bn}_nohead.vcf
done
```

Concatenate all vcfs and add a header (header was just grepped earlier)
```
cat *.vcf > highcov_tigers_phased.vcf
cat header.vcf highcov_tigers_phased_header.vcf
```

Rerun the chrom renaming and pull major chroms so it is consistent with file we are using
```
sed -e 's/-HRSCAF_/HRSCAF/g' highcov_tigers_phased_header.vcf >highcov_tigers_phased_headerv1.vcf
java -jar picard.jar CreateSequenceDictionary -R panthera_tigris_11Jul2019_VGYpX.fasta -O panthera_tigris_11Jul2019_VGYpX.dict
java -jar picard.jar SortVcf I=highcov_tigers_phased_headerv1.vcf O=highcov_tigers_phased_headerv1.sorted.vcf R=panthera_tigris_11Jul2019_VGYpX.fasta SEQUENCE_DICTIONARY=panthera_tigris_11Jul2019_VGYpX.dict
bgzip highcov_tigers_phased_headerv1.sorted.vcf 
tabix -p vcf highcov_tigers_phased_headerv1.sorted.vcf.gz
zcat highcov_tigers_phased_eheader.vcf.gz | sed -e 's/-HRSCAF_//g' | bgzip -c > highcov_tigers_phased_eheader.vcf.gz
zcat genotypes.highcov.biallelic2.vcf.gz | sed -e 's/Scaffold_120HRSCAF207/chrA1/g' | sed -e 's/Scaffold_13HRSCAF59/chrA2/g' | sed -e 's/Scaffold_11HRSCAF30/chrA3/g' | sed -e 's/Scaffold_1HRSCAF1/chrB1/g' | sed -e 's/Scaffold_7HRSCAF14/chrB2/g' | sed -e 's/Scaffold_3HRSCAF9/chrB3/g' | sed -e 's/Scaffold_430HRSCAF546/chrB4/g' | sed -e 's/Scaffold_22HRSCAF95/chrC1/g' | sed -e 's/Scaffold_6HRSCAF13/chrC2/g' | sed -e 's/Scaffold_9HRSCAF19/chrD1/g' | sed -e 's/Scaffold_8HRSCAF18/chrD2/g' | sed -e 's/Scaffold_14HRSCAF68/chrD3/g' | sed -e 's/Scaffold_10HRSCAF21/chrD4/g' | sed -e 's/Scaffold_119HRSCAF206/chrE1/g' | sed -e 's/Scaffold_90HRSCAF172/chrE2/g' | sed -e 's/Scaffold_18HRSCAF89/chrE3/g' | sed -e 's/Scaffold_15HRSCAF71/chrF1/g' | sed -e 's/Scaffold_2HRSCAF5/chrF2/g' > highcov_tigers_phased_eheader_renamechroms.vcf.gz

```
