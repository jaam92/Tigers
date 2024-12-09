ml biology vcftools
ml biology plink/1.90b5.3

#Covert vcf to plink
#vcftools --gzvcf 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.vcf.gz --chrom_map chrom_map.txt --plink --out 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.vcf.gz 

#Convert ped and map to binaries
# plink --ped 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.vcf.gz.ped --map 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.vcf.gz.map --allow-extra-chr --make-bed --out 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap

#update family id
#plink --bfile 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap --update-ids recodeFamID.txt --allow-extra-chr --make-bed --out 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.updatedFID

#recode to tped and tfam per family
#for fam in $(awk '{print $1}' 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.updatedFID.fam | sort | uniq); do echo $fam | plink --bfile 20FMISSDP5.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.updatedFID --keep-fam /dev/stdin --recode transpose --out $fam; done

