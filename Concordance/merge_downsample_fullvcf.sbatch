#!/bin/bash
#SBATCH --job-name=merge-related
#SBATCH --time=2-00:00:00
#SBATCH --mem=80000
#SBATCH --partition=dpetrov
#SBATCH -n 12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology bcftools
module load biology vcftools
for file in *-gencove-pcc-filter.vcf.gz;do
        bn=`basename $file -gencove-pcc-filter.vcf.gz`
        bcftools merge generic.highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.lowcov86.filter.nodups.vcf.gz ${bn}-gencove-pcc-filter.vcf.gz -O z -o ${bn}_genmerge_relate.vcf.gz
        tabix -p vcf ${bn}_genmerge_relate.vcf.gz
        vcftools --gzvcf ${bn}_genmerge_relate.vcf.gz --relatedness --out ${bn}
done
