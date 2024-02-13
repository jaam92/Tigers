#!/bin/bash
#SBATCH --job-name=GATK-cats
#SBATCH --time=6-00:00:00
#SBATCH --mem=240000
#SBATCH --partition=dpetrov
#SBATCH -n 36
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ellieearmstrong@gmail.com

module load biology bwa
module load biology samtools
module load biology gatk
module load biology bcftools

REFERENCE=/oak/stanford/groups/dpetrov/armstrong/captive_tigers/new-2022/concordance/snp_calling/panthera_tigris_11Jul2019_VGYpX.fasta
THREADS_BWA=46

bwa index -a bwtsw $REFERENCE
samtools faidx $REFERENCE
java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar CreateSequenceDictionary R=$REFERENCE

for file in *_1.fastq.gz;do
       bn=`basename $file _1.fastq.gz`
       bwa mem -t $THREADS_BWA $REFERENCE ${bn}_1.fastq.gz ${bn}_2.fastq.gz | samtools view -bS - | samtools sort - > ${bn}.sorted.bam
       samtools index ${bn}.sorted.bam
done

for file in *.sorted.bam;do
       bn=`basename $file .sorted.bam`
       java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar AddOrReplaceReadGroups I=${bn}.sorted.bam O=${bn}.sorted.rg.bam SORT_ORDER=coordinate RGID=${bn} RGLB=${bn} RGPL=illumina RGSM=${bn} RGPU=${bn} CREATE_INDEX=True
       java -jar /oak/stanford/groups/dpetrov/armstrong/applications/picard.jar MarkDuplicates I=${bn}.sorted.rg.bam O=${bn}.sorted.rg.md.bam M=${bn}.metrics.md.txt 
       samtools index ${bn}.sorted.rg.md.bam
       gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.haplotypecaller.g.vcf -ERC GVCF -R $REFERENCE
done

for file in *.sorted.rg.md.bam;do
       bn=`basename $file .sorted.rg.md.bam`
       gatk HaplotypeCaller -I ${bn}.sorted.rg.md.bam -O ${bn}.rg.md.haplotypecaller.g.vcf -ERC GVCF -R $REFERENCE
done

gatk CombineGVCFs -R $REFERENCE --variant SRR7651465.rg.md.haplotypecaller.g.vcf --variant SRR7651468.rg.md.haplotypecaller.g.vcf --variant SRR836354.rg.md.haplotypecaller.g.vcf -O tiger-impute.g.vcf.gz
gatk --java-options "-Xmx4g" GenotypeGVCFs -R $REFERENCE -V tiger-impute.g.vcf.gz -O tiger-impute.vcf.gz
bcftools view -R ../../relatedness/all_merge/highcov.nodups.highcov.biallelic.AN.pCC.DPQUALtest.genmap.bed -O z -o tiger-impute-ec-pcc-filterhigh.vcf.gz tiger-impute-ec-pcc.vcf.gz
