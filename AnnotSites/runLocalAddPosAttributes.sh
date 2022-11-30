#!/bin/sh
#SBATCH --job-name=annotate
#SBATCH --output=/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/annot.out
#SBATCH --error=/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/annot.err
#SBATCH --time=10:00:00
#SBATCH -p normal,hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=jaam@stanford.edu

#Load python
module load python/2.7.13 

#set variables
for CHROM in {A1,A2,A3,B1,B2,B3,B4,C1,C2,D1,D2,D3,D4,E1,E2,E3,F1,F2}
do

for i in {18032FL-56-01-02,18032FL-56-02-14,18032FL-56-02-16,18032FL-61-01V1-02,18032FL-61-01V1-16,5594-DP-0001_S3,5594-DP-0002_S4,5594-DP-0003_S5,AMU1,AMU11,AMU15,AMU17,AMU18,AMU2,AMU20,AMU21,AMU22,AMU24,AMU4,AMU5,AMU7,AMU8,AMU9,BEN_CI2,BEN_CI3,BEN_CI4,BEN_CI5,BEN_CI6,BEN_NE1,BEN_NE2,BEN_NE3,BEN_NOR2,BEN_SA3,BEN_SA4,BEN_SI1,BEN_SI2,BEN_SI3,BEN_SI4,BEN_SI5,BEN_SI6,GEN1,GEN2,GEN22,GEN5,GEN8,MAL1,MAL10,MAL11,MAL12,MAL13,MAL14,MAL15,MAL2,MAL26,MAL3,MAL4,MAL7,MAL9,SRR5591010,SRR7152379,SRR7152380,SRR7152381,SRR7152382,SRR7152383,SRR7152384,SRR7152385,SRR7152386,SRR7152387,SRR7152388,SRR7152390,SRR7152391,SRR7152393,SRR7152394,SRR7152395,SRR7152396,SRR7152397,SRR7152398,SRR7152399,SRR7152400,SRR7152401,SRR7152402,SRR7152403,SRR7152404,SRR7152405,SRR7152406,SRR7152407,SRR7152408,SRR7152409,SRR7152410,SRR836306,SUM1,SUM13,SUM14,SUM2,SUM3,SUM6,SUM8,SUM9,T14,T17,T18}
do

#give files
inputAnnotFile='/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/VEPAnnots/'$CHROM'_annotatedGTwithVEP_'$i'.txt'
outputFile1='/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/'$CHROM'_annotatedGTwithVEP_'$i'_addROHAnnot.txt'

ROHrange1='/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/ROHs/chr'$CHROM'_typeA_'$i'.txt'
ROHrange2='/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/ROHs/chr'$CHROM'_typeB_'$i'.txt'
ROHrange3='/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/AnnotatedVCF/ROHs/chr'$CHROM'_typeC_'$i'.txt'

#Step 1: Intersect ROH Ranges
python addPositionAttributes.py ${inputAnnotFile} 2 ${outputFile1} ranges ${ROHrange1} TypeA ${ROHrange2} TypeB ${ROHrange3} TypeC 

done


done
