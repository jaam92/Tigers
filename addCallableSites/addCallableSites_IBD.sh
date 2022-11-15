#!/bin/sh

#takes like 10 minutes to run 

#load software
ml biology bedtools
ml biology R/4.0.2

#convert truffle index to vcf positions for intersecting with callable sites
Rscript convertTruffle2VCFpos.R

echo "done converting truffle vcf positions"

#call ibd and concatenate
for f in {1..18} 
do 
coverageBed -a IBD/"$f"_truffle_allSubSpecies_calledPerSpecies.segments.convert2VCFpos.bed -b "$f"_fixed_hom_18SFSvar.sorted.fmt.merged.bed > IBD/"$f"_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed 

done

echo "done with bedtools"

for f in {1..18}
do

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11}' IBD/"$f"_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed >> IBD/allChroms_truffle_allSubSpecies_calledPerSpecies.segments.coverage.bed 

done

echo "done contenating"
