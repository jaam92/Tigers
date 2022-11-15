#!/bin/sh

#load software
ml biology bedtools

for f in {1..18} 
do 

coverageBed -a ROH/"$f"_allIndivs_garlicROH_700bpWindow_sorted.bed -b "$f"_fixed_hom_18SFSvar.sorted.fmt.merged.bed > ROH/"$f"_allIndivs_garlicROH_700bpWindow_coverage.bed 

done

echo "done annotating"

for f in {1..18} 
do 

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}' ROH/"$f"_allIndivs_garlicROH_700bpWindow_coverage.bed  >> ROH/allChroms_allIndivs_garlicROH_700bpWindow_CoverageCallableSites.bed 

done

echo "done concatenating results"



