
###bypass centromeres
##awk '{print $2"\t"0"\t"0}' /scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/plink/chrom_map.txt > fakeCentromeres.txt


for p in {Amur,Bengal,Generic,Indochinese,Malayan,South_China,Sumatran} 
do 

/scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/software/garlic/bin/linux/garlic --tped /scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/plink/"$p".tped --tfam /scratch/users/elliea/jazlyn-ellie/oct2022-captives-usethese/plink/"$p".tfam --centromere fakeCentromeres.txt --error 0.001 --winsize 700 --out "$p"

done 
