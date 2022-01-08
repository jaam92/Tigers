First install admixture with conda
```
conda create --name admixture
conda activate admixture
conda install -c bioconda admixture
```

Convert files to plink
```
plink2 --vcf highcov-nofilter-biallelic-AN-MM-pcc.vcf.gz --allow-extra-chr --make-bed --out highcov-nofilter-biallelic-AN-MM-pcc
```

Admixture only accepts numerical chromosomes, so edit bim file, e.g.
```
sed -i 's/chrA1/1/g' highcov-nofilter-biallelic-AN-MM-pcc.bim
```

For supervised method, also need to create a .pop file and name it with the same prefix as the plink files

To do this, you can get a list of the individuals with 

```
bcftools query -l *.vcf.gz > individuals.txt
```

Then, make a file that has one line per indiviudal with the subspecies group. Group that you want to be called (in this case the Generics) should just be blank lines.
```
admixture highcov-nofilter-biallelic-AN-MM-pcc.bed 6 --supervised
```
