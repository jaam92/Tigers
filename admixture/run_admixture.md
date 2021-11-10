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
```
admixture highcov-nofilter-biallelic-AN-MM-pcc.bed 6 --supervised
```
