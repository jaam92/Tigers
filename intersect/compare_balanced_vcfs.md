First make subsets of various (N=3,6,10) file sets

```
bcftools view -S N6.txt ../admixture_N3/high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz -O z -o N6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz
```

Next split these according to population
```
bcftools view -S N6_amur.txt ../admixture_N3/high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz --min-ac=1 -O z -o amurN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz
```


Run intersection on all the vcffiles to find what is shared and singleton b/c we dont give a shit about anything else
```
bcftools isec -n -2 amurN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz bengalN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz /
indochineseN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz malayanN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz /
sumatranN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz genericN6-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz /
-p test_N6_compare
```
