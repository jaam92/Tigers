We want to know how variation is shared across populations, for this, use the balanced datasets (N=3,6,10)

First create subsets of each population, include --min-ac or else it will just reprint the variable sites from all populations and the overlap will screw up
```
bcftools view high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz -S amur_N3.txt --min-ac=1 -O z -o amurN3-highcorr-nodups-biallelic-AN-MM-pcc.vcf.gz
```

Next, run bcftools isec across the datasets to compare them, like so:
```
bcftools isec -n +2 amurN10-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz bengalN10-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz \
malayanN10-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz sumatranN10-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz \
genericN10-high-corr-nodups-biallelic-AN-MM-pcc.vcf.gz -p N10_compare
```

Then plot, yay.
