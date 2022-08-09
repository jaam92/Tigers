Here, we are combining the phased reference files so that we can try out local ancestry calling.

```
for file in $(find . -type f -name "*.vcf.gz");do
  bn = `basename $file .vcf.gz`
  zcat $file | grep -v '^#' > ${bn}_nohead.vcf
done

```
