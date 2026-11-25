# Genomic diversity and admixture history of captive tigers in the United States

Code and analysis notes accompanying:

> Armstrong, E. E., Mooney, J. A., et al. (2024). Unraveling the genomic diversity and admixture history of captive tigers in the United States. *Proceedings of the National Academy of Sciences*, 121(37), e2402924121. [Full text on PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11441546/)

This repository contains the scripts used to filter genotype data and to estimate population structure, ancestry, relatedness, inbreeding, and the distribution of deleterious variation in captive U.S. tigers compared with reference individuals from recognized subspecies (Amur, Bengal, Indochinese, Malayan, South China, and Sumatran). Captive tigers of unknown or mixed ancestry are referred to as "Generic" throughout the code.

## Data availability

| Resource | Location |
|---|---|
| Code (archived release of this repository) | Zenodo: https://doi.org/10.5281/zenodo.13540809 |
| Whole-genome sequence data | NCBI BioProject [PRJNA976043](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA976043) |
| Raw VCFs | Dryad: https://doi.org/10.5061/dryad.k0p2ngff1 |

Large input files (VCFs, BAMs, and the reference genome) are not stored in this repository; download them from the sources above.

## Repository structure

Directories are listed roughly in the order the analyses were run.

| Directory | Contents |
|---|---|
| `IndivFiles/` | Sample metadata and the lists of unrelated individuals used downstream. `AnnotateAndSummarizeMetadata.R` summarizes the sample set. |
| `filtering/` | `prepare_vcfs.md`: selecting biallelic SNPs, renaming scaffolds to domestic-cat chromosome names (chrA1–chrF2), call-rate/depth/quality filters, removal of low-mappability sites, and merging high- and low-coverage samples. |
| `Relatedness/` | Pairwise relatedness with VCFtools and SNPRelate/GENESIS (`running_SNPRelate.md`, `vcftools_truffle.sh`). |
| `PCA/` | PCA with PLINK 2 (`pca_plink.md`), plotting (`plot_plink_PCAs.md`, `PCA_finalData.R`), identifying duplicates and unrelated individuals (`ID_dupsAndunrelateds.R`), and private/shared variation across subspecies (`private_variation.md`). |
| `IBS_Clustering.R` | IBS-based clustering of individuals. |
| `admixture/` | Unsupervised and supervised ADMIXTURE runs (`run_admixture.md`). |
| `Concordance/` | Downsampling high-coverage genomes and comparing ancestry estimates across coverage levels (full, 5x, 2x, 1x, 0.5x, 0.25x). |
| `local_ancestry/` | Preparing phased reference panels and collapsing local ancestry calls. |
| `intersect/`, `SharedVariation/` | Comparing variation shared among subspecies using sample-size-balanced subsets (N = 3, 6, 10). |
| `FSNP_Het/` | Genome-wide heterozygosity and plotting. |
| `ROH/` | Runs of homozygosity called with GARLIC (`step1`–`step3` shell scripts), comparison with bcftools/roh, and inbreeding coefficients (F<sub>ROH</sub>) in `InbreedingAnalyses.R`. |
| `IBD/` | IBD segments (TRUFFLE) and IBD-based comparisons within and between subspecies. |
| `addCallableSites/` | Correcting ROH and IBD estimates for callable sites. |
| `AnnotSites/` | Annotating polarized variants with VEP, SIFT, and GERP scores, and counting deleterious genotypes inside and outside ROH. |
| `sfs/` | Site frequency spectra from polarized VCFs, including SFS split by SIFT annotation. |

## Summary files

Several directories include small summary tables (`.csv` and `.txt`) with results or inputs that can be used without rerunning the full pipeline.

| File | Description |
|---|---|
| `IndivFiles/individual_ids.csv` | Metadata for every sample: subspecies assignment (two schemes, `Subspecies` and `Subspecies2`), sequencing coverage (high/low), coat phenotype, and sex. |
| `IndivFiles/Unrelateds_basedOnSubspecies1_highCov.csv`, `IndivFiles/Unrelateds_basedOnSubspecies2_highCov.csv` | IDs of the unrelated high-coverage individuals used in downstream analyses, under each subspecies assignment scheme. |
| `IBD/IBDScores_subspecies1Unrelateds.txt`, `IBD/IBDScores_subspecies2Unrelateds.txt` | IBD scores per subspecies: sample size, normalization constant, group IBD score, score normalized per Mb, and score relative to Generic tigers. |
| `ROH/TrueROH_propCoveredwithin1SDMean_gr100kb_allChroms_highCov_runSpeciesSep_garlic.txt` | Filtered ROH calls (>100 kb) from GARLIC for high-coverage individuals: ROH size class, start, end, length, individual, chromosome, proportion of the ROH covered by callable sites, and whether the ROH was retained. |
| `ROH/chrom_map.txt`, `ROH/fakeCentromeres.txt` | Chromosome name mapping and placeholder centromere positions used as GARLIC inputs. |
| `AnnotSites/GTAnnotationCountResults_Nov2022_Tigers_addSIFT_noID.txt` | Counts of deleterious variants per individual. |
| `Concordance/ancestry_concordance_all.csv` | Estimated ancestry fractions for each source subspecies, per individual and per downsampled coverage level. |

## Software

The analyses used the following tools; please see the Methods and SI Appendix of the paper for versions and parameters.

Command line: bcftools, tabix/bgzip, bedtools, VCFtools, PLINK 1.9 and PLINK 2, ADMIXTURE, GARLIC, TRUFFLE, GATK, Picard, Ensembl VEP, and SIFT.

R (v4.0+): tidyverse, data.table, ggpubr, patchwork, SNPRelate, GENESIS, GWASTools, GenomicRanges, cluster, seriation, corrplot, ggdendro, dendextend, and optparse.

## Usage notes

Most steps are documented as annotated commands in Markdown files rather than as a single automated pipeline, so they are best read and run step by step. Many scripts contain hard-coded paths from the authors' local machines and the Stanford Sherlock cluster (for example `setwd()` calls and `/scratch/...` paths); update these before running. Several shell scripts use `ml`/`module load` commands specific to Sherlock.

The reference genome is the *Panthera tigris* assembly `panthera_tigris_11Jul2019_VGYpX`, with major scaffolds renamed to their putative domestic-cat chromosome equivalents (the mapping is given in `filtering/prepare_vcfs.md` and `ROH/chrom_map.txt`).

## Citation

If you use this code or data, please cite the paper:

```bibtex
@article{armstrong2024captivetigers,
  title   = {Unraveling the genomic diversity and admixture history of captive tigers in the United States},
  author  = {Armstrong, Ellie E. and Mooney, Jazlyn A. and others},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {121},
  number  = {37},
  pages   = {e2402924121},
  year    = {2024},
  doi     = {10.1073/pnas.2402924121}
}
```

The archived code can also be cited directly via Zenodo: https://doi.org/10.5281/zenodo.13540809

## License

This work is licensed under the [Creative Commons Attribution-ShareAlike 4.0 International License](LICENSE).

## Contact

For questions, please open an issue on this repository or contact the corresponding authors listed in the paper.
