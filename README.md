![Build Status](https://gitlab.com/pages/plain-html/badges/master/build.svg)

---

## Short Read RNA-seq data analysis tutorial

---

This repository contains scripts and documentation associated with the analyses of Short Read RNA-seq data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Quality Control using FastQC and MultiQC](#Quality Control using FastQC and MultiQC)
- [Generate genome index file using STAR](#Generate genome index file using using STAR)
- [Aligning Reads using STAR](#Aligning Reads using STAR)
- [Differential Expression Genes analysis](#Differential Expression Genes analysis)
- [Scripts](#Scripts)


<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Quality Control using FastQC and MultiQC

A step-by-step description of how to analyze short read RNA-seq data using **FastQC**, **MultiQC** and **STAR** tools.

- **STEP 1: (Quality Control using FastQC)**

  We used FastQC and MultiQC to check the quality of the raw sequence data coming from high throughput sequencing. FastQC provides a set of analyses which gives a quick impression of whether the data has any problems.

  `./01_FastQC.sh`
  
  `./02_MultiQC.sh`


## Generate genome index file using STAR

- **STEP 2:(Generate genome index file using STAR)**

  The first step of mapping sequencing data is to build a genome index. We used the human genome GRCh38 sequence fasta file and an annotation GTF file to generate genome index.

  `./03_Create_STAR_index.sh`

## Aligning Reads using STAR

- **STEP 3:(Aligning Reads using STAR)**

  We used the human genome GRCh38 to map the reads. The reads were mapped to the GRCh38 reference genome using STAR. 

  `./04_STAR.sh`

  By running STAR in the 2-pass mode, novel junctions will be discovered. The 2-pass mode does not increase the number of detected novel junctions, but allows to detect more splices reads mapping to novel junctions. The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the first pass, and use them as ”annotated” junctions for the 2nd pass mapping. STAR will perform the 1st pass mapping, then it will automatically extract junctions, insert them into the genome index, and, finally, re-map all reads in the 2nd mapping pass. This option can be used with annotations, which can be included either at the run-time (see #1), or at the genome generation step.

  `./05_STAR_2nd_pass_option.sh`

- **STEP 4: (Create Gene Counts)**

  Generate gene count file for each samples.

  `./06_STAR_gene_counts`

## Differential Expression Genes analysis using EdgeR

- **STEP 5: (Build Eset Raw file)**

First step of DEGs (Differntial Expression Genes) analysis is to load count of the genes, sample annotations and gene annotations files. Build a eset with raw counts.

`mkdir counts`
`mkdir save`
`mkdir data`

`07_buildEsetRaw.Rmd`

- **STEP 6: (Build Contrasts between samples/conditions)**

In the next step, normalize the data and set up the contrasts between samples/conditions. Build a eset with normalized data that includes contrasts information.

`08_buildContrasts.Rmd`

- **STEP 7: (Analysis Report)**

In the final step, create a analysis report that includes heatmap of DEGs and pathway analysis. The report also includes table of DEGs and pathway information with p.value and adj p.value.

`09_Analysis_Report.Rmd`

## Scripts

- FastQC (Quality Control)
  - 01_FastQC.sh

- MultiQC (Quality Control)
  - 02_MultiQC.sh
      
- STAR (Aligning the Reads)
  - 03_Create_STAR_index.sh
  - 04_STAR.sh
  - 05_STAR_2nd_pass_option.sh
  - 06_STAR_gene_counts.sh

- EdgeR (Differential Expression Analysis)
  - 07_buildEsetRaw.Rmd
  - 08_buildContrasts.Rmd
  - 09_Analysis_Report.Rmd
