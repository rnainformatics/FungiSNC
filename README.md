# FungiSNC

## Introduction
FungiSNC is presented to discovery, profile and functional annotate diverse kind of sncRNAs (including miRNA, piRNA, tRFs, siRNA, snRNA, snoRNA, rRNA, phasiRNA and natsiRNA) which can not only be used in sRNA-Seq, but it is can also be used in single cell sRNA-seq.
FungiSNC have these main features:
		• sRNA can be detected and profiled for as many as 29 species.
		• Unbiased classify the sncRNA into different categories.
		• Differential expression analysis of sncRNA with paired cases small RNAs transcriptome.
		• Differential expression analysis with group case small RNAs transcriptome.
		• Gene targets and function prediction for sncRNA.
		• sncRNA modification calling using cleavage-based method
		• Very user-friendly web interfaces and convenient data analysis queue system.

The webserver can be accessed: https://bioinformatics.sc.cn/FungiSNC/.
## Installation
**If you want to use this tool on a Linux system, you can install the required libraries and dependencies using the links below.**
```
git clone https://github.com/rnainformatics/FungiSNC.git
```

### Required Languages

- **Perl >= 5.16** (with multithread support)
- **Python 3**
- **R**

### Required Perl CPAN Modules

- Scalar::Util  
- Data::Dumper  
- Parallel::ForkManager  
- Getopt::Long  
- experimental  
- SVG  
- File::Spec  
- List::Util  
- Math::CDF  
- Try::Tiny  
- JSON  

### Required R Packages

- seqinr  
- XML  
- RCurl  
- data.table  
- ggpubr  
- jsonlite  
- GOstats  
- GenomicAlignments  
- GenomicFeatures  

### Required External Tools

- [Bowtie (v1 & v2)](http://bowtie-bio.sourceforge.net/)
- [Samtools v1.9](https://github.com/samtools/samtools)
- [Tapir](http://bioinformatics.psb.ugent.be/webtools/tapir/)
- [BEDTools](https://bedtools.readthedocs.io/)
- [ViennaRNA >= v2](https://www.tbi.univie.ac.at/RNA/)


## Bundled Tools (already included in the distribution)

- megablast  
- ShortStack  
- mireap  
- mirtrace  
- phasiRNAClassifier  
- PHASIS  
- TargetFinder_1.6  

All are available under `sRNAtools/program/` and don't require separate installation.


## Installation Resources

| Tool         | Download URL |
|--------------|--------------|
| Perl         | https://www.perl.org/ |
| Python 3     | https://www.python.org/ |
| Bowtie2      | [Download](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip) |
| Bowtie1      | [Download](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip) |
| Samtools     | https://github.com/samtools/samtools/releases/tag/1.9 |
| Tapir        | http://bioinformatics.psb.ugent.be/webtools/tapir/tapir-1.2.tar.gz |
| BEDTools     | https://github.com/arq5x/bedtools2


## Configuration
	
	* db path in DBCONFIG.conf
	
	* ViennaRNA perl lib path in /program/mireap/bin/mireap.pl

	* ViennaRNA perl lib path in /program/miRDeep_v2

	* ViennaRNA perl lib path in /program/isomiR2Function_rev.pl

	* The program included (such as 'megablast' and 'ShortStack') excutable

## USAGE

 REQUIRED ARGUMENTS
```
 perl run.pl <parameters>
  -mismatch     <num>   The mismatch allowed in the tags mapping [1]
  -infile       <str>   The path of input file (fasta) [required *]
  -species      <str>   The reference species [ath]
  -minlen       <num>   The minimum value of length interval [16]
  -maxlen       <num>   The maximum value of length interval [32]
  -mindepth     <num>   The least number of alignments (Read Per Million) to call tRFs [20]
  -pvalue       <num>   P-value inferred based on Binomial method to distinguish random fragments [0.01]
  -percent      <num>   Minimum percentage of tRF tags in according tRF-producing region on tRNA [0.03]
  -ncpu <num>   Number of threads used in the analysis [1]
  -outdir       <str>   The output directory of results [tRFtools_time]
```
Authors:

Qi Liu,Xiaoqiang Lang,Guoxian	Liu
	
 Email:
 
  rnainfor@gmail.com
  
  langxiaoqiang@foxmail.com
  
  guoxianliu0205@163.com
