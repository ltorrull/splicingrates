# SplicingRates

**[Citation](https://elifesciences.org/articles/32537): Pai AA, Henriques T, McCue K, Burkholder A, Adelman K, and Burge CB. (2017). The kinetics of pre-mRNA splicing in the Drosophila genome and the influence of gene architecture. *eLife.* 6:e32537**

The splicing rates pipeline is designed to estimate rates of intron splicing from a timecourse of genomewide nascent RNA 4sU-seq data. Rates are measured as splicing half-lives, which represent the time it takes for half of the transcripts from the gene to splice out the intron. The current model is designed to measure rates of consititutively spliced introns, since the model assumes that splicing of the intron goes to completion among the population of steady-state mature mRNA molecules. Modifications to estimate rates of splicing for alternatively spliced introns are in progress.

The pipeline consists of 3 steps, described in detail below:
1. ```splicingrates_introns```: Identify and get coordinates for specific introns for which to calculate half-lives.
2. ```splicingrates_reads```: Get exon-exon and intron-exon junction reads for the designated introns.
3. ```splicingrates_model```: Estimate half-lives for the designated introns from 4sU-seq data.

Other folders in this repository include:

*eLife_manuscript/* Files to recreate the analyses and figures shown in Pai *et al.* eLife 2017
- ```SplicingRates.Rmd``` : R markdown describing analysis of mRNA splicing rates from Drosophila S2 4sU-seq data
- ```SplicingRates_figures.Rmd``` : R markdown to recreate figures in Pai *et al.* eLife 2017

*simulations/* Files for the simulations conducted in Pai *et al.* eLife 2017
- ```simulation.R``` : simulate nascent 4sU-seq reads across timepoints, genomic, and library parameters
- ```miso_modeling.R``` : simulation of PSI decreate method, where inputs are MISO summary PSI values computed with reads from ```simulation.R```
- ```junction_modeling.R``` : simulation of junction dynamics method, where inputs are reads output from ```simulation.R```

### Requirements (versions used for development)
- samtools (v1.3)
- bedtools (v2.26.0)
- python (v3.8.2)
- R (v4.0.2)

#### Python Dependencies
- pysam (v0.16)
- pandas (v0.25.3)

## Table of contents
[Overview of Splicing Rates pipeline](#overview-of-splicing-rates-pipeline)

[Detailed Tutorial](#detailed-tutorial-to-run-the-splicing-rates-pipeline)


## Overview of Splicing Rates pipeline

The splicing rates pipeline was designed to be run on bam files containing mapped reads from a timecourse of 4sU-seq experiments. Here, we describe the usage of the three main steps of the splicing rates pipeline. Below, we provide sections that discuss alternative parameter usage for classification and quantification, as well as a tutorial walking through all the steps necessary to run the splicing rates pipeline.

### Splicing Rates Introns
Get intron coordinates from a gtf file, parsing by type of intron. This step includes (a) identifying introns that are constitutive or alternative (or not restricting by either), and (b) calculating the distance from the 3' splice site to the end of the gene (3' distance).

```
usage: python3 splicingrates_introns.py [-h] --gtf x.gtf 
                                            [--introns {allintrons,alternative,constitutive}] 
                                             --outname

optional arguments:
  -h, --help            show this help message and exit

Input:
  --gtf X.gtf           gtf file from which to get annotations (full path) (default: None)

Parameters:
  --introns {allintrons, alternative, constitutive}     
                        comma-separated list of types of annotations to output
                        into a file (default: constitutive)

Output:
  --outname             basename of output file (full path) (default: None)
```

### Splicing Rates Reads
Identify and quantify exon-exon and intron-exon junction reads for designated introns. This step includes (a) specifying the intron regions from which to isolate reads, accounting for read length, with ```--intronRegions``` and (b) quantifying junction reads for specific introns with ```--junctionReads```.

```
usage: python3 splicingrates_reads.py [-h] --introns x.introntype.bed --readlength readlength
                                          [--intronRegions] [--overlap]
                                          [--junctionReads] [--bam x.bam] [--outdir]
                                          [--readtype {single,paired}]
                                          [--readstrand {fr-unstrand,fr-firststrand,fr-secondstrand}]

Get relevant intron-exon and exon-exon junction reads for splicing rate
estimation.

optional arguments:
  -h, --help            show this help message and exit
  --introns x.introntype.bed
                        bed file with specific introns for which to estimate
                        splicing rates. required if --intronRegions (default:
                        None)
  --readlength readlength
                        length of reads (default: None)

region information:
  --intronRegions       Calculate specific intronic regions to use for
                        grabbing reads (default: False)
  --overlap             overlap of junction reads with region of interest (nt)
                        (default: 10)                        

read quantification:
  --junctionReads       Quantify junction reads from specified introns
                        (default: False)
  --bam x.bam           bam from which to extract junction reads. required if
                        --junctionReads (default: None)
  --outdir              directory for temporary and final output files. required if
                        --junctionReads (default: None)
  --readtype {single,paired}
                        type of read (default: paired)
  --readstrand {fr-unstrand,fr-firststrand,fr-secondstrand}
                        directionality of RNA-seq data 
                        (default: fr-firststrand)
```

### Splicing Rates Model
Estimate splicing rates (half-lives) for designated introns. This step includes estimating half-lives for either/or (a) each replicate individually and (b) junction read counts summed across replicates.

```
usage: python3 splicingrates_model.py [-h] --basename x_Xm_repX --timepoints t1,t2,...tn 
                                          [--replicates] [--summed]
                                           --outname OUTNAME [--txnrate]
                                           --threedist x_threedist.txt

Estimate splicing rates (half-lives) using ie and ee junction reads.

optional arguments:
  -h, --help            show this help message and exit

read information:
  --basename BASENAME   basename of combo files, including full path. 
                        expected to have '_Xm_repX' prefix
                        (default: None)
  --timepoints t1,t2,...tn
                        comma separated list of 4sU labeling timepoints
                        (minutes) (default: None)
  --replicates          number of replicates for each timepoint (default: 1)
  --summed              Estimate splicing rates with summed junction reads
                        (across replicates) (default: False)

model assumptions:
  --outname OUTNAME     name for final splicing rate output file (full path)
                        (default: None)
  --txnrate             transcription elongation rate in nt/min (default =
                        1500 nt/min) (default: 1500)
  --threedist x_threedist.txt
                        file with 3' distances for each intron from
                        splicinrates_introns.py (full path) (default: None)
```


## Detailed Tutorial to run the splicing rates pipeline

In this tutorial, we walk through all the steps to run the splicing rates pipeline, in the minimum number of command lines and each step individually. For each step, we discuss the possible parameters that can be changed, how to do so, and the considerations involved in each of the parameters. Finally, we show example inputs and outputs of each step (with column explanations) so the user knows what to expect and can make custom files as needed.

### Jump to a step:
[Step 0: Genome Alignment](#step-0-genome-alignment)

[Step 1: Intron Delineation](#step-1-intron-delineation)

[Step 2: Get Overlapping Junction Reads](#step-2-get-overlapping-junction-reads)

[Step 2A: Specifying intron regions](#step-2A-specifying-intron-regions)

[Step 2B: Extracting Junction Reads](#step-2B-extracting-junction-reads)

[Step 3: Estimating Splicing Rates](#step-3-estimating-splicing-rates)

### Step 0: Genome Alignment

Align raw reads in fastq format to the genome with your favorite splicing-aware mapper (ie. STAR | hisat2) to obtain a sorted, indexed bam file. When building a STAR index or running hisat2, we recommend using the same gtf annotation that you will use for downstream steps.

For instance, to map with STAR (using ENCODE parameters) and index the bam:

```
STAR --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate

samtools index [bamfile].bam
```

**SAMPLE NAMING CONVENTION** <br>
***Downstream steps assume that the bam file is named with the following convention: ```[samplename]_[T]m_rep[N].bam```, where T is the 4sU labeling timepoint and N is the replicate. T can be a string (i.e. "totalRNA") as long as the "m" is included. Even if there is just 1 replicate, the "rep1" is expected to be included in the file name. THIS IS IMPORTANT FOR SAMPLES TO BE PROPERLY RECOGNIZED AND PARSED IN LATER STEPS! ***

### Step 1: Intron Delineation

***NOTE: ONLY NEEDS TO BE PERFORMED ONCE PER GENOME BUILD***

This step takes in an annotation file (gtf file) and outputs a bed file of intron coordinates for the type of introns specified.

Example usage:
```
python3 splicingrates_introns.py --gtf [annotations].gtf --introns constitutive 
                                 --outname [genomename]
```

**Types of introns**

Users can choose which type of intron to output. Choices: allintrons, alternative, constitutive (default). The user provided basename will have either ```allintrons```, ```alternative_introns```, or ```constitutive_introns``` added before the extensions for each output file.

<p align="center">
<img src="./readme/introns.png" width="50%" height="50%">
</p>

**Output Files**

Two files are output: </br>
1. Bed file with precise intron coordinates and ID with gene/intron coordinates [bed file, ```[genomename].[introntype].bed```] <br>
2. Tab-separated file with distance from 3' splice site to the end of the annotated gene as shown in the figure above, which we call the 3' distance [txt file, ```[genomename].[introntype].threedist```]

Example intron coordinate bed file:
```
1	12227	12613	ENSG00000223972:1:12227-12613:+	intron	+
1	14501	15005	ENSG00000227232:1:14501-15005:-	intron	-
1	15038	15796	ENSG00000227232:1:15038-15796:-	intron	-
1	15947	16607	ENSG00000227232:1:15947-16607:-	intron	-
1	16765	16858	ENSG00000227232:1:16765-16858:-	intron	-
1	17055	17233	ENSG00000227232:1:17055-17233:-	intron	-
```

Example 3' distance file, where column 2 is the 3' distance:
```
ENSG00000223972:1:12227-12613:+	  1796
ENSG00000227232:1:14501-15005:-	  97
ENSG00000227232:1:15038-15796:-	  634
ENSG00000227232:1:15947-16607:-	  1543
ENSG00000227232:1:16765-16858:-	  2361
ENSG00000227232:1:17055-17233:-	  2651
```

### Step 2: Get Overlapping Junction Reads

Identify and quantify exon-exon and intron-exon junction reads for designated introns. This step includes (a) specifying the intron regions from which to isolate reads, accounting for read length, with ```--intronRegions``` and (b) quantifying junction reads for specific introns with ```--junctionReads``` (for each sample that will be used in the model).

Example usage using default parameters to run both steps in tandem:

```
python3 splicingrates_reads.py --introns [genome].[introntype].bed --readlength 
                               --intronRegions --overlap 10
                               --junctionReads --bam [sample]_[T]m_rep[N].bam --outdir [/path/for/output]
                               --readtype paired --readstrand fr-firststrand
```

### Step 2A: Specifying Intron Regions

***NOTE: This step accounts for the specific read length of the samples being analyzed. It only needs to be performed once per read length for each genome build.***

Before proceeding, you must identify intron regions which junction reads must overlap, using the bed file generated by ```splicingrates_introns.py```:

```
python3 splicingrates_reads.py --intronRegions --introns [genome].[introntype].bed --readlength
```

Three **output** files will be generated (as indicated in figure below):

1. Intron-exon boundary regions: ```[genome].[introntype]_[readlength]nt_intron_ie.bed``` <br>
2. Upstream exon region: ```[genome].[introntype]_[readlength]nt_exonup_ee.bed``` <br>
3. Downstream exon region: ```[genome].[introntype]_[readlength]nt_exondown_ee.bed```

<p align="center">
<img src="./readme/junctionregions.png" width="50%" height="50%">
</p>

### Step 2B: Extracting Junction Reads

***NOTE: This step needs to be performed independently for each sample that will be used for modeling in Step 3.***

If you already have intron region files (from running Step 2A), then you can proceed directly to extracting junction reads:

```
python3 splicingrates_reads.py --introns [genome].[introntype].bed --readlength RL
                               --junctionReads --bam [sample]_[T]m_rep[N].bam --outdir [/path/for/output/]
                               --readtype paired --readstrand fr-firststrand 
```

Junction reads are extracted by parsing the CIGAR strings of mapped reads. To correctly assign junction reads the user needs to provide information about read type and strandedness of the reads:
1. Read type can be changed with ```--readtype``` with option {single or paired}, default: paired </br>
2. Strandedness of the reads can be changed with ```--readstrand``` with options {fr-firststrand, fr-secondstrand, fr-unstrand}, default: fr-firststrand 

Strandedness is determined by the type of library preparation protocol. We borrow the library strandedness naming convention from Tophat/Bowtie:

<p align="center">
<img src="./readme/readStrand.png" width="50%" height="50%">
</p>

Junction reads are assigned to introns based on an matches to either the 3' intron-exon boundary (for intron-exon junction reads) or evidence for splicing together of two specific exon regions corresponding to the 5'ss and 3'ss of the designed intron (for exon-exon split junction reads). Split reads that arise from splicing using only one splice site of the pair will not be counted as a junction read for that intron. The minimum overlap length for junction reads is determined by ```--overlap```, with a default of 10nt as specified when running ```splicingrates_reads.py --intronRegions```. 

<p align="center">
<img src="./readme/junctionreads.png" width="90%" height="90%">
</p>


**Outputs** <br>
This step results in 5 types of files, grouped into three categories:

*intron-exon (IE) reads* <br>
1. bed file(s) of read start sites for non-split reads, named using ```[/path/for/output/][sample]_[T]m_rep[N]_startsites_read[1/2].bed.gz```, with two files corresponding to each mate if using paired end reads <br>
2. bed file(s) of ie reads based on read start falling into the intron-exon boundary region identified in Step 1, named using ```[/path/for/output/][sample]_[T]m_rep[N]_iejunc.bed.gz```.<br>

*exon-exon (EE) reads* <br>
1. bam files containing only the split exon-exon junction reads, named using ```[/path/for/output/][sample]_[T]m_rep[N]_junctions.bam```. <br>
2. bed files containing only exon-exon junction reads overlapping the corresponding exon-exon boundary region identified in Step 1, named using ```[/path/for/output/][sample]_[T]m_rep[N]_ee[up/down]junc.bed.gz```

*combined junction reads* <br>
- file with final ie and ee junction read counts for each intron, named using ```[/path/for/output/][sample]_[T]m_rep[N]_junctionCombo.coverage``` and with the following format:

```
intron	ie_count	ee_count
ENSG00000223972:1:12227-12613:+	  1	0
ENSG00000227232:1:14501-15005:-	  2	36
ENSG00000227232:1:15038-15796:-	  4	11
ENSG00000227232:1:15947-16607:-	  8	29
ENSG00000227232:1:16765-16858:-	  60	110
ENSG00000227232:1:17055-17233:-	  74	19
```

### Step 3: Estimating Splicing Rates

To calculate splicing half-lives:
```
python3 splicingrates_model.py --basename [sample] --timepoints t1,t2,tn
                               --replicates N --outname [sample] 
                               --threedist [genome].[introntype].threedist
                               --summed --txnrate 1500
```

This script takes in all samples from the study, across timepoints and replicates, and calculates splicing half-lives in one of two ways:
1. For each replicate individually, using samples from each timepoint within the replicate (assumes timepoints from the same replicate have been conducted using the same cells or in a single experimental workflow).
2. Across all the samples, using the sum of junction reads across each replicate to create a single sample per timepoint. To run the script in this mode, use ```--summed```.

**Inputs**

For a 4sU-seq dataset in Fly S2 cells with 3 timepoints and 3 replicates, using the desired naming convention will result in the following files, each of which should have been processed using ```splicingrates_reads.py``` (and the output saved in /path/for/output/):

replicate 1: *FlyS2_5m_rep1.bam*, *FlyS2_10m_rep1.bam*, *FlyS2_20m_rep1.bam* <br>
replicate 2: *FlyS2_5m_rep2.bam*, *FlyS2_10m_rep2.bam*, *FlyS2_20m_rep2.bam* <br>
replicate 3: *FlyS2_5m_rep3.bam*, *FlyS2_10m_rep3.bam*, *FlyS2_20m_rep3.bam*

*Example Inputs:*
- <pre><code>--basename /path/for/output/<i>FlyS2</i></code></pre>
- <pre><code>--timepoints <i>5,10,20</i></code></pre>
- <pre><code>--replicates <i>3</i></code></pre>
- <pre><code>--outname <i>FlyS2_splicingrates</i></code></pre>
- <pre><code>--threedist <i>drosmel6_constitutive_introns.threedist</i></code></pre>

**Output**

This step results in ```.halflives``` file(s) with the following columns (one per each replicates with filename ```FlyS2_splicingrates_rep1.halflives```, and one for summed reads if specified with filename ```FlyS2_splicingrates_summed.halflives```):

| Column Name | Description |
| ----------- | ----------- |
| intron | intron name, with coordinates and gene name |
| Tm_ie | count of intron-exon junction reads in timepoint T (repeated for each timepoint) |
| Tm_ee | count of exon-exon junction reads in timepoint T (repeated for each timepoint) |
| Tm_ratio | ratio of ie/ee reads in timepoint T (repeated for each timepoint) |
| Tm_Dprime | D' value, short-form simplification for model fitting |
| Tm_Rprime | R' value, short-form simplification for model fitting |
| half.life | splicing half-life (minutes) |
| fit.error | error around the fit from the R optim function |
| result | tag to interpret missing results |

The final column of the ```.halflives``` file contains a tag specifying the interpretation of missing values in the results, with the following options:
1. *success*: Successful fit with half-lives that can be measured using the experimental time points
2. *no_fit*: Introns for which there are either no junction reads (not expressed) or the optimization algorithm failed to find a local minimum.
3. *noEE_verySlow*: Introns without exon-exon junction reads in the last timepoint. This is evidence that the intron is either not spliced or spliced very slowly, since there is no evidence for splicing by the last timepoint in the time course.
4. *noIE_veryFast*: Introns without intron-exon junctions in the first timepoint. This is evidence that the intron is spliced very quickly, since there is no evidence for unspliced transcripts at the first timepoint in the time course.

