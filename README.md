---
layout: tutorial_page
permalink: /Bioinformatics-for-genomics-medecine_2017_module2_lab
title: BFX-Genomics-medecine Lab 2
header1: Workshop Pages for Students
header2: Bioinformatics for Genomics Medicine Module 2 Lab
image: /site_images/CBW_population_icon.jpg
home: https://bioinformaticsdotca.github.io/
---

-----------------------

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

-----------------------

# CBW  BFX-Genomics-medicine Module 2 - Genome Alignment   

 
by Mathieu Bourgey, _Ph.D_



## Introduction

The workshop is focused on the application of bioinformatics analysis to provide genomics oriented medecine. In most of the case this type of approach focus on exome sequencing. In more speecific case, for example cancer, the analysis will focus on whole genome. In both the cases the analysis is very similar and the main differences are the use of restricted area of  search for the different steps. In this practical we will show you how to launch individual first steps of a DNA-Seq pipeline which apply for both exome and whole genome data.

We will be working on a 1000 genome sample, NA12878. You can find the whole raw data on the 1000 genome website:
<http://www.1000genomes.org/data>

NA12878 is the child of the trio while NA12891 and NA12892 are her parents.

|Mother|Father|Child|
|---|---|---|
|NA12892|NA12891|NA12878|


If you finish early, feel free to perform the same steps on the other two individuals: NA12891 & NA12892. 


For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
We're going to focus on the reads extracted from a 300 kbp stretch of chromosome 1

|Chromosome|Start|End|
|---|---|---|
|1|17704860|18004860|



## Original Setup

### Compute Canada HPC server

Read these [directions](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/scripts/logging_into_the_Compute_Canada_HPC.md) for information on how to log into your assigned working node. 

### Software requirements

These are all already installed, but here are the original links.

  * [BVATools](http://bitbucket.org/mugqic/bvatools/downloads/)
  * [SAMTools](http://sourceforge.net/projects/samtools/)
  * [BWA](http://bio-bwa.sourceforge.net/)
  * [Genome Analysis Toolkit](http://www.broadinstitute.org/gatk/)
  * [Picard](http://broadinstitute.github.io/picard/)
  * [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### Environment setup

```
#set up
echo "export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6" >> ~/.bashrc
echo "module use $MUGQIC_INSTALL_HOME/modulefiles"  >> ~/.bashrc

source ${HOME}/.bashrc

export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/
export WORK_DIR="${HOME}/bfx_genomic_medecine/module2"

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/samtools/1.4 
module load mugqic/bwa/0.7.12 mugqic/GenomeAnalysisTK/3.7 mugqic/picard/1.123 
module load mugqic/trimmomatic/0.36 mugqic/R_Bioconductor/3.3.2_3.4


rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
ln -s /home/partage/genomic_medecine/Module2/* .
```

### Data files

The initial structure of your folders should look like this:   


```
ROOT
|-- raw_reads/               # fastqs from the center (down sampled)
    `-- NA12878/             # Child sample directory
    `-- NA12891/             # Father sample directory
    `-- NA12892/             # Mother sample directory
`-- scripts/                 # command lines scripts
`-- saved_results/           # precomputed final files
```

### Cheat sheets

* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)
* [commands file of this module](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/scripts/commands.sh)


## First data glance

So you've just received an email saying that your data is ready for download from the sequencing center of your choice.

**What should you do?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_data.md)


### Fastq files

Let's first explore the fastq file.

Try these commands

```
zless -S raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz

```

These are fastq file. 

**Could you descride the fastq format?** [Solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_fastq1.md)

```
zcat raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz | head -n4
zcat raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz | head -n4
```

**What was special about the output and why was it like that?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_fastq2.md)


You could also count the reads

```
zgrep -c "^@SN1114" raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz
```

We found  56512 reads

**Why shouldn't you just do?**

```
zgrep -c "^@" raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz
```

[solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_fastq3.md)


### Quality

We can't look at all the reads. Especially when working with whole genome 30x data. You could easilly have Billions of reads.

Tools like FastQC and BVATools readsqc can be used to plot many metrics from these data sets.

Let's look at the data:

```
mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/
```

TOCHANGE scp or ssh -X open a web browser on your laptop, and navigate to `http://cbwXX.dyndns.info/`, where `XX` is the id of your node. You should be able to find there the directory hierarchy under `~/workspace/` on your node. open ```originalQC``` folder and open the images.


**What stands out in the graphs?**

[solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_fastqQC1.md)


All the generated graphics have their uses. This being said 2 of them are particularly useful to get an overal picture of how good or bad a run went. 


These are the Quality box plots 

<img src="https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/img/QualityBoxPlot.png?raw=true" alt="Quality box plots" width="750" />  


and the nucleotide content graphs.

<img src="https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/img/nucleotide_content.png?raw=true" alt="nucleotide_content" width="750" />  

The Box plot shows the quality distribution of your data. The Graph goes > 100 because both ends are appended one after the other.


The quality of a base is computated using the Phread quality score.

<img src="https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/img/phred_formula.png?raw=true" alt="phred_formula" width="750" />  

The formula outputs an integer that is encoded using an [ASCII](http://en.wikipedia.org/wiki/ASCII) table. 

<img src="https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/img/ascii_table.png?raw=true" alt="ascii_table" width="750" />  

The way the lookup is done is by taking the the phred score adding 33 and using this number as a lookup in the table. The Wikipedia entry for the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format) has a summary of the varying values.

Older illumina runs were using phred+64 instead of phred+33 to encode their fastq files.



### Trimming

After this careful analysis of the raw data we see that

- Some reads have bad 3' ends.
- no read has adapter sequences in it.

Although nowadays this doesn't happen often, it does still happen. In some cases, miRNA, it is expected to have adapters. Since they are not part of the genome of interest they should be removed if enough reads have them.


To be able to remove adapters and low qualtity beses we will use Trimmomatic. 

The adapter file is already in your reference  folder.

We can look at the adapters

```
cat adapters.fa
```

**Why are there 2 different ones ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_trim1.md)


Let's try removing them and see what happens.

```
mkdir -p reads/NA12878/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12878/NA12878.trim.out

cat reads/NA12878/NA12878.trim.out
```

**What does Trimmomatic says it did ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_trim2.md)

Let's look at the graphs now

```
mkdir -p postTrimQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/
```

**How does it look now?** 
[solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_trim3.md)

**Could we have done a better job?** 
[solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_trim4.md)


## Alignment

The raw reads are now cleaned up of artefacts we can align the read to the reference.

In case you have multiple readsets or library you should align them separatly !

**Why should this be done separatly?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_aln1.md)

```
mkdir -p alignment/NA12878/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPU:runNA12878_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  reads/NA12878/NA12878_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12878/NA12878_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_HOME}/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12878/NA12878.sorted.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000
```

**Why is it important to set Read Group information?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_aln2.md)

The details of the fields can be found in the SAM/BAM specifications [Here](http://samtools.sourceforge.net/SAM1.pdf)
For most cases, only the sample name, platform unit and library one are important. 

**Why did we pipe the output of one to the other? Could we have done it differently?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_aln3.md)



### Lane merging (optional)

In case we generate multiple lane of sequencing or mutliple library. It is not practical to keep the data splited and all the reads should be merge into one massive file. 

Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.


### SAM/BAM

Let's spend some time to explore bam files.

try

```
samtools view alignment/NA12878/NA12878.sorted.bam | head -n4
```

Here you have examples of alignment results.
A full description of the flags can be found in the [SAM specification](http://samtools.sourceforge.net/SAM1.pdf)

Try using [picards explain flag site](https://broadinstitute.github.io/picard/explain-flags.html) to understand what is going on with your reads

The flag is the 2nd column.

**What do the flags of the first 4 reads mean?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_sambam1.md)

Let's take the 3nd one and find it's pair.

try

```
samtools view alignment/NA12878/NA12878.sorted.bam | grep "1313:19317:61840"

```

**Why did searching one name find both reads?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_sambam4.md)


You can use samtools to filter reads as well.

```
# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12878/NA12878.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12878/NA12878.sorted.bam

```

**How many reads mapped and unmapped were there?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_sambam2.md)


Another useful bit of information in the SAM is the CIGAR string.
It's the 6th column in the file. This column explains how the alignment was achieved.

 * M == base aligns *but doesn't have to be a match*. A SNP will have an M even if it disagrees with the reference.
 * I == Insertion
 * D == Deletion
 * S == soft-clips. These are handy to find un removed adapters, viral insertions, etc.

An in depth explanation of the CIGAR can be found [here](http://genome.sph.umich.edu/wiki/SAM)
The exact details of the cigar string can be found in the SAM spec as well.
Another good site

## Cleaning up alignments

We started by cleaning up the raw reads. Now we need to fix and clean some alignments.

### Indel realignment

The first step for this is to realign around indels and snp dense regions.
The Genome Analysis toolkit has a tool for this called IndelRealigner.

It basically runs in 2 steps

1- Find the targets
2- Realign them.

```
java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/realign.intervals \
  -I alignment/NA12878/NA12878.sorted.bam \
  -L 1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/NA12878/realign.intervals \
  -o alignment/NA12878/NA12878.realigned.sorted.bam \
  -I alignment/NA12878/NA12878.sorted.bam

```

**How could we make this go faster?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_realign1.md)

**How many regions did it think needed cleaning?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_realign2.md)

### FixMates (optional)

This step shouldn't be necessary...But it is some time.

This goes through the BAM file and find entries which don't have their mate information written properly.

This used to be a problem in the GATKs realigner, but they fixed it. It shouldn't be a problem with aligners like BWA, but there are always corner cases that create one-off corrdinates and such.

This happened a lot with bwa backtrack. This happens less with bwa mem and recent GATK so we will skip today.

```
#java -Xmx8G -jar ${PICARD_HOME}/FixMateInformation.jar \
#VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000 \
#INPUT=alignment/NA12878/NA12878.realigned.sorted.bam \
#OUTPUT=alignment/NA12878/NA12878.matefixed.sorted.bam
```

### Mark duplicates

As the step says, this is to mark duplicate reads.

**What are duplicate reads ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_markdup1.md)

**What are they caused by ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_markdup2.md)


**What are the ways to detect them ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_markdup3.md)

Here we will use picards approach:

```
java -Xmx8G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/NA12878/NA12878.realigned.sorted.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.bam \
  METRICS_FILE=alignment/NA12878/NA12878.sorted.dup.metrics
```

We can look in the metrics output to see what happened.

```
less alignment/NA12878/NA12878.sorted.dup.metrics
```
**How many duplicates were there ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_markdup4.md)

This is very low, we expect in general <2%.

We can see that it computed seperate measures for each library.

**Why is this important to do and not combine everything ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_markdup5.md)


### Recalibration

This is the last BAM cleaning up step.

The goal for this step is to try to recalibrate base quality scores. The vendors tend to inflate the values of the bases in the reads.
Also, this step tries to lower the scores of some biased motifs for some technologies.

It runs in 2 steps, 
1- Build covariates based on context and known snp sites
2- Correct the reads based on these metrics

```
java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -knownSites ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  -L 1:17700000-18100000 \
  -o alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -I alignment/NA12878/NA12878.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -BQSR alignment/NA12878/NA12878.sorted.dup.recalibration_report.grp \
  -o alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -I alignment/NA12878/NA12878.sorted.dup.bam
```


## Extract Metrics

Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.

### Compute coverage

If you have data from a capture kit, you should see how well your targets worked

Both GATK and BVATools have depth of coverage tools. We wrote our own in BVAtools because
- GATK was deprecating theirs, but they changed their mind
- GATK's is very slow
- We were missing some output that we wanted from the GATK's one (GC per interval, valid pairs, etc)

Here we'll use the GATK one

```
java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12878/NA12878.sorted.dup.recal.coverage \
  -I alignment/NA12878/NA12878.sorted.dup.recal.bam \
  -L  1:17700000-18100000

#### Look at the coverage
less -S alignment/NA12878/NA12878.sorted.dup.recal.coverage.sample_interval_summary
```

Coverage is the expected ~30x. 

summaryCoverageThreshold is a usefull function to see if your coverage is uniform.
Another way is to compare the mean to the median. If both are almost equal, your coverage is pretty flat. If both are quite different, that means something is wrong in your coverage. A mix of WGS and WES would show very different mean and median values.

### Insert Size

```
java -Xmx8G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#look at the output
less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.insertSize.tsv
```

**What is the insert size and the corresponding standard deviation ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_insert1.md)

**Is the insert-size important ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_insert2.md)

### Alignment metrics

For the alignment metrics, we used to use ```samtools flagstat``` but with bwa mem since some reads get broken into pieces, the numbers are a bit confusing.
You can try it if you want.

We prefer the Picard way of computing metrics

```
java -Xmx8G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12878/NA12878.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#### explore the results

less -S alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv

```

**What is the percent of aligned reads ?** [solution](https://github.com/mbourgey/CBW_BFX_GENOMIC_MEDECINE_module2/blob/master/solutions/_alnMetrics1.md)

## Summary

In this lab, we aligned reads from the sample NA12878 to the reference genome `hg19`:

- We became familiar with FASTQ and SAM/BAM formats. 
- We checked read QC with BVAtools.
- We trimmed unreliable bases from the read ends using Trimmomatic. 
- We aligned the reads to the reference using BWA. 
- We sorted the alignments by chromosome position using PICARD. 
- We realigned short indels using GATK. 
- We fixed mate issues using PICARD.
- We recalibrate the Base Quality using GATK.
- We generate alignment metrics using GATK and PICARD.  


## Acknowledgments

I would like to thank and acknowledge Louis Letourneau for this help and for sharing his material. The format of the tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI.
