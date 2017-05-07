#set up
export REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/
export WORK_DIR="${HOME}/bfx_genomic_medecine/module2"

module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/samtools/1.4 
module load mugqic/bwa/0.7.12 mugqic/GenomeAnalysisTK/3.7 mugqic/picard/1.123 
module load mugqic/trimmomatic/0.36 mugqic/R_Bioconductor/3.3.2_3.4


cd $WORK_DIR

#NA12891

mkdir -p originalQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/


mkdir -p reads/NA12891/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12891/NA12891_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12891/NA12891_CBW_chr1_R2.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12891/NA12891.trim.out


mkdir -p postTrimQC/
java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/

mkdir -p alignment/NA12891/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12891\tSM:NA12891\tLB:NA12891\tPU:runNA12891_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  reads/NA12891/NA12891_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12891/NA12891_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_HOME}/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12891/NA12891.sorted.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000


# Say you want to count the *un-aligned* reads, you can use
samtools view -c -f4 alignment/NA12891/NA12891.sorted.bam

# Or you want to count the *aligned* reads you, can use
samtools view -c -F4 alignment/NA12891/NA12891.sorted.bam

java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12891/realign.intervals \
  -I alignment/NA12891/NA12891.sorted.bam \
  -L 1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/NA12891/realign.intervals \
  -o alignment/NA12891/NA12891.realigned.sorted.bam \
  -I alignment/NA12891/NA12891.sorted.bam


java -Xmx8G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/NA12891/NA12891.realigned.sorted.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.bam \
  METRICS_FILE=alignment/NA12891/NA12891.sorted.dup.metrics


java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -knownSites ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  -L 1:17700000-18100000 \
  -o alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -I alignment/NA12891/NA12891.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -BQSR alignment/NA12891/NA12891.sorted.dup.recalibration_report.grp \
  -o alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -I alignment/NA12891/NA12891.sorted.dup.bam
  

java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12891/NA12891.sorted.dup.recal.coverage \
  -I alignment/NA12891/NA12891.sorted.dup.recal.bam \
  -L  1:17700000-18100000


java -Xmx8G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12891/NA12891.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12891/NA12891.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY


java -Xmx8G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12891/NA12891.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12891/NA12891.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY

#NA12892

java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  --read2 raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  --threads 1 --regionName ACTL8 --output originalQC/


mkdir -p reads/NA12892/

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 \
  raw_reads/NA12892/NA12892_CBW_chr1_R1.fastq.gz \
  raw_reads/NA12892/NA12892_CBW_chr1_R2.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_S2.t20l32.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 \
  2> reads/NA12892/NA12892.trim.out



java -Xmx8G -jar ${BVATOOLS_JAR} readsqc \
  --read1 reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  --read2 reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  --threads 1 --regionName ACTL8 --output postTrimQC/

  
mkdir -p alignment/NA12892/

bwa mem -M -t 2 \
  -R '@RG\tID:NA12892\tSM:NA12892\tLB:NA12892\tPU:runNA12892_1\tCN:Broad Institute\tPL:ILLUMINA' \
  ${REF}/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  reads/NA12892/NA12892_CBW_chr1_R1.t20l32.fastq.gz \
  reads/NA12892/NA12892_CBW_chr1_R2.t20l32.fastq.gz \
  | java -Xmx8G -jar ${PICARD_HOME}/SortSam.jar \
  INPUT=/dev/stdin \
  OUTPUT=alignment/NA12892/NA12892.sorted.bam \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1500000

java -Xmx8G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12892/realign.intervals \
  -I alignment/NA12892/NA12892.sorted.bam \
  -L 1

java -Xmx8G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -targetIntervals alignment/NA12892/realign.intervals \
  -o alignment/NA12892/NA12892.realigned.sorted.bam \
  -I alignment/NA12892/NA12892.sorted.bam

java -Xmx8G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/NA12892/NA12892.realigned.sorted.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.bam \
  METRICS_FILE=alignment/NA12892/NA12892.sorted.dup.metrics

java -Xmx8G -jar ${GATK_JAR} \
  -T BaseRecalibrator \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -knownSites ${REF}/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  -L 1:17700000-18100000 \
  -o alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -I alignment/NA12892/NA12892.sorted.dup.bam

java -Xmx8G -jar ${GATK_JAR} \
  -T PrintReads \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -BQSR alignment/NA12892/NA12892.sorted.dup.recalibration_report.grp \
  -o alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -I alignment/NA12892/NA12892.sorted.dup.bam

java  -Xmx8G -jar ${GATK_JAR} \
  -T DepthOfCoverage \
  --omitDepthOutputAtEachBase \
  --summaryCoverageThreshold 10 \
  --summaryCoverageThreshold 25 \
  --summaryCoverageThreshold 50 \
  --summaryCoverageThreshold 100 \
  --start 1 --stop 500 --nBins 499 -dt NONE \
  -R ${REF}/genome/Homo_sapiens.GRCh37.fa \
  -o alignment/NA12892/NA12892.sorted.dup.recal.coverage \
  -I alignment/NA12892/NA12892.sorted.dup.recal.bam \
  -L  1:17700000-18100000


java -Xmx8G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12892/NA12892.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.tsv \
  HISTOGRAM_FILE=alignment/NA12892/NA12892.sorted.dup.recal.metric.insertSize.histo.pdf \
  METRIC_ACCUMULATION_LEVEL=LIBRARY


java -Xmx8G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
  VALIDATION_STRINGENCY=SILENT \
  REFERENCE_SEQUENCE=${REF}/genome/Homo_sapiens.GRCh37.fa \
  INPUT=alignment/NA12892/NA12892.sorted.dup.recal.bam \
  OUTPUT=alignment/NA12892/NA12892.sorted.dup.recal.metric.alignment.tsv \
  METRIC_ACCUMULATION_LEVEL=LIBRARY



