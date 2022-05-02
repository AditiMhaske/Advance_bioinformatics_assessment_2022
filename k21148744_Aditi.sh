
#!/bin/bash #



## Making directories for the NGS pipeline
## following are the directories we are making to save the raw data, trimmed data and results. We will also make log directory to save the commands or errors which
## we might come across

mkdir ngs_course
mkdir ngs_course/dnaseq
cd ngs_course/dnaseq
mkdir data meta results log
ls -lF
cd ~/ngs_course/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

## Downloading the raw data from the provided links using wget

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

## downloading the annotation file we will need for the alignment

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

## Moving the raw fastq files to the untrimmed fastq directory and the annotation.bed file to the data directory
mv *fastq.qz ~/ngs_course/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/ngs_course/dnaseq/data

## Downloading the reference file hg19 and moving the downloaded file to the data directory

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/ngs_course/dnaseq/data/

## Downloading the Miniconda installer for doing the NGS pipeline and tools within the miniconda package.
cd ~/
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./ Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib

## checking the contents of the untrimmed fastq
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
ls -lart

## using zcat to decompress and concatenate the .fastq.qz to .fastq files
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq

## doing the quality assessment for the raw files usimg the fastqc tool
fastqc NGS0001.R1.fastq NGS0001.R2.fastq

## making directory to save our fastqc results
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

cd ~/ngs_course/dnaseq/data/untrimmed_fastq

## trimmming the adapters from the raw data

trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq \
  -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R \
 ILLUMINACLIP:/home/ubuntu//miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
 TRAILING:25 MINLEN:50




## doing quality assessment again on the trimmed files
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P

## making directory for the reference file

mkdir -p ~/ngs_course/dnaseq/data/reference
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/

## indexing the low divergent seqeunces in fasta format
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
ls ~/ngs_course/dnaseq/data/reference

## making space to run the remaining pipeline

rm -r ~/ngs_course/dnaseq/data/untrimmed_fastq

## making directory for saving the aligned data
mkdir ~/ngs_course/dnaseq/data/aligned_data

## The BWA-MEM algorithm performs local alignment

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

cd ~/ngs_course/dnaseq/data/aligned_data

## using samtools to convert the sam file to binary version .bam
samtools view -h -b NGS0001.sam > NGS0001.bam

## sortimg the .bam file
samtools sort NGS0001.bam > NGS0001_sorted.bam

## indexing the sorted bam file
samtools index NGS0001_sorted.bam

cd ~/ngs_course/dnaseq/data/aligned_data

## alignment statistics using samtools
samtools flagstat NGS0001_sorted.bam
samtools idxstats NGS0001_sorted.bam

## deleting the sam file to more space for the analysis
cd ~/ngs_course/dnaseq/data/aligned_data
rm NGS0001.sam

## Using picard to locate and tag duplicate reads in the sorted BAM file
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

## indexing the sorted marked .bam file
samtools index NGS0001_sorted_marked.bam

## Setting Minimum MAPQ quality score to 20 and bitwise flag to 1796
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

## indexing the soreted.filtered.bam file
samtools index NGS0001_sorted_filtered.bam

## decompressing the hg19 genome file
zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa

## converting the reference which is in the text format, index it with samtools faidx
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

## calling variants with freebayes

freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf
~/ngs_course/dnaseq/results/NGS0001.vcf

##compressing the vcf file
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

## filtering the vcf file based on the quality scores
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

cd ~/ngs_course/dnaseq/data/aligned_data

## Using bedtools intersect to screen the overlaps between the vcf file and the provided annotation bed file
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf -b ../annotation.bed > ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## compressing the vcf file
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## Indexig the vcf with tabix
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz

## making space by deleting the trimmed fastq
rm -r ~/ngs_course/dnaseq/data/trimmed_fastq

## We downloaded the annovar by registering on the annovar website and downloaded the the .tar.gz file from the link provided by email. This file was then ## transferred to the terminal using filezilla
## Using annovar for annotation of the vcf file
cd ~/

## unzipping the annovar file
tar -zxvf annovar.latest.tar.gz

## downloading the databases required for annotation through annovar
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/

## convert VCF file to ANNOVAR inputfile without losing any VCF-specific information

./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered_R_vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_R.avinput

## table_annovar takes the annovar input file and generate a tab-delimited output file with many columns, each representing one set of annotation

./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_R.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

