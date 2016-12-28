# RNASeq_MM

## Install packages

### STAR installation on EMACLINUX

Do all the command as SuperUser

```
cd /opt
git clone https://github.com/alexdobin/STAR.git

# Compile it
cd STAR/source
sudo make

# Add /opt to PATH
cd ~
vi .bashrc
# Add this line to the bottom of the files
PATH=$PATH:/opt
```

## Project Procedures

### Unify the reads chunks deriving from PTP sequencing

```
# Entering in the FASTQ dir
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/FASTQ.files

# Aggregate READ 1
for i in `ls`; do
  echo "Entering $i dir";
  cd $i;
  cat *R1*.fastq.gz > "$i"_R1.fastq.gz ;
  echo "$i finished"  ;
  cd .. ;
done

# Aggregate READ 2
for i in `ls`; do
  echo "Entering $i dir";
  cd $i;
  cat *R2*.fastq.gz > "$i"_R2.fastq.gz ;
  echo "$i finished"  ;
  cd .. ;
done

# Remove previous fastq files
for i in `ls`; do
  echo "Entering $i dir";
  cd $i;
  rm *R1_* ; rm *R2_* ;
  echo "$i finished"  ;
  cd .. ;
done
```

## FASTQC

## STAR

### Genome index

```
# run on xlabserver6
wget ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# STAR --runThreadN 1 --runMode genomeGenerate --genomeDir Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /media/emaglinux/0DBF12730DBF1273/DATA/Genome/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf --sjdbOverhang 88

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/fastq_luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/fastq_luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /storage/fastq_luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88 >> stderr_STAR.log &

#Fatal INPUT FILE error, no valid exon lines in the GTF file: /storage/fastq_luca/genome/gencode.v25.annotation.gtf
#Solution: check the formatting of the GTF file. Most likely cause is the difference in chromosome naming between GTF and FASTA file.

# Rerun after changing fast indexing from "1" to "chr1", etc
sed '/>[1-9][XYM]/s/>/>chr/g' Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa > Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa
sed 's/chrMT/chrM/g' Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa > temp.fa; mv temp.fa Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa 

# run on xlabserver4

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /storage/luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88 >> stderr_STAR.log &
## Fatal INPUT FILE error, no valid exon lines in the GTF file: /storage/luca/genome/gencode.v25.annotation.gtf
```

### Mapping

```
STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/homo_sapiens.releas.INGM/ --runThreadN 8 --readFilesIn ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R1.fastq.gz ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix Sample_MM-431_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMunmapped Within --outBAMsortingThreadN 8 > Sample_MM-431.STAR_mapping.log &
```
