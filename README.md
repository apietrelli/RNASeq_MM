# RNASeq_MM

## Install packages

### STAR installation on EMAGLINUX

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

### MULTIQC

```
# summarize analysis FASTQC / STAR / DESEQ-HTSEQ etc
multiqc .
```


### Genome index

```
# run on xlabserver6
wget ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# STAR --runThreadN 1 --runMode genomeGenerate --genomeDir Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /media/emaglinux/0DBF12730DBF1273/DATA/Genome/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf --sjdbOverhang 88

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/fastq_luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/fastq_luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /storage/fastq_luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88 >> stderr_STAR.log &

#Fatal INPUT FILE error, no valid exon lines in the GTF file: /storage/fastq_luca/genome/gencode.v25.annotation.gtf
#Solution: check the formatting of the GTF file. Most likely cause is the difference in chromosome naming between GTF and FASTA file.

# Rerun after changing fasta indexing from "1" to "chr1", etc
sed '/>[1-9][XYM]/s/>/>chr/g' Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa > Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa
sed 's/chrMT/chrM/g' Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa > temp.fa; mv temp.fa Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/fastq_luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/fastq_luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa --sjdbGTFfile /storage/fastq_luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88 &

# run on xlabserver4

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /storage/luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88 >> stderr_STAR.log &
## Fatal INPUT FILE error, no valid exon lines in the GTF file: /storage/luca/genome/gencode.v25.annotation.gtf
## change fasta indexing file

./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /storage/luca/genome/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /storage/luca/genome/Homo_sapiens.GRCh38.dna.primary_assembly.renamed.fa --sjdbGTFfile /storage/luca/genome/gencode.v25.annotation.gtf --sjdbOverhang 88
```
### the pre-compiled binary on server is too much memory consuming >>> try compiling from source using a stand-alone gcc installation

```
wget http://ftp.gnu.org/gnu/gsrc/gsrc-2014.10.11.tar.gz
tar -zxvf gsrc-2014.10.11.tar.gz
cd gsrc-2014.10.11/
./bootstrap                       # to create the configure script
./configure --prefix=$HOME/gnu    # --prefix is directory to install the packages
                                  # Pick your --prefix by your wishes.
. ./setup.sh                      # This just sets some ENV variables and appends to PATH
                                 # and other variables to allow GSRC to work seamlessly.
                                     # Put this line in your .bashrc.


cd gsrc
make -C gnu/gcc MAKE_ARGS_PARALLEL="-j12"
#(or make -C gnu/gcc MAKE_ARGS_PARALLEL="-jN" to speed up for a N-core system)
### error extracting /storage/luca/gsrc-2014.10.11/gnu/gcc/download/gcc-4.9.1.tar.gz
### do it manually. create "work" dir
tar -zxvf /storage/luca/gsrc-2014.10.11/gnu/gcc/download/gcc-4.9.1.tar.gz  --no-same-owner --no-same-permissions -C work

### doesn't work!

make -C gnu/gcc install
```

### Mapping

```

### sample command lines
#  STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/homo_sapiens.release.INGM/ --runThreadN 8 --readFilesIn ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R1.fastq.gz ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix Sample_MM-431_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMunmapped Within --outBAMsortingThreadN 8 --outSAMtype BAM > Sample_MM-431.STAR_mapping.log &

# cycle for all fastq
# cd ../../FASTQ.files/
for i in `ls`; do
  echo "Entering $i dir";
  cd $i;
  STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/homo_sapiens.release.GRCh38.29122016/ --runThreadN 8 --readFilesIn "$i"_R1.fastq.gz "$i"_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix "$i"_ --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMunmapped Within --outBAMsortingThreadN 8 --outSAMtype BAM > "$i".STAR_mapping.log ;
  echo "$i mapped"  ;
  cd .. ;
done &
```

### IGV

```
# check if remove unannotated
# test whether STAR create BWig or similar formats
# or: convert BAM to BED
```

### HTSEq

### DESEQ
