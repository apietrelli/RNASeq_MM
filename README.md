# RNASeq_MM

## Sample Description

Samples sequenced on **DATE**

The Rna-Seq Illumina sequencing Library used was:

** Total RNA Stranded (read2 strand) **

With this library the **Read 2** follows the sense of the transcript instead of
the first read

*NB The strand of the library does not follow the sense of the transcript*


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

# NOT WORKING:
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

# ATTEMPT2
```
# http://stackoverflow.com/questions/9450394/how-to-install-gcc-piece-by-piece-with-gmp-mpfr-mpc-elf-without-shared-libra

# gcc infrastructure
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/

# install GMP
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-4.3.2.tar.bz2
bunzip2 gmp-4.3.2.tar.bz2
tar xvf gmp-4.3.2.tar
cd gmp-4.3.2
./configure --prefix=/storage/luca/tmp/gcc
make && make check && make install

# install MPFR
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-2.4.2.tar.bz2
bunzip2 mpfr-2.4.2.tar.bz2
tar xvf mpfr-2.4.2.tar
cd mpfr-2.4.2
./configure --prefix=/storage/luca/tmp/gcc --with-gmp=/storage/luca/tmp/gcc
make && make check && make install

# install MPC
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-0.8.1.tar.gz
tar zxvf mpc-0.8.1.tar.gz
cd mpc-0.8.1
./configure --prefix=/storage/luca/tmp/gcc --with-gmp=/storage/luca/tmp/gcc --with-mpfr=/storage/luca/tmp/gcc
make && make check && make install

# install ELF
wget http://www.mr511.de/software/libelf-0.8.13.tar.gz
tar zxvf libelf-0.8.13.tar.gz
cd libelf-0.8.13
./configure --prefix=/storage/luca/tmp/gcc
make && make check && make install

# install GCC - different scratch dir in/storage/luca   
mkdir tmp2
wget ftp://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-4.8.5/gcc-4.8.5.tar.gz
tar zxvf gcc-4.8.5.tar.gz
mkdir -p tmp2/gcc-4.8.5-scratch
cd tmp2/gcc-4.8.5-scratch

 ./../gcc-4.8.5/configure --disable-libstdcxx-pch --enable-languages=all --enable-libgomp --enable-lto --enable-threads=posix --enable-tls --with-gmp=/storage/luca/tmp/gcc --with-mpfr=/storage/luca/tmp/gcc --with-mpc=/storage/luca/tmp/gcc --with-libelf=/storage/luca/tmp/gcc --with-fpmath=sse
make && make install

# error >> checking for suffix of object files... configure: error: in `/home/manu/gcc/gcc/i686-pc-linux-gnu/libgcc':
# see https://gcc.gnu.org/wiki/FAQ#configure_suffix
# try suggested solution
./contrib/download_prerequisites
mkdir gcc-build
# Then run the configure either by fully qualified path or by relative path while in the the gcc-build current working directory.
# A makefile will be created in the gcc-build directory. Run make in the gcc-build current working directory to begin the build of GCC.

# In the meanwhile, Alex Dobin indicized genome for us! So kind.
http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/GENCODE/GRCh38_Gencode25/

```


### Mapping

```
# dir on EMAGLINUX

cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/STAR_mapping

# load files from Alex Dobin's genome index directory
STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeLoad LoadAndRemove --runThreadN 8 --readFilesIn ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R1.fastq.gz ../../FASTQ.files/Sample_MM-431/Sample_MM-431_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix Sample_MM-431_ --outSAMunmapped Within --outBAMsortingThreadN 8 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outWigType bedGraph > Sample_MM-431.STAR_mapping.log &

# cycle for all fastq WITH bedGraph output
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/STAR_mapping
ls ../../FASTQ.files/ > elenco
head -10 elenco > elenco.part1
for i in `cat elenco.part1` ; do
  echo "Mapping $i" ;
  mkdir "$i";
  STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeLoad LoadAndKeep --runThreadN 8 --readFilesIn ../../FASTQ.files/"$i"/"$i"_R1.fastq.gz ../../FASTQ.files/"$i"/"$i"_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix "$i"_  --outSAMunmapped Within --outBAMsortingThreadN 8 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outWigType bedGraph > "$i".STAR_mapping.log ;
  echo "$i mapped" ;
done > RNA-SEQ_30MM.STARmapping.13012016.log &

# cycle for all fastq WITHOUT bedGraph output
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/STAR_mapping
ls ../../FASTQ.files/ > elenco
head -10 elenco > elenco.part1
for i in `cat elenco.attempt3` ; do
  echo "Mapping $i" ;
  date ;
  mkdir -p "$i" ;
  cd "$i" ;
  STAR --genomeDir /media/emaglinux/0DBF12730DBF1273/DATA/STAR/Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeLoad LoadAndKeep --runThreadN 8 --readFilesIn ../../../FASTQ.files/"$i"/"$i"_R1.fastq.gz ../../../FASTQ.files/"$i"/"$i"_R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix "$i"_  --outSAMunmapped Within --outBAMsortingThreadN 8 --outSAMtype BAM Unsorted > "$i".STAR_mapping.log ;
  echo "$i mapped" ;
  cd .. ;
  date ;
done > RNA-SEQ_30MM.STARmapping.13012016.log &

# create elenco.part.2
ls ../../FASTQ.files/ > elenco
tail -22 elenco > elenco.part2
sed '/Sample_MM-431/d' elenco.part2 > temp ; mv temp elenco.part2
# relaunch previuos w elenco.part2
### error >> could not shut down?
cat elenco.part1 elenco.part2 > elenco.attempt3

```

### Sorting and indexing BAM files

---

#### Example with Sample_KMS-11

- Sorting Procedure using `samtools`

Limiting memory size to 200Mb to better handle the RAM limit because the sorting
process is RAM consuming

```
samtools sort -@ 5 -m 200M Sample_KMS-11_Aligned.out.bam Sample_KMS-11_Aligned.out.sorted
[bam_sort_core] merging from 140 files...
```
Time elapsed **~ 33min**

The BAM SortedByCoordinate is largely reduced in size respect to the Unsorted one

Sample_KMS-11_Aligned.out.bam | Sample_KMS-11_Aligned.out.sorted.bam
--- | ---
9.2G | **6.1G**

NB : *Think about removing the Unsorted after sorting completed*

- Indexing always using `samtools`

```
samtools index Sample_KMS-11_Aligned.out.sorted.bam
```

- Create BigWig using `deeptools` in particular with the `bedCoverage` command
following the [manual] (http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)
for the usage

No Normalization

```
# Reverse Strand
bamCoverage -p 3 -b Sample_KMS-11_Aligned.out.sorted.bam --filterRNAstrand reverse -o Sample_KMS-11.rev.bw &

# Forward Strand
bamCoverage -p 3 -b Sample_KMS-11_Aligned.out.sorted.bam --filterRNAstrand forward -o Sample_KMS-11.fwd.bw
```
With RPKM normalization allowed by bamCoverage
```
# Reverse Strand
bamCoverage -p 7 -b Sample_KMS-11_Aligned.out.sorted.bam --normalizeUsingRPKM --filterRNAstrand reverse -o Sample_KMS-11.RPKM.rev.bw

# Forward Strand
bamCoverage -p 7 -b Sample_KMS-11_Aligned.out.sorted.bam --normalizeUsingRPKM --filterRNAstrand forward -o Sample_KMS-11.RPKM.fwd.bw
```

**-p** is the number of processor used

#### Cycle for all the samples

```
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/STAR_mapping/
for i in `ls`; do
  echo "Start $i" `date`
  echo "Sorting BAM file.."
  samtools sort -@ 5 -m 3GB "$i"/"$i"_Aligned.out.bam "$i"/"$i".out.sorted
  echo "Indexing sorted file.."
  samtools index "$i"/"$i".out.sorted.bam
  echo "Creating bigWig.."
  bamCoverage -p 3 -b "$i"/"$i".out.sorted.bam --normalizeUsingRPKM --filterRNAstrand forward -o "$i"/"$i".RPKM.fwd.bw &
  bamCoverage -p 3 -b "$i"/"$i".out.sorted.bam --normalizeUsingRPKM --filterRNAstrand reverse -o "$i"/"$i".RPKM.rev.bw
  echo "End sample"
done

```

### IGV

```
# check if remove unannotated
# test whether STAR create BWig or similar formats
# or: convert BAM to BED

```

## Counts per gene

### HTSEq

```
# bam files sorted by pos. stranded

htseq-count -r pos -f bam ../STAR_mapping/Sample_KMS-11/Sample_KMS-11_Aligned.out.sorted.bam /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf 2> Sample_KMS-11.count.stats > Sample_KMS-11.count &

htseq-count -r pos -s reverse -f bam ../STAR_mapping/Sample_KMS-11/Sample_KMS-11.out.sorted.bam /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf 2> Sample_KMS-11.REV.count.stats > Sample_KMS-11.REV.count &

cut -f1,7 Sample_KMS-11.union.featureCounts > Sample_KMS-11.featureCounts

### check warning:
### Warning: Mate pairing was ambiguous for 24853 records; mate key for first such record: ('HWI-1KL157:151:C4H67ACXX:5:1104:14343:23511', 'first', 'chr1', 28530, 'chr1', 28741, 312).
### 44056929 SAM alignment pairs processed.

```

### featureCounts

SuperSpeed!!!
No cycle needed, because it takes multiple input file on the command-line

```
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/
mkdir featureCounts

# Create sample index for featureCounts
ls STAR_mapping/*/*.bam > Samples_list.featureCounts

# Run featureCounts multiple samples
featureCounts -p -a /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf -s 2 -T 5 -O -o featureCounts/RNA-SEQ_30MM.paired.union.featureCounts `cat Samples_list.featureCounts` &> RNA-SEQ_30MM.featureCounts.log &

# Run modified GRCh38_Gencode25 (gencode.v25.annotation.v2.gtf) WITH WHSC2-2
cd /media/emaglinux/0DBF12730DBF1273/Rshared/RNA-SEQ_30MM/Analisi/
mkdir featureCounts.v2
cp /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.v2.gtf
### added lncWHSC-2 entry
### REMEMBER: in gtf file mandatory structure is: gene > transcript > exon entries

featureCounts -p -a /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.v2.gtf -s 2 -T 5 -O -o featureCounts.v2/RNA-SEQ_30MM.paired.union.featureCounts.v2 `cat Samples_list.featureCounts` &> RNA-SEQ_30MM.featureCounts.v2.log &
### try on single sample unstranded
featureCounts -p -a /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.v2.gtf -T 5 -O -o featureCounts.v2/Sample_MM-263.featureCounts.v2 STAR_mapping/Sample_MM-263/Sample_MM-263.out.sorted.bam &
```
### DESEQ

### normalization

Check mean fragment size
```
samtools view sample.bam | awk '{print $9}' | sort | uniq -c
```
### RSEM

Install RSEM and prepare reference genome. By default RSEM uses Bowtie. Use star to customize alignments.  
```
rsem-prepare-reference --gtf ./../../../../DATA/Genome/Annotation/gencode.v25.annotation.gtf -p 6 --star ./../../../../DATA/Genome/FASTA/GRCh38.primary_assembly.genome.fa ref/Gencode.v25
### need running star only if re-run alignments on fastq files

```
Calculate expression based on STAR_mapping
```
rsem-calculate-expression --paired-end --alignments -p 6 --bam Sample_KMS-11/Sample_KMS-11.out.sorted.bam ref/Gencode.v25 Sample_KMS-11.out.sorted  
### doesn't work > should re-run star with '--quantMode TranscriptomeSAM' option?
```                         


### FASTQ BACKUP ON NAS


Move FASTQ with "mv" command

*NB Do as SuperUser*

```
sudo su

for i in `ls` ; do
   echo "$i"_R1;
   split -b2G --numeric-suffixes=1 "$i"/"$i"_R1.fastq.gz "$i"/"$i"_R1.fastq.gz.part- ;
   echo "$i"_R2 ;
   split -b2G --numeric-suffixes=1 "$i"/"$i"_R2.fastq.gz "$i"/"$i"_R2.fastq.gz.part- ;
   mkdir /media/NAS/RNAseq.30MM/FASTQ/"$i" ;
   mv "$i"/*.part* /media/NAS/RNAseq.30MM/FASTQ/"$i" ;
done
```

# Note and others

### rsync syntax

```
#   rsync -ah --progress --exclude "*.gz" "$i"/ /media/NAS/RNAseq.30MM/FASTQ/ ;
```
