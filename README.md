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
```

## FASTQC

## STAR

### Genome index

```
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir Homo_sapiens.GRCh38.release.87_GENECODE.v25/ --genomeFastaFiles /media/emaglinux/0DBF12730DBF1273/DATA/Genome/FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /media/emaglinux/0DBF12730DBF1273/DATA/Genome/Annotation/gencode.v25.annotation.gtf --sjdbOverhang 88
```
