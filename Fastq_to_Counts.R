#---------Bulk RNAseq Analysis

#Set up the working directory----
setwd("/home/marcelo/Downloads/Raw_files-20230220T104319Z-001/Raw_files")

#Loading required packages----
#install.packages("R.utils")
library(R.utils)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)

#install.packages("data.table")
library(data.table)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RUVSeq")
library(RUVSeq)
install.packages("XML")
#(if the above installation of RUVseq didn't work, try this)
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("DESeq2")
library(DESeq2)

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

install.packages("pheatmap")
library(pheatmap)

#install.packages("RColorBrewer")
library(RColorBrewer)

#install.packages("ggplot2")
library(ggplot2)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("Rqc")
library(Rqc)

#Downloading FastQ files----
url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_1.fastq.gz"
destination<-"SRR5924196_1.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_2.fastq.gz"
destination<-"SRR5924196_2.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_1.fastq.gz"
destination<-"SRR5924198_1.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_2.fastq.gz"
destination<-"SRR5924198_2.fastq.gz"
download.file(url,destination)


#Downloading Genome file----
url<-"ftp://ftp.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
download.file(url,destination)
gunzip(destination)

#Downloading GTF file----
url<-"ftp://ftp.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
download.file(url,destination)
gunzip(destination)

install.packages("fastqcr")
#Your system should also have JAVA installed.
#visit www.java.com for installation
library(fastqcr)
fastqc_install()
fastqc()
qc <- qc_aggregate(getwd())
qc


#Building Index----
buildindex("Sc_full_index_rsubread",
           "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
           indexSplit=F)

reads1 <- list.files(pattern = "_1.fastq.gz$" )
reads2 <- list.files(pattern = "_2.fastq.gz$" )
all.equal(length(reads1),length(reads2))

#Performing Alignment----
align(index="Sc_full_index_rsubread",
      readfile1=reads1,
      readfile2=reads2,
      input_format="gzFASTQ",
      output_format="BAM",
      nthreads=10)

#Checking the BAM files generated----
bam.files <- list.files(pattern = ".BAM$", full.names = TRUE)
bam.files

#Checking the mapping quality----
props <- propmapped(files=bam.files)
props

#Generating Feature counts----
fcLim <- featureCounts(files = bam.files,
                       annot.ext="Saccharomyces_cerevisiae.R64-1-1.96.gtf",
                       GTF.featureType="exon",
                       GTF.attrType="gene_id",
                       isGTFAnnotationFile=TRUE,
                       isPairedEnd = TRUE)

fc <- data.frame(fcLim[["counts"]])
colnames(fc) <- c("Normal", "Tumor")


