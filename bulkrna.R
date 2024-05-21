#----library calls ----
library(recount3)

#----brain data----
brain <- readRDS("C:/Users/meemo/Downloads/Genomics and Transcriptomics/Transcriptomics Exam/rse_brain.RDS")
head(brain)

#transform counts
assays(brain)$counts <- transform_counts(brain)
brain

#look at what information is contained within colData
names(colData(brain))

#observe the RIN of each sample
head(colData(brain)$gtex.smrin)
#sample sources
table(colData(brain)$gtex.smtsd)
#sample sexes
table(colData(brain)$gtex.sex)
#sample ages
table(colData(brain)$gtex.age)

#the number of reads in each sample
head(colData(brain)$"recount_qc.star.number_of_input_reads_both") #both because the sequencing was paired end
#the percentage of uniquely mapped reads
head(colData(brain)$'recount_qc.star.uniquely_mapped_reads_%_both')

#plotting the distribution of reads per sample
inputreads <- colData(brain)$"recount_qc.star.number_of_input_reads_both"
boxplot(inputreads)

#plotting the percentage of mapped reads
mapped <- colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped)

#estimated percent of reads coming from ribosomal RNA
boxplot(colData(brain)$gtex.smrrnart)
#estimated percent of reads coming from mitochondrial genes
boxplot(colData(brain)$"recount_qc.aligned_reads%.chrm")
#rRNA and mitochondrial genes
plot(colData(brain)$gtex.smrrnart, colData(brain)$"recount_qc.aligned_reads%.chrm")

#look at what information is contained within rowData
names(rowData(brain))

#types of genes
table(rowData(brain)$gbkey)
#genomic coordinates of each gene
rowRanges(brain)
#how many genes come from each chromosome
table(rowRanges(brain)@seqnames) #indicated that non-canonical chromosomes were included

#---- check assigned columns ----
#todo: make it a function