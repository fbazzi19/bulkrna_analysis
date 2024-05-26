#remember to set the working directory accordingly!

#----library calls ----
library(recount)
library(recount3)
library(edgeR)

#----load in and transform counts data----
brain <- readRDS("rse_brain.RDS")
head(brain)
blood <- readRDS("rse_blood.RDS")
head(blood)
liver <- readRDS("rse_liver.RDS")
head(liver)

#transform counts
assays(brain)$counts <- transform_counts(brain)
assays(blood)$counts <- transform_counts(blood)
assays(liver)$counts <- transform_counts(liver)


#---- Data Overview (Brain only) ----
# Compute TPMs
assays(brain)$TPM <- recount::getTPM(brain)

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
table(rowRanges(brain)@seqnames) #indicates that non-canonical chromosomes were included



#---- Overall Quality Plots ----
#TODO: make pretty
#brain
boxplot(colData(brain)$gtex.smrin)
boxplot(colData(brain)$gtex.smrrnart)
boxplot(colData(brain)$"recount_qc.star.uniquely_mapped_reads_%_both")
#blood
boxplot(colData(blood)$gtex.smrin)
boxplot(colData(blood)$gtex.smrrnart)
boxplot(colData(blood)$"recount_qc.star.uniquely_mapped_reads_%_both")
#liver
boxplot(colData(liver)$gtex.smrin)
boxplot(colData(liver)$gtex.smrrnart)
boxplot(colData(liver)$"recount_qc.star.uniquely_mapped_reads_%_both")

#----Clean up factors causing bias----
#function to clean up factors potentially causing bias
filter_data <- function(rse){
  #challenge 1a and 2: remove genes annotated on mitochondrial DNA
  #and keep genes only on canonical chromosomes
  canonchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY")
  rse<- rse[rowRanges(rse)@seqnames %in% canonchr]
  
  #challenge 1b and 1c: remove pseudogenes and rRNA genes
  rse<- rse[((rowData(rse)$gbkey!="Gene" & rowData(rse)$gbkey!="rRNA") | is.na(rowData(rse)$gbkey))]
  
  #challenge 1d: remove genes with annotated length <200
  rse<- rse[rowData(rse)$bp_length>=200]
  
  return(rse)
}

brain <- filter_data(brain)
blood <- filter_data(blood)
liver <- filter_data(liver)


#---- check assigned columns ----
#This functions starts at a given column and checks if that sample meets
#the required thresholds. It continues iteratively until it finds 3 samples.
get_my_samples <- function(rse, start){
  cols <- c()
  while (length(cols)<3){
    if (colData(rse)$gtex.smrin[start]>6 & colData(rse)$gtex.smrrnart[start]<0.1 & colData(rse)$"recount_qc.star.uniquely_mapped_reads_%_both"[start]>85){
      cols <- c(cols, start)
      print(paste0("Sample ", start, " has RIN ", colData(rse)$gtex.smrin[start]))
      print(paste0("Sample ", start, " has rRNA fraction ", colData(rse)$gtex.smrrnart[start]))
      print(paste0("Sample ", start, " has % mapped reads ", colData(rse)$"recount_qc.star.uniquely_mapped_reads_%_both"[start]))
    }
    start <- start+1
  }
  return (cols)
}

#call the function for each tissue type
brain_cols <- get_my_samples(brain, 30)
blood_cols <- get_my_samples(blood, 30)
liver_cols <- get_my_samples(liver, 30)


#---- Set up count table for edger ----
#count tables of selected columns
counts_brain_selected <- assays(brain[,brain_cols])$counts
counts_blood_selected <- assays(blood[,blood_cols])$counts
counts_liver_selected <- assays(liver[,liver_cols])$counts

#define the row names as the official gene names
rownames(counts_brain_selected) <- rowData(brain)$gene_name
rownames(counts_blood_selected) <- rowData(blood)$gene_name
rownames(counts_liver_selected) <- rowData(liver)$gene_name

#bind the count tables together
x <- cbind(counts_brain_selected,counts_blood_selected,counts_liver_selected)
#assign the column names
colnames(x) <- c("Brain30", "Brain31","Brain33",
                 "Blood30", "Blood31","Blood32",
                 "Liver31","Liver32","Liver33")

#----DGE Object----
y <- DGEList(counts = x)
