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
# Compute TPMs
assays(brain)$TPM <- recount::getTPM(brain)
assays(blood)$TPM <- recount::getTPM(blood)
assays(liver)$TPM <- recount::getTPM(liver)

#---- Data Overview (Brain only) ----

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

#----Clean up factors causing bias (OPTIONAL)----
#function to clean up factors potentially causing bias
filter_data <- function(rse){
  #challenge 1a: remove genes annotated on mitochondrial DNA
  rse<- rse[rowRanges(rse)@seqnames !="chrM"]
  
  #challenge 1b and 1c: remove pseudogenes and rRNA genes
  rse<- rse[((rowData(rse)$gbkey!="Gene" & rowData(rse)$gbkey!="rRNA") | is.na(rowData(rse)$gbkey))]
  
  #challenge 1d: remove genes with annotated length <200
  rse<- rse[rowData(rse)$bp_length>=200]
  
  return(rse)
}

canonchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
             "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
             "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
             "chrX", "chrY", "chrM")
brain_filtered<- brain[rowRanges(brain)@seqnames %in% canonchr]
blood_filtered<- blood[rowRanges(blood)@seqnames %in% canonchr]
liver_filtered<- liver[rowRanges(liver)@seqnames %in% canonchr]

brain_filtered <- filter_data(brain)
blood_filtered <- filter_data(blood)
liver_filtered <- filter_data(liver)

#TODO: add some details to the challenge, see forum

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
#define how the replicates are grouped
group <- as.factor(c("Brain","Brain","Brain",
                     "Blood","Blood","Blood",
                     "Liver","Liver","Liver"))

y$samples$group <- group

#add information to the samples regarding the quality
y$samples$rin <- as.factor(c(colData(brain[,brain_cols])$gtex.smrin,
                             colData(blood[,blood_cols])$gtex.smrin,
                             colData(liver[,liver_cols])$gtex.smrin))

y$samples$slice <- as.factor(c(colData(brain[,brain_cols])$gtex.smtsd,
                               colData(blood[,blood_cols])$gtex.smtsd,
                               colData(liver[,liver_cols])$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(brain[,brain_cols])$gtex.sex,
                             colData(blood[,blood_cols])$gtex.sex,
                             colData(liver[,liver_cols])$gtex.sex))

y$samples$age <- as.factor(c(colData(brain[,brain_cols])$gtex.age,
                             colData(blood[,blood_cols])$gtex.age,
                             colData(liver[,liver_cols])$gtex.age))

y$samples$rRNA <- as.factor(c(colData(brain[,brain_cols])$gtex.smrrnart,
                              colData(blood[,blood_cols])$gtex.smrrnart,
                              colData(liver[,liver_cols])$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(brain[,brain_cols])$"recount_qc.star.uniquely_mapped_reads_%_both", 
                                colData(blood[,blood_cols])$"recount_qc.star.uniquely_mapped_reads_%_both",
                                colData(liver[,liver_cols])$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(brain[,brain_cols])$"recount_qc.aligned_reads%.chrm", 
                              colData(blood[,blood_cols])$"recount_qc.aligned_reads%.chrm",
                              colData(liver[,liver_cols])$"recount_qc.aligned_reads%.chrm"))
y

#genes with zero counts
table(rowSums(y$counts==0)==9)

#remove genes with zero or low expression
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

#----Normalization----
#log cpm before normalization
logcpm_before <- cpm(y, log=TRUE)
#normalization
y <- calcNormFactors(y, method = "TMM")
#log cpm after normalization
logcpm_after <- cpm(y, log=TRUE)
#boxplots
#TODO: make pretty
boxplot(logcpm_before, notch=T)
boxplot(logcpm_after, notch=T)

#----Linear Model----
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

#visualize how samples cluster together
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)
#label the samples by rRNA
plotMDS(logcpm, labels=y$samples$rRNA)
#label by percent of aligned reads
plotMDS(logcpm, labels=y$samples$chrm)
#label by age
plotMDS(logcpm, labels=y$samples$age)

#plot biological coefficient of variation
y <- estimateDisp(y, design)
plotBCV(y)

#fit linear model
fit <- glmQLFit(y, design)
#brain (top) vs blood (bottom)
qlfBBL <- glmQLFTest(fit, contrast=c(-1,1,0))
#liver (top) vs blood (bottom)
qlfLBL <- glmQLFTest(fit, contrast=c(-1,0,1))
#liver (top) vs brain (bottom)
qlfLB <- glmQLFTest(fit, contrast=c(0,-1,1))

qlfBBL

#DE genes, sorted by p-value, with the FDR added
topTags(qlfBBL, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfLBL, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfLB, n=10,adjust.method = "BH", sort.by = "PValue")

#get the whole table and write to a file
resultsBBL <- topTags(qlfBBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsLBL <- topTags(qlfLBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsLB <- topTags(qlfLB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsBBL, "resultsBBL.txt")
write.table(resultsLBL, "resultsLBL.txt")
write.table(resultsLB, "resultsLB.txt")

#summary of DE genes
summary(decideTests(qlfBBL, p.value=0.05, adjust.method = "BH", lfc=0)) #can adjust p value and log fold change
summary(decideTests(qlfLBL, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfLB, p.value=0.05, adjust.method = "BH", lfc=0))

#find genes that are upregulated 1v3, ie, upregulated in one compared to both
#first find row names where it is upregulated against one 
brain_up_idx <- which(topTags(qlfBBL, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC>0)
brain_up_genes <- row.names(topTags(qlfBBL, n = 10000000, adjust.method = "BH", 
                                    sort.by = "PValue", p.value = 0.01)[["table"]][brain_up_idx,])
#find the rows where it is upregulated against the other
brain_up_idx <- which(topTags(qlfLB, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC<0)
#perform the intersection of the two
brain_up_genes <- intersect(brain_up_genes, 
                            row.names(topTags(qlfLB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)[["table"]][brain_up_idx,]))

#repeat for blood and liver
blood_up_idx <- which(topTags(qlfBBL, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC<0)
blood_up_genes <- row.names(topTags(qlfBBL, n = 10000000, adjust.method = "BH", 
                                    sort.by = "PValue", p.value = 0.01)[["table"]][blood_up_idx,])
#find the rows where it is upregulated against the other
blood_up_idx <- which(topTags(qlfLBL, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC<0)
#perform the intersection of the two
blood_up_genes <- intersect(blood_up_genes, 
                            row.names(topTags(qlfLBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)[["table"]][blood_up_idx,]))

liver_up_idx <- which(topTags(qlfLBL, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC>0)
liver_up_genes <- row.names(topTags(qlfLBL, n = 10000000, adjust.method = "BH", 
                                    sort.by = "PValue", p.value = 0.01)[["table"]][liver_up_idx,])
#find the rows where it is upregulated against the other
liver_up_idx <- which(topTags(qlfLB, n = 10000000, adjust.method = "BH", 
                              sort.by = "PValue", p.value = 0.01)[["table"]]$logFC>0)
#perform the intersection of the two
liver_up_genes <- intersect(liver_up_genes, 
                            row.names(topTags(qlfLB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)[["table"]][liver_up_idx,]))

#check the expression across all samples, not just my assigned ones
which(rowData(brain)$gene_name == "UGT1A1") #gene found in liver cells
boxplot(assays(brain)$TPM[20205,],assays(blood)$TPM[20205,], assays(liver)$TPM[20205,], outline=F )
#use a statistical test to prove this expression across all samples is significant
wilcox.test(assays(brain)$TPM[20205,], assays(blood)$TPM[20205,], alternative='two.sided',exact = F,correct = T)
wilcox.test(assays(brain)$TPM[20205,], assays(liver)$TPM[20205,], alternative='two.sided',exact = F,correct = T)
wilcox.test(assays(blood)$TPM[20205,], assays(liver)$TPM[20205,], alternative='two.sided',exact = F,correct = T)
#TODO: adjust so that the comparison is liver vs blood + brain