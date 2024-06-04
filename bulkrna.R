#remember to set the working directory accordingly!

#----library calls ----
library(recount)
library(recount3)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(limma)
library(ggrepel)
library(openxlsx)

#----Colors----
#color pallete to use for figure generation throughout
colorpallete <- c("#DA9C7E", "#66B386", "#E6D939", 
                  "#EDC2C2", "#B1C2FF", "#BCB0CE",
                  "#EB7B70", "#B06CC8", "#626CB2")

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

#removing genes with little/no annotation
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "LOC")]), ]
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "LINC")]), ]
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "MIR")]), ]
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "SNORD")]), ]
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "RPL")]), ] 
x <- x[!(row.names(x) %in%  row.names(x)[startsWith(row.names(x), "RPS")]), ]

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
datlcb <- stack(as.data.frame(logcpm_before)) #dataframe for ggplot
beforeplt <- ggplot(datlcb, aes(x=ind, y=values, fill=ind)) + 
  geom_boxplot(alpha=0.3, notch=TRUE) +
  scale_fill_manual(values = colorpallete)+
  ylab("Log(CPM Reads Mapped)")+
  xlab("Sample")+
  ggtitle("CPM Reads Mapped before Normalization")+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size=18))

datlca <- stack(as.data.frame(logcpm_after)) #dataframe for ggplot
afterplt <- ggplot(datlca, aes(x=ind, y=values, fill=ind)) + 
  geom_boxplot(alpha=0.3, notch=TRUE) +
  scale_fill_manual(values = colorpallete)+
  ylab("Log(CPM Reads Mapped)")+
  xlab("Sample")+
  ggtitle("CPM Reads Mapped after Normalization")+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 18))

ggarrange(beforeplt, afterplt,
          ncol = 1, nrow = 2)
#look at the normalization factors directly
y$samples['norm.factors']

#----Linear Model----
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

#visualize how samples cluster together
logcpm <- cpm(y, log=TRUE)
par(bg = "#F9FCF3")
mds <-plotMDS(logcpm)
mdsgroup <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = group)
mdsgroupplot <- ggplot(mdsgroup, aes(Dim1, Dim2, colour = group))+ 
  geom_point()+ 
  geom_label_repel(aes(label = Group),
                   box.padding   = 0.35, 
                   point.padding = 0.5) +
  scale_color_manual(values = c(colorpallete[2], colorpallete[7], colorpallete[8])) +
  ylab("Leading logFC Dimension 2")+
  xlab("Leading logFC Dimension 1")+
  ggtitle("Tissue")+
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 18))
#label the samples by rRNA
mdsrrna <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = y$samples$rRNA)
mdsrrnaplot <- ggplot(mdsrrna, aes(Dim1, Dim2, colour = group))+ 
  geom_point()+ 
  geom_label_repel(aes(label = Group),
                   box.padding   = 0.35, 
                   point.padding = 0.5) +
  scale_color_manual(values = c(colorpallete[2], colorpallete[7], colorpallete[8])) +
  ylab("Leading logFC Dimension 2")+
  xlab("Leading logFC Dimension 1")+
  ggtitle("% rRNA")+
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 18))
#label by percent of mitochondrial reads
mdschrm <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = y$samples$chrm)
mdschrmplot <- ggplot(mdschrm, aes(Dim1, Dim2, colour = group))+ 
  geom_point()+ 
  geom_label_repel(aes(label = Group),
                   box.padding   = 0.35, 
                   point.padding = 0.5) +
  scale_color_manual(values = c(colorpallete[2], colorpallete[7], colorpallete[8])) +
  ylab("Leading logFC Dimension 2")+
  xlab("Leading logFC Dimension 1")+
  ggtitle("% Mitochondrial Genes")+
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 18))
#label by slice
mdsslice <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = y$samples$slice)
mdssliceplot <- ggplot(mdsslice, aes(Dim1, Dim2, colour = group))+ 
  geom_point()+ 
  geom_label_repel(aes(label = Group),
                   box.padding   = 0.35, 
                   point.padding = 0.5) +
  scale_color_manual(values = c(colorpallete[2], colorpallete[7], colorpallete[8])) +
  ylab("Leading logFC Dimension 2")+
  xlab("Leading logFC Dimension 1")+
  ggtitle("Slice")+
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 18))

ggarrange(mdsgroupplot, mdsrrnaplot, mdschrmplot, mdssliceplot,
          ncol = 4, nrow = 1)

#plot biological coefficient of variation
y <- estimateDisp(y, design)
plotBCV(y, col.common = colorpallete[9], col.trend = colorpallete[7], col.tagwise = colorpallete[6])

#fit linear model
fit <- glmQLFit(y, design)
#brain (top) vs blood (bottom)
qlfBBL <- glmQLFTest(fit, contrast=c(-1,1,0))
#liver (top) vs blood (bottom)
qlfLBL <- glmQLFTest(fit, contrast=c(-1,0,1))
#liver (top) vs brain (bottom)
qlfLB <- glmQLFTest(fit, contrast=c(0,-1,1))

#genes that are differentially expressed in one compared to both
#blood vs both
qlfBL <- glmQLFTest(fit, contrast=c(1,-0.5,-0.5))
#liver vs both 
qlfL <- glmQLFTest(fit, contrast=c(-0.5,-0.5,1))
#brain vs both
qlfB <- glmQLFTest(fit, contrast = c(-0.5,1,-0.5))

qlfBBL

#trim lists to exclude logcpm<0 (cpm<1)
qlfBBL <- qlfBBL[which(qlfBBL$table$logCPM>=0),]
qlfLBL <- qlfLBL[which(qlfLBL$table$logCPM>=0),]
qlfLB <- qlfLB[which(qlfLB$table$logCPM>=0),]
qlfBL <- qlfBL[which(qlfBL$table$logCPM>=0),]
qlfL <- qlfL[which(qlfL$table$logCPM>=0),]
qlfB <- qlfB[which(qlfB$table$logCPM>=0),]

#DE genes, sorted by p-value, with the FDR added
topTags(qlfBBL, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfLBL, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfLB, n=10,adjust.method = "BH", sort.by = "PValue")
#1 v 2
topTags(qlfB, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfBL, n=10,adjust.method = "BH", sort.by = "PValue")
topTags(qlfL, n=10,adjust.method = "BH", sort.by = "PValue")

#get the whole table and write to a file
resultsBBL <- topTags(qlfBBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsLBL <- topTags(qlfLBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsLB <- topTags(qlfLB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsBBL, "resultsBBL.txt")
write.table(resultsLBL, "resultsLBL.txt")
write.table(resultsLB, "resultsLB.txt")
#1 v 2
resultsB <- topTags(qlfB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsBL <- topTags(qlfBL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
resultsL <- topTags(qlfL, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsB, "resultsB.txt")
write.table(resultsBL, "resultsBL.txt")
write.table(resultsL, "resultsL.txt")
write.xlsx(resultsB$table, "resultsB.xlsx", rowNames=T, colNames=T)
write.xlsx(resultsBL$table, "resultsBL.xlsx", rowNames=T, colNames=T)
write.xlsx(resultsL$table, "resultsL.xlsx", rowNames=T, colNames=T)

#summary of DE genes
summary(decideTests(qlfBBL, p.value=0.01, adjust.method = "BH", lfc=1)) #can adjust p value and log fold change
summary(decideTests(qlfLBL, p.value=0.01, adjust.method = "BH", lfc=1))
summary(decideTests(qlfLB, p.value=0.01, adjust.method = "BH", lfc=1))
#summary 1 v 2
summary(decideTests(qlfB, p.value=0.01, adjust.method = "BH", lfc=1)) #can adjust p value and log fold change
summary(decideTests(qlfBL, p.value=0.01, adjust.method = "BH", lfc=1))
summary(decideTests(qlfL, p.value=0.01, adjust.method = "BH", lfc=1))


#check the expression across all samples, not just my assigned ones
which(rowData(brain)$gene_name == "UGT1A1") #gene found in liver cells
# Make individual data frames
a <- data.frame(group = "Brain", value = assays(brain)$TPM[32683,])
b <- data.frame(group = "Blood", value = assays(blood)$TPM[32683,])
c <- data.frame(group = "Liver", value = assays(liver)$TPM[32683,])

# Combine into one long data frame
plot.data <- rbind(a, b, c)

ggplot(plot.data, aes(x=group, y=value, fill=group)) + 
  geom_boxplot(alpha=0.3, notch=TRUE) +
  scale_fill_manual(values = c(colorpallete[1], colorpallete[2], colorpallete[8]))+
  ylab("Transcript per Million (TPM)")+
  xlab("Tissue")+
  ggtitle("TPM Levels for UGT1A1")+
  theme(legend.position="none",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size=18))

#use a statistical test to prove this expression across all samples is significant
wilcox.test(assays(liver)$TPM[32683,], append(assays(brain)$TPM[32683,],assays(blood)$TPM[32683,]), alternative='two.sided',exact = F)
