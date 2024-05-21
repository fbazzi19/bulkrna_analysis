#library calls
library(recount3)

#brain data
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
