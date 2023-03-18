###################
#Set work directory
###################
setwd("~/18071_data_analysis/GDCdata")

###############
#Load libraries
###############
suppressWarnings(suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr); 
  library("DESeq2"); library(biomaRt); library(survival);
  library(stringr); library(edgeR); 
  library(survminer); library(data.table);
  library(tidyverse)}))

####################################################
#Download gene expression data of melanoma from TCGA
####################################################
#download link
melanoma.exp_mat.link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.htseq_counts.tsv.gz"

#retrieve gene expression data
download.file(melanoma.exp_mat.link, destfile= paste(getwd(), "mel_dataset.tsv.gz", sep ="/"))
melanoma.exp_mat <- fread(paste0(getwd(), "/mel_dataset.tsv.gz"))

#View(melanoma.exp_mat)

#extract Ensembl_ID
Ensembl_ID <- as.data.frame(melanoma.exp_mat$Ensembl_ID)
colnames(Ensembl_ID) <- "Ensembl_ID"

##remove Ensembl_ID from columns
melanoma.exp_mat <- melanoma.exp_mat[,-1]

#transform gene expression data from pseudo-count to raw count matrix
pseudocnt_count <- function(count){
  raw_count <- as.integer(2^(count)-1)
  return(raw_count)
}

#apply function on dataframe to get raw count
melanoma.exp_mat <- apply(melanoma.exp_mat, 2, FUN = pseudocnt_count)

#View(melanoma.exp_mat)

#cbind Ensembl_ID with melanoma.exp.mat
melanoma.exp_mat  <- cbind(Ensembl_ID, melanoma.exp_mat) %>% as.data.frame() %>% column_to_rownames(var = "Ensembl_ID")

#View(melanoma.exp_mat)

#mapping gene symbol to Ensembl_ID and feature reduction
#remove version ID from rownames in count matrix
rownames(melanoma.exp_mat) <- gsub("\\..*","",rownames(melanoma.exp_mat))

#create object for annotating features (i.e Ensembl_ID)
mart <- biomaRt::useMart(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl")
transcript.info <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id","entrezgene_id",
                                                 "transcript_biotype",  "transcript_length"), mart = mart)

#save transcript information
saveRDS(transcript.info, paste0(getwd(), "/transcript.info.rds"))

#select protein coding genes only from biomaRt annotation data
transcript.info <- filter(transcript.info, transcript_biotype %in% c("protein_coding"))

#remove duplicates in transcript.info based on "hgnc_symbol", 
transcript.info <- transcript.info %>% distinct(hgnc_symbol, .keep_all = TRUE) #in this case I kept only the first occurrence of duplicates

#annotate "ensembl_ID" with "hgnc_symbol"
annotation_table <- transcript.info[,1:2]

#create a column containing ensembl_ID in brca raw_count 
melanoma.exp_mat <- rownames_to_column(melanoma.exp_mat)

#rename new column as "ensembl_gene_id"
colnames(melanoma.exp_mat)[1] <- "ensembl_gene_id"

#map melanoma.exp_mat to created annotation_table by "ensembl_gene_id"
melanoma.exp_mat <- left_join(melanoma.exp_mat, annotation_table, 
                              by='ensembl_gene_id') #note: this create a column called "hgnc_symbol" at the end of the raw_count df

#rm rows containing NA
melanoma.exp_mat <- na.omit(melanoma.exp_mat)

#make "hgnc_symbol" the new rownames and delete "ensembl_gene_id" column
rownames(melanoma.exp_mat) <- melanoma.exp_mat$hgnc_symbol

#delete the column containing "ensembl_gene_id" & "hgnc_symbol"
melanoma.exp_mat <- melanoma.exp_mat[, !names(melanoma.exp_mat) %in% c("ensembl_gene_id", "hgnc_symbol")]

#download link for phenotype data
phenoData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.GDC_phenotype.tsv.gz"

#retrieve phenodata
download.file(phenoData_link, destfile= paste(getwd(), "mel_phenoData.tsv.gz", sep ="/"))
mel_phenodata <- fread(paste0(getwd(), "/mel_phenoData.tsv.gz"))

#subset only samples found in melanoma.exp_mat
mel_phenodata <- mel_phenodata %>% filter(submitter_id.samples %in% colnames(melanoma.exp_mat))

#subset "Metastatic" $ "primary Tumor"
mel_phenodata <- mel_phenodata %>% 
  filter(sample_type.samples %in% c("Metastatic", "Primary Tumor")) %>% as.data.frame() %>% 
  column_to_rownames(var = "submitter_id.samples")

#subset only samples that are "Metastatic" & "Primary Tumor" in melanoma.exp_mat
melanoma.exp_mat <- melanoma.exp_mat[, colnames(melanoma.exp_mat) %in% rownames(mel_phenodata)] %>% as.data.frame()

#check if desired sample type is appropriately retrieved
#table(mel_phenodata$sample_type.samples)

#View phenodata 
#head(mel_phenodata)

#retrieve survival data 
survData_link <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.survival.tsv"

#retrieve SurvData
download.file(survData_link, destfile= paste(getwd(), "mel_survData.tsv.gz", sep ="/"))
mel_survdata <- fread(paste0(getwd(), "/mel_survData.tsv.gz"))

#convert overall survival time from days to month
mel_survdata <- mel_survdata %>% mutate(OS.time.year = OS.time / 365) %>%
  as.data.frame() %>% column_to_rownames(var = "sample")

#subset samples in "melanoma.exp_mat" and "mel_phenodata" that is contained in "mel_survdata" (i.e samples with survival details)
melanoma.exp_mat <- melanoma.exp_mat[, colnames(melanoma.exp_mat) %in% rownames(mel_survdata)] %>% as.data.frame()
mel_phenodata <- mel_phenodata %>% filter(rownames(mel_phenodata) %in% rownames(mel_survdata)) %>% as.data.frame()

#subset samples in "mel_survdata" that are contained in "melanoma.exp_mat"
mel_survdata <- mel_survdata %>% filter(rownames(mel_survdata) %in% colnames(melanoma.exp_mat)) %>% as.data.frame()

#check dimension of "melanoma.exp_mat", "mel_phenodata" & "mel_survdata"; note that n_samples must be the same in all three dataframe
dim(melanoma.exp_mat)
dim(mel_phenodata)
dim(mel_survdata)

#add a reference column for joining mel_phenodata and mel_survdata
mel_phenodata[["sample"]] <- rownames(mel_phenodata)
mel_survdata[["sample"]] <- rownames(mel_survdata)

#create a new dataframe containing phenotype & survival information
new_mel_phenodata <- left_join(mel_phenodata, mel_survdata, 
                               by= "sample") %>% column_to_rownames(var = "sample") 

#Reordering rownames of mel_phenodata to match melanoma.exp_mat; NB: this should be done to prevent error creating a summarized experiment
reorder_ndx <- match(colnames(melanoma.exp_mat), rownames(new_mel_phenodata))
new_mel_phenodata <- new_mel_phenodata[reorder_ndx, ]

#View new phenodata
#head(new_mel_phenodata)

#retrieve transcript_length from biomart
transcript <- data.frame(transcript.info[, c("hgnc_symbol", "transcript_length")])

#make a dataframe of genes present in melanoma.exp_mat
hgnc_symbol <- as.data.frame(rownames(melanoma.exp_mat), row.names = rownames(melanoma.exp_mat))
colnames(hgnc_symbol) <- "hgnc_symbol"

#join hgnc_symbol dataframe $ transcript dataframe to get 
transcript <- left_join(hgnc_symbol, transcript, 
                        by='hgnc_symbol') #%>% column_to_rownames(var = "hgnc_symbol")

rownames(transcript) <- transcript$hgnc_symbol

#View transcript
#head(transcript)

###########################################
#Normalization by TPM using custom function
###########################################
#normalize count matrix by TPM
#define TPM function
tpm_transform <- function(count, transcript_length){
  rate <- count/transcript_length
  rate/sum(rate) * 1e6
}

#perform TPM normalization on melanoma.exp_mat using count and transcript_length
melanoma.exp_mat_tpm <- as.data.frame(apply(melanoma.exp_mat, 2, function(x) tpm_transform(x, transcript$transcript_length)))
#View(melanoma.exp_mat_tpm)

#quick check to confirm data normalization
# colMeans(melanoma.exp_mat_tpm)
# colSums(melanoma.exp_mat_tpm)

##############################
#Normalization by TMM in edgeR
##############################
#Organize data in a summarized experiment (se) for downstream
mel_se <- SummarizedExperiment(assays = melanoma.exp_mat,
                               rowData = transcript,
                               colData = new_mel_phenodata)

#save summarized experiment as RDS for posterity
saveRDS(mel_se, "~/18071_data_analysis/GDCdata/mel_se.rds")

#read summarized experiment
mel_se <- readRDS("mel_se.rds")
#assay(mel_se)
#colData(mel_se)
#rowData(mel_se)

#normalization by trimmed mean of M (TMM)
#create factor for all samples
group <- factor(mel_se$sample_type.samples)

#create dgeList
dge <- DGEList(counts=assay(mel_se), group = group,
               samples=colData(mel_se),
               genes=as.data.frame(rowData(mel_se)))

#re-filter data using edgeR
keep <- rowSums(cpm(dge)>100) >= 2 #we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
dge <- dge[keep, , keep.lib.sizes=FALSE]

#free memory
rm(keep)

#reset library size
dge$samples$lib.size <- colSums(dge$counts)
head(dge$samples)

# Normalization (by TMM)
dge <- calcNormFactors(dge, method="TMM")

#save dge as RDS
saveRDS(object = dge,
        file = "mel_dge.RDS",
        compress = FALSE)

#read dge
dge <- readRDS("mel_dge.RDS")

#get tmm normalized count
melanomaData_tmm_normalize <- cpm(dge, log = FALSE) %>% as.matrix.default()

#View normalized count matrix
#View(melanomaData_tmm_normalize)

#########################################
#Exploratory analysis of Normalized count
#########################################
#set seed
set.seed(213455324)

#subset random samples from raw count matrix
random_samples <- sample(seq_len(ncol(melanoma.exp_mat)),42L)
sample_cnt  <- melanoma.exp_mat[,random_samples]
sample_tpm  <- melanoma.exp_mat_tpm[,random_samples]
sample_tmm  <- melanomaData_tmm_normalize[,random_samples]

#define plot labels
x1 <- "log(1+TPM)"
y1 <- "log(1+count)"
y2 <- "log(1+TPM)"
y3 <- "log(1+count)"

#define plot parameters
par(mai=c(1.2,.7,.1,.05))
par(mfrow=c(2,2))

#plot sampling
#un-normalized count mtx
boxplot(log(1+sample_cnt), las=2, ylab = y1)

#tpm-normalized count mtx
boxplot(log(1+sample_tpm), las=2, ylab = y1)

#tmm-normalized count mtx
boxplot(log(1+sample_tmm), las=2, ylab = y2)

###############################################
#subset gene of interests for survival analysis
###############################################
#subset TP53, MDM2, & BCL6
survival_genes <- melanomaData_tmm_normalize[rownames(melanomaData_tmm_normalize) %in% c("TP53", "MDM2", "BCL6"),]

#transpose survival_genes dataframe
survival_genes <- as.data.frame(t(survival_genes))

#subset overall survival (OS) and overall survival time (OS.time.year) from summarized experiment
new_survival_info <- as.data.frame(colData(mel_se)[, c("OS","OS.time.year")], row.names = rownames(colData(mel_se)))

#rename column name in new_survival_info to fit survminer package
colnames(new_survival_info) <- c(OS = "event", OS.time.year = "time")

#create the dataframe for survival analysis
survival_df <- cbind(survival_genes, new_survival_info)

#save survival_df as rds 
saveRDS(survival_df, "~/18071_data_analysis/GDCdata/survival_df.rds")

#read survival_df
survival_df <- readRDS("survival_df.rds")

#View survival genes
head(survival_df)

#survival analysis
#determine the optimal cut-point of genes
surv.cut <- surv_cutpoint(survival_df, time = "time", event = "event",
                          variables = c("TP53", "MDM2", "BCL6"))

#view cut-point summary
summary(surv.cut)

#categorize expression level of genes based on the optimal cut-point
surv.cat <- surv_categorize(surv.cut)

#view categories
surv.cat

#fit survival curves for "TP53" and visualize
surv.fit <- survfit(Surv(time, event) ~TP53, data = surv.cat)
ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T)

#fit survival curves for "MDM2" and visualize
surv.fit <- survfit(Surv(time, event) ~MDM2, data = surv.cat)
ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T)

#fit survival curves for "BCL6" and visualize
surv.fit <- survfit(Surv(time, event) ~BCL6, data = surv.cat)
ggsurvplot(surv.fit, data = surv.cat, risk.table = F, conf.int = T, surv.median.line = "hv",
           pval = T)
