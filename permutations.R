##compares predictiveness of isoforms for a known biomarker (gene) among ranked stability indicies (are the most "stable" isoforms, the most predictive?)

library(PharmacoGx)
library(survcomp)
library(gimme)
library(dplyr)
library(forestplot)

setwd({"~/Desktop/permutation"})

#load biological replicates
load("inter.RData")

#read in PSets
gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")
CCLE <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CCLE.rds")
GDSC <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC2-8.2.rds")

#modify annotations so that the molecular profiles can be read in the permutation function
S4Vectors::metadata(gCSI@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"
S4Vectors::metadata(CCLE@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"
S4Vectors::metadata(GDSC@molecularProfiles[["Kallisto_0.46.1.isoforms.counts"]])$annotation <- "isoform"

#function that creates data.frame of every single combination between two vectors (e.g. drug and gene combinations)
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#read in known biomarkers (curated by BHKLAB), and filter for expression biomarkers only with drugs that exist across gCSI,CCLE,GDSC
#biomarkers selected are: 1) ERBB2 + Lapatinib; 2) ALK + Crizotinib; 3) PHB + Paclitaxel; 4) ESR2 + Erlotinib; 5) EGFR & Lapatinib
intersected_drugs <- Reduce(intersect, list(gCSI@sensitivity$info$drugid, CCLE@sensitivity$info$drugid, GDSC@sensitivity$info$drugid))
biomarkers <- as.data.frame(readxl::read_xlsx("biomarkers.xlsx"))
biomarkers <- biomarkers[grep("EXPR", biomarkers$Alteration.type),]
biomarkers <- biomarkers[which(biomarkers$compound %in% intersected_drugs),]
genes_keep <- c("ERBB2","ALK","PHB","ESR2","EGFR")
biomarkers <- biomarkers[which(biomarkers$gene %in% genes_keep),]
biomarkers <- biomarkers[c(1,2,4,5,6),]

################################
#Compute Concordance Index (CI)#
################################

#function for computing CI
computeCI <- function(PSet, mData, features, cells){
  
  rnaseq <- summarizeMolecularProfiles(PSet, mDataType = mData)
  rnaseq <- as.data.frame(rnaseq@assays@data$exprs)
  
  sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(PSet, sensitivity.measure = "aac_recomputed",  fill.missing = F))
  intersected_rnacells <- cells
  drugs <- unique(rownames(sensitivity_data))
  commonSamples <- intersected_rnacells
  commonSamples <- intersected_rnacells[which(intersected_rnacells %in% colnames(sensitivity_data))]
  combinations <- as.data.frame(expand.grid.unique(features, drugs, include.equals = TRUE))
  combinations$ci <- NA
  combinations$pvalue <- NA
  combinations$se <- NA
  combinations$upper <- NA
  combinations$lower <- NA
  colnames(combinations) <- c("gene","drug","ci","pvalue","se","upper","lower")
  
  for (i in 1:nrow(combinations)){
    print(paste0(i, " out of ", nrow(combinations), " complete"))
    tt <- sensitivity_data[combinations[,2][i], commonSamples]
    ci <- survcomp::concordance.index(as.numeric(tt), surv.time = as.numeric(unlist(-rnaseq[combinations[,1][i], commonSamples])), 
                                      surv.event = rep(1,length(sensitivity_data[commonSamples])), 
                                      outx = TRUE, method="noether", na.rm = TRUE)
    
    combinations$pvalue[i] <- ci$p.value
    combinations$ci[i] <- ci$c.index
    combinations$se[i] <- ci$se
    combinations$upper[i] <- ci$upper
    combinations$lower[i] <- ci$lower
  }
  
  return(combinations)
  
}

#function computes CI for a data.frame of biomarkers and organizes results in a data.frame
computeCIBiomarker <- function(biomarker_matrix, mData, cells) {
  
  biomarkers_CI <- as.data.frame(matrix(ncol = 17, nrow = nrow(biomarker_matrix)))
  colnames(biomarkers_CI) <- c("gene","drug", "gCSI_CI", "CCLE_CI","GDSC_CI", "gCSI_pvalue", "CCLE_pvalue", "GDSC_pvalue", "gCSI_se", "CCLE_se","GDSC_se", "gCSI_upper",
                               "CCLE_upper","GDSC_upper","gCSI_lower","CCLE_lower","GDSC_lower")
  for (i in 1:nrow(biomarker_matrix)){
    
    drug <- biomarker_matrix$compound[i]
    gene <- biomarker_matrix$gene[i]
    
    if (!mData == "Kallisto_0.46.1.isoforms.counts"){
      gene_select <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name == gene)]
    }else{
      gene_select = gene
    }
    
    if (length(cells) == 48){
      gcsi_CI <- computeCI(PSet=gCSI, mData = mData, features = gene_select, cells = cells)
      ccle_CI <- computeCI(PSet=CCLE, mData = mData, features = gene_select, cells = cells)
      gdsc_CI <- computeCI(PSet=GDSC, mData = mData, features = gene_select, cells = cells)
    }else{
      gcsi_CI <- computeCI(PSet=gCSI, mData = mData, features = gene_select, cells = gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq$cellid)
      ccle_CI <- computeCI(PSet=CCLE, mData = mData, features = gene_select, cells = CCLE@molecularProfiles$Kallisto_0.46.1.rnaseq$cellid)
      gdsc_CI <- computeCI(PSet=GDSC, mData = mData, features = gene_select, cells = GDSC@molecularProfiles$Kallisto_0.46.1.rnaseq$cellid)
      
    }
    
    biomarkers_CI$drug[i] <- drug
    biomarkers_CI$gene[i] <- gene
    
    biomarkers_CI$gCSI_CI[i] <- gcsi_CI$ci[which(gcsi_CI$drug == drug)]
    biomarkers_CI$gCSI_pvalue[i] <- gcsi_CI$pvalue[which(gcsi_CI$drug == drug)]
    biomarkers_CI$gCSI_se[i] <- gcsi_CI$se[which(gcsi_CI$drug == drug)]
    biomarkers_CI$gCSI_upper[i] <- gcsi_CI$upper[which(gcsi_CI$drug == drug)]
    biomarkers_CI$gCSI_lower[i] <- gcsi_CI$lower[which(gcsi_CI$drug == drug)]
    
    biomarkers_CI$CCLE_CI[i] <- ccle_CI$ci[which(ccle_CI$drug == drug)]
    biomarkers_CI$CCLE_pvalue[i] <- ccle_CI$pvalue[which(ccle_CI$drug == drug)]
    biomarkers_CI$CCLE_se[i] <- ccle_CI$se[which(ccle_CI$drug == drug)]
    biomarkers_CI$CCLE_upper[i] <- ccle_CI$upper[which(ccle_CI$drug == drug)]
    biomarkers_CI$CCLE_lower[i] <- ccle_CI$lower[which(ccle_CI$drug == drug)]
    
    biomarkers_CI$GDSC_CI[i] <- gdsc_CI$ci[which(gdsc_CI$drug == drug)]
    biomarkers_CI$GDSC_pvalue[i] <- gdsc_CI$pvalue[which(gdsc_CI$drug == drug)]
    biomarkers_CI$GDSC_se[i] <- gdsc_CI$se[which(gdsc_CI$drug == drug)]
    biomarkers_CI$GDSC_upper[i] <- gdsc_CI$upper[which(gdsc_CI$drug == drug)]
    biomarkers_CI$GDSC_lower[i] <- gdsc_CI$lower[which(gdsc_CI$drug == drug)]
    
  }
  
  return(biomarkers_CI)
  
}

############
#Gene-Level#
############

#gene_CI <- computeCIBiomarker(biomarker_matrix = biomarkers, mData = "Kallisto_0.46.1.rnaseq", cells = intersected_rnacells)
#save(gene_CI, file="gene_CI.RData")

load("gene_CI.RData")



###############
#Isoform-Level#
###############

##import transcript stability data generated by transcript_stability.R
transcript_stability <- read.csv("transcript_stability.csv")

##get isoforms ID's for each gene biomarker (only isoforms present in generated by transcript_stability.R)
erbb2_isoforms <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == "ERBB2")]
erbb2_isoforms <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% erbb2_isoforms)]
erbb2_isoforms <- as.data.frame(expand.grid.unique(erbb2_isoforms, "Lapatinib", include.equals = TRUE))
colnames(erbb2_isoforms) <- c("gene","compound")

phb_isoforms <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == "PHB")]
phb_isoforms <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% phb_isoforms)]
phb_isoforms  <- as.data.frame(expand.grid.unique(phb_isoforms, "Paclitaxel", include.equals = TRUE))
colnames(phb_isoforms) <- c("gene","compound")

alk_isoforms <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == "ALK")]
alk_isoforms <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% alk_isoforms)]
alk_isoforms  <- as.data.frame(expand.grid.unique(alk_isoforms, "Crizotinib", include.equals = TRUE))
colnames(alk_isoforms) <- c("gene","compound")

egfr_isoforms <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == "EGFR")]
egfr_isoforms <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% egfr_isoforms)]
egfr_isoforms  <- as.data.frame(expand.grid.unique(egfr_isoforms, "Lapatinib", include.equals = TRUE))
colnames(egfr_isoforms) <- c("gene","compound")

esr2_isoforms <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$transcript_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata$gene_name == "ESR2")]
esr2_isoforms <- transcript_stability$transcript_id[which(transcript_stability$transcript_id %in% esr2_isoforms)]
esr2_isoforms  <- as.data.frame(expand.grid.unique(esr2_isoforms, "Erlotinib", include.equals = TRUE))
colnames(esr2_isoforms) <- c("gene","compound")

#erbb2_CI <- computeCIBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#phb_CI <- computeCIBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#alk_CI <- computeCIBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#egfr_CI <- computeCIBiomarker(biomarker_matrix = egfr_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#esr2_CI <- computeCIBiomarker(biomarker_matrix = esr2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)

#erbb2_CI$gene_name <- "ERBB2"
#phb_CI$gene_name <- "PHB"
#alk_CI$gene_name <- "ALK"
#egfr_CI$gene_name <- "EGFR"
#esr2_CI$gene_name <- "ESR2"

#isoform_CI <- rbind(erbb2_CI, phb_CI, alk_CI, egfr_CI, esr2_CI)

#save(isoform_CI, file="isoform_CI.RData")

load("isoform_CI.RData")


######################
#Compute Permutations#
######################

#function to call drugSensitivitySig for a given PSet, cells, and drugs
PermutationSig <- function(PSet, mData, features, cells, drugs){
  
  signature <- drugSensitivitySig(PSet, mData, drugs=drugs, 
                                  features=features, sensitivity.measure = "aac_recomputed", modeling.method="pearson", 
                                  inference.method="resampling", cells=cells, nthread=2, parallel.on = "gene")
  return(signature)
}


#function to compute permutations for a biomarker data.frame
computePerBiomarker <- function(biomarker_matrix, mData, cells) {
  
  biomarkers_perm <- as.data.frame(matrix(ncol = 17, nrow = nrow(biomarker_matrix)))
  colnames(biomarkers_perm) <- c("gene","drug", "gCSI_pearson", "CCLE_pearson","GDSC_pearson", "gCSI_pvalue", "CCLE_pvalue", "GDSC_pvalue", "gCSI_sig", "CCLE_sig","GDSC_sig", "gCSI_upper", "gCSI_lower","CCLE_upper","CCLE_lower","GDSC_upper","GDSC_lower")
  
  for (i in 1:nrow(biomarker_matrix)){
    drug <- biomarker_matrix$compound[i]
    gene <- biomarker_matrix$gene[i]
    
    if (!mData == "Kallisto_0.46.1.isoforms.counts"){
      gene_select <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[match(gene, gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name)]
    }else{
      gene_select = gene
    }
    
    if (length(cells) == 48){
      gcsi_per <- PermutationSig(PSet = gCSI, mData = mData, features = gene_select, 
                                 cells = cells, drugs = drug)
      
      ccle_per <- PermutationSig(PSet = CCLE, mData = mData, features = gene_select, 
                                 cells = cells, drugs = drug)
      
      gdsc_per <- PermutationSig(PSet = GDSC, mData = mData, features = gene_select, 
                                 cells = cells, drugs = drug)
    }else{
      gcsi_per <- PermutationSig(PSet = gCSI, mData = mData, features = gene_select, 
                                 cells = gCSI@molecularProfiles$Kallisto_0.46.1.isoforms$cellid, drugs = drug)
      
      ccle_per <- PermutationSig(PSet = CCLE, mData = mData, features = gene_select, 
                                 cells = CCLE@molecularProfiles$Kallisto_0.46.1.isoforms$cellid, drugs = drug)
      
      gdsc_per <- PermutationSig(PSet = GDSC, mData = mData, features = gene_select, 
                                 cells = GDSC@molecularProfiles$Kallisto_0.46.1.isoforms$cellid, drugs = drug)
      
    }
    biomarkers_perm$drug[i] <- drug
    biomarkers_perm$gene[i] <- gene
    
    biomarkers_perm$gCSI_pearson[i] <- gcsi_per[,,"estimate"]
    biomarkers_perm$gCSI_pvalue[i] <- gcsi_per[,,"pvalue"]
    biomarkers_perm$gCSI_sig[i] <- gcsi_per[,,"significant"]
    biomarkers_perm$gCSI_upper[i] <- gcsi_per[,,"upper"]
    biomarkers_perm$gCSI_lower[i] <- gcsi_per[,,"lower"]
    
    biomarkers_perm$CCLE_pearson[i] <- ccle_per[,,"estimate"]
    biomarkers_perm$CCLE_pvalue[i] <- ccle_per[,,"pvalue"]
    biomarkers_perm$CCLE_sig[i] <- ccle_per[,,"significant"]
    biomarkers_perm$CCLE_upper[i] <- ccle_per[,,"upper"]
    biomarkers_perm$CCLE_lower[i] <-ccle_per[,,"lower"]
    
    biomarkers_perm$GDSC_pearson[i] <- gdsc_per[,,"estimate"]
    biomarkers_perm$GDSC_pvalue[i] <- gdsc_per[,,"pvalue"]
    biomarkers_perm$GDSC_sig[i] <- gdsc_per[,,"significant"]
    biomarkers_perm$GDSC_upper[i] <- gdsc_per[,,"upper"]
    biomarkers_perm$GDSC_lower[i] <- gdsc_per[,,"lower"]
    
    
  }
  return(biomarkers_perm)
  
}

############
#Gene-Level#
############

#gene_permut <- computePerBiomarker(biomarker_matrix = biomarkers, mData = "Kallisto_0.46.1.rnaseq", cells = intersected_rnacells)
#save(gene_permut, file="gene_permut.RData")

load("gene_permut.RData")

#convert significance from 1/0 to YES/NO
gene_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] <- ifelse(gene_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] == 0 | is.na(gene_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")]), "NO", "YES")

#add standard error (using 95% confidence interval - 3.92) (upper-lower/3.92)
gene_permut$gCSI_se <- (gene_permut$gCSI_upper - (gene_permut$gCSI_lower))/3.92
gene_permut$CCLE_se <- (gene_permut$CCLE_upper - (gene_permut$CCLE_lower))/3.92
gene_permut$GDSC_se <- (gene_permut$GDSC_upper - (gene_permut$GDSC_lower))/3.92

###############
#Isoform-Level#
###############

#erbb2_permut <- computePerBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#phb_permut <- computePerBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#alk_permut <- computePerBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#egfr_permut <- computePerBiomarker(biomarker_matrix = egfr_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#esr2_permut <- computePerBiomarker(biomarker_matrix = esr2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)

#erbb2_permut$gene_name <- "ERBB2"
#phb_permut$gene_name <- "PHB"
#alk_permut$gene_name <- "ALK"
#egfr_permut$gene_name <- "EGFR"
#esr2_permut$gene_name <- "ESR2"

#isoform_permut <- rbind(erbb2_permut, phb_permut, alk_permut, egfr_permut, esr2_permut)
#save(isoform_permut, file="isoform_permut.RData")

load("isoform_permut.RData")

isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] <- ifelse(isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] == 0 | is.na(isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")]), "NO", "YES")

isoform_permut$gCSI_se <- (isoform_permut$gCSI_upper - (isoform_permut$gCSI_lower))/3.92
isoform_permut$CCLE_se <- (isoform_permut$CCLE_upper - (isoform_permut$CCLE_lower))/3.92
isoform_permut$GDSC_se <- (isoform_permut$GDSC_upper - (isoform_permut$GDSC_lower))/3.92


#######################
######Forest_Plot######
#######################

#add stabilities (gcsi/gdsc, gcsi/ccle, gdsc/ccle) to isoform permutations
isoform_permut$gcsi_ccle_stab <- transcript_stability$gcsi_ccle_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$gdsc_ccle_stab <- transcript_stability$gdsc_ccle_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$gcsi_gdsc_stab <- transcript_stability$gcsi_gdsc_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$mean_stability <- rowMeans(isoform_permut[,c(-1:-21)])


isoform_CI$gcsi_ccle_stab <- transcript_stability$gcsi_ccle_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$gdsc_ccle_stab <- transcript_stability$gdsc_ccle_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$gcsi_gdsc_stab <- transcript_stability$gcsi_gdsc_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$mean_stability <- rowMeans(isoform_CI[,c(-1:-18)])


#genes-drug associations to be included in forest plot
genes <- c("PHB", "ALK","ERBB2","EGFR","ESR2")
drugs <- c("Paclitaxel", "Crizotinib", "Lapatinib","Lapatinib","Erlotinib")

for (i in 1:length(genes)){
  gene <- genes[i]
  drug <- drugs[i]
  
  gene.CI <- gene_CI[which(gene_CI$gene == gene & gene_CI$drug == drug),]
  gene_perm <- gene_permut[which(gene_permut$gene == gene & gene_permut$drug == drug),]
  
  transcript_CI <- isoform_CI[which(isoform_CI$gene_name == gene & isoform_CI$drug == drug),]  
  transcript_perm <- isoform_permut[which(isoform_permut$gene_name == gene & isoform_permut$drug == drug),]
  
  #order transcripts based on stability (mean stability)
  transcript_CI <- transcript_CI[order(transcript_CI$mean_stability, decreasing = TRUE),] 
  transcript_perm <- transcript_perm[order(transcript_perm$mean_stability, decreasing = TRUE),] 
  
  
  plot_combined_meta <- data.frame(matrix(nrow = nrow(transcript_perm), ncol=8))
  colnames(plot_combined_meta) <- c("ci_estimate", "ci_se", "ci_upper","ci_lower", "perm_estimate","perm_se","perm_upper","perm_lower")
  rownames(plot_combined_meta) <- transcript_perm$gene
  
  #combined gene CI
  combined_ci <- combine.est(
    c(
      gene.CI$gCSI_CI, gene.CI$CCLE_CI, gene.CI$GDSC_CI
    ),
    c(
      gene.CI$gCSI_se, gene.CI$CCLE_se, gene.CI$GDSC_se
    ),na.rm = TRUE,hetero = TRUE)
  
  #combined gene pearson correlation from permutations
  combined_pear <- combine.est(
    c(
      gene_perm$gCSI_pearson, gene_perm$CCLE_pearson, gene_perm$GDSC_pearson
    ),
    c(
      gene_perm$gCSI_se, gene_perm$CCLE_se, gene_perm$GDSC_se
    ),na.rm = TRUE,hetero = TRUE)
  
  #loop for combined isoform pearson correlation from permutations
  for (i in 1:nrow(transcript_perm)){
    plot_combined_meta$ci_estimate[i] <- combine.est(
      c(
        transcript_CI$gCSI_CI[i], transcript_CI$CCLE_CI[i], transcript_CI$GDSC_CI[i]
      ),
      c(
        transcript_CI$gCSI_se[i], transcript_CI$CCLE_se[i], transcript_CI$GDSC_se[i]
      ),na.rm = TRUE,hetero = TRUE)$estimate
    
    plot_combined_meta$ci_se[i] <- combine.est(
      c(
        transcript_CI$gCSI_CI[i], transcript_CI$CCLE_CI[i], transcript_CI$GDSC_CI[i]
      ),
      c(
        transcript_CI$gCSI_se[i], transcript_CI$CCLE_se[i], transcript_CI$GDSC_se[i]
      ),na.rm = TRUE,hetero = TRUE)$se
    
    
    plot_combined_meta$perm_estimate[i] <- combine.est(
      c(
        transcript_perm$gCSI_pearson[i], transcript_perm$CCLE_pearson[i], transcript_perm$GDSC_pearson[i]
      ),
      c(
        transcript_perm$gCSI_se[i], transcript_perm$CCLE_se[i], transcript_perm$GDSC_se[i]
      ),na.rm = TRUE,hetero = TRUE)$estimate
    
    
    plot_combined_meta$perm_se[i] <- combine.est(
      c(
        transcript_perm$gCSI_pearson[i], transcript_perm$CCLE_pearson[i], transcript_perm$GDSC_pearson[i]
      ),
      c(
        transcript_perm$gCSI_se[i], transcript_perm$CCLE_se[i], transcript_perm$GDSC_se[i]
      ),na.rm = TRUE,hetero = TRUE)$se
    
    plot_combined_meta$ci_lower[i] <- plot_combined_meta$ci_estimate[i] + qnorm(0.025, lower.tail=TRUE) *  plot_combined_meta$ci_se[i]
    plot_combined_meta$ci_upper[i] <- plot_combined_meta$ci_estimate[i] + qnorm(0.025, lower.tail=FALSE) *  plot_combined_meta$ci_se[i]
    
    plot_combined_meta$perm_lower[i] <- plot_combined_meta$perm_estimate[i] + qnorm(0.025, lower.tail=TRUE) *  plot_combined_meta$perm_se[i]
    plot_combined_meta$perm_upper[i] <- plot_combined_meta$perm_estimate[i] + qnorm(0.025, lower.tail=FALSE) *  plot_combined_meta$perm_se[i]
    
  }
  
  
  combined_ci_lower <- combined_ci$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_ci$se
  combined_ci_upper <- combined_ci$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_ci$se
  
  combined_pear_lower <- combined_pear$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_pear$se
  combined_pear_upper <- combined_pear$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_pear$se
  
  #combined_ci_p <- pnorm((combined_ci$estimate - 0.5)/combined_ci$se, lower.tail = combined_ci$estimate < 0.5) * 2
  
  rows_l <- length(c(NA, combined_ci$estimate,
                     plot_combined_meta$ci_estimate))
  
  c_indices <- structure(
    list(
      mean  = c(NA, combined_ci$estimate,
                plot_combined_meta$ci_estimate),
      
      lower = c(NA, combined_ci_lower,
                plot_combined_meta$ci_lower),
      
      upper = c(NA, combined_ci_upper,
                plot_combined_meta$ci_upper)
    ),
    .Names = c("C-index    ", "lower", "upper"),
    row.names = c(NA, -(rows_l)), 
    class = "data.frame"
  )
  
  
  c_tabletext <- cbind(
    c("Gene/Isoform", paste0(gene, " (Gene)"), c(transcript_CI$gene)),
    c("N",rep(144,rows_l-1)),
    
    c("C-index", formatC(combined_ci$estimate, format = "e", digits = 2),
      formatC(plot_combined_meta$ci_estimate, format = "e", digits = 2)),

    c("Pearson \n(permutation)", formatC(combined_pear$estimate, format = "e", digits = 2),
      formatC(plot_combined_meta$perm_estimate, format = "e", digits = 2)),
    
    
    c("Mean stability", NA, formatC(transcript_perm$mean_stability, format = "e", digits = 2))
  )
  
  
  fn <- local({
    i = 0
    
    b_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))))
    l_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))))
    function(..., clr.line, clr.marker){
      i <<- i + 1
      fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
      #fpDrawSummaryCI(...,col=s_clrs[i])
    }
  })
  
  fn1 <- local({
    i = 0
    
    s_clrs = c(c("darkgreen",rep("red", nrow(transcript_CI))))
    function(..., col){
      i <<- i + 1
      fpDrawSummaryCI(...,col=s_clrs[i])
    }
  })
  
  fileName = paste0("figures/",gene,"_",drug,".pdf")
  pdf(fileName, width=9 , height=8, onefile=FALSE)
  
  forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T), xlab="Concordance Index", 
             title=paste0(gene,"_",drug), zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:5, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                            ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                            xlab=gpar(fontfamily = "", cex=0.8, fontface=2), 
                            legend=gpar(fontfamily = "", cex = 0.5, fontface=1)),
             fn.ci_sum = fn1,
             col = fpColors(text="black"),
             xticks= c(.3, .4, .5, .6, .7, .8)
  )
  
  dev.off()
  
}
