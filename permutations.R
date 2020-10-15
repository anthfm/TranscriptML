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
#biomarkers selected are: 1) ERBB2 + Lapatinib; 2) ALK + Crizotinib; 3) PHB + Paclitaxel
intersected_drugs <- Reduce(intersect, list(gCSI@sensitivity$info$drugid, CCLE@sensitivity$info$drugid, GDSC@sensitivity$info$drugid))
biomarkers <- as.data.frame(readxl::read_xlsx("biomarkers.xlsx"))
biomarkers <- biomarkers[grep("EXPR", biomarkers$Alteration.type),]
biomarkers <- biomarkers[which(biomarkers$compound %in% intersected_drugs),]
genes_keep <- c("ERBB2","ALK","PHB")
biomarkers <- biomarkers[which(biomarkers$gene %in% genes_keep),]
biomarkers <- biomarkers[c(1,2,4),]

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

#erbb2_CI <- computeCIBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#phb_CI <- computeCIBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#alk_CI <- computeCIBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)

#erbb2_CI$gene_name <- "ERBB2"
#phb_CI$gene_name <- "PHB"
#alk_CI$gene_name <- "ALK"

#isoform_CI <- rbind(erbb2_CI, phb_CI, alk_CI)

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
  
  biomarkers_perm <- as.data.frame(matrix(ncol = 11, nrow = nrow(biomarker_matrix)))
  colnames(biomarkers_perm) <- c("gene","drug", "gCSI_pearson", "CCLE_pearson","GDSC_pearson", "gCSI_pvalue", "CCLE_pvalue", "GDSC_pvalue", "gCSI_sig", "CCLE_sig","GDSC_sig")
  
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
    
    biomarkers_perm$CCLE_pearson[i] <- ccle_per[,,"estimate"]
    biomarkers_perm$CCLE_pvalue[i] <- ccle_per[,,"pvalue"]
    biomarkers_perm$CCLE_sig[i] <- ccle_per[,,"significant"]
    
    biomarkers_perm$GDSC_pearson[i] <- gdsc_per[,,"estimate"]
    biomarkers_perm$GDSC_pvalue[i] <- gdsc_per[,,"pvalue"]
    biomarkers_perm$GDSC_sig[i] <- gdsc_per[,,"significant"]
    
    
    
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


###############
#Isoform-Level#
###############

#erbb2_permut <- computePerBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#phb_permut <- computePerBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)
#alk_permut <- computePerBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms.counts", cells = intersected_rnacells)

#erbb2_permut$gene_name <- "ERBB2"
#phb_permut$gene_name <- "PHB"
#alk_permut$gene_name <- "ALK"

#isoform_permut <- rbind(erbb2_permut, phb_permut, alk_permut)
#save(isoform_permut, file="isoform_permut.RData")

load("isoform_permut.RData")

isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] <- ifelse(isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")] == 0 | is.na(isoform_permut[,c("gCSI_sig","CCLE_sig","GDSC_sig")]), "NO", "YES")



#######################
######Forest_Plot######
#######################

#add stabilities (gcsi/gdsc, gcsi/ccle, gdsc/ccle) to isoform permutations
isoform_permut$gcsi_ccle_stab <- transcript_stability$gcsi_ccle_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$gdsc_ccle_stab <- transcript_stability$gdsc_ccle_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$gcsi_gdsc_stab <- transcript_stability$gcsi_gdsc_spearman[match(isoform_permut$gene, transcript_stability$transcript_id)]
isoform_permut$mean_stability <- rowMeans(isoform_permut[,c(-1:-12)])


isoform_CI$gcsi_ccle_stab <- transcript_stability$gcsi_ccle_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$gdsc_ccle_stab <- transcript_stability$gdsc_ccle_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$gcsi_gdsc_stab <- transcript_stability$gcsi_gdsc_spearman[match(isoform_CI$gene, transcript_stability$transcript_id)]
isoform_CI$mean_stability <- rowMeans(isoform_CI[,c(-1:-18)])


#genes-drug associations to be included in forest plot
genes <- c("PHB", "ALK","ERBB2")
drugs <- c("Paclitaxel", "Crizotinib", "Lapatinib")

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
  
  rows_l <- length(c(NA, gene.CI$gCSI_CI,
                     c(transcript_CI$gCSI_CI),
                     gene.CI$CCLE_CI,
                     c(transcript_CI$CCLE_CI),
                     gene.CI$GDSC_CI,
                     c(transcript_CI$GDSC_CI)))
  
  c_indices <- structure(
    list(
      mean  = c(NA, gene.CI$gCSI_CI,
                c(transcript_CI$gCSI_CI),
                gene.CI$CCLE_CI,
                c(transcript_CI$CCLE_CI),
                gene.CI$GDSC_CI,
                c(transcript_CI$GDSC_CI)),
      
      lower = c(NA, gene.CI$gCSI_lower,
                c(transcript_CI$gCSI_lower),
                gene.CI$CCLE_lower,
                c(transcript_CI$CCLE_lower),
                gene.CI$GDSC_lower,
                c(transcript_CI$GDSC_lower)),
      
      upper = c(NA, gene.CI$gCSI_upper,
                c(transcript_CI$gCSI_upper),
                gene.CI$CCLE_upper,
                c(transcript_CI$CCLE_upper),
                gene.CI$GDSC_upper,
                c(transcript_CI$GDSC_upper))
    ),
    .Names = c("C-index    ", "lower", "upper"),
    row.names = c(NA, -rows_l), 
    class = "data.frame"
  )
  
  
  c_tabletext <- cbind(
    c("PSet", "gCSI(gene)", c(transcript_CI$gene), "CCLE(gene)", c(transcript_CI$gene), "GDSC2(gene)", c(transcript_CI$gene)),
    
    c("C-index", formatC(gene.CI$gCSI_CI, format = "e", digits = 2),
      formatC(c(transcript_CI$gCSI_CI), format = "e", digits = 2),
      formatC(gene.CI$CCLE_CI, format = "e", digits = 2),
      formatC(c(transcript_CI$CCLE_CI), format = "e", digits = 2),
      formatC(gene.CI$GDSC_CI, format = "e", digits = 2),
      formatC(c(transcript_CI$GDSC_CI), format = "e", digits = 2)),
    
    c("P-value \n(permutation)", formatC(gene_perm$gCSI_pvalue, format = "e", digits = 2),
      formatC(transcript_perm$gCSI_pvalue, format = "e", digits = 2),
      formatC(gene_perm$CCLE_pvalue, format = "e", digits = 2),
      formatC(transcript_perm$CCLE_pvalue, format = "e", digits = 2),
      formatC(gene_perm$GDSC_pvalue, format = "e", digits = 2),
      formatC(transcript_perm$GDSC_pvalue, format = "e", digits = 2)),
    
    c("Pearson \n(permutation)", formatC(gene_perm$gCSI_pearson, format = "e", digits = 2),
      formatC(transcript_perm$gCSI_pearson, format = "e", digits = 2),
      formatC(gene_perm$CCLE_pearson, format = "e", digits = 2),
      formatC(transcript_perm$CCLE_pearson, format = "e", digits = 2),
      formatC(gene_perm$GDSC_pearson, format = "e", digits = 2),
      formatC(transcript_perm$GDSC_pearson, format = "e", digits = 2)),
    
    c("Significant \n(permutation)", gene_perm$gCSI_sig,
      transcript_perm$gCSI_sig,
      gene_perm$CCLE_sig,
      transcript_perm$CCLE_sig,
      gene_perm$GDSC_sig,
      transcript_perm$GDSC_sig),
    
    c("Mean stability", NA, formatC(transcript_perm$mean_stability, format = "e", digits = 3),
      NA,
      formatC(transcript_perm$mean_stability, format = "e", digits = 3),
      NA,
      formatC(transcript_perm$mean_stability, format = "e", digits = 3))
  )
  
  
  fn <- local({
    i = 0
    
    b_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))),c("darkgreen",rep("red", nrow(transcript_CI))), c("darkgreen",rep("red", nrow(transcript_CI))))
    l_clrs =  c(c("darkgreen",rep("red", nrow(transcript_CI))),c("darkgreen",rep("red", nrow(transcript_CI))), c("darkgreen",rep("red", nrow(transcript_CI))))
    function(..., clr.line, clr.marker){
      i <<- i + 1
      fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
      #fpDrawSummaryCI(...,col=s_clrs[i])
    }
  })
  
  fileName = paste0("figures/",gene,"_",drug,".pdf")
  pdf(fileName, width=12 , height=15, onefile=FALSE)
  
  forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,rep(F,rows_l)), xlab="Concordance Index", 
             title=paste0(gene,"_",drug), zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:6, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                            ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                            xlab=gpar(fontfamily = "", cex=0.8, fontface=2),
                            legend=gpar(fontfamily = "", cex = 0.5, fontface=1)),
             fn.ci_norm = fn,
             col=fpColors(summary="blue"), 
             xticks= c(.3, .4, .5, .6, .7, .8)
  )
  
  dev.off()
  
}
