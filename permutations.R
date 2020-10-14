library(PharmacoGx)
library(survcomp)
library(gimme)
library(dplyr)
library(forestplot)

#goal: is most stable isoforms most predictive of drug response?

setwd({"~/Desktop/permutation"})

load("inter.RData")

gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")
CCLE <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CCLE.rds")
GDSC <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC2-8.2.rds")

S4Vectors::metadata(gCSI@molecularProfiles[["Kallisto_0.46.1.isoforms"]])$annotation <- "isoform"
S4Vectors::metadata(CCLE@molecularProfiles[["Kallisto_0.46.1.isoforms"]])$annotation <- "isoform"
S4Vectors::metadata(GDSC@molecularProfiles[["Kallisto_0.46.1.isoforms"]])$annotation <- "isoform"

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


#common drugs across all 3 datasets
intersected_drugs <- Reduce(intersect, list(gCSI@sensitivity$info$drugid, CCLE@sensitivity$info$drugid, GDSC@sensitivity$info$drugid))

#top known EXPRESSION biomarkers with intersected drugs (n = 15)
biomarkers <- as.data.frame(readxl::read_xlsx("~/Desktop/permutation/biomarkers.xlsx"))
biomarkers <- biomarkers[grep("EXPR", biomarkers$Alteration.type),]
biomarkers <- biomarkers[which(biomarkers$compound %in% intersected_drugs),]


###################
#Concordance Index#
###################

###########
###GENES###
###########

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


computeCIBiomarker <- function(biomarker_matrix, mData, cells) {

biomarkers_CI <- as.data.frame(matrix(ncol = 17, nrow = nrow(biomarker_matrix)))
colnames(biomarkers_CI) <- c("gene","drug", "gCSI_CI", "CCLE_CI","GDSC_CI", "gCSI_pvalue", "CCLE_pvalue", "GDSC_pvalue", "gCSI_se", "CCLE_se","GDSC_se", "gCSI_upper",
                             "CCLE_upper","GDSC_upper","gCSI_lower","CCLE_lower","GDSC_lower")
for (i in 1:nrow(biomarker_matrix)){
  
  drug <- biomarker_matrix$compound[i]
  gene <- biomarker_matrix$gene[i]
  
  if (!mData == "Kallisto_0.46.1.isoforms"){
    gene_select <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[which(gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name == gene)]
  }else{
    gene_select = gene
  }
  
  gcsi_CI <- computeCI(PSet=gCSI, mData = mData, features = gene_select, cells = cells)
  ccle_CI <- computeCI(PSet=CCLE, mData = mData, features = gene_select, cells = cells)
  gdsc_CI <- computeCI(PSet=GDSC, mData = mData, features = gene_select, cells = cells)
  
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


#biomarkers_CI_inter <- computeCIBiomarker(biomarker_matrix = biomarkers, mData = "Kallisto_0.46.1.rnaseq", cells = intersected_rnacells)
#save(biomarkers_CI_inter, file="gene_compute_CI.RData")

load("gene_compute_CI.RData")

##############
###Isoforms###
##############


##import transcript stability data generated by transcript_stability.R
transcript_stability <- read.csv("transcript_stability.csv")

##get isoforms for each biomarker
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

#erbb2_isoforms_CI <- computeCIBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)
#phb_isoforms_CI <- computeCIBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)
#alk_isoforms_CI <- computeCIBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)

#save(erbb2_isoforms_CI, phb_isoforms_CI, alk_isoforms_CI, file="isoforms_compute_CI.RData")

load("isoforms_compute_CI.RData")

#ERBB2

erbb2 <- transcript_stability[which(transcript_stability$transcript_id %in% erbb2_isoforms$gene),c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman","transcript_id")]
erbb2$mean_stability <- rowMeans(erbb2[,-4])
erbb2 <- erbb2[order(erbb2$mean_stability, decreasing = T),]
#erbb2 <- erbb2[c(1,6,11),] #get range of stabilities (high to low)



######################
####Permutations######
######################


PermutationSig <- function(PSet, mData, features, cells, drugs){
  
  signature <- drugSensitivitySig(PSet, mData, drugs=drugs, 
                                  features=features, sensitivity.measure = "aac_recomputed", modeling.method="pearson", 
                                  inference.method="resampling", cells=cells, nthread=2, parallel.on = "gene")
  return(signature)
}



computePerBiomarker <- function(biomarker_matrix, mData, cells) {
  
  biomarkers_perm <- as.data.frame(matrix(ncol = 11, nrow = nrow(biomarker_matrix)))
  colnames(biomarkers_perm) <- c("gene","drug", "gCSI_pearson", "CCLE_pearson","GDSC_pearson", "gCSI_pvalue", "CCLE_pvalue", "GDSC_pvalue", "gCSI_sig", "CCLE_sig","GDSC_sig")
  
  for (i in 1:nrow(biomarker_matrix)){
    drug <- biomarker_matrix$compound[i]
    gene <- biomarker_matrix$gene[i]
    
    if (!mData == "Kallisto_0.46.1.isoforms"){
      gene_select <- gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_id[match(gene, gCSI@molecularProfiles$Kallisto_0.46.1.rnaseq@elementMetadata$gene_name)]
    }else{
      gene_select = gene
    }
    

    gcsi_per <- PermutationSig(PSet = gCSI, mData = mData, features = gene_select, 
                               cells = cells, drugs = drug)
    
    ccle_per <- PermutationSig(PSet = CCLE, mData = mData, features = gene_select, 
                               cells = cells, drugs = drug)

    gdsc_per <- PermutationSig(PSet = GDSC, mData = mData, features = gene_select, 
                               cells = cells, drugs = drug)
    
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

###########
###GENES###
###########

#biomarkers_perm_inter <- computePerBiomarker(biomarker_matrix = biomarkers, mData = "Kallisto_0.46.1.rnaseq", cells = intersected_rnacells)
#save(biomarkers_perm_inter, file="gene_compute_permutation.RData")

load("gene_compute_permutation.RData")

##############
###Isoforms###
##############

#erbb2_perm_inter <- computePerBiomarker(biomarker_matrix = erbb2_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)
#phb_perm_inter <- computePerBiomarker(biomarker_matrix = phb_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)
#alk_perm_inter <- computePerBiomarker(biomarker_matrix = alk_isoforms, mData = "Kallisto_0.46.1.isoforms", cells = intersected_rnacells)
#save(erbb2_perm_inter, phb_perm_inter, alk_perm_inter, file = "biomarker_compute_permutation_iso.RData")

load("biomarker_compute_permutation_iso.RData")



#############################
######Forest_Plot (Gene)#####
#############################

genes <- c("PHB", "ALK","ERBB2")
drugs <- c("Paclitaxel", "Crizotinib", "Lapatinib")

for (i in 1:length(genes)){
  gene <- genes[i]
  drug <- drugs[i]
  
  gene_CI <- biomarkers_CI_inter[which(biomarkers_CI_inter$gene == gene & biomarkers_CI_inter$drug == drug),]
  gene_perm <- biomarkers_perm_inter[which(biomarkers_perm_inter$gene == gene & biomarkers_perm_inter$drug == drug),]
  
  
  c_indices <- structure(
    list(
      mean  = c(NA, gene_CI$gCSI_CI, gene_CI$CCLE_CI, gene_CI$GDSC_CI),
      lower = c(NA, gene_CI$gCSI_lower, gene_CI$CCLE_lower, gene_CI$GDSC_lower),
      upper = c(NA, gene_CI$gCSI_upper, gene_CI$CCLE_upper, gene_CI$GDSC_upper)
    ),
    .Names = c("C-index    ", "lower", "upper"),
    row.names = c(NA, -4L), 
    class = "data.frame"
  )
  
  
  c_tabletext <- cbind(
    c("PSet", "gCSI", "CCLE", "GDSC2"),
    c("N", 48, 48, 48), #common samples for each dataset
    c("C-index", formatC(gene_CI$gCSI_CI, format = "e", digits = 2), formatC(gene_CI$CCLE_CI, format = "e", digits = 2),formatC(gene_CI$GDSC_CI, format = "e", digits = 2)),
    c("P-value \n(permutation)", formatC(gene_perm$gCSI_pvalue, format = "e", digits = 2), formatC(gene_perm$CCLE_pvalue, format = "e", digits = 2), formatC(gene_perm$GDSC_pvalue, format = "e", digits = 2)),
    c("Pearson \n(permutation)", formatC(gene_perm$gCSI_pearson, format = "e", digits = 2),formatC(gene_perm$CCLE_pearson, format = "e", digits = 2), formatC(gene_perm$GDSC_pearson, format = "e", digits = 2)),
    c("Significant \n(permutation)", formatC(gene_perm$gCSI_sig),formatC(gene_perm$CCLE_sig), formatC(gene_perm$GDSC_sig))
    
  )
  
  
  fileName = paste("figures/",gene,"_",drug,".pdf")
  pdf(fileName, width=9, height=3, onefile=FALSE)
  
  forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,F,F,F), xlab="Concordance Index", 
             title="", zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:6, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                            ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                            xlab=gpar(fontfamily = "", cex=0.8, fontface=2),
                            legend=gpar(fontfamily = "", cex = 1, fontface=1)),
             col=fpColors(box=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          line=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          summary="blue"),
             xticks= c(.4, .5, .6, .7, .8)
  )
  
  dev.off()
  
}




















######################
#####Forest_Plot######
######################

#selecting biomarkers with gene/drug significance computed through permutations

#ERBB2 + Lapatinib
biomarkers_per_inter_erbb2 <- biomarkers_per_inter[7,] 
biomarkers_ci_inter_erbb2 <- biomarkers_CI_inter[7,] 

#ALK + Crizotinib
biomarkers_per_inter_alk <- biomarkers_per_inter[5,] 
biomarkers_ci_inter_alk <- biomarkers_CI_inter[5,] 

#PHB + Paclitaxel
biomarkers_per_inter_phb <- biomarkers_per_inter[3,] 
biomarkers_ci_inter_phb <- biomarkers_CI_inter[3,] 


#ERBB2 + Lapatinib (PLOT)

c_indices <- structure(
  list(
    mean  = c(NA, biomarkers_ci_inter_erbb2$gCSI_CI, biomarkers_ci_inter_erbb2$CCLE_CI, biomarkers_ci_inter_erbb2$GDSC_CI),
    lower = c(NA, biomarkers_ci_inter_erbb2$gCSI_lower, biomarkers_ci_inter_erbb2$CCLE_lower, biomarkers_ci_inter_erbb2$GDSC_lower),
    upper = c(NA,biomarkers_ci_inter_erbb2$gCSI_upper, biomarkers_ci_inter_erbb2$CCLE_upper, biomarkers_ci_inter_erbb2$GDSC_upper)
  ),
  .Names = c("C-index    ", "lower", "upper"),
  row.names = c(NA, -4L), 
  class = "data.frame"
)


c_tabletext <- cbind(
  c("PSet", "gCSI", "CCLE", "GDSC2"),
  c("N", 48, 48, 48), #common samples for each dataset
  c("C-index", formatC(biomarkers_ci_inter_erbb2$gCSI_CI, format = "e", digits = 2), formatC(biomarkers_ci_inter_erbb2$CCLE_CI, format = "e", digits = 2),formatC(biomarkers_ci_inter_erbb2$GDSC_CI, format = "e", digits = 2)),
  c("P-value \n(permutation)", formatC(biomarkers_per_inter_erbb2$gCSI_pvalue, format = "e", digits = 2), formatC(biomarkers_per_inter_erbb2$CCLE_pvalue, format = "e", digits = 2), formatC(biomarkers_per_inter_erbb2$GDSC_pvalue, format = "e", digits = 2)),
  c("Pearson \n(permutation)", formatC(biomarkers_per_inter_erbb2$gCSI_pearson, format = "e", digits = 2),formatC(biomarkers_per_inter_erbb2$CCLE_pearson, format = "e", digits = 2), formatC(biomarkers_per_inter_erbb2$GDSC_pearson, format = "e", digits = 2))
)


fileName = "erbb2_lapatinib.pdf"
pdf(fileName, width=7, height=3, onefile=FALSE)

forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,F,F,F), xlab="Concordance Index", 
           title="", zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:5, col = "#000044")),
           txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                          ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                          xlab=gpar(fontfamily = "", cex=0.8, fontface=2),
                          legend=gpar(fontfamily = "", cex = 1, fontface=1)),
           col=fpColors(box=RColorBrewer::brewer.pal(n=4, name="Set2"),
                        line=RColorBrewer::brewer.pal(n=4, name="Set2"),
                        summary="blue"),
           xticks= c(.4, .5, .6, .7, .8)
)

dev.off()





