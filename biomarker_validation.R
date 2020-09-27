library(survcomp)
library(gimme)
library(dplyr)
library(forestplot)

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


gcsi_ccle <- read.csv("~/Desktop/transcript_stability_gcsi_ccle.csv")
gcsi_gdsc<- read.csv("~/Desktop/transcript_stability_gcsi_gdsc.csv")
gdsc_ccle<- read.csv("~/Desktop/transcript_stability_gdsc_ccle.csv")

gcsi_ccle <- gcsi_ccle[which(gcsi_ccle$stability=="YES"),]
gcsi_gdsc <- gcsi_gdsc[which(gcsi_gdsc$stability=="YES"),]
gdsc_ccle <- gdsc_ccle[which(gdsc_ccle$stability=="YES"),]


gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")
CCLE <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CCLE.rds")
GDSC <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC2-8.2.rds")
GDSC8.0 <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC2-8.0.rds")
GDSC1_8.2 <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/GDSC1-8.2.rds")
CTRPv2 <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/CTRPv2.rds")

#intersect stable transcripts across all dataset pairs

intersected_transcripts <- Reduce(intersect, list(gcsi_ccle$transcript_id_ver, gcsi_gdsc$transcript_id_ver, gdsc_ccle$transcript_id_ver))


gCSI_sens <- as.data.frame(summarizeSensitivityProfiles(gCSI, sensitivity.measure = "aac_recomputed",  fill.missing = F))
CCLE_sens <- as.data.frame(summarizeSensitivityProfiles(CCLE, sensitivity.measure = "aac_recomputed",  fill.missing = F))
GDSC_sens <- as.data.frame(summarizeSensitivityProfiles(GDSC, sensitivity.measure = "aac_recomputed",  fill.missing = F))
CTRPv2_sens <- as.data.frame(summarizeSensitivityProfiles(CTRPv2, sensitivity.measure = "aac_recomputed",  fill.missing = F))


gCSI_rna <- summarizeMolecularProfiles(gCSI, mDataType = "Kallisto_0.46.1.isoforms.counts")
CCLE_rna <- summarizeMolecularProfiles(CCLE, mDataType = "Kallisto_0.46.1.isoforms.counts")
GDSC_rna <- summarizeMolecularProfiles(GDSC, mDataType = "Kallisto_0.46.1.isoforms.counts")


gCSI_rna <- as.data.frame(gCSI_rna@assays@data$exprs)
CCLE_rna <- as.data.frame(CCLE_rna@assays@data$exprs)
GDSC_rna <- as.data.frame(GDSC_rna@assays@data$exprs)


computeCI <- function(rnaseq, sensitivity_data){
drugs <- unique(rownames(sensitivity_data))
commonSamples <- intersected_rnacells
commonSamples <- intersected_rnacells[which(intersected_rnacells %in% colnames(sensitivity_data))]
combinations <- as.data.frame(expand.grid.unique(intersected_transcripts, drugs, include.equals = TRUE))
combinations$ci <- NA
combinations$pvalue <- NA
combinations$se <- NA
combinations$upper <- NA
combinations$lower <- NA
colnames(combinations) <- c("transcript","drug","ci","pvalue","se","upper","lower")

for (i in 1:nrow(combinations)){
  print(paste0(i, " out of ", nrow(combinations), " complete"))
  tt <- sensitivity_data[combinations[,2][i], commonSamples]
  tt[which(is.na(tt))] <- 0 #some sensitivities are NA due to filterNoisyCurve function, which causes error when running CI with survcomp
  ci <- survcomp::concordance.index(as.numeric(tt), surv.time = as.numeric(unlist(-rnaseq[combinations[,1][i], commonSamples])), surv.event = rep(1,length(sensitivity_data[commonSamples])),outx = F, method="noether")
  combinations$pvalue[i] <- ci$p.value
  combinations$ci[i] <- ci$c.index
  combinations$se[i] <- ci$se
  combinations$upper[i] <- ci$upper
  combinations$lower[i] <- ci$lower
}

return(combinations)

}


gcsi_CI <- computeCI(gCSI_rna, gCSI_sens)
ccle_CI <- computeCI(CCLE_rna, CCLE_sens)
gdsc_CI <- computeCI(GDSC_rna, GDSC_sens)
#ctrpv2_CI <-  computeCI(CCLE_rna, CTRPv2_sens)

#save(gcsi_CI, ccle_CI, gdsc_CI, file="~/Desktop/stable_CI.RData")
load("~/Desktop/stable_CI.RData")

#filter for CI > 0.65 & p-value < 0.05 
gcsi_CI_fil <- gcsi_CI[which(gcsi_CI$ci > 0.6 & gcsi_CI$pvalue < 0.05),]
ccle_CI_fil <- ccle_CI[which(ccle_CI$ci > 0.6 & ccle_CI$pvalue < 0.05),]
gdsc_CI_fil <- gdsc_CI[which(gdsc_CI$ci > 0.6 & gdsc_CI$pvalue < 0.05),]

rownames(gcsi_CI_fil) <- paste0(gcsi_CI_fil$transcript,"_",gcsi_CI_fil$drug)
rownames(ccle_CI_fil) <- paste0(ccle_CI_fil$transcript,"_",ccle_CI_fil$drug)
rownames(gdsc_CI_fil) <- paste0(gdsc_CI_fil$transcript,"_",gdsc_CI_fil$drug)

intersected_ids <- Reduce(intersect, list(rownames(gcsi_CI_fil), rownames(ccle_CI_fil), rownames(gdsc_CI_fil)))


gcsi_xx <- gcsi_CI_fil[intersected_ids,]
ccle_xx <- ccle_CI_fil[intersected_ids,]
gdsc_xx <- gdsc_CI_fil[intersected_ids,]


gcsi_xx <- arrange(gcsi_xx, pvalue, desc(ci))
ccle_xx <- arrange(ccle_xx, pvalue, desc(ci))
gdsc_xx <- arrange(gdsc_xx, pvalue, desc(ci))



####ENST00000275493.7_Erlotinib validation
gcsi_ci <- gcsi_xx[grep("ENST00000275493.7_Erlotinib", rownames(gcsi_xx)),"ci"]
ccle_ci <- ccle_xx[grep("ENST00000275493.7_Erlotinib", rownames(ccle_xx)),"ci"]
gdsc_ci <- gdsc_xx[grep("ENST00000275493.7_Erlotinib", rownames(gdsc_xx)),"ci"]

gcsi_serr <- gcsi_xx[grep("ENST00000275493.7_Erlotinib", rownames(gcsi_xx)),"se"]
ccle_serr <- ccle_xx[grep("ENST00000275493.7_Erlotinib", rownames(ccle_xx)),"se"]
gdsc_serr <- gdsc_xx[grep("ENST00000275493.7_Erlotinib", rownames(gdsc_xx)),"se"]

gcsi_lower <- gcsi_xx[grep("ENST00000275493.7_Erlotinib", rownames(gcsi_xx)),"lower"]
ccle_lower <- ccle_xx[grep("ENST00000275493.7_Erlotinib", rownames(ccle_xx)),"lower"]
gdsc_lower <- gdsc_xx[grep("ENST00000275493.7_Erlotinib", rownames(gdsc_xx)),"lower"]

gcsi_upper <- gcsi_xx[grep("ENST00000275493.7_Erlotinib", rownames(gcsi_xx)),"upper"]
ccle_upper <- ccle_xx[grep("ENST00000275493.7_Erlotinib", rownames(ccle_xx)),"upper"]
gdsc_upper <- gdsc_xx[grep("ENST00000275493.7_Erlotinib", rownames(gdsc_xx)),"upper"]

gcsi_pvalue <- gcsi_xx[grep("ENST00000275493.7_Erlotinib", rownames(gcsi_xx)),"pvalue"]
ccle_pvalue <- ccle_xx[grep("ENST00000275493.7_Erlotinib", rownames(ccle_xx)),"pvalue"]
gdsc_pvalue <- gdsc_xx[grep("ENST00000275493.7_Erlotinib", rownames(gdsc_xx)),"pvalue"]


combined_ci <- combine.est(
  c(
    gcsi_ci, ccle_ci, gdsc_ci
  ),
  c(
    gcsi_serr, ccle_serr, gdsc_serr
  ),na.rm = TRUE,hetero = TRUE)


combined_ci_lower <- combined_ci$estimate + qnorm(0.025, lower.tail=TRUE) *  combined_ci$se
combined_ci_upper <- combined_ci$estimate + qnorm(0.025, lower.tail=FALSE) *  combined_ci$se
combined_ci_p <- pnorm((combined_ci$estimate - 0.5)/combined_ci$se, lower.tail = combined_ci$estimate < 0.5) * 2


######CREATE VALIDATION EGFR FOREST PLOT######

c_indices <- structure(
  list(
    mean  = c(NA, gcsi_ci, ccle_ci, gdsc_ci, combined_ci$estimate),
    lower = c(NA, gcsi_lower, ccle_lower, gdsc_lower, combined_ci_lower),
    upper = c(NA, gcsi_upper, ccle_upper, gdsc_upper, combined_ci_upper)
  ),
  .Names = c("C-index    ", "lower", "upper"),
  row.names = c(NA, -5L), 
  class = "data.frame"
)

c_tabletext <- cbind(
  c("PSet", "gCSI", "CCLE", "GDSC2","Meta analysis"),
  c("N", 48, 48, 48, 48), #common samples for each dataset
  c("C-index", formatC(gcsi_ci, format = "e", digits = 2), formatC(ccle_ci, format = "e", digits = 2),formatC(gdsc_ci, format = "e", digits = 2) ,formatC(combined_ci$estimate, format = "e", digits = 2)),
  c("P-value", formatC(gcsi_pvalue, format = "e", digits = 2), formatC(ccle_pvalue, format = "e", digits = 2), formatC(gdsc_pvalue, format = "e", digits = 2),formatC(combined_ci_p, format = "e", digits = 2))
)



fileName = "~/Desktop/ci/CI.pdf"
pdf(fileName, width=7, height=3, onefile=FALSE)

forestplot(c_tabletext, c_indices, new_page = TRUE, boxsize = 0.3, is.summary=c(T,F,F,F), xlab="Concordance Index", 
           title="", zero=c(.49, .51),hrzl_lines=list("2"=gpar(lty=2, columns=1:4, col = "#000044")),
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


#ENST00000452441.5    Erlotinib (gCSI) #DDR1

#ENST00000275493.7_Erlotinib 




