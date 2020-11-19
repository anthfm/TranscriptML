#The coverage value shown in the output of StringTie (and in other genomics programs) is an average of 
#all per-base coverages across the length of genomic segment (exon) or set of segments (transcript)
options(stringsAsFactors = FALSE)
#need to get first(5')/last(3') 100bp regions of transcript to calculate 5' and 3' bias

gCSI <- readRDS("/Users/anthmam/Desktop/Projects/BHKLAB/PSets/gCSI.rds")
feat <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms@elementMetadata
feat <- feat[,c("seqnames","start","end","strand","transcript_id")]

#Isolate + strand transcripts (5' begins at "start")
feat_pos <- as.data.frame(feat[which(feat$strand=="+"),])
feat_pos$five.prime.start <- feat_pos$start
feat_pos$five.prime.end <- feat_pos$start + 99

feat_pos$three.prime.start <- feat_pos$end - 99
feat_pos$three.prime.end <- feat_pos$end


#Isolate - strand transcripts (5' begins at "end")
feat_neg <- as.data.frame(feat[which(feat$strand=="-"),])
feat_neg$five.prime.start <- feat_neg$end - 99
feat_neg$five.prime.end <- feat_neg$end

feat_neg$three.prime.start <- feat_neg$start
feat_neg$three.prime.end <- feat_neg$start + 99


combined_five_prime <- rbind(feat_pos[,c("seqnames","five.prime.start","five.prime.end","transcript_id")],
                             feat_neg[,c("seqnames","five.prime.start","five.prime.end","transcript_id")])

colnames(combined_five_prime) <- c("chr","start","end","name")
combined_five_prime$start <- as.integer(combined_five_prime$start)
combined_five_prime$end <- as.integer(combined_five_prime$end)

combined_three_prime <- rbind(feat_pos[,c("seqnames","three.prime.start","three.prime.end","transcript_id")],
                             feat_neg[,c("seqnames","three.prime.start","three.prime.end","transcript_id")])

colnames(combined_three_prime) <- c("chr","start","end","name")
combined_three_prime$start <- as.integer(combined_three_prime$start)
combined_three_prime$end <- as.integer(combined_three_prime$end)

write.table(combined_five_prime, file="~/Desktop/bed/five-prime.bed", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(combined_three_prime, file="~/Desktop/bed/three-prime.bed", col.names = F, row.names = F, quote = F, sep = "\t")


#./mosdepth -b five-prime.bed test587007 /cluster/projects/bhklab/procdata/circRNA_human/rnaseq/202007/Bias/HISAT2/gCSI/test587007/test587007_sorted.bam
