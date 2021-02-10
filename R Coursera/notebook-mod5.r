
# no library to load

# cell #1
cnvLogs <- read.table("readonly/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", 
            header = T, fill = T)
dim(cnvLogs)
# 284458 6
head(cnvlogs)
#             Sample Chromosome     Start       End Num_Probes Segment_Mean
#1 TCGA-3C-AAAU-10A-01D-A41E-01          1   3218610  95674710      53225       0.0055

# cell #2
summary(cnvLogs)
# mean = -0.1132 median = 0


# cell #3
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-01A-11D-A41E-01"))

# cell #4
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-01A-11D-A41E-01" & Segment_Mean > 0))

# cell #5
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-01A-11D-A41E-01" & Segment_Mean < 0))

# cell #6
mean(subset(cnvLogs, Sample == "TCGA-3C-AAAU-01A-11D-A41E-01")[["Segment_Mean"]])
