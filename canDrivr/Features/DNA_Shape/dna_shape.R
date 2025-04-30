# Get DNA shapes of 10 BP regions overlapping with variant
# Gets the shape of the wild-type regions with the ref allele

# Get DNA shapes of 10 BP regions overlapping with variant

# options(repos = c(CRAN = "https://cran.rstudio.com/"))
# install.packages(c("bios2mds", "BSgenome.Hsapiens.UCSC.hg38", "DNAshapeR", "usethis", "dplyr", "withr", "ggplot2", "tzdb", "readr", "tidyverse", "data.table", "ps", "devtools", "foreach", "iterators", "doParallel", "doSNOW", "Peptides", "stringr", "tidyr"))

library(bios2mds)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DNAshapeR)
library(usethis)
library(dplyr)
library(withr)
library(ggplot2)
library(tzdb)
library(readr)
library(tidyverse)
library(data.table)
library(ps)
library(devtools)
library(foreach)
library(iterators)
library(doParallel)
library(doSNOW)
library(Peptides)
library(stringr)
library(tidyr)

args <- commandArgs()
# print(args)
# Reads in variant file in the format: "chrom", "pos", "end", "ref", "alt"
featureDir=args[6]
featureOutputDir=args[8]
variants = paste(featureDir, args[7], sep = "") # name of the variants file

# Import variants for shaping
variants <- read.csv(variants, header = TRUE)
# print(colnames(variants))

colnames(variants) <-  c("chrom", "start", "stop", "ref_allele", "alt_allele", "driver_stat")
varTable <- variants
# print(head(varTable))


# Get the desired base pair range for DNA shape (range before was 10bp)
variants[2] = variants[2] - 10
variants[3] = variants[3] + 10


# # Make a GRRanges object
# # print(variants)
variants <- makeGRangesFromDataFrame(variants, start.field="start", end.field="stop")
variants <- trim(variants)
# Ensure row consistency
# valid_idx <- which(varTable$start %in% start(variants))
invalid_idx <- which((varTable$start %in% start(variants)))
varTable <- varTable[setdiff(1:nrow(varTable), invalid_idx), , drop = FALSE]
print(head(variants))

# # Get the 10bp fasta for each variant
# variants$wildtype_seq = getSeq(Hsapiens, names=variants$chrom, start=variants$start_pos, end=variants$end_pos)
# print(class(variants))  # Check its type
# print(str(variants))    # View structure
# print(head(variants))   # See first few rows


getFasta(variants, BSgenome = Hsapiens, width = 20, filename = "hg38.fa")
fn <- "hg38.fa"

# # Get the shape properties of each position for each variant
pred <- getShape(fn)


# pred <- pred[!sapply(pred, is.null)]

# # Reduce all of the properties into a single matrix
dnaShape = Reduce("cbind", pred)

if (nrow(varTable) != nrow(dnaShape)) {
    print(nrow(varTable))
    print(nrow(dnaShape))
    stop("Final row mismatch! Check filtered variants.")
}

dnaShape = cbind(varTable, dnaShape)

colnames(dnaShape) = c(colnames(dnaShape)[1:6], paste(1:20, "MGW", sep = "_"), paste(1:19, "HelT", sep = "_"), paste(1:20, "ProT", sep = "_"),
                 paste(1:19, "Roll", sep = "_"), paste(1:20, "EP", sep = "_"))

dnaShape = dnaShape[-3]
# print(length(dnaShape))

# Write DNA shape properties to CSV
name = paste(featureOutputDir, "dnaShape_indtest.txt", sep = "")

# Remove columns where all values are NA
dnaShape = dnaShape[,colSums(is.na(dnaShape))<nrow(dnaShape)]

write.table(dnaShape, name, quote = FALSE, row.names = FALSE, sep = ",")

unlink(paste(featureDir, "*.fa*", sep = ""))
