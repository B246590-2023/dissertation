
setwd("\\\\cmvm.datastore.ed.ac.uk\\cmvm\\eb\\groups\\clark_grp2\\Yuxin\\usefulScriptsAndResources\\")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rtracklayer")
install.packages("dplyr")
install.packages("tximport")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("readr")
install.packages("tidyr")
install.packages("dplyr")
install.packages("purrr")
install.packages("skimr")
install.packages("writexl")
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}


library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(DESeq2)
library(dplyr)
library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(skimr)
library(writexl)
library(readxl)

list.files()

gtf <- import('./GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.gtf.gz')

gtf_df <- as.data.frame(gtf)
gtfdf2 <- filter(gtf_df, type == "transcript")
glimpse(gtfdf2)


tx2gene <- gtfdf2[c("transcript_id","gene")]
head(tx2gene)

skimr::skim(tx2gene)
write.table(tx2gene, file = "tx2gene.txt")
getwd()

head(tx2gene)

dir <- "./kallisto_output"
subdirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
files <- file.path(subdirs, "abundance.tsv") 

samples <- data.frame(
  sample = basename(subdirs),
  files = files,
  stringsAsFactors = FALSE
)


example_quant_file <- read.table(files[1], header = TRUE)
head(example_quant_file)




FIXXER <- function(FILE, ADDRESS){
  tmp <- read_tsv(FILE)
  tmp <- tidyr::separate(tmp, col = target_id, into = c("col1", "col2", "col3", "col4", "col5", "col6"), sep = "_", remove = TRUE, extra = "merge", fill = "right")
  missing_parts <- tmp %>% filter(is.na(col4) | is.na(col5))
  if (nrow(missing_parts) > 0) {
    warning("There are rows with missing parts in target_id. These rows will be skipped.")
    print(missing_parts)  
    tmp <- tmp %>% filter(!is.na(col4) & !is.na(col5))
  }
  tmp <- dplyr::select(tmp, -col1, -col2, -col3, -col6)
  tmp$target_id <- paste(tmp$col4, tmp$col5, sep = "_")
  tmp <- dplyr::select(tmp, -col4, -col5)
  tmp <- dplyr::relocate(tmp, target_id)
  write_tsv(tmp, paste0(ADDRESS, "/abundance_new.tsv"))
}


write_dir <- subdirs
purrr::map2(files, write_dir, FIXXER)
samples$new_files <- file.path(write_dir, "abundance_new.tsv")


file_check <- file.exists(samples$new_files)
print(file_check)


missing_files <- samples$new_files[!file.exists(samples$new_files)]
if (length(missing_files) > 0) {
  print("Do not exist：")
  print(missing_files)
} else {

  txi <- tximport(samples$new_files, type = "kallisto", txOut = FALSE, tx2gene = tx2gene)
}
summary(txi)
head(txi$counts)
head(txi$abundance)
head(txi$length)
txi$counts[1:2, 1:4]

str(txi)
names(txi)
txi$abundance[1:6, 1:4]

saveRDS(txi, file = "txi.rds")

txi_counts_df <- as.data.frame(txi$counts)
txi_counts_df$Gene <- rownames(txi_counts_df)


txi_counts_df <- txi_counts_df[, c("Gene", colnames(txi_counts_df)[1:ncol(txi$counts)])]

write_xlsx(txi_counts_df, path = "C:/Users/name/Desktop/txi_counts_with_gene_names.xlsx")



getwd()

saveRDS(txi, file = "txi.rds")

txi_list <- list(
  counts = as.data.frame(txi$counts),
  length = as.data.frame(txi$length),
  abundance = as.data.frame(txi$abundance),
  countsFromAbundance = as.data.frame(txi$countsFromAbundance)
)
write_xlsx(txi_list, path = "txi_data.xlsx")

write_xlsx(txi_list, path = "C:/Users/name/Desktop/txi_data.xlsx")
sample_names <- colnames(txi$counts)

data <- read_excel("C:/Users/name/Desktop/stron.xlsx")
data <- read_excel("C:/Users/name/Desktop/nem.xlsx")

str(data)

head(data)



















