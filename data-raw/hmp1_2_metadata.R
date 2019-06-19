library(usethis)

hmp1_2_metadata <- read.csv2("data-raw/hmp1_2_metadata.tsv",
                          sep = "\t", header = T)
hmp1_2_metadata <- hmp1_2_metadata[which(hmp1_2_metadata$SRS != "#N/A"),]

use_data(hmp1_2_metadata, overwrite = TRUE)