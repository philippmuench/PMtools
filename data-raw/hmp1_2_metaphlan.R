hmp1_2_metaphlan <- read.table("data-raw/hmp1_2_metaphlan.txt", header=T, sep="\t")

usethis::use_data(hmp1_2_metaphlan , overwrite = TRUE)