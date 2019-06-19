humann2_table <- as.data.frame(read.table("data-raw/humann2_table.tsv", header = T, sep = '\t', stringsAsFactors = F))
usethis::use_data(humann2_table, overwrite = TRUE)