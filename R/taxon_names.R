#' @title shortenTaxon
#'
#' @description shortens taxon names
#' @param dat character vector holding taxonomic names
#' @param level column that holds taxonomic level of interest
#' @param delim delimiter that seperates levels (e.g. |)
#'   information
#' @export
shortenTaxons <- function(dat, level = 2, delim = "|") {
  seperated <-
    as.data.frame(do.call(rbind, strsplit(dat, paste0("\\", delim))))
  stopifnot(ncol(seperated) > 1)
  seperated <- seperated[[level]]
  if (length(grep("__", seperated) > 0)) {
    taxa <-
      stringr::str_sub(seperated, 4, stringr::str_length(seperated))
  }
  if (length(grep("_", taxa) > 0)) {
    taxa <- gsub("_", " ", taxa)
  }
  return(taxa)
}