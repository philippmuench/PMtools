% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humann2_utils.R
\name{plotHummanMetadata}
\alias{plotHummanMetadata}
\title{plotHummanMetadata}
\usage{
plotHummanMetadata(
  metaphlan,
  distance.method = "bray",
  cluster.method = "single",
  metadata = hmp1_2_metadata
)
}
\arguments{
\item{metaphlan}{MetaPhlAn2 table}

\item{distance.method}{Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".}

\item{cluster.method}{The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).}
}
\description{
outputs order of samples by bray curtis
}
