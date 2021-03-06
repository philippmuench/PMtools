% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humann2_tools.R
\name{humann2Barplot}
\alias{humann2Barplot}
\title{humann2Barplot}
\usage{
humann2Barplot(
  humann2.table,
  num.bugs = "auto",
  num.bugs.explained.fraction = 0.25,
  feature = "Cas1",
  featue.column = 1,
  taxa.column = 2,
  metadata,
  metadata.id = 9,
  metadata.factor = 4,
  use.custom.column.stratification = T,
  order.by = "bc",
  column.stratification.order = c("Oral", "Skin", "Vaginal", "Gut"),
  custom.order,
  last.plot.colors = NULL
)
}
\arguments{
\item{humann2.table}{a dataframe containing the HUMAnN2 output in tsv format}

\item{num.bugs}{number of most abundant bugs to show as statification, set to
to "auto" if this should be ajusted based on num.bugs.explained.fraction}

\item{num.bugs.explained.fraction}{fraction that needs to be explained by the
number of bugs shown}

\item{feature}{name of gene in humann2_table that will be plotted}

\item{featue.column}{index of column in humann2 table that holds the feature}

\item{taxa.column}{index of column in humann2 table that holds the taxon infomation}

\item{metadata}{a dataframe containing the metadata}

\item{metadata.id}{column in metadata file that holds the smaple id (colum names in humann2 table)}

\item{metadata.factor}{column in metadata file that holds the annotation of interest}

\item{use.custom.column.stratification}{custom orderin of row statification}

\item{order.by}{how to order bars within one statification}

\item{column.stratification.order}{oder of entries in metadata.factor}

\item{custom.order}{custom sample order}

\item{last.plot.colors}{dataframe of plot colors for additional cross-lookup of taxa
information}
}
\description{
Generates barplots statified by taxon and metadata for one
HUMAnN2 feature, similar to `humann2_barplot` Optimized for viewability of
many features in one figure by limiting the number of statifications.
}
