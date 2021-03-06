% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/humann2_tools.R
\name{makeHumann2Barplot}
\alias{makeHumann2Barplot}
\title{makeHumann2Barplot}
\usage{
makeHumann2Barplot(
  dat,
  last.plot.colors,
  scale = "proportional-log",
  use.random.colors = T,
  hide.legend = T,
  space = "fixed",
  show.all.taxa = T,
  fixed.floor = NULL,
  fixed.ymax = NULL,
  sample.threshold = 10,
  hide.strata.legend = F
)
}
\arguments{
\item{dat}{table holding preprocessed humann2 information using `humann2Barplot`}

\item{last.plot.colors}{dataframe of plot colors}

\item{scale}{how to scale the height of bars, on default proportional-log}

\item{use.random.colors}{use randomcoloR instead of RColorBrewer}

\item{hide.legend}{boolean information if ledgend should be included}

\item{space}{free or fixed (x scale)}

\item{show.all.taxa}{boolean}

\item{fixed.floor}{set ymin to a fixed value to prevent focus on minor bugs}

\item{fixed.ymax}{set xmax to a fixed value keep y axis on multiple plots comparable}

\item{sample.threshold}{minimum number of samples per strata required for plotting}

\item{hide.strata.legend}{hide color legend}
}
\description{
Generates barplots statified by taxon and metadata
}
