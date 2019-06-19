#' @title humann2Barplot
#'
#' @description Generates barplots statified by taxon and metadata for one
#' HUMAnN2 feature, similar to `humann2_barplot` Optimized for viewability of
#' many features in one figure by limiting the number of statifications.
#'
#' @param humann2.table a dataframe containing the HUMAnN2 output in tsv format
#' @param num.bugs number of most abundant bugs to show as statification
#' @param feature name of gene in humann2_table that will be plotted
#' @param featue.column index of column in humann2 table that holds the feature
#' @param taxa.column index of column in humann2 table that holds the taxon infomation
#' @param metadata a dataframe containing the metadata
#' @param metadata.id column in metadata file that holds the smaple id (colum names in humann2 table)
#' @param metadata.factor column in metadata file that holds the annotation of interest
#' @param use.custom.column.stratification custom orderin of row statification
#' @param column.stratification.order oder of entries in metadata.factor
#'   information
#' @param order.by how to order bars within one statification
#' @export
humann2Barplot <- function(humann2.table,
                           num.bugs = 3,
                           feature = "Cas1",
                           featue.column = 1,
                           taxa.column = 2,
                           metadata,
                           metadata.id = 9,
                           metadata.factor = 4,
                           use.custom.column.stratification = T,
                           order.by = "bc",
                           column.stratification.order =  c("Oral",  "Skin", "Vaginal", "Gut")) {
  stopifnot(num.bugs > 0)
  stopifnot(any(humann2.table[, 1] == feature)) # feature is not in table
  # reduce table to relevant featue
  humann2.table <-
    humann2.table[which(humann2.table[, featue.column] == feature),]
  # get total abundance for feature for bugs
  humann2.table$abundance <-
    rowSums(humann2.table[, 3:ncol(humann2.table)])

  humann2.unclassified <-
    humann2.table[which(humann2.table[, taxa.column] == "unclassified"),]
  # replace name of "unclassified" to "other"
  humann2.unclassified$taxa <- "other"
  humann2.classified <-
    humann2.table[which(humann2.table[, taxa.column] != "unclassified"),]

  # get the top $bugs number of taxa in classified subset
  lst <-
    sort(humann2.classified$abundance,
         index.return = TRUE,
         decreasing = TRUE)
  top.index <-
    lapply(lst, "[", lst$x %in% head(unique(lst$x), num.bugs))

  # set taxa description of all non-top taxa to "other"
  if (nrow(humann2.classified) > num.bugs) {
    humann2.classified[-top.index$ix,][, taxa.column] <- "known"
  }

  humann.top.bugs <- rbind(humann2.unclassified, humann2.classified)
  humann.top.bugs$abundance <- NULL
  humann.top.bugs.m <- reshape2::melt(humann.top.bugs)
  # add metadata relevant for stratification
  humann.top.bugs.m$meta <-
    metadata[match(humann.top.bugs.m$variable, metadata[, metadata.id]), metadata.factor]

  if (use.custom.column.stratification)
    humann.top.bugs.m$meta <-
    factor(humann.top.bugs.m$meta, levels = column.stratification.order)

  if (order.by == "bc") {
    datalist <- list()
    for (meta in unique(humann.top.bugs.m$meta)) {
      if (is.na(meta))
        next

      meta.samples <-
        metadata[which(metadata[, metadata.factor] == meta), metadata.id]
      meta.community <-
        humann2.table[which(humann2.table[, featue.column] == feature), ]
      rownames(meta.community) <- meta.community[, taxa.column]
      # limit to e.g. one body site
      meta.community <-
        meta.community[, which(names(meta.community) %in% meta.samples)]

      meta.community <- as.matrix(t(meta.community))
      bc <-
        as.matrix(vegan::vegdist(meta.community, method = "bray"))
      bc[is.na(bc)] <- 0
      bc.clusters <- stats::hclust(as.dist(bc), method = "single")
      bc.order.index <-
        stats::order.dendrogram(stats::as.dendrogram(bc.clusters))

      datalist[[length(datalist) + 1]] <-
        data.frame(
          feature = feature,
          meta = meta,
          samples = rownames(meta.community)[bc.order.index]
        )
    }
    humann.top.bugs.bc <- do.call(rbind, datalist)
    # change the facort order
    humann.top.bugs.m$variable <- factor(humann.top.bugs.m$variable,
                                         levels = humann.top.bugs.bc$samples)
    message("Finished sorting by BC.")
  }
  return(humann.top.bugs.m)
}

#' @title makeHumann2Barplot
#'
#' @description Generates barplots statified by taxon and metadata
#' @param dat table holding preprocessed humann2 information using `humann2Barplot`
#' @param scale how to scale the height of bars, on default sqrt
#' @export
makeHumann2Barplot <-
  function(dat,
           scale = "sqrt",
           palette = "Set1",
           hide.legend = T) {
    y_limit <-
      max(stats::aggregate(value ~ variable, data = dat, FUN = sum)$value)
    p <-
      ggplot2::ggplot(dat = dat, ggplot2::aes(x = variable, y = value, fill = taxa))
    p <- p + ggplot2::geom_bar(stat = "identity")
    p <- p + ggplot2::ggtitle(dat[1, 1])
    if (scale == "sqrt")
      p <-
      p + ggplot2::scale_y_sqrt(expand = c(0, 0), limits = c(0, y_limit))
    p <-
      p +  ggplot2::facet_grid(. ~ meta, space = "free_x", scales = "free_x")
    p <-
      p +  PMtools::themePM() + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    # make facet_grid nicer
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "black"))
    p <- p + ggplot2::scale_fill_brewer(palette = palette)
    if (hide.legend){
      p <- p + ggplot2::theme(legend.position = "none")
    } else {
      p <- p + ggplot2::theme(legend.position = "bottom")
    }
    return(p)
  }