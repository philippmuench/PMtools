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
    lapply(lst, "[", lst$x %in% utils::head(unique(lst$x), num.bugs))

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
      print(meta)
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
      meta.community <- as.data.frame(meta.community)
      #meta.community <- as.matrix(t(meta.community))
      meta.community.t <- data.table::transpose(meta.community)
      colnames(meta.community.t) <- rownames(meta.community)
      rownames(meta.community.t) <- colnames(meta.community)
      zero.count.ids <- which(rowSums(meta.community.t) == 0)
      pos.count.ids <- which(rowSums(meta.community.t) != 0)
      if (length(zero.count.ids) >= nrow(meta.community.t)) next # all are zero
      if (length(zero.count.ids) > 0) {
         meta.community.t <- meta.community.t[pos.count.ids,]
         message(paste("removed", length(zero.count.ids), "samples with empty rows"))
      }
      if (nrow(meta.community.t) < 3) {
        # pseudo sort since clustering would require more samples
        datalist[[length(datalist) + 1]] <-
          data.frame(
            feature = feature,
            meta = meta,
            samples = rownames(meta.community.t)
          )
        next
      }
      bc <-
        as.matrix(vegan::vegdist(meta.community.t, method = "bray"))
      bc[is.na(bc)] <- 0
      bc.clusters <- stats::hclust(stats::as.dist(bc), method = "single")
      bc.order.index <-
        stats::order.dendrogram(stats::as.dendrogram(bc.clusters))
      datalist[[length(datalist) + 1]] <-
        data.frame(
          feature = feature,
          meta = meta,
          samples = rownames(meta.community.t)[bc.order.index]
        )
    }
    humann.top.bugs.bc <- do.call(rbind, datalist)
    # change the facort order
    humann.top.bugs.m$variable <- factor(humann.top.bugs.m$variable,
                                         levels = humann.top.bugs.bc$samples)
    message("Finished sorting by BC.")
  }
  humann.top.bugs.m <- humann.top.bugs.m[which(humann.top.bugs.m$value != 0),]
  # sum up all known taxa per stratum
  humann.top.bugs.m.agg <- stats::aggregate(value ~ SRS + taxa + variable + meta, data = humann.top.bugs.m, FUN = sum)
  return(humann.top.bugs.m.agg)
}

#' @title makeHumann2Barplot
#'
#' @description Generates barplots statified by taxon and metadata
#' @param dat table holding preprocessed humann2 information using `humann2Barplot`
#' @param scale how to scale the height of bars, on default sqrt
#' @param bugs.colors html color codes of bar plots, must be same length as num.bugs
#' @param hide.legend boolean information if ledgend should be included
#' @export
makeHumann2Barplot <-
  function(dat,
           scale = "sqrt",
           bugs.colors = c("#1b9e77","#d95f02", "#7570b3"),
           hide.legend = T) {
    if (scale == "log10+1") {
      dat$value <- log10(dat$value + 1)
    }
    if (scale == "pseudolog") {
      dat$value <- pseudoLog10(dat$value)
    }
    p <-
      ggplot2::ggplot(dat = dat, ggplot2::aes(x = variable, y = value, fill = taxa))
    p <- p + ggplot2::geom_bar(stat = "identity")
    p <- p + ggplot2::ggtitle(dat[1, 1])
    if (scale == "sqrt") {
      p <- p + ggplot2::scale_y_sqrt(expand = c(0, 0))
      p <- p + ggplot2::ylab("abundance (sqrt)")
    } else if (scale == "pseudolog") {
      p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
   #   p <- p +  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), expand = c(0, 0))
      p <- p + ggplot2::ylab("abundance (pseudolog)")
    } else if (scale == "log10+1"){
      p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
      p <- p + ggplot2::ylab("abundance (log10+1)")
    } else {
      p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
      p <- p + ggplot2::ylab("abundance")
    }

    p <- p +  PMtools::themePM()
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank())
    p <-
      p +  ggplot2::facet_grid(. ~ meta, space = "free_x", scales = "free_x")
    if (hide.legend) {
      p <- p + ggplot2::theme(legend.position = "none")
    } else {
      p <- p + ggplot2::theme(legend.position = "bottom")
    }
    # coloring
    if (length(bugs.colors) < length(unique(dat$taxa)) - 2) {
      message("Not enough colors provided, using RColorBrewer")
      bugs.colors <- RColorBrewer::brewer.pal(length(unique(dat$taxa)) - 2, "Set1")
    }
    p <- p + ggplot2::scale_fill_manual(values =  c(bugs.colors, "grey80", "grey60"))
    p <- p + ggplot2::guides(fill = ggplot2::guide_legend(title = "", ncol = 2))
    return(p)

  }
