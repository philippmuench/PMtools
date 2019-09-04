#' @title humann2Barplot
#'
#' @description Generates barplots statified by taxon and metadata for one
#' HUMAnN2 feature, similar to `humann2_barplot` Optimized for viewability of
#' many features in one figure by limiting the number of statifications.
#'
#' @param humann2.table a dataframe containing the HUMAnN2 output in tsv format
#' @param num.bugs number of most abundant bugs to show as statification, set to
#' to "auto" if this should be ajusted based on num.bugs.explained.fraction
#' @param num.bugs.explained.fraction fraction that needs to be explained by the
#' number of bugs shown
#' @param feature name of gene in humann2_table that will be plotted
#' @param featue.column index of column in humann2 table that holds the feature
#' @param taxa.column index of column in humann2 table that holds the taxon infomation
#' @param metadata a dataframe containing the metadata
#' @param metadata.id column in metadata file that holds the smaple id (colum names in humann2 table)
#' @param metadata.factor column in metadata file that holds the annotation of interest
#' @param use.custom.column.stratification custom orderin of row statification
#' @param order.by how to order bars within one statification
#' @param column.stratification.order oder of entries in metadata.factor
#'   information
#' @param custom.order custom sample order
#' @export
humann2Barplot <- function(humann2.table,
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
                           column.stratification.order =  c("Oral",  "Skin", "Vaginal", "Gut"),
                           custom.order) {
  stopifnot(num.bugs > 0)
  stopifnot(any(humann2.table[, 1] == feature)) # feature is not in table

  # reduce table to relevant featue
  humann2.table <-
    humann2.table[which(humann2.table[, featue.column] == feature), ]
  # get total abundance for feature for bugs
  humann2.table$abundance <-
    rowSums(humann2.table[, 3:ncol(humann2.table)])

  humann2.unclassified <-
    humann2.table[which(humann2.table[, taxa.column] == "unclassified"), ]
  # replace name of "unclassified" to "other"
  humann2.unclassified$taxa <- "Unclassified"
  humann2.classified <-
    humann2.table[which(humann2.table[, taxa.column] != "unclassified"), ]

  # get the top $bugs number of taxa in classified subset
  lst <-
    sort(humann2.classified$abundance,
         index.return = TRUE,
         decreasing = TRUE)

  # if num.bugs set to "auto", get the number of bugs that explain 50% of variance
  if (num.bugs == "auto") {
    index <-
      lapply(lst, "[", lst$x %in% unique(lst$x))
    sorted <- humann2.classified[index$ix, ]
    sorted$cum_abundance <- cumsum(sorted$abundance)
    breakpoint <-
      sum(sorted$abundance) * num.bugs.explained.fraction
    num.bugs <- min(which(sorted$cum_abundance > breakpoint))
  }

  top.index <-
    lapply(lst, "[", lst$x %in% utils::head(unique(lst$x), num.bugs))

  # set taxa description of all non-top taxa to "known"
  if (nrow(humann2.classified) > num.bugs) {
    humann2.classified[-top.index$ix, ][, taxa.column] <- "Other"
  }

  # shorten taxa name
  taxa.names <- humann2.classified[top.index$ix, ][, taxa.column]
  humann2.classified[top.index$ix, ][, taxa.column] <-
    PMtools::shortenTaxons(taxa.names)

  humann.top.bugs <- rbind(humann2.unclassified, humann2.classified)
  humann.top.bugs$abundance <- NULL
  humann.top.bugs.m <- reshape2::melt(humann.top.bugs)
  # add metadata relevant for stratification
  humann.top.bugs.m$meta <-
    metadata[match(humann.top.bugs.m$variable, metadata[, metadata.id]), metadata.factor]

  if (use.custom.column.stratification)
    humann.top.bugs.m$meta <-
    factor(humann.top.bugs.m$meta, levels = column.stratification.order)

  if (order.by == "custom") {
    # change the facort order
    humann.top.bugs.m$variable <- factor(humann.top.bugs.m$variable,
                                         levels = custom.order)
    message("Finished sorting by custom order.")
  }
  if (order.by == "bc") {
    datalist <- list()
    for (meta in unique(humann.top.bugs.m$meta)) {
      print(meta)
      if (is.na(meta))
        next
      meta.samples <-
        metadata[which(metadata[, metadata.factor] == meta), metadata.id]
      meta.community <-
        humann2.table[which(humann2.table[, featue.column] == feature),]
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
      if (length(zero.count.ids) >= nrow(meta.community.t))
        next # all are zero
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
      bc.clusters <-
        stats::hclust(stats::as.dist(bc), method = "single")
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
  #ordering <-  humann.top.bugs.bc$samples
  #humann.top.bugs.m <- humann.top.bugs.m[which(humann.top.bugs.m$value != 0),]

  # sum up all known taxa per stratum
  humann.top.bugs.m.agg <-
    stats::aggregate(value ~ SRS + taxa + variable + meta,
                     data = humann.top.bugs.m,
                     FUN = sum)

  # if we have no other category, add dummy
  if (length(which(humann.top.bugs.m.agg$taxa == "Other")) == 0) {
    dummy <-
      data.frame(
        SRS = feature,
        taxa = "Other",
        variable = humann.top.bugs.m.agg$variable[1],
        meta = unique(humann.top.bugs.m$meta)[1],
        value = 0
      )
    humann.top.bugs.m.agg <- rbind(humann.top.bugs.m.agg, dummy)
  }
  return(humann.top.bugs.m.agg)
}

#' @title makeHumann2Barplot
#'
#' @description Generates barplots statified by taxon and metadata
#' @param dat table holding preprocessed humann2 information using `humann2Barplot`
#' @param last.plot.colors dataframe of plot colors
#' @param scale how to scale the height of bars, on default proportional-log
#' @param use.random.colors use randomcoloR instead of RColorBrewer
#' @param hide.legend boolean information if ledgend should be included
#' @param space free or fixed (x scale)
#' @export
makeHumann2Barplot <-
  function(dat,
           last.plot.colors,
           scale = "proportional-log",
           use.random.colors = T,
           hide.legend = T,
           space = "free") {
    unclassified.name <- "Unclassified"
    other.name <- "Other"
    # get taxon names for coloring
    taxon.names <- unique(dat$taxa)
    if (length(grep(other.name, taxon.names)) > 0)
      taxon.names <- taxon.names[which(taxon.names != other.name)]
    if (length(grep(unclassified.name, taxon.names)) > 0)
      taxon.names <-
      taxon.names[which(taxon.names != unclassified.name)]

    # scaling value transformation
    if (scale == "log10+1") {
      message("using log10(dat$value + 1)")
      dat$value <- log10(dat$value + 1)
    }
    if (scale == "pseudolog") {
      message("using PMtools::pseudoLog10(dat$value)")
      dat$value <- PMtools::pseudoLog10(dat$value)
    }
    if (scale == "proportional-log") {
      # remove stratification
      rescale <- TRUE
      if (rescale) {
        c_epsilon <- 1e-10

        dat$match <- paste0(dat$variable, dat$meta)

        # aggregate by sample
        dat.no.taxon.strata <-  stats::aggregate(value ~ SRS + variable + meta + match,
                                                 data = dat,
                                                 FUN = sum)

        table.colsums <-  dat.no.taxon.strata$value
        ymin <- min(table.colsums[which(table.colsums > 0)])
        floor <- floor(log10(ymin))

        if (log10(ymin) - floor < c_epsilon) {
          floor <- floor - 1
        }
        floors <-  rep(floor , length(table.colsums))

        dat.no.taxon.strata$crests <- dat.no.taxon.strata$value

        # log10 values that are > 0.01
        dat.no.taxon.strata$crests[which(dat.no.taxon.strata$crests > 10**floor)] <- log10(dat.no.taxon.strata$crests[which(dat.no.taxon.strata$crests > 10**floor)])
        # floor values that are < 0.01
        dat.no.taxon.strata$crests[which( dat.no.taxon.strata$crests <= 10**floor)] <- floor

        dat.no.taxon.strata$heights <- dat.no.taxon.strata$crests - floor
        dat$agg <- dat.no.taxon.strata[match(dat$match, dat.no.taxon.strata$match),]$value
        dat$agg_heights <- dat.no.taxon.strata[match(dat$match, dat.no.taxon.strata$match),]$heights

        dat$value <- (dat$value / dat$agg) * dat$agg_heights
        ymax <- ceiling(log10(max(dat[which(dat$agg > 0),]$agg)))

      } else {
        dat$match <- paste0(dat$variable, dat$meta)

        # aggregate by sample
        dat.no.taxon.strata <-  stats::aggregate(value ~ SRS + variable + meta + match,
                                                 data = dat,
                                                 FUN = sum)
        dat.no.taxon.strata$log10 <- log10(dat.no.taxon.strata$value)
        ## do log10 on non-stratified table
        #dat.no.taxon.strata <- dat.no.taxon.strata[which(dat.no.taxon.strata$value > 0),]

        dat$agg <- dat.no.taxon.strata[match(dat$match, dat.no.taxon.strata$match),]$value
        dat$agg_log <- dat.no.taxon.strata[match(dat$match, dat.no.taxon.strata$match),]$log10

        # calculate back the proportion
        dat$value <- (dat$value / dat$agg) * dat$agg_log

      }

    }

    order.levels <- c(taxon.names, other.name, unclassified.name)
    p <-
      ggplot2::ggplot(dat = dat, ggplot2::aes(
        x = variable,
        y = value,
        fill = factor(taxa, levels = order.levels)
      ))
    p <- p + ggplot2::geom_bar(stat = "identity")

    # scaling axis manipulation
    if (scale == "sqrt") {
      message("using ggplot2::scale_y_sqrt()")
      p <- p + ggplot2::scale_y_sqrt(
        expand = c(0, 0),
        breaks = function(x)
          unique(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))))
      )

      p <- p + ggplot2::ylab("abundance (sqrt)")
    } else if (scale == "pseudolog") {
      p <- p + ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = function(x)
          unique(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))))
      )
      #   p <- p +  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), expand = c(0, 0))
      p <- p + ggplot2::ylab(dat[1, 1])
    } else if (scale == "log10+1") {
      p <- p + ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = function(x)
          unique(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))))
      )
      p <- p + ggplot2::ylab("abundance (log10+1)")
    } else if (scale == "none") {
      message("using no scaling")
      p <- p + ggplot2::ylab("abundance (no scaling)")
    } else if (scale == "proportional-log") {
      message("using propotional log")
      if (rescale) {
        print(floor)

        breaks_fun <-  function(x){
          unique(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))))
        }

        labels_fun <-  function(x){
          print(x + floor)
        }

        p <- p + ggplot2::scale_y_continuous(
          expand = c(0, 0),
          labels = labels_fun,
          breaks = breaks_fun)


        #  # Change y tick mark labels
        #   p <- p + ggplot2::scale_y_discrete(breaks = c(0:ymax),
        #                       labels = c(rep("aaa",length(0:ymax) )))
      }
      p <- p + ggplot2::ylab("log10")
    } else if (scale == "ggplot2-log10") {
      message("using ggplot2 log10 axis scaling")
      p <- p + ggplot2::scale_y_log10()
      p <- p + ggplot2::ylab("abundance (ggplot2 scaling)")
    } else {
      message("using no scaling but adjust breaks that it looks pretty")
      p <- p + ggplot2::scale_y_continuous(
        expand = c(0, 0),
        breaks = function(x)
          unique(floor(pretty(seq(
            0, (max(x) + 1) * 1.1
          ))))
      )
      p <- p + ggplot2::ylab("abundance")
    }


    p <- p + PMtools::themePM(base.size = 7, axis.family = "mono")
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )
    if (space == "free") {
      p <-
        p +  ggplot2::facet_grid(. ~ meta, space = "free_x", scales = "free_x")
    } else {
      p <-
        p +  ggplot2::facet_grid(. ~ meta,
                                 space = "free_x",
                                 scales = "free_x",
                                 shrink = T)
    }

    if (hide.legend) {
      p <- p + ggplot2::theme(legend.position = "none")
    } else {
      p <- p + ggplot2::theme(legend.position =  c(0.5, 1))
    }

    # coloring
    taxa_list <- unique(dat$taxa)
    taxa_list <- taxa_list[taxa_list != other.name &
                             taxa_list != unclassified.name]
    if (use.random.colors) {
      colors.df <-
        data.frame(
          taxa = taxa_list,
          color = randomcoloR::randomColor(count = length(taxa_list)),
          stringsAsFactors = F
        )
    } else {
      colors.df <-
        data.frame(taxa = taxa_list,
                   color = RColorBrewer::brewer.pal(taxa_list, "Set1"))
    }

    # replace colors if needed
    if (!missing(last.plot.colors) & !is.null(last.plot.colors)) {
      colors.df$new <- last.plot.colors[match(colors.df$taxa, last.plot.colors$taxa),]$color
      colors.df[which(is.na(colors.df$new)),]$new  <- colors.df[which(is.na(colors.df$new)),]$color
      colors.df$color <- NULL
      colnames(colors.df)[2] <- "color"
    }

    # add other and unclassified to colors
    colors.df.extended <-
      rbind(colors.df, data.frame(
        taxa = c(other.name, unclassified.name),
        color = c("grey80", "grey60")
      ))
    #write.table(file =  "log.txt", colors.df.extended, sep = "\t", quote = F, append = T)
    p <-
      p + ggplot2::scale_fill_manual(values =  colors.df.extended$color, breaks = colors.df.extended$taxa)
    p <-
      p + ggplot2::guides(fill = ggplot2::guide_legend(
        title = "",
        ncol = length(colors.df.extended$taxa)
      ))

    # remove facet_grid legend
    p <- p + ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )
    # legend size
    p <- p + ggplot2::scale_size(range = c(5, 20), guide = "none")

    # reduce legend point size
    p <-
      p + ggplot2::theme(legend.key.size = ggplot2::unit(0.2, "line"))
    return(list("gplot" = p, "colors" = colors.df))
  }