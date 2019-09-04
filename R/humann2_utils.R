
#' @title orderHumannBySimilarity
#'
#' @description outputs order of samples by bray curtis
#' @param metaphlan MetaPhlAn2 table
#' @param distance.method Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
#' @param cluster.method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @export
orderHumannBySimilarity <-
  function(metaphlan,
           distance.method = "bray",
           cluster.method = "single") {
    # generate community matrix

    meta.community <- as.data.frame(metaphlan)
    rownames(meta.community) <- meta.community$taxa
    meta.community$taxa <- NULL

    meta.community.t <- data.table::transpose(meta.community)
    colnames(meta.community.t) <- meta.community$taxa
    rownames(meta.community.t) <- colnames(meta.community)

    bc <-
      as.matrix(vegan::vegdist(meta.community.t, method = distance.method))

    bc[is.na(bc)] <- 0
    bc.clusters <-
      stats::hclust(stats::as.dist(bc), method = cluster.method)
    bc.order.index <-
      stats::order.dendrogram(stats::as.dendrogram(bc.clusters))

    ordering <- rownames(meta.community.t)[bc.order.index]
    return(ordering)
  }


#' @title plotHummanMetadata
#'
#' @description outputs order of samples by bray curtis
#' @param metaphlan MetaPhlAn2 table
#' @param distance.method Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
#' @param cluster.method The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @export
plotHummanMetadata <-
  function(metaphlan,
           distance.method = "bray",
           cluster.method = "single",
           metadata = hmp1_2_metadata) {
    # generate community matrix
    meta.community <- as.data.frame(metaphlan)
    rownames(meta.community) <- meta.community$taxa
    meta.community$taxa <- NULL

    meta.community.t <- data.table::transpose(meta.community)
    colnames(meta.community.t) <- meta.community$taxa
    rownames(meta.community.t) <- colnames(meta.community)

    bc <-
      as.matrix(vegan::vegdist(meta.community.t, method = distance.method))

    bc[is.na(bc)] <- 0
    bc.clusters <-
      stats::hclust(stats::as.dist(bc), method = cluster.method)
    bc.order.index <-
      stats::order.dendrogram(stats::as.dendrogram(bc.clusters))

    ordering <- rownames(meta.community.t)[bc.order.index]
    area <- metadata[match(ordering, metadata$SRS),]$STSite
    df <- data.frame(site = ordering, area = area)
    return(df)
  }

humann2BarplotDummy <- function(custom.order.area,
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

  humann2.table <- humann2.table[1,]
  humann2.table[,c(3:ncol(humann2.table))] <- 1

  humann2.classified <-
    humann2.table[which(humann2.table[, taxa.column] != "unclassified"), ]

  # get the top $bugs number of taxa in classified subset
  lst <-
    sort(humann2.classified$abundance,
         index.return = TRUE,
         decreasing = TRUE)

  top.index <-
    lapply(lst, "[", lst$x %in% utils::head(unique(lst$x), 1))

  humann.top.bugs <- humann2.classified
  humann.top.bugs$abundance <- NULL
  humann.top.bugs.m <- reshape2::melt(humann.top.bugs)
  # add metadata relevant for stratification
  humann.top.bugs.m$meta <-
    metadata[match(humann.top.bugs.m$variable, metadata[, metadata.id]), metadata.factor]

  humann.top.bugs.m$meta <-
    factor(humann.top.bugs.m$meta, levels = column.stratification.order)

  # change the facort order
  humann.top.bugs.m$variable <- factor(humann.top.bugs.m$variable,
                                       levels = custom.order)
  humann.top.bugs.m$meta2 <- custom.order.area[match(humann.top.bugs.m$variable, custom.order.area$site),]$area


  # sum up all known taxa per stratum
  humann.top.bugs.m.agg <-
    stats::aggregate(value ~ SRS + taxa + variable + meta + meta2,
                     data = humann.top.bugs.m,
                     FUN = sum)

  return(humann.top.bugs.m.agg)
}


makeHumann2BarplotDummy <-
  function(dat,
           last.plot.colors = NULL,
           scale = "sqrt",
           use.random.colors = T,
           hide.legend = T,
           space = "free") {
    unclassified.name <- "Unclassified"
    other.name <- "Other"
    message(paste(nrow(last.plot.colors$colors), "colors provided"))
    # get taxon names for coloring
    taxon.names <- unique(dat$taxa)
    if (length(grep(other.name, taxon.names)) > 0)
      taxon.names <- taxon.names[which(taxon.names != other.name)]
    if (length(grep(unclassified.name, taxon.names)) > 0)
      taxon.names <-
      taxon.names[which(taxon.names != unclassified.name)]

    if (scale == "log10+1") {
      dat$value <- log10(dat$value + 1)
      message("use log10(dat$value + 1) for scaling")
    }
    if (scale == "pseudolog") {
      message("use PMtools::pseudoLog10 for scaling")
      dat$value <- PMtools::pseudoLog10(dat$value)
    }

    order.levels <- c(taxon.names, other.name, unclassified.name)
    p <-
      ggplot2::ggplot(dat = dat, ggplot2::aes(
        x = variable,
        y = value,
        fill = factor(meta2)
      ))
    p <- p + ggplot2::geom_bar(stat = "identity")
    # p <- p + ggplot2::ggtitle(dat[1, 1])
    if (scale == "sqrt") {
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
    } else {
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

    if (scale == "ggplot2") {
      message("use ggplot2::scale_y_log10() scaling")
      p <- p + ggplot2::scale_y_log10()
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

    # remove facet_grid legend
    p <- p + ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )
    # legend size
    p <- p + ggplot2::scale_size(range = c(5, 20), guide = "none")

    annot_col = list(area = c("Gut" = "brown", "Skin" = "purple", "Oral" = "blue", "Vaginal" = "Orange"),
                     site = c("Buccal_mucosa" = "darkblue",
                              "Hard_palate" = "lightcyan",
                              "Keratinized_gingiva" = "mediumblue",
                              "Palatine_Tonsils" = "darkorchid",
                              "Saliva" = "cyan",
                              "Subgingival_plaque" = "dodgerblue",
                              "Supragingival_plaque" = "grey20",
                              "Throat" = "cornflowerblue",
                              "Tongue_dorsum" = "grey50",
                              "Stool" = "saddlebrown",
                              "Anterior_nares" = "orange",
                              "L_Retroauricular_crease" = "red",
                              "R_Antecubital_fossa" = "yellow",
                              "R_Retroauricular_crease" = "tomato",
                              "Mid_vagina" = "green",
                              "Posterior_fornix" = "seagreen",
                              "Vaginal_introitus" = "olivedrab"))

    p <- p + ggplot2::scale_fill_manual(values = annot_col$site)
    # reduce legend point size
    p <-
      p + ggplot2::theme(legend.key.size = ggplot2::unit(0.2, "line"))
    return(list("gplot" = p, "colors" = colors.df))
  }