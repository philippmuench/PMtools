[![Build Status](https://travis-ci.org/philippmuench/PMtools.svg?branch=master)](https://travis-ci.org/philippmuench/PMtools)
[![codecov](https://codecov.io/gh/philippmuench/PMtools/branch/master/graph/badge.svg)](https://codecov.io/gh/philippmuench/PMtools)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# PMtools

Miscellaneous helper functions mostly for comparative genomics, metagenomics and package development from P. Münch at LMU Munich & Helmholtz Centre for Infection Research.

```r
install.packages("devtools")
devtools::install_github("philippmuench/PMtools")
```

## usage

### HUMAnN2 tools

### generation of barplots

``` r
# load example datasets
data(humann2_table)
data(hmp1_2_metadata)
data(hmp1_2_metaphlan)
```

![](man/figures/README-example-1.png)

Command to generate a figure of a single feature. Set `num.bugs = "auto"` to automatically adjust the number of bugs needed to show 25% of RA

``` r
# generate the sample order
custom.order <-
  orderHumannBySimilarity(hmp1_2_metaphlan, distance.method = "bray")

# generate the data used for plotting
dat <-
  humann2Barplot(
    humann2_table,
    metadata = hmp1_2_metadata,
    feature = "Cas2",
    num.bugs = "auto",
    order.by = "custom",
    custom.order = custom.order
  )
  
# generate the plot
p <-
  makeHumann2Barplot(
    dat,
    NULL,
    hide.legend = F,
    scale = "pseudolog",
    space = "fixed"
  )

# show figure
print(p$gplot)
```

re-use colors

``` r
# generate the data used for plotting
dat_plot1 <-
  humann2Barplot(
    humann2_table,
    metadata = hmp1_2_metadata,
    feature = "Cas1",
    num.bugs = "auto",
    order.by = "custom",
    custom.order = custom.order
  )
  
# generate the data used for plotting
dat_plot2 <-
  humann2Barplot(
    humann2_table,
    metadata = hmp1_2_metadata,
    feature = "Cas2",
    num.bugs = "auto",
    order.by = "custom",
    custom.order = custom.order
  )
  
# generate the plot
p1 <-
  makeHumann2Barplot(
    dat_plot1,
    NULL,
    hide.legend = F,
    scale = "pseudolog",
    space = "fixed"
  )

# generate the plot
p2 <-
  makeHumann2Barplot(
    dat_plot2,
    p1$colors,
    hide.legend = F,
    scale = "pseudolog",
    space = "fixed"
  )
# same taxa now have the same color
print(p1$gplot)
print(p2$gplot)
```

you can also plot multiple features into one figure

``` r
cas_plots <- vector('list', 10)
plot.colors <- NULL
for (cas in paste0("Cas", 1:10)) {
  cas_plots[[cas]]  <- local({
    dat <-
      humann2Barplot(
        humann2_table,
        metadata = hmp1_2_metadata,
        feature = cas,
        num.bugs = "auto",
        num.bugs.explained.fraction = 0.35,
        order.by = "custom",
        custom.order = custom.order
      )
    p <-
      makeHumann2Barplot(
        dat,
        plot.colors,
        hide.legend = F,
        scale = "pseudolog",
        space = "fixed"
        )
    plot.colors <<- rbind(plot.colors, p$colors)
    write.table(p$colors, file="log.txt", append=T, sep ="\t", row.names=F,
    col.names=F, quote=F)
    print(p$gplot)
  })
}

# you can plot these together
pdf("all_cas.pdf", width = 6, height = 7)
print(multiplot(plotlist = cas_plots, cols = 2))
dev.off()

# or access each plot individually 
print(cas_plots[["Cas1"]])
# and the data underlying the plot
cas_plots[["Cas1"]]$dat
```


## License and copyright
Copyright 2019 Philipp Münch

Source code to PMtools is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). PMtools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

