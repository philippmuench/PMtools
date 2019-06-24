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
library(PMtools)
data(humann2_table)
data(hmp1_2_metadata)
dat <- humann2Barplot(humann2_table, metadata = hmp1_2_metadata, feature = "Cas2")
p <- makeHumann2Barplot(dat)
pdf("plot.pdf", width = 5, height = 5)
print(p)
dev.off()
```

## License and copyright
Copyright 2019 Philipp Münch

Source code to PMtools is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). PMtools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

