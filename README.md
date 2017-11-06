# Sewage pattern generator (SPG)

R-package that provides functions to generate flow and substance patterns
  in complex sewer systems.

The main aim is to efficiently model realistic short-term variations
  of flows and substances in sewers and to assess their effect on
  sampling.  In addition to complex sewer networks, the package is
  also particularly useful to optimize sampling at the effluent of
  individual premises (e.g. hospital, school, prison) where short-term
  variations are typically highest and challenge the reliable
  collection of representative average samples over a period of time.

  Different scenarios (e.g. distribution of substance across
  sub-catchments) and sampling setups can be evaluated at relatively
  low compuational cost. The processes advection and dispersion are
  considered.

## Installation

1. Install [R](https://cloud.r-project.org/) and [R-Studio](https://www.rstudio.com/products/RStudio/) or any other editor.

2. Install `devtools` (type in the R command line)
```
install.packages("devtools")
```

3. Install `SPG` (type in the R command line)
```
library(devtools)
install_github("scheidan/SPG/SPG")
