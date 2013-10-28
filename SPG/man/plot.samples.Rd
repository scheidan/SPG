\name{plot.samples}
\alias{plot.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot function for \code{samples} objects
}
\description{
  Plot function for \code{samples} objects.
}
\usage{
\method{plot}{samples}(x, type='relative', title = NULL, ...)
}

\arguments{
  \item{x}{
    object of class \code{samples}
  }
  \item{type}{
    type of error to plot (unit on y-axis). Character string, either \code{"relative"},
    or \code{"absolute"} or \code{"factor"}.
}
  \item{title}{
    optional character vector for titel.
  }
  \item{\dots}{
    further arguments passed to \code{plot} see \code{\link{par}}.
  }
}

\details{
  If \code{type="factor"} boxplots for logarithm of sampled concentration divided by
  true concentration are produced. This representation has the advantage, that for 
  example an underestimation of 50\% generates the same visual impression as an 
  overestimation by 200\%. \code{type="relative"} should be used if a scale in \% is 
  desired and \code{type="absolute"} is used to visualize differences in 
}

\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
}

\seealso{
  \code{\link{SPG-package}}, \code{\link{def.grav}},  \code{\link{def.pump}}, \code{\link{take.samples}}, \code{\link{plot.flow}}
}