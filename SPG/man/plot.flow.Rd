\name{plot.flow}
\alias{plot.flow}
\title{
  Plot function for \code{flow} objects
}
\description{
  This function allows you to plot the entire modeled pattern or sections thereof.
}
\usage{
\method{plot}{flow}(x, days=NULL, from=0, title=NULL, ...)
}

\arguments{
  \item{x}{
    object of class \code{flow}
  }
  \item{days}{
    number of days to be plotted. If \code{NULL}
    all days are plotted.
  }
  \item{from}{
    start point in days.
  }
  \item{title}{
    character vector for main title.
  }
  \item{\dots}{
    further arguments passed to \code{plot} see \code{\link{par}}.
  }
}

\details{
  Substance loads, wastewater flows and concentration patterns are plotted in three pannels.
  Weekends are marked with a green bar, weekdays are labelled blue. 
  
  For an example see the documentation of \code{\link{SPG-package}}.
}

\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
}


\seealso{
  \code{\link{SPG-package}}, \code{\link{def.grav}},  \code{\link{def.pump}}, \code{\link{take.samples}}, \code{\link{plot.samples}}
}
