\name{take.samples}
\alias{take.samples}

\title{
  Take composite samples
}
\description{
  Time-, volume- and flow-proportional sampling modes are implemented to simulate the 
  collection of a composite sample over a defined period of time. Average concentrations
  from these composite samples are compared to the "true" average concentration calculated 
  from the highest temporal resolution of the modeled time series.
}
\usage{
take.samples(flow, composite.duration = 24, sampling.interval = NULL,
  v.prop.vol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flow}{
    an object of class \code{flow} calculated based on the functions \code{\link{def.grav}} 
	and \code{\link{def.pump}}; for an example see \code{\link{SPG-package}}.
  }
  \item{composite.duration}{
    duration of one composite sample [h]
  }
  \item{sampling.interval}{
    time between samples for time- and flow-proportional sampling mode [min].
	If \code{NULL} the sampling interval is calculated to collect approximately 
	the same amount of samples per composite sample as with the volume-proportional 
	sampling mode. Note: Either \code{sampling.interval} or \code{v.prop.vol} (or both)
	must be specified.
  }
  \item{v.prop.vol}{
    volume increment to trigger sampling in the volume-proportional sampling 
	mode [m^3]. 
	If \code{NULL} the increment is calculated to approximately take 
	the same amount of samples as with the time- and flow-proportional sampling
	mode. Note: Either \code{sampling.interval} or \code{v.prop.vol} (or both)
	must be specified.
  }
}

\details{
  Modeled time series of wastewater flows and mass fluxes of the substance
  of interest at high temporal resolution are the basis to assess the performance of 
  different sampling modes and intervals. Typically, 24-hour composite samples are 
  collected to obtain daily average concentrations. The following variables must be 
  specified: i) the duration over which a representative average concentration should be 
  obtained (e.g. 24h), ii) a sampling interval in minutes (for time- and flow-proportional 
  sampling mode, usually between 5min to 1h) and/or iii) a volume increment in cubicmeters
  (for the volume-proportional sampling mode) must be specified. All three sampling modes 
  are calculated simultaneously to facilitate comparison. If either \code{sampling.interval} or 
  \code{v.prop.vol} are not specified the missing argument will be calculated such that the number 
  of samples to form the composite samples is approximately the same. If both \code{sampling.interval} 
  and \code{v.prop.vol} are specified, average sampling may differ. In any case, the intervals are 
  displayed when plotting the result. Typically, the samples based on longer sampling 
  intervals will lead to larger relative errors.
    

  For an example see the documentation of \code{\link{SPG-package}}.
}

\value{
  A list of class \code{samples} with the following components:
  \item{C.true}{vector containg the "true" average concentrations based on 
    the data series highest temporal resolution for the durations of each 
	composite sample}
  \item{C.time.prop}{vector containing the concentrations obtained with
    time-proportional sampling}
  \item{C.Vol.prop}{vector containing the concentrations obtained with
    volume-proportional sampling}
  \item{C.flow.prop}{vector containing the concentrations obtained with
    flow-proportional sampling}
  \item{Vol}{vector containing the water volumes for each period [m^3]}

}

\references{
  Ort C., Lawrence M.G., Rieckermann J. and Joss A. (2010) Monitoring Pharmaceuticals
  and Personal Care Products (PPCPs) and Illicit Drugs in Wastewater Systems: Are Your 
  Conclusions Valid? A Critical Review. Environmental Science & Technology 44 (16), 6024-6035. 
}
\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
}

\seealso{
  \code{\link{SPG-package}}, \code{\link{def.grav}},  \code{\link{def.pump}}, \code{\link{plot.flow}}, \code{\link{plot.samples}}
}
