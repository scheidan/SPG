\name{def.pump}
\alias{def.pump}

\title{
  Define a pump
}
\description{
  This function allows you to define a pump that is operated in an intermittent mode of operation. 
}

\usage{
def.pump(V.max, V.min = 0, pump.rate, distance = 0, v.flow = 1, Disp = 0.16)
}

\arguments{
  \item{V.max}{
    volume of the pump sump [m^3], i.e. the volume at the trigger level at which the 
	pump starts emptying the pump sump.
  }
  \item{V.min}{
    minimum volume of water that remains in the pump sump when the pump turns off [m^3]
  }
  \item{pump.rate}{
    pump rate at which the pump sump is emptied [l/s]
  }
  \item{distance}{
    distance of pump to sampling point or until it joins the next gravity system [m/s]
  }
  \item{v.flow}{
    flow velocity in the sewer (open channel gravity flow) after the pump [m/s] 
  }
  \item{Disp}{
    dispersion coefficient in the gravity sewer after the pump [m^2/s]
  }
}

\details{
  This function computes the flow and substance patterns after a pump, and accounts for the processes 
  advection and dispersion in the subsequent gravity sewer pipe until this sewer joins another gravity
  system or the sampling point.
  
  The wastewater flows Q defined in \code{def.grav} are the driving factor for the 
  operation of the pumps. The pump sump is filled until the usual storage capacity 
  is reached and then the pump starts emptying the at the defined rate \code{pump.rate}.  
  
  The distance that needs to be specified is the length of the sewer pipe after the pump until the
  wastewater is collected either in the next pump sump or until the wastewater reaches the sampling 
  point. This also must be specified if the pump is located within another gravity system (see graph
  in VIGNETTE.PDF (provide link).

  For definitions and values of velocity \code{v.flow} and dispersion \code{Disp} see \code{\link{def.grav}}   

  For an example see documentation of \code{\link{SPG-package}}.

}

\value{ A function that computes the flow and substance patterns after a
  pump as well as the dispersion in the consequent sewer. It expects a
  \code{flow} object as argument.


}

\references{
  Ort et al. SEWPROF D1.1 (in prep., this software tool)
}
\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
}

\seealso{

  \code{\link{SPG-package}}, \code{\link{def.grav}}, \code{\link{take.samples}}, \code{\link{plot.flow}}, \code{\link{plot.samples}}
}
% \examples{
% ## example here
% }
