\name{SPG-internals}
\alias{disp.trans}
\alias{gen.Q}
\alias{gen.S}
\alias{print.grav}
\alias{print.pump}
\alias{print.samples}
\alias{sampling}
\alias{pump.trans}
\alias{+.flow}

\title{
  Internal functions
}

\description{
  This functions are only internally used and are not visible for the user.
}

\usage{

gen.Q(t.sim, population, temp.res.sim, pattern = 1, Q.c.day = 250,
    Exp.pulses = 100, flow.distances, pulse.dur, v.flow = 1,
    Q.base = 0, Disp = 0.16)

gen.S(t.sim, exp.n.pulses.WD, pattern.WD = 1,
    exp.n.pulses.WE, pattern.WE,
    pulse.masses, pulse.dur.min, pulse.dur.max, temp.res.sim,
    v.flow, Disp, flow.distances)

sampling(flow, composite.duration, sampling.interval, v.prop.vol = NULL)

pump.trans(flow, V.max, V.min = 0, pump.rate)

disp.trans(flow, distance, v.flow = 1, Disp = 0.16)


\method{+}{flow}(e1, e2 = NULL)

\method{print}{grav}(x, ...)
\method{print}{pump}(x, ...)
\method{print}{samples}(x, ...)

}

\arguments{
  \item{t.sim}{
    length of imulated patters [sec]
}
  \item{population}{
    population in catchment
}
  \item{temp.res.sim}{
    temporal resolution of generated pattern [sec]
}
  \item{pattern}{
    vector with rate of diurnal drug excretion pattern [-].
}
  \item{Q.c.day}{
    daily water consumption per day per person [l/c/day]
  }
  \item{pulse.dur}{
    vector of flow distances to sample from [m]
    }
 \item{Exp.pulses}{
    expected number of water pulses per person per day
}
  \item{v.flow}{
    flow velocity [m/s]
}
  \item{Disp}{
    dispersion coefficient
}
  \item{flow.distances}{
    vector with flow distances [m]
}
  \item{Q.base}{
    base flow that is added to the pattern [l/s]
  }
  \item{exp.n.pulses.WD}{
    expected number of pulses per day and user on weekdays
  }
 \item{pattern.WD}{
    vector with rate of diurnal drug excretion pattern for weekdays
    [-]
  }
 \item{exp.n.pulses.WE}{
    expected number of pulses per day and user on weekends
  }
 \item{pattern.WE}{
    vector with rate of diurnal drug excretion pattern for weekends
    [-]
}
  \item{pulse.masses}{
    vector with pulse masses to sample from [mg]
  }
  \item{pulse.dur.min}{
    min width of a pulse [sec]
  }

  \item{pulse.dur.max}{
    max width of a pulse [sec]
  }

  \item{composite.duration}{
    duration of one composit sample [sec]
}
  \item{sampling.interval}{
    time between samples for time and flow proportional sampling [sec]
  }
  \item{v.prop.vol}{
    volume after that a sample is take for volume proportional sampling
    [l]. If \code{NULL} the volume is chosen so that roughly same amount of
    samples are take as with time-proportional and flow-proportional sampling.
  }
   \item{V.max}{
    volume of the pump sump [l], i.e. the pump is turned on if the
    volume is reached.
}
  \item{V.min}{
    volume where the pump is turned off [l]
}
  \item{pump.rate}{
    pump capacity [l/s]
  }
  \item{distance}{
    flow distance [m]
  }
 \item{flow, e1, e2}{
    object of class \code{flow}
}
  \item{x}{
    object of class \code{grav}, \code{pump}, \code{samples},
    respectively.
 }
  \item{...}{
    further arguments passed the generic functions.
}
}

\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
}
