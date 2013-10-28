\name{SPG-package}
\alias{SPG-package}
\alias{SPG}
\docType{package}
\title{
  Sewage pattern generator
}
\description{
  The package provides functions to generate flow and substance patterns
  in complex sewer systems. 
  
}

\details{
  \tabular{ll}{
    Package: \tab SPG\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2013-08-15\cr
    License: \tab GPL-2\cr
  }

  Gravity systems (open channel flow) and pumps can be linked to build models 
  of complex sewer systems. The main aim is to efficiently model realistic short-term 
  variations of flows and substances in sewers and to assess their effect on sampling.
  In addition to complex sewer networks, the package is also particularly useful to 
  optimize sampling at the effluent of individual premises (e.g. hospital, school, 
  prison) where short-term variations are typically highest and challenge the reliable 
  collection of representative average samples over a period of time.

  Different scenarios (e.g. distribution of substance across sub-catchments) and sampling 
  setups can be evaluated at relatively low compuational cost. The processes advection 
  and dispersion are considered.   
  
  To getting started see SPG \href{../doc/SPG.pdf}{Vignette}. It contains an example on how to i) define
  (sub-)catchments and pumps, ii) setup a sewer system (i.e. link sub-
  catchments, iii) run simulations and visualize flows and iv) compare
  different sampling strategies.
  
  For details see links to individual functions below.
  

}
\author{
  Andreas Scheidegger \email{andreas.scheidegger@eawag.ch} (for R-specific coding questions)
  
  Christoph Ort \email{christoph.ort@eawag.ch} (for input data and application-related questions)
  
  
}
\references{
  Ort C., Schaffner C., Giger W. and Gujer W. (2005) Modeling Stochastic Load Variations in Sewer 
  Systems. Water Science & Technology 52 (5), 113-122. 
  
  Ort C. and Gujer W. (2006) Sampling for representative micropollutant loads in sewer systems. 
  Water Science & Technology 54 (6-7), 169-176. 
  
  Ort C., Lawrence M.G., Rieckermann J. and Joss A. (2010) Monitoring Pharmaceuticals
  and Personal Care Products (PPCPs) and Illicit Drugs in Wastewater Systems: Are Your 
  Conclusions Valid? A Critical Review. Environmental Science & Technology 44 (16), 6024-6035.
  
  Ort C., Lawrence M.G., Reungoat J. and Mueller J.F. (2010) Sampling for PPCPs in sewers: 
  Comparison of Different Sampling Modes and Optimization Strategies. Environmental Science & 
  Technology 44 (16), 6289-6296.
  
  Ort et al. SEWPROF D1.1 (in prep., this software tool)
}

\seealso{
  \code{\link{def.grav}},  \code{\link{def.pump}}, \code{\link{take.samples}}, \code{\link{plot.flow}}, \code{\link{plot.samples}}
}


\examples{
## =========================
## 1) define gravity systems
## =========================

# properties that are valid for all catchments can be specified in global variables here
# and can be reused when defining the individual gravitiy systems.

# the following vector describes the diurnal variation in 2-h steps, in this
# example it is the same for the substance and for the flow
# --------------------------------------------------------------------------

diurnal.variation.Q.S <- c(4, 8, 15, 20, 12, 10, 11, 7, 6, 10, 8, 5)


# in a first scenario we assume that the substance of interest is consumed by 
# 0.5% of the population and that the consumers are distributed homogeneously
# over all sub-catchments
# ---------------------------------------------------------------------------

G1.S1 <- def.grav(population=40000, 
                     frac.consumers=0.005, 
                     pulse.mass=5,
                     diurnal.S.WD=diurnal.variation.Q.S,
                     diurnal.S.WE=diurnal.variation.Q.S,
                     Q.c.day=160, 
                     diurnal.Q=diurnal.variation.Q.S,
                     Q.base=0, 
                     flow.distance=3000)

G2.S1 <- def.grav(population=9000, 
                     frac.consumers=0.005, 
                     pulse.mass=5,
                     diurnal.S.WD=diurnal.variation.Q.S,
                     diurnal.S.WE=diurnal.variation.Q.S,
                     Q.c.day=160, 
                     diurnal.Q=diurnal.variation.Q.S,
                     Q.base=0, 
                     flow.distance=1000)

G3.S1 <- def.grav(population=1000, 
                     frac.consumers=0.005, 
                     pulse.mass=5,
                     diurnal.S.WD=diurnal.variation.Q.S,
                     diurnal.S.WE=diurnal.variation.Q.S,
                     Q.c.day=160, 
                     diurnal.Q=diurnal.variation.Q.S,
                     Q.base=0, 
                     flow.distance=250)

G4.S1 <- def.grav(population=500, 
                     frac.consumers=0.005,
                     pulse.mass=5,
                     diurnal.S.WD=diurnal.variation.Q.S,
                     diurnal.S.WE=diurnal.variation.Q.S,
                     Q.c.day=160, 
                     diurnal.Q=diurnal.variation.Q.S,
                     Q.base=0, 
                     flow.distance=100)


# print definitions
# -----------------

G1.S1
G2.S1
G3.S1
G4.S1



## ===============
## 2) define pumps
## ===============

P2 <- def.pump(V.max=200, 
                V.min=1, 
                pump.rate=100, 
                distance=3200,
                v.flow=1, 
                Disp=0.16)

P3 <- def.pump(V.max=10, 
                V.min=1, 
                pump.rate=25, 
                distance=1200,
                v.flow=1, 
                Disp=0.16)

P4 <- def.pump(V.max=2, 
               V.min=1, 
               pump.rate=5, 
               distance=100,
               v.flow=1, 
               Disp=0.16)

# print definitions
# -----------------

P2
P3
P4



## =============================================
## 3) build the system and generate the patterns
## =============================================

sim.dur <- 10   # simulate 10 days to minimize computational time. 
                # NOTE: for stable results about 100 days must be simulated.

temp.res <- 2   # a temporal resolution of 2 minutes is recommended


# calling the gravity system functions defined above (Gx.Sy) calculates wastewater 
# flows and substance patterns for the number of days (sim.dur) and temporal
# resolution (temp.res). Compose the drainage system by simply adding the 
# different flows. If the flows from a gravity system are collected in a pump sump and
# then pumped towards the point where the pattern should be evaluated (i.e. the 
# sampling point) or towards the next pump station you simply embrace the gravity system
# function with a pump function:
# --------------------------------------------------------------------------------------

flow.S1 = G1.S1(sim.dur, temp.res) + 
          P2(G2.S1(sim.dur, temp.res) + P3(G3.S1(sim.dur, temp.res))) + 
          P4(G4.S1(sim.dur, temp.res))


# plot flows
# ----------

plot(flow.S1)                     # plots all modeled days
plot(flow.S1, days=2, from=4)     # plots two days (here days 5 and 6)



## ========================
## 4) take composit samples
## ========================

# take 24-hour composite samples with 15 minutes sampling interval
# ----------------------------------------------------------------

samples.S1.15min <- take.samples(flow.S1, composite.duration=24, sampling.interval=15)


# compare results, plot and print sampling errors
# -----------------------------------------------

plot(samples.S1.15min, type="relative")    # the errors will be displayed as %
plot(samples.S1.15min, type="factor")      # the errors will be displayed as factor
plot(samples.S1.15min, type="absolute")    # the errors will be displayed as absolute values

samples.S1.15min
}
