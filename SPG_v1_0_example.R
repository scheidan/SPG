## ================
## SPG: Get started
## ================

# Installing R, RStudio and the package SPG
# -----------------------------------------

# Please follow the instrucitons in the file "SPG_install.PDF" if you have not previously
# worked with R and RStudio and to facilitate the installation of the SPG package. For more 
# details about this example and a graphical representation of the case study the reader is 
# referred to the file "SPG_vignette.PDF".



## ==============
## 0) preparation
## ==============

# define where your results should be written to (the folder must exist already):
# -------------------------------------------------------------------------------

path <- "C:/temp/results"
setwd(path)

# load the SPG package that you have previously installed
# -------------------------------------------------------

library(SPG)

# start the help for the SPG package
# ----------------------------------

?SPG 



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