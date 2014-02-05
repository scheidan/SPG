## =======================================================
## Project: pattern generator
##
## Description: file contains all functions for the patter generator
##              not all function are intended to be used uy the enduser
##
## File: functions.R
## Path: q:/Abteilungsprojekte/Eng/SWWData/Illicit_Drugs_Verantwortung_C_Ort/Pumps/Package/
##
## February  5, 2014 -- Andreas Scheidegger, Christoph Ort
##
## andreas.scheidegger@eawag.ch
## =======================================================


## =======================================================
##
## Internal functions, not directly accessed by the user
##
## All internal functions work with the following units only:
##  length:  meter
##  volumes: liter
##  time:    seconds
##  mass:    mg
## =======================================================

## =================================
## Function to generate flow pattern
## =================================

## t.sim          time to simulate [sec]
## population     number of person
## temp.res.sim   time resolution [sec]
## pattern        vector with rate of diurnal usage pattern [-]
## Q.c.day        daily water consumption per person and day [!! l/day !!]
## Exp.pulses     number of expected pulses for one person per day
## flow.distances vector of flow distances to sample from [m]
## pulse.dur    initial pulse width [sec]
## v.flow         flow velocity [m/s]
## Q.base         baseflow [m/s]
## Disp           dipersation coefficient [m^2/s]
##
## -> play with Exp.pulses and the initial sd of the pulses


gen.Q <- function(t.sim, population, temp.res.sim, pattern=1, Q.c.day=250, Exp.pulses=100,
                  flow.distances, pulse.dur,  v.flow=1, Q.base=0, Disp=0.16) {

    ## ---------------------------------
    ## define diurnal pattern for one household

    ## define rate of pulses (rate of a NHPP)
    lambda <- rep(pattern, each=24*60*60/temp.res.sim/length(pattern))
    lambda <- lambda/sum(lambda)*Exp.pulses/temp.res.sim

    ## ---------------------------------
    ## make parts with different flow distances

    flow.distances <- quantile(flow.distances, prob=seq(0.1,1,length=10)) # get quantiles of flow dist

    Q <- matrix(NA, ncol=length(unique(flow.distances)), nrow=t.sim/temp.res.sim)

    ## loop over different flow distances
    for(k in 1:length(unique(flow.distances))) {
        ## ---------------------------------
        ## sample pulses

        ## calculate pulse rate for every time step
        m <- rep(lambda*temp.res.sim*population, length.out=t.sim/temp.res.sim)/length(unique(flow.distances))
        
        ## sample from poisson distributions
        n.pulses <- rpois(t.sim/temp.res.sim, lambda=m)
        
        ## ---------------------------------
        ## convolute with distribution

        ## define impulse response function in relation of the sewer length
        ## parametrized as gamma distribution

        ## define delay and sd -> shape is fix
        delay <- flow.distances[k]/v.flow    # [sec]
        sd <- sqrt(pulse.dur^2 + 2*Disp*flow.distances[k]/v.flow^3) # [sec.]
        ## calculate scale
        shape <- 3
        scale <- sd/sqrt(shape)

        Resp <- dgamma(seq(0, t.sim-1, temp.res.sim)-delay, scale=scale, shape=shape)

        ## pad zeros to ensure a 'good' length for the FFT
        ln <- length(n.pulses)
        nn <- nextn(ln, 2)               #find good length
        n.pulses <- c(n.pulses, rep(0, nn-ln))
        Resp <- c(Resp, rep(0, nn-ln))
        ## convolute n.pulses with Resp
        Q[,k] <- convolve(n.pulses, rev(Resp), type="open")[1:ln] *
            Q.c.day/Exp.pulses

    }

    return(rowSums(Q) + Q.base)
}



## Q <- gen.Q(t.sim=1*24*3600, population=10, temp.res.sim=180, pattern=c(1,1,2,3,1,6,6,1,1,1),
##       Q.c.day=250, Exp.pulses=100,
##       flow.distances=c(100, 400, 600, 2000, 4000), pulse.dur=300,  v.flow=1, Q.base=0, Disp=0.16)
## plot(Q, type='l')



## =================================
## Function to generate substance load pattern
## =================================

## t.sim             time to simulate [sec]
## exp.n.pulses.WD   expected number of pulses per day and user
## pattern.WD        vector with rate of diurnal usage pattern for weekday [-]
## exp.n.pulses.WE   expected number of pulses par day and user
## pattern.WE        vector with rate of diurnal usage pattern, for weekend [-]
## pulse.masses       vector with pulse masses to sample from [mg]
## pulse.dur.min   min width of a pulse [sec]
## pulse.dur.max   max of puls with [sec]
## temp.res.sim      time resolution [sec]
## flow.distances    vector with flow distances [m]
## v.flow            flow velocity [m/s]
## Disp              dipersation coefficient [m^2/s]
##
##
## OUT:
## pattern of substance [mg/sec]


gen.S <- function(t.sim, exp.n.pulses.WD, pattern.WD=1,
                  exp.n.pulses.WE,
                  pattern.WE,
                  pulse.masses, pulse.dur.min, pulse.dur.max, temp.res.sim,
                  v.flow, Disp, flow.distances) {

    ## ---------------------------------
    ## define diurnal and weekly pattern one user

    ## define rate of pulses (rate of a NHPP)
    lambda.WD <- rep(pattern.WD, each=24*60*60/temp.res.sim/length(pattern.WD))
    lambda.WE <- rep(pattern.WE, each=24*60*60/temp.res.sim/length(pattern.WE))
    lambda.WD <- lambda.WD/sum(lambda.WD) * exp.n.pulses.WD
    lambda.WE <- lambda.WE/sum(lambda.WE) * exp.n.pulses.WE

    ## calculate pulse rate for every time step
    m <- c(lambda.WD, rep(c( rep(lambda.WD, 5), rep(lambda.WE, 2)), length.out=t.sim/temp.res.sim-length(lambda.WE)))
    ## (the first day will be removed later in function def.grav())

    ## sample from poisson distributions
    n.pulses <- rpois(t.sim/temp.res.sim, lambda=m)

    ## ----------------------
    ## model dispersion and pulse mass

    ## approximate pnorm(x,0,1) for faster execution
    pnorm.approx <- approxfun(seq(-6, 6, length=1000), pnorm(seq(-6, 6, length=1000)),
                              yleft=0, yright=1)

    pattern <- rep(0, t.sim/temp.res.sim)
    times <- seq(temp.res.sim, t.sim, length=t.sim/temp.res.sim)
    for(t in which(n.pulses>0)) {
        ## sample initial pulse width
        pulse.dur <- runif(n.pulses[t], pulse.dur.min, pulse.dur.max)
        ## sample pulse mass
        pulse.mass <- sample(pulse.masses, n.pulses[t], replace=TRUE)
        ## sample flow dist
        flow.dist <- sample(flow.distances, n.pulses[t], replace=TRUE)

        ## downstream sd of pulse
        pulse.sd <- sqrt(pulse.dur^2 + 2*Disp*flow.dist/v.flow^3)

        ## time of peak downstream
        t.peak <- t*temp.res.sim + flow.dist/v.flow # [sec]

        for(i in 1:n.pulses[t]) {
            pattern <- pattern +
              ## c(0, diff(pnorm(times, t.peak[i], pulse.sd[i]))*pulse.mass[i])/temp.res.sim
              c(0, diff(pnorm.approx((times-t.peak[i])/pulse.sd[i]))*pulse.mass[i])/temp.res.sim
          }
    }
    return(pattern)
}




## =================================
## Function to calculate a transformation due to a pump
## =================================

## flow = matrix with flows (column 1: Q[l/sec], column 2: S[mg/sec])
## V.max = Volume where pump is switched on [l]
## V.min = Volume where pump is switched off [l]
## pump.rate = pump capacity [l/sec]

pump.trans <- function(flow, V.max, V.min=0, pump.rate) {

    ## read temporal resolution [sec]
    temp.res.sim <- attr(flow, "temp.res.sim")

    ## V = Volume of water in 'tank'
    V <- rep(0, nrow(flow))
    ## S = mass of substance in 'tank'
    S <- rep(0, nrow(flow))

    ## water flow out of pump
    Q.out <- rep(0, nrow(flow))
    ## substance flow out of pump
    S.out <- rep(0, nrow(flow))

    ## state of pump
    state <- 'off'
    ## mass balances
    for(t in 2:nrow(flow)) {
        ## switch pump 'on' if volume is larger than V.max
        if(V[t-1] > V.max) {
            state <- 'on'
        }
        ## switch pump 'off' if V < V.min
        if(state=='on' & V[t-1] <= V.min) {
            state <- 'off'
        }

        ## set Q.out
        if(state=='on') {
            Q.out[t-1] <- min(pump.rate, V[t-1]/temp.res.sim + flow[t-1, 1]) #
            ## pump not more out than what is in the tank
        } else {
            Q.out[t-1] <- 0
        }

        ## mass balance for water
        V[t] <- V[t-1] + (flow[t-1, 1] - Q.out[t-1])*temp.res.sim

        ## mass balance for SUBSTANCE
        S.out[t-1] <- S[t-1]/(V[t-1]+flow[t-1, 1]*temp.res.sim)*Q.out[t-1]
        if(is.nan(S.out[t-1])) S.out[t-1] <- 0 # occurs only when V == 0
        S[t] <- S[t-1] + flow[t-1, 2]*temp.res.sim - S.out[t-1]*temp.res.sim

        ## prevent rounding errors (due to subtraction of similar numbers)
        if(S[t]<0) S[t] <- 0

    }

    ## return matrix with flows
    flow.out <- matrix(c(Q.out, S.out), ncol=2)
    dimnames(flow.out) <- dimnames(flow)
    class(flow.out) <- "flow"
    attr(flow.out, "temp.res.sim") <- temp.res.sim # store resolution as attribute

    return(flow.out)

}



## =================================
## Function to calculate sewer dispersation
## =================================

## flow:  flow object
## distance: sewer length [meters]
## flow.v: flow velocity [meter/sec]

disp.trans <- function(flow, distance, v.flow=1, Disp=0.16) {

    ## read temporal resolution [sec]
    temp.res.sim <- attr(flow, "temp.res.sim")
    t.sim <- nrow(flow)*temp.res.sim

    ## Calculate dispersion in function of the distance
    ## define delay and sd -> shape is fix
    delay <- distance/v.flow    # [sec]
    sd <- sqrt(1^2 + 2*Disp*distance/v.flow^3) # [sec.]
    ## calculate scale
    shape <- 3
    scale <- sd/sqrt(shape)

    Resp <- dgamma(seq(0, t.sim-1, temp.res.sim)-delay, scale=scale, shape=shape)*temp.res.sim

    ## convolute inputs with Resp

    Q.in <- flow[,1]
    S.in <- flow[,2]
    ln <- length(Q.in)
    
    ## pad zeros to ensure a 'good' length for the FFT
    nn <- nextn(ln, 2)               #find good length
    Q.in <- c(Q.in, rep(0, nn-ln))
    S.in <- c(S.in, rep(0, nn-ln))
    Resp <- c(Resp, rep(0, nn-ln))
    
    Q.out <- convolve(Q.in, rev(Resp), type="open")[1:ln]
    S.out <- convolve(S.in, rev(Resp), type="open")[1:ln]

    ## return matrix with flow
    flow.out <- matrix(c(Q.out[1:nrow(flow)], S.out[1:nrow(flow)]), ncol=2)
    dimnames(flow.out) <- dimnames(flow)
    ## set very small values (due to numerical errors) to zero
    flow.out[flow.out < sqrt(.Machine$double.eps)] <- 0

    class(flow.out) <- "flow"
    attr(flow.out, "temp.res.sim") <- temp.res.sim # store resolution as attribute

    return(flow.out)

}



## =================================
## Function to 'grap' samples
## =================================

## flow                flow object
## composite.duration  duration of one composit sample [sec]
## sampling.interval   how frequent samples are taken sampling interval [sec]
## v.prop.vol          after how many litre volume propoprtional samples are taken? [l]
##                     if 'NULL' it is calculated so that ~same number of samples is taken as
##                     with the other sampling approaches.


sampling <- function(flow, composite.duration, sampling.interval=NULL, v.prop.vol=NULL) {

    Q <- flow[,1]
    S <- flow[,2]
    
    composite.duration <- composite.duration/60 #[min]
    if(!is.null(sampling.interval)) sampling.interval <- sampling.interval/60 #[min]
    
    ## read out time and time.step
    temporal.resolution.sim <- attr(flow, "temp.res.sim") # [sec]
    time <- seq(temporal.resolution.sim, by=temporal.resolution.sim, length=length(Q))/60 # [min]
    time.step <- temporal.resolution.sim/60 # [min]
    
    ## determine number of available days
    nr.days <- max(time)/24/60
    
    
    # if(!is.null(sampling.interval)) sampling.interval <- sampling.interval*60 # change unit to seconds??????
    
    if(is.null(sampling.interval)) {
      average.daily.ww.volume <- sum(Q)*temporal.resolution.sim/nr.days
      t.prop.n <- average.daily.ww.volume/v.prop.vol
      sampling.interval <- 24*3600/t.prop.n
      sampling.interval <- round(sampling.interval/60,1)   #[min]
    }
    
    
    ## determine number of composite samples within the simulation period
    nr.composite.samples <- floor(max(time)/composite.duration)

    ## determine start and end points of composite samples
    start.times <- seq(from=min(time),by=composite.duration,length.out=nr.composite.samples) # [min]
    start.index <- start.times/time.step

    length.composite.sample.index <- diff(start.index)[1]

    ## define indeces for each sample
    indeces.composite.samples <- rep(c(1:nr.composite.samples),each=length.composite.sample.index)

    ## Determine true concentration for each time step
    C.true <- S/Q

    ## function to interpolate values
    get.C <- approxfun(time, C.true, rule=2)
    get.Q <- approxfun(time, Q, rule=2)
    
    ## Determine true full scale loads per composite sample
    true.flows <- tapply(Q, indeces.composite.samples, function(x) sum(x) * temporal.resolution.sim)

    ## Determine true full scale loads per composite sample
    true.loads <- tapply(S, indeces.composite.samples, function(x) sum(x) * temporal.resolution.sim)

    
    ## Time-proportional sampling
    ## ==========================

    ## determine points in time when samples are taken for time-proportional sampling
    ## ------------------------------------------------------------------------------
    t.prop.sample.time <- seq(sampling.interval, max(time), sampling.interval)

    ## take time-prop samples and calculate loads based on mean composite concentration
    ## -------------------------------------------------------------------------------
    conc.t.samples <- get.C(t.prop.sample.time)
    sample.number <- as.numeric(cut(t.prop.sample.time, c(start.times, Inf), include.lowest=TRUE, right=FALSE))
    t.sampled.loads <- tapply(conc.t.samples, sample.number, mean) * true.flows

    ## Flow-proportional sampling
    ## ==========================

    ## determine points in time when samples are taken for flow-proportional sampling
    ## ------------------------------------------------------------------------------
    f.prop.sample.time <- t.prop.sample.time
                       
    ## take flow-prop samples and calculate loads based on weighted concentrations
    ## ---------------------------------------------------------------------------
    load.f.samples.weighted <- get.C(f.prop.sample.time) * get.Q(f.prop.sample.time)
    sample.number <- as.numeric(cut(t.prop.sample.time, c(start.times, Inf), include.lowest=TRUE, right=FALSE))
    volumes.f.samples <- tapply(get.Q(f.prop.sample.time), sample.number, sum) #total Q over each composite period
    f.sampled.loads <- tapply(load.f.samples.weighted, sample.number, sum) / volumes.f.samples * true.flows


    ## Volume-proportional sampling
    ## ============================

    ## if not given, determine sample increment "v.prop.vol" for volume proportional sampling
    ## ------------------------------------------------------------------------

    average.daily.ww.volume <- sum(Q)*temporal.resolution.sim/nr.days
    
    if(is.null(v.prop.vol)) {
      v.prop.n <- composite.duration/sampling.interval # determine number of samples equivalent to time- or flow-prop. sampling
      v.prop.vol <- average.daily.ww.volume/24*composite.duration/60/v.prop.n # increment "v.prop.vol"
    }

    ## waring if t.res is to large for volume prop sampling

    if(max(Q)*temporal.resolution.sim > v.prop.vol) {
      warning("Volume proportional sampling is unreliable because 't.res' is to large!", call.=FALSE)
    }
   
    ## determine points in time when samples are taken with volume-proportional sampling mode
    ## --------------------------------------------------------------------------------------
    v.prop.sample.index <- rep(FALSE, length(Q))
    cum.vol <- 0

    for(i in 1:length(Q))
    {
        cum.vol <- cum.vol+Q[i]*temporal.resolution.sim
        if(cum.vol>v.prop.vol)
        {
            v.prop.sample.index[i] <- TRUE
            cum.vol <- 0
        }
    }
    v.prop.sample.times <- time[v.prop.sample.index]

    ## take volume-prop samples and calculate loads based on mean composite concetration
    ## ---------------------------------------------------------------------------------

    conc.v.samples <- get.C(v.prop.sample.times) #OK
    sample.number <- as.numeric(cut(v.prop.sample.times, c(start.times, Inf), include.lowest=TRUE, right=FALSE))
    periods.with.sample.v.prop <- unique(sample.number)
    v.sampled.loads <- rep(NA, length(true.flows))
    v.sampled.loads[periods.with.sample.v.prop] <- tapply(conc.v.samples, sample.number, mean) *
    true.flows[periods.with.sample.v.prop]


    ## calculate concentration
    ## ------------------------

    C.true <- true.loads/true.flows
    C.time.prop <- t.sampled.loads/true.flows
    C.flow.prop <- f.sampled.loads/true.flows
    C.Vol.prop <- v.sampled.loads/true.flows

    
     ## ----------------------
    ## return concentrations
    res <- list(C.true=C.true, C.time.prop=C.time.prop, C.flow.prop=C.flow.prop,
                C.Vol.prop=C.Vol.prop, Vol=true.flows/1000)
    class(res) <- "samples"
    attr(res, "sampling.interval.t.f") <- sampling.interval
    attr(res, "sampling.interval.v") <- round(24*60/(average.daily.ww.volume/v.prop.vol),1)
    attr(res, "sampling.interval.v.min") <- round(v.prop.vol/min(Q)/60,0)
    attr(res, "sampling.interval.v.max") <- round(v.prop.vol/max(Q)/60,0)
    attr(res, "v.prop.vol") <- round(v.prop.vol/1000,1)
    return(res)
}




## =======================================================
##
## Enduser functions
##
## !! non-SI units can be used here as inputs !!
## =======================================================

## =================================
## Generate a function for a gravity system
## =================================

## the resulting function has only two arguments: sim.dur, temp.res

def.grav <- function(population, 
                      frac.consumers=NULL, 
                      exp.n.pulses=NULL,
                      WE.factor=1,
                      pulse.mass,
                      diurnal.S.WD=1,
                      diurnal.S.WE=diurnal.S.WD,
                      Q.c.day=250, 
                      diurnal.Q=1,
                      Q.base=0, 
                      flow.distance=2000,
                      v.flow=1, 
                      Disp=0.16,
                      pulse.dur.min=5, 
                      pulse.dur.max=95) {
    ## checks
    if(is.null(frac.consumers) & is.null(exp.n.pulses)){
        stop("Either 'frac.consumers' or 'exp.n.pulses' must be specified!")
    }
    if(!is.null(frac.consumers) & !is.null(exp.n.pulses)){
        stop("Only 'frac.consumers' or 'exp.n.pulses' can be specified!")
    }

    ## calcualte some default values
    if(!is.null(frac.consumers)){
        exp.n.pulses.WD <- population*frac.consumers
        exp.n.pulses.WE <- population*frac.consumers*WE.factor
    }
    if(!is.null(exp.n.pulses)){
        exp.n.pulses.WD <- exp.n.pulses
        exp.n.pulses.WE <- exp.n.pulses*WE.factor
    }

    ## generate puls masses if only mean is given
    if(length(pulse.mass)==1) pulse.masses <- pmax(0, rnorm(1000, pulse.mass, pulse.mass/5))
    if(length(pulse.mass)>1) pulse.masses <- pulse.mass

    ## generate flow distances if only mean is given
    if(length(flow.distance)==1) flow.distances <- pmax(0, rnorm(1000, flow.distance, flow.distance/4))
    if(length(flow.distance)>1) flow.distances <- flow.distance


    ff <- function(sim.dur, temp.res=5) {
   
        temp.res.sim <- floor(temp.res*60)

        ## generate substance pattern
        S <- gen.S(t.sim=(sim.dur+1)*86400,
                   exp.n.pulses.WD=exp.n.pulses.WD,
                   exp.n.pulses.WE=exp.n.pulses.WE,
                   pattern.WD=diurnal.S.WD,
                   pattern.WE=diurnal.S.WE,
                   pulse.masses=pulse.masses,
                   pulse.dur.min=pulse.dur.min, pulse.dur.max=pulse.dur.max,
                   temp.res.sim=temp.res.sim,
                   v.flow=v.flow, Disp=Disp, flow.distances=flow.distances)

        ## generate Q pattern
        Q <- gen.Q(t.sim=(sim.dur+1)*86400, population=population, temp.res.sim=temp.res.sim,
                   pattern=diurnal.Q, Q.c.day=Q.c.day, Exp.pulses=100,
                   flow.distances=flow.distances, pulse.dur=300,  v.flow=v.flow,
                   Q.base=Q.base, Disp=Disp)

        ## compile to a 'flow' object
        flows <- matrix(c(Q,S), ncol=2)
        flows <- flows[-(1:(86400/temp.res.sim)),] # remove the 'burn-in'-day
        ## vector for weekend and work days
        day <- rep(c(rep("WD", 86400/temp.res.sim*5), rep("WE", 86400/temp.res.sim*2)), length=nrow(flows))
        rownames(flows) <- paste(seq(temp.res, nrow(flows)*temp.res, length=nrow(flows)), day, sep='_')

        colnames(flows) <- c("Q", "S")
        class(flows) <- "flow"
        attr(flows, "temp.res.sim") <- temp.res.sim # store resolution as attribute
        return(flows)
    }

    ## add attributes for print function
    attr(ff, "population") <- population
    attr(ff, "frac.consumers") <- frac.consumers
    attr(ff, "exp.n.pulses") <- exp.n.pulses
    attr(ff, "WE.factor") <- WE.factor
    attr(ff, "pulse.mass") <- pulse.mass
    attr(ff, "diurnal.S.WD") <- diurnal.S.WD
    attr(ff, "diurnal.S.WE") <- diurnal.S.WE
    attr(ff, "diurnal.Q") <- diurnal.Q
    attr(ff, "v.flow") <- v.flow
    attr(ff, "Disp") <- Disp
    attr(ff, "flow.distance") <- flow.distance
    attr(ff, "Q.base") <- Q.base

    class(ff) <- "grav"
    return(ff)
}



## =================================
## generate pump function
## =================================

## combines pump and dispersion
##
## flows = flow object
## V.max = Volume where pump is switched on [m^3]
## V.min = Volume where pump is switched off [m^3]
## pump.rate = pump capacity [l/sec]
## distance: sewer distance [meters]
## flow.v: flow velocity [meter/sec]
## Disp    dispersion

## return a function that thake a 'flow'-object as input

def.pump <- function(V.max, V.min=0, pump.rate, distance=0, v.flow=1, Disp=0.16) {

    ff <- function(flow) {
        flows <- pump.trans(flow, V.max=V.max*1e3, V.min=V.min*1e3, pump.rate=pump.rate)
        ## compute dispersion of distance > 0
        if(distance>0) {
            flows <- disp.trans(flows, distance=distance, v.flow=v.flow, Disp=Disp)
        }
        return(flows)
    }

    ## add attributes for print function
    attr(ff, "V.max") <- V.max
    attr(ff, "V.min") <- V.min
    attr(ff, "pump.rate") <- pump.rate
    attr(ff, "distance") <- distance
    attr(ff, "v.flow") <- v.flow
    attr(ff, "Disp") <- Disp
    class(ff) <- "pump"

    return(ff)
}



## =================================
## Sampling function
## =================================

## flows               flow object
## composite.duration  duration of one composit sample [h]
## sampling.interval   how frequent samples are taken sampling interval (for time and flow proportional) [min]
## v.prop.vol          volume after that a sample is take for volume proportional sampling [m^3]

take.samples <- function(flow, composite.duration=24, sampling.interval=NULL, v.prop.vol=NULL) {

  if(is.null(sampling.interval) & is.null(v.prop.vol)) {
    stop("Please provide either a sampling interval (t-/f-prop sampling) and/or sampling increment (v-prop. sampling)", call.=FALSE)
  }
  
  if(!is.null(v.prop.vol)) v.prop.vol <- v.prop.vol*1000 # change unit to litre
  if(!is.null(sampling.interval)) sampling.interval <- sampling.interval*60
  
  C <- sampling(flow, composite.duration=composite.duration*3600,
                sampling.interval=sampling.interval,
                v.prop.vol=v.prop.vol)
  attr(C, "composite.duration") <- composite.duration
  if(!is.null(sampling.interval)) attr(C, "sampling.interval.t.f") <- sampling.interval/60
  if(!is.null(v.prop.vol)) attr(C, "v.prop.vol") <- v.prop.vol/1000   # in m3!
  return(C)
}


## =======================================================
##
## Plot and print functions for objects of the classes
## "pump", "grav", "flow", "samples"
##
## =======================================================


## =================================
## Functions for 'flow' objects
## =================================

## ----------------------
## plot function

plot.flow <- function(x, days=NULL, from=0, title=NULL, ...) {

    temp.res.sim <- attr(x, "temp.res.sim")

    day.max <- nrow(x)*temp.res.sim/(24*3600) # number of simulated days
    if(from >= day.max) stop(paste("Only", day.max, "days are simulated! Choose 'from' smaller."))
    if(is.null(days)) days <- day.max

    t.min <- from*24*3600
    t.max <- (days+from)*24*3600
    tt <- (1:nrow(x))*temp.res.sim
    select <- tt>t.min & tt<t.max       # vector to select values to be plotted

    ## position of end of days
    t.days <- (0:(days+from+1))*24*3600/temp.res.sim
    ## position of begin and end of weekends
    t.WE.beginn <- (seq(5, max(7, (days+from+1)), 7)-from)  * 24*3600/temp.res.sim
    ##t.WE.beginn <- (7*(1:(days+from+1))-from-1)*24*3600/temp.res.sim
    t.WE.end <- t.WE.beginn+2*24*3600/temp.res.sim

    ## color for WD/WE box
    col.WD <- rgb(0.60, 0.80, 0.85)
    col.WE <- rgb(0.45, 0.65, 0.45)

    ## save par settings
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    par(mfrow=c(3,1))
    ## space for main title
    if(!is.null(title)) par(oma=c(0,0,2.5,0))

    if(days>1) {     
      xlab <-'days'
      max.day.plot <- min(from+days, nrow(x)*temp.res.sim/(24*3600))
      xlabels <- seq(from, max.day.plot, by=2)
      at.x <- t.days[xlabels-from+1]      #plotting position
      
      ## change unit is hours
    } else  {
      xlab <-'hours'
      xlabels <- pretty(c(0, min(days, nrow(x)*temp.res.sim/(24*3600))*24), n=12, min.n=6)
      at.x <- seq(0, sum(select), length=length(xlabels))
    }

    
    ## plot S
    plot(x[select,2], type='n',
         main="substance pattern", xlab=xlab, ylab="load [mg/sec]",
         ylim=c(0, max(x[select,2])*1.1),
         las=1, xaxt='n', ...)
    rect(1, -1e6, sum(select), 0, border=NA, col=col.WD)
    rect(t.WE.beginn, -1e6, t.WE.end, 0, border=NA, col=col.WE)
    abline(v=t.days, lty=3, col=4)
    lines(x[select,2])
    axis(side=1, at=at.x, labels=xlabels)
    box()

    ## plot Q
    plot(x[select,1], type='n',
         main="Q pattern", xlab=xlab, ylab="Q [l/sec]", ylim=c(0, max(x[select,1])*1.1),
         las=1, xaxt='n', ...)
    rect(1, -1e6, sum(select), 0, border=NA, col=col.WD)
    rect(t.WE.beginn, -1e6, t.WE.end, 0, border=NA, col=col.WE)
    abline(v=t.days, lty=3, col=4)
    lines(x[select,1])
    axis(side=1, at=at.x, labels=xlabels)
    box()
    
    ## plot concentration
    conc <- x[select,2]/x[select,1]*1000
    plot(conc, type='n', main="concentration", xlab=xlab,
         ylab="c [microgram/l]",
         ylim= c(0, max(conc[is.finite(conc)])*1.1),
         las=1, xaxt='n', ...)
    rect(1, -1e6, sum(select), 0, border=NA, col=col.WD)
    rect(t.WE.beginn, -1e6, t.WE.end, 0, border=NA, col=col.WE)
    abline(v=t.days, lty=3, col=4)
    lines(conc)
    axis(side=1, at=at.x, labels=xlabels)
    box()
    
    ## add main title
    if(!is.null(title)) mtext(title, outer=TRUE)
}


## ----------------------
## function to add two "flow" objects

## checks temp.resolution and length before adding the flows

"+.flow" <- function(e1, e2 = NULL) {
    temp.res.sim1 <- attr(e1, "temp.res.sim")
    temp.res.sim2 <- attr(e2, "temp.res.sim")
    if(temp.res.sim1 !=  temp.res.sim2) stop("The two flow objects must have the same 'temp.res'!")
    if(nrow(e1) != nrow(e2)) stop("The two flow objects must have the same simulation duration!")

    NextMethod()                        # call method to add matrices
}




## =================================
## Functions for 'samples' objects
## =================================

## ----------------------
## plot function

plot.samples <- function(x, type='relative', title=NULL, ...) {

    if(type!='absolute' & type!='relative' & type!='factor')
        stop("Argument 'type' must either be 'relative', 'absolute', or 'factor'!")

    composite.duration <- attr(x, "composite.duration")
    sampling.interval.t.f <- attr(x, "sampling.interval.t.f")
    sampling.interval.v <- attr(x, "sampling.interval.v")
    sampling.max.dt <- attr(x, "sampling.interval.v.min") 
    sampling.min.dt <- attr(x, "sampling.interval.v.max")
    v.prop.vol <- attr(x, "v.prop.vol")
    
    n.samples <- length(x[[1]])
   
    ## compute relative errors
    if(type=='relative') {
        Errors <- data.frame(time.prop=(x$C.time.prop-x$C.true)/x$C.true,
                             vol.prop=(x$C.Vol.prop-x$C.true)/x$C.true,
                             flow.prop=(x$C.flow.prop-x$C.true)/x$C.true)
        
        text.main <- paste("Relative errors based on ", n.samples, " ", composite.duration,
                           "h-composite samples\n sampling interval(t/f-prop.): dt=",
                           sampling.interval.t.f, "min\n avrg. sampling interval(vol-prop.*): dt=",
                           sampling.interval.v, "min", sep="")
        
        text.sub <- paste("*dV=", v.prop.vol,"m3 (dT=", sampling.min.dt, " - ", sampling.max.dt, "min)", sep="")
        
        text.y <- "relative error [%]"
        boxplot(Errors*100, main=text.main, sub=text.sub, ylab=text.y, ...)
        abline(h=c(-100,0), lty=3)
    }
    if(type=='factor') {
        Errors <- data.frame(time.prop=x$C.time.prop/x$C.true,
                             vol.prop=x$C.Vol.prop/x$C.true,
                             flow.prop=x$C.flow.prop/x$C.true)
        Errors[Errors<=0] <- NA

        text.main <- paste("Errors as factor based on ", n.samples, " ", composite.duration,
                           "h-composite samples\n sampling interval(t/f-prop.): dt=",
                           sampling.interval.t.f, "min\n avrg. sampling interval(vol-prop.*): dt=",
                           sampling.interval.v, "min", sep="")
        
        text.sub <- paste("*dV=", v.prop.vol,"m3 (dT=", sampling.min.dt, " - ", sampling.max.dt, "min)", sep="")
        
        text.y <- "log10 (measured C / true C)"

        boxplot(log10(Errors), main=text.main, sub=text.sub, ylab=text.y, ...)
        abline(h=0, lty=3)
    }
     if(type=='absolute') {
         Errors <- data.frame(time.prop=(x$C.time.prop-x$C.true),
                             vol.prop=(x$C.Vol.prop-x$C.true),
                             flow.prop=(x$C.flow.prop-x$C.true))

        text.main <- paste("Absolute errors based on ", n.samples, " ", composite.duration,
                           "h-composite samples\n sampling interval(t/f-prop.): dt=",
                           sampling.interval.t.f, "min\n avrg. sampling interval(vol-prop.*): dt=",
                           sampling.interval.v, "min", sep="")
         
        text.sub <- paste("*dV=", v.prop.vol,"m3 (dT=", sampling.min.dt, " - ", sampling.max.dt, "min)", sep="")
         
        text.y <- "absolute error [microgram/l]"
                boxplot(1000*Errors, main=text.main, sub=text.sub, ylab=text.y, ...)
                abline(h=0, lty=3)
    }

}

## ----------------------
## print function

print.samples <- function(x, ...) {
    composite.duration <- attr(x, "composite.duration")
    sampling.interval <- attr(x, "sampling.interval")
    sampling.interval <- attr(x, "sampling.interval")
    v.prop.vol <- attr(C, "v.prop.vol") # in litre!
    
    n.samples <- length(x[[1]])

    ## compute relative errors
    error.time.prop <- (x$C.time.prop-x$C.true)/x$C.true
    error.flow.prop <- (x$C.flow.prop-x$C.true)/x$C.true
    error.Vol.prop <- (x$C.Vol.prop-x$C.true)/x$C.true

    cat("Samples: \n", ...)
    cat("Sampling errors are based on ", n.samples, " ",
        composite.duration, "h composit samples, sampling interval: ",
        sampling.interval, "min\n", sep='', ...)
    cat("\n relative error with time proportional sampling:\n")
    print(summary(error.time.prop[is.finite(error.time.prop)]), ...)
    cat("\n relative error with volume proportional sampling:\n")
    print(summary(error.Vol.prop[is.finite(error.Vol.prop)]), ...)
    cat("\n relative error with flow proportional sampling:\n")
    print(summary(error.flow.prop[is.finite(error.flow.prop)]), ...)
}

## =================================
## Functions for 'pump' objects
## =================================


## ----------------------
## print function

print.pump <- function(x, ...) {

    ## read attributes
    V.max <- attr(x, "V.max")
    V.min <- attr(x, "V.min")
    pump.rate <- attr(x, "pump.rate")
    distance <- attr(x, "distance")
    v.flow <- attr(x, "v.flow")
    Disp <- attr(x, "Disp")

    ## print
    print("Pump:", ...)
    cat("V.max [m^3]:", V.max,
        "\nV.min [m^3]:", V.min,
        "\npump.rate [l/s]:", pump.rate,
        "\ndistance [m]:", distance,
        "\nDispersion [m^2/s]:", Disp,
        "\nv.flow [m/s]:", v.flow,'\n', ...)
}

## =================================
## Functions for 'grav' objects
## =================================

## ----------------------
## print function

print.grav <- function(x, ...) {

    ## read attributes
    population <- attr(x, "population")
    exp.n.pulses <- attr(x, "exp.n.pulses")
    WE.factor <- attr(x, "WE.factor")
    frac.consumers <- attr(x, "frac.consumers")
    pulse.mass <- attr(x, "pulse.mass")
    diurnal.S.WD <- attr(x, "diurnal.S.WD")
    diurnal.S.WE <- attr(x, "diurnal.S.WE")
    diurnal.Q <- attr(x, "diurnal.Q")
    flow.distance <- attr(x, "flow.distance")
    Q.base <- attr(x, "Q.base")
    v.flow <- attr(x, "v.flow")
    Disp <- attr(x, "Disp")

    ## print
    print("Gravity system:")
    if(!is.null(frac.consumers)) {
        cat("population [-]:", population,
            "\nconsumer fraction [-]:", frac.consumers,
            "\nweekend factor [-]:", WE.factor,
            "\nmean pulse mass [mg]:", mean(pulse.mass),
            "\nbase flow [l/s]:", Q.base,
            "\nDispersion [m^2/s]:", Disp,
            "\nflow velocity [m/s]:", v.flow, '\n', ...)
    }
    if(is.null(frac.consumers)) {
        cat("population [-]:", population,
            "\nexpected number of pulses per day [-]:", exp.n.pulses,
            "\nweekend factor [-]:", WE.factor,
            "\nmean pulse mass [mg]:", mean(pulse.mass),
            "\nbase flow [l/s]:", Q.base,
            "\nDispersion [m^2/s]:", Disp,
            "\nflow velocity [m/s]:", v.flow, '\n', ...)
    }
}


## =======================================================
