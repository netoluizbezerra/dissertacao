rm(list = ls())
# Libraries #
library(xlsx)
library(gsubfn)
require(pracma)
require(numDeriv)
library(fGarch)
library(miscTools)
library(quantmod)
library(dplyr)

###############################################
# Loglikelihood Function SP - GAS (normal)    #
# Inputs: vp = vector of initial parameters   #
#         vy = sample                         #
# Output: Loglikelihood                       #
###############################################

LogLikelihoodGasVolaUniv = function(vp, vy){
  vy = as.numeric(vy)
  cT = length(vy)
  dloglik = 0
  vf = rep(0, cT); vscore = rep(0, cT)
  domega = vp[1]; vA = vp[2]; vB = vp[3]; dmu = vp[4]
  vf[1] = domega/(1 - sum(vB))
  dist.param <- c(0,1)
  q = rep(0,cT); nabla = rep(0,cT)
  for(t in 1:cT){
    arg <- (vy[t] - dmu)/exp(vf[t]/2)
    q[t] = dnorm(arg, dist.param[1], dist.param[2])
    nabla[t]<- ddnorm(arg, dist.param[1], dist.param[2])/q[t]
    vscore[t] <-  - 1/2 * (nabla[t]) * arg  - 1/2     #Poderia usar direto: - 1/2*arg^2 - 1/2
    vf[t+1] <- domega + 2*vscore[t]*vA + vf[t]*vB
  }
  dloglik = -1/cT*(-1/2*(sum(vf)) + sum(log(q)))
  if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE" | sum(vB) > 1){
    dloglik <- 100
  }
  return(dloglik)
}

###############################################
# Loglikelihood Function SP - GAS (t-student) #
# Inputs: vp = vector of initial parameters   #
#         vy = sample                         #
# Output: Loglikelihood                       #
###############################################

llhGasVolaUniv.t = function(vp, vy){
  cT = length(vy)
  vy = as.numeric(vy)
  dloglik = 0
  vf = rep(0, cT); vscore = rep(0, cT)
  domega = vp[1]; vA = vp[2]; vB = vp[3]; dmu = vp[4]; ddf = vp[5]
  #Setup for SP-GAS:
  q = as.matrix(rep(0,(cT)))
  ft = as.matrix(rep(0,(cT)))
  nabla = as.matrix(rep(0,(cT)))
  vf[1] = domega/(1 - sum(vB))
  for(t in 1:cT){
    arg <- (vy[t] - dmu)/exp(vf[t]/2)
    q[t] <- dt(arg, ddf, log = TRUE)
    nabla.q = -(arg*(ddf + 1))/(ddf + arg^2)
    vscore[t] <- -1/2*(nabla.q*arg + 1)
    vf[t+1] <-  domega + (2*(ddf+3)/ddf)*vscore[(t)]*vA + vf[(t)]*vB
  }
  ft = vf[1:(cT-1)]
  dloglik = -1/cT*(-1/2*(sum(ft)) + sum(q))
  if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE" | sum(vB) > 1-1e-6){
    dloglik <- 100
  }
  return(dloglik)
}

###############################################
# Optimizer - SP GAS                          #
# Input : vp = vector of initial parameters   #
# Output: list of optimized parameters        #
###############################################

param <- function(optimiser){
  if(is.null(optimiser$par[5]) == TRUE){
    domega = optimiser$par[1]
    vA = optimiser$par[2]
    vB = optimiser$par[3]
    dmu = optimiser$par[4]
    list(domega = domega, vA = vA, vB = vB, dmu = dmu)
  } else{
    domega = optimiser$par[1]
    vA = optimiser$par[2]
    vB = optimiser$par[3]
    dmu = optimiser$par[4]
    ddf = optimiser$par[5]
    list(domega = domega, vA = vA, vB = vB, dmu = dmu, ddf = ddf)}
}


###############################################
# Kernel Estimator  - SP GAS                  #
# Inputs: vp = vector of initial parameters   #
#         vy = sample                         #
# Output: Loglikelihood                       #
###############################################

normal <- function(x){
  (1/sqrt(2*pi))*exp((-x^2)/2)
}
d.normal <- function(x){
  (-x/sqrt(2*pi))*exp((-x^2)/2)
}
dens.np <- function(x, x.i, h, k){
  n = length(x)
  result = rep(0,n)
  result = sum(k((x-x.i)/h))/(n*h)
  return(result)
}
ddens.np <- function(x, x.i, h,  k){
  n = length(x)
  result = 0
  result = sum(k((x-x.i)/h))/(n*h^2)
  return(result)
}

#############################################################
# Optimizer - SP GAS                                        #
# Input : vp = vector of parameters from parametric model   #
#         vy = sample                                       #
#         std.res = Stadardized residuals                   #
#         h = Bandwidth kernel density estimator            #
#         h2 = Bandwidth 1der. density kernel estimator     #
# Output: Loglikelihood SP - model                          #
#############################################################

sp.gas.likelihood = function(vp, vy, std.res, h, h2){
  vy <- as.numeric(vy)
  domega = vp[1]; vA = vp[2]; vB = vp[3]; dmu = vp[4]
  cT = length(vy)
  dloglik = 0
  vf = rep(0, cT); vscore = rep(0, cT)
  q.store = as.matrix(rep(0,cT))
  q.line = 0
  ft = as.matrix(rep(0,cT))
  nabla = as.matrix(rep(0,cT)) 
  vf[1] = domega/(1-vB)  
  h = h
  h2 = h2
  for(t in 1:length(vy)){
    arg <- (vy[t] - dmu)/exp(vf[t]/2)
    q.store[t] <- dens.np(x = std.res, x.i = arg, h = h, k = normal)
    q.line <- -ddens.np(x = std.res, x.i = arg, h = h2, k = d.normal)
    nabla[t] <- q.line/q.store[t]
    vscore[t] <-  -1/2 * (nabla[t]) * arg  - 1/2
    vf[t+1] <- domega + vscore[t]*vA + vf[t]*vB
  }
  dloglik = -1/cT*(-1/2*(sum(vf)) + sum(log(q.store)))
  
  if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE"){
    dloglik <- 100
  }
  return(dloglik)
}

get.ft = function(vy, horizon, gas){
  list[domega, vA, vB, dmu, ddf] <- param(gas$optimizer.new)
  f_t = rep(0,(horizon+1))
  score = rep(0,(horizon))
  domega = domega/(1-vB)
  f_t[1] = domega
  if(is.na(ddf) == TRUE){
    for(i in 1:horizon){
      arg <- (vy[i] - dmu)/exp(f_t[i]/2)
      nabla <- ddnorm(arg, 0, 1)/dnorm(arg, 0, 1)
      score[i] <- -1/2*(arg*nabla + 1)
      f_t[i+1] <- (1-vB)*domega + 2*vA*score[i] + vB*f_t[i]  
    }} else{
      for(i in 1:horizon){
        arg <- (vy[i] - dmu)/exp(f_t[i]/2)
        nabla <- -(arg*(ddf + 1))/(ddf + arg^2)
        score[i] <- -1/2*(arg*nabla + 1)
        f_t[i+1] <- domega*(1-vB) + (2*(ddf+3)/ddf)*score[i]*vA + vB*f_t[i]   
      } 
    }
  sigma = exp(f_t)
  list(f_t = f_t, sigma = sigma)
}

get.ft.sp = function(vy, df, h, h2, horizon, optimizer.sp, std.res){
  list[domega, vA, vB, dmu, ddf] <- param(optimizer.sp)
  f_t = rep(0,(horizon+1));
  vscore = rep(0,(horizon))
  q.store= rep(0,(horizon))
  domega = domega/(1-vB)
  f_t[1] = domega
  for(i in 1:horizon){
    arg <- (vy[i] - dmu)/exp(f_t[i]/2)
    q.store[i] <- dens.np(x = std.res, x.i = arg, h = h, normal)
    q.line <- -ddens.np(x = std.res, x.i = arg, h = h2, d.normal)
    nabla <- q.line/q.store[i]
    vscore[i] <-  -1/2*(arg*nabla + 1)
    f_t[i+1] <- (1-vB)*domega + 2*vscore[i]*vA + f_t[i]*vB
  }
  sigma = exp(f_t)
  list(f_t = f_t, sigma = sigma, q.store = q.store)
}

f.vp0 <- function(domega, vA, vB, dmu, ddf){
  if(ddf == 0){
    vp0 <<- c(domega, vA, vB, dmu)    
  } else{
    vp0 <<- c(domega, vA, vB, dmu, ddf)   
  }
}

optim.gas <- function(vy, df){
  if(df == 0){
    domega = 1; vA = 0.10; vB = 0.5; dmu = 0; ddf = 0
    vp0 <- f.vp0(domega, vA, vB, dmu, ddf)
    optimizer.new = nlminb(start = vp0, objective = LogLikelihoodGasVolaUniv, 
                           vy = vy, control = list(trace = 1), 
                           lower = -Inf, upper = Inf)
  } else {
    domega = 1; vA = 0.10; vB = 0.5; dmu = 0; ddf = 8
    vp0 <- f.vp0(domega, vA, vB, dmu, ddf)
    optimizer.new = nlminb(start = vp0, objective = llhGasVolaUniv.t, 
                           vy = vy, control = list(trace = 1), 
                           lower = -Inf, upper = Inf)
  }
  list(optimizer.new = optimizer.new)
}

sp.gas <- function(vy, h, h2, horizon, gas){
  start_time <- Sys.time()
  list[vf.temp, sigma] = get.ft(vy, horizon, gas)
  list[domega, vA, vB, dmu, ddf] = param(gas$optimizer.new)
  vf.temp = vf.temp[1:(length(vf.temp)-1)]
  std.res <- (vy - dmu)/exp(vf.temp/2)
  std.res <- scale(std.res, scale = TRUE)
  former.log.like = gas$optimizer.new$objective + 1
  new.log.like = gas$optimizer.new$objective
  check.convergence = rep(0,10)
  n = 0
  i = 0
  while(new.log.like < former.log.like & i < 4){
    i = i + 1
    print(paste("SP-GAS", i))
    if(i == 1){
      if(is.na(ddf) == TRUE){
        vP <- c(domega, vA, vB, dmu)  
      } else {
        vP <- c(domega, vA, vB, dmu, ddf)  
      }
    } else { 
      if(is.na(ddf) == TRUE){
        list[domega, vA, vB, dmu, ddf] = param(optimizer.sp)
        vP <- c(domega, vA, vB, dmu)  
      } else {
        list[domega, vA, vB, dmu, ddf] = param(optimizer.sp)
        vP <- c(domega, vA, vB, dmu, ddf)  
      }
    }
    former.log.like = new.log.like
    optimizer.sp = nlminb(start = vP, objective = sp.gas.likelihood, 
                          gradient = NULL, hessian = NULL, 
                          vy = vy, std.res = std.res,  h, h2,
                          control = list(trace = 1), 
                          lower = -Inf, upper = Inf)
    
    new.log.like <- optimizer.sp$objective
    check.convergence[i] <- optimizer.sp$convergence
    list[vf.temp, sigma] = get.ft.sp(vy, h = 0.5, h2 = 0.5, horizon = length(vy), optimizer.sp = optimizer.sp, std.res = std.res)
    vf.temp = vf.temp[1:(length(vf.temp)-1)]
    std.res <- (vy - dmu)/exp(vf.temp/2)
    std.res <- scale(std.res, scale = TRUE)
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  print(paste("Procedimento finalizado na iteração", i))
  list(optimizer.sp = optimizer.sp, check.convergence = check.convergence, i = i)
}

