###############################################
# Densities and Kernel Estimator  - SP GAS    #
# Inputs: vp = vector of initial parameters   #
#         vy = sample                         #
# Output: Loglikelihood                       #
###############################################

kernel.t <- function(x, ddf = 5) gamma((ddf + 1)/2)/(sqrt(ddf*pi)*gamma(ddf/2)) * (1 + x^2/ddf)^(-(ddf + 1)/2)

d.t <- function(x, ddf = 5) gamma((ddf + 1)/2)/(sqrt(ddf*pi)*gamma(ddf/2)) * (-ddf - 1) * x *((ddf + (x)^2)/ddf)^(-(ddf + 3)/2)

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
    vf[t+1] <- domega + 2*vscore[t]*vA + vf[t]*vB
  }
  dloglik = -1/cT*(-1/2*(sum(vf)) + sum(log(q.store)))
  
  if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE"){
    dloglik <- 100
  }
  return(dloglik)
}

sp.gas.likelihood_2 = function(vp, vy, std.res, h, h2, ddf = 5){
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
    q.store[t] <- dens.np(x = std.res, x.i = arg, h = h, k = kernel.t)
    q.line <- -ddens.np(x = std.res, x.i = arg, h = h2, k = d.t)
    nabla[t] <- q.line/q.store[t]
    vscore[t] <-  -1/2 * (nabla[t]) * arg  - 1/2
    vf[t+1] <- domega + 1/2*vscore[t]*vA + vf[t]*vB
  }
  dloglik = -1/cT*(-1/2*(sum(vf)) + sum(log(q.store)))
  if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE"){
    dloglik <- 100
  }
  return(dloglik)
}



sp.gas_2 <- function(vy, h, h2, horizon, gas){
  start_time <- Sys.time()
  list = get.ft(vy, horizon, gas)
  vf.temp = list$f_t; sigma = list$sigma
  par <- param(gas$optimizer.new)
  domega  = par$domega; vA = par$vA; vB = par$vB; dmu = par$dmu; ddf = par$ddf
  vf.temp = vf.temp[1:(length(vf.temp)-1)]
  std.res <- (vy - dmu)/exp(vf.temp/2)
  std.res <- scale(std.res, scale = TRUE)
  former.log.like = gas$optimizer.new$objective + 1
  new.log.like = gas$optimizer.new$objective
  check.convergence = rep(0,10)
  n = 0
  i = 0
  while(new.log.like < former.log.like & i < 2){
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
        par <- param(optimizer.sp) ####VER ESTA LINHA 
        domega  = par$domega; vA = par$vA; vB = par$vB; dmu = par$dmu; ddf = par$ddf
        vP <- c(domega, vA, vB, dmu)  
      } else {
        par <- param(optimizer.sp)
        domega  = par$domega; vA = par$vA; vB = par$vB; dmu = par$dmu; ddf = par$ddf
        vP <- c(domega, vA, vB, dmu, ddf)  
      }
    }
    former.log.like = new.log.like
    optimizer.sp = nlminb(start = vP, objective = sp.gas.likelihood_2, 
                          gradient = NULL, hessian = NULL, 
                          vy = vy, std.res = std.res,  h, h2,
                          control = list(trace = 1), 
                          lower = -Inf, upper = Inf)
    
    new.log.like <- optimizer.sp$objective
    check.convergence[i] <- optimizer.sp$convergence
    list_new = get.ft.sp(vy, h = h, h2 = h2, horizon = length(vy), optimizer.sp = optimizer.sp, std.res = std.res)
    vf.temp = list_new$f_t; sigma = list_new$sigma
    vf.temp = vf.temp[1:(length(vf.temp)-1)]
    std.res <- (vy - dmu)/exp(vf.temp/2)
    std.res <- scale(std.res, scale = TRUE)
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  print(paste("Procedimento finalizado na iteração", i))
  list(optimizer.sp = optimizer.sp, check.convergence = check.convergence, i = i)
}
