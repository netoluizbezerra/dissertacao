####################
# Bandwidth issues #
####################
rm(list = ls())
library(kedd)
ndraws=100
horizon.simul=2000
alpha1=0.10
beta1 =0.85
mu=0.00
omega= 2
burnin = 100
model = garchSpec(model = list(mu = mu ,omega = omega, alpha = alpha1, beta = beta1, shape = 5), cond.dist = "std")
R <- ndraws #Numero de repeticoes do estudo 
paramEstimados_GAS <-  t(data.frame(omega = rep(0,R) ,alpha = rep(0,R), beta = rep(0,R), mu = rep(0,R), loglike = rep(0,R)))
paramEstimados_GAS_band <-  t(data.frame(omega.band = rep(0,R) ,alpha.band = rep(0,R), beta.band = rep(0,R), mu.band = rep(0,R), loglike.band = rep(0,R)))
paramEstimados_GAS.prior <-  t(data.frame(omega.prior = rep(0,R) ,alpha.prior = rep(0,R), beta.prior = rep(0,R), mu.prior = rep(0,R), loglike.prior = rep(0,R)))
Bandwidth <- t(data.frame(h = rep(0,R) , h2 = rep(0,R)))

for(n in 1:ndraws){
  vy <- as.numeric(garchSim(model, n = (horizon.simul)))
  vy = vy[burnin:horizon.simul]
  
  #Gas Paramétrico - Gaussiano
  df=0
  gas <- optim.gas(vy, df)
  beta.til <- gas$optimizer.new$par[3] - gas$optimizer.new$par[2]
  omega <- gas$optimizer.new$par[1]/(1 - beta.til)  
  
  #bandwidth selection
  vP = c(0.5,0.5) 
  list = get.ft(vy, horizon = (horizon.simul - burnin), gas)
  vf.temp = list$f_t; sigma = list$sigma
  dmu = gas$optimizer.new$par[4]
  std.res <- (vy - dmu)/exp(vf.temp/2)
  std.res <- scale(std.res, scale = TRUE)
  
  optimizer.sp.band = nlminb(start = vP, objective = sp.gas.likelihood.bandwidth, 
                             gradient = NULL, hessian = NULL, 
                             vy = vy, std.res = std.res, gas,
                             control = list(trace = 1), 
                             lower = 0.01, upper = 0.9)
  
  #Gas Semi-Paramétrico
  band = optimizer.sp.band$par
  h = band[1]
  h2 = band[2]
  
  sp_gas = sp.gas(vy, h, h2, horizon = length(vy), gas)
  beta.sp.til <- sp_gas$optimizer.sp$par[3] - sp_gas$optimizer.sp$par[2]
  omega.sp <- sp_gas$optimizer.sp$par[1]/(1 - beta.sp.til)
  
  #Gas Semi-Paramétrico - t-student
  h = 0.5
  h2 = 0.5
  sp_gas_prior = sp.gas(vy, h, h2, horizon = length(vy), gas)
  beta.sp.til_prior <- sp_gas_prior$optimizer.sp$par[3] - sp_gas_prior$optimizer.sp$par[2]
  omega.sp_prior <- sp_gas_prior$optimizer.sp$par[1]/(1 - beta.sp.til)
  
  paramEstimados_GAS[,n] <- unlist(c(omega, gas$optimizer.new$par[2], beta.til, gas$optimizer.new$par[4],-length(vy)*gas$optimizer.new$objective))
  paramEstimados_GAS_band[,n] <- unlist(c(omega.sp, sp_gas$optimizer.sp$par[2], beta.sp.til, sp_gas$optimizer.sp$par[4],-length(vy)*sp_gas$optimizer.sp$objective))
  paramEstimados_GAS.prior[,n] <- unlist(c(omega.sp_prior, sp_gas_prior$optimizer.sp$par[2], beta.sp.til_prior, sp_gas_prior$optimizer.sp$par[4] ,-length(vy)*sp_gas_prior$optimizer.sp$objective))
  Bandwidth[,n] <- unlist(c(band[1],band[2]))
}
  
