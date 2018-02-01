####################
# DGP - GARCH - T ##
####################
rm(list=ls())

ndraws=500
horizon.simul=2000
alpha1=0.10
beta1 =0.85
mu=0.00
omega= 2
burnin = 101
model = garchSpec(model = list(omega = omega, alpha = alpha1, beta = beta1, shape = 5), cond.dist = "std")
R <- ndraws #Numero de repeticoes do estudo 
paramEstimados_GARCH <-  t(data.frame(mu = rep(0,R) ,omega = rep(0,R), alpha = rep(0,R), beta = rep(0,R), loglike = rep(0,R)))
paramEstimados_GARCH.t <-  t(data.frame(mu = rep(0,R) ,omega = rep(0,R), alpha = rep(0,R), beta = rep(0,R),  ddf = rep(0,R),loglike = rep(0,R)))
paramEstimados_GAS <-  t(data.frame(omega = rep(0,R) ,alpha = rep(0,R), beta = rep(0,R), mu = rep(0,R), loglike = rep(0,R)))
paramEstimados_GAS.t <-  t(data.frame(omega = rep(0,R) ,alpha = rep(0,R), beta = rep(0,R), mu = rep(0,R), ddf = rep(0,R), loglike = rep(0,R)))
paramEstimados_SPGAS <- t(data.frame(omega.gas.sp = rep(0,R) ,alpha.gas.sp = rep(0,R), beta.gas.sp = rep(0,R), mu.gas.sp = rep(0,R), loglike.gas.sp = rep(0,R), "iteração.gas.sp" = rep(0,R), h.gas.sp = rep(0,R), h2.gas.sp = rep(0,R)))
paramEstimados_SPGAS_t <- t(data.frame(omega.gas.st = rep(0,R) ,alpha.gas.st = rep(0,R), beta.gas.st = rep(0,R), mu.gas.st = rep(0,R), loglike.gas.st = rep(0,R), "iteração" = rep(0,R), h = rep(0,R), h2 = rep(0,R)))


for(t in 1:ndraws){
  vy <- as.numeric(garchSim(model, n = horizon.simul))
  
  #Gas Paramétrico - Normal
  df=0
  gas <- optim.gas(vy, df)
  beta.til <- gas$optimizer.new$par[3] - gas$optimizer.new$par[2]
  omega <- gas$optimizer.new$par[1]/(1 - gas$optimizer.new$par[3])
  
  #Gas Semi-Paramétrico
  h = 0.5
  h2 = 0.5
  sp_gas = sp.gas(vy, h, h2, horizon = length(vy), gas)
  beta.sp.til <- sp_gas$optimizer.sp$par[3] - sp_gas$optimizer.sp$par[2]
  omega.sp <- sp_gas$optimizer.sp$par[1]/(1 - sp_gas$optimizer.sp$par[3])
  
  #Gas Semi-Paramétrico - t-student
  h = 0.5
  h2 = 0.5
  sp_gas_t = sp.gas_2(vy, h, h2, horizon = length(vy), gas)
  beta.sp.til_t <- sp_gas_t$optimizer.sp$par[3] - sp_gas_t$optimizer.sp$par[2]
  omega.sp_t <- sp_gas_t$optimizer.sp$par[1]/(1 - sp_gas_t$optimizer.sp$par[3])
  
  #GARCH  - N
  garch <- garchFit(data=vy, cond.dist="norm", trace=F)
  
  #GARCH  - T
  garch.t <- garchFit(data=vy, cond.dist="std", trace=F)
  
  #Gas Paramétrico - T-Student
  df=1
  gas.t <- optim.gas(vy, df)
  beta.t.til <-  gas.t$optimizer.new$par[3] - gas.t$optimizer.new$par[2]
  omega.t <- gas.t$optimizer.new$par[1]/(1 - gas.t$optimizer.new$par[3])
  
  paramEstimados_GARCH[,t] <- unlist(c(garch@fit$par,-garch@fit$llh))
  paramEstimados_GARCH.t[,t] <- unlist(c(garch.t@fit$par,-garch.t@fit$llh))
  paramEstimados_GAS[,t] <- unlist(c(omega, gas$optimizer.new$par[2], beta.til, gas$optimizer.new$par[4],-length(vy)*gas$optimizer.new$objective))
  paramEstimados_GAS.t[,t] <- unlist(c(omega.t, gas.t$optimizer.new$par[2], beta.t.til, gas.t$optimizer.new$par[4], gas.t$optimizer.new$par[5],-length(vy)*gas.t$optimizer.new$objective))
  paramEstimados_SPGAS[,t] <- unlist(c(omega.sp, sp_gas$optimizer.sp$par[2], beta.sp.til ,sp_gas$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas$optimizer.sp[2]), sp_gas$i, h, h2))
  paramEstimados_SPGAS_t[,t] <- unlist(c(omega.sp_t, sp_gas_t$optimizer.sp$par[2], beta.sp.til_t ,sp_gas_t$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas_t$optimizer.sp[2]), sp_gas_t$i, h, h2))
}

