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
model = garchSpec(model = list(omega = omega, alpha = alpha1, beta = beta1, skew = 2, shape = 5), cond.dist = "sstd")
R <- ndraws #Numero de repeticoes do estudo 

paramEstimados_GARCH.t <-  t(data.frame(mu.garch.t = rep(0,R) ,omega.garch.t = rep(0,R), alpha.garch.t = rep(0,R), beta.garch.t = rep(0,R),  ddf.garch.t = rep(0,R),loglike.garch.t = rep(0,R)))
paramEstimados_GARCH.st <-  t(data.frame(mu.garch.st = rep(0,R) ,omega.garch.st = rep(0,R), alpha.garch.st = rep(0,R), beta.garch.st = rep(0,R), skew.garch.st = rep(0,R), ddf.garch.st = rep(0,R),loglike.garch.st = rep(0,R)))
paramEstimados_GAS <-  t(data.frame(omega.gas = rep(0,R) ,alpha.gas = rep(0,R), beta.gas = rep(0,R), mu.gas = rep(0,R), loglike.gas = rep(0,R)))
paramEstimados_GAS.t <-  t(data.frame(omega.gas.t = rep(0,R) ,alpha.gas.t = rep(0,R), beta.gas.t = rep(0,R), mu.gas.t = rep(0,R), ddf.gas.t = rep(0,R), loglike.gas.t = rep(0,R)))
paramEstimados_SPGAS <- t(data.frame(omega.gas.st = rep(0,R) ,alpha.gas.st = rep(0,R), beta.gas.st = rep(0,R), mu.gas.st = rep(0,R), loglike.gas.st = rep(0,R), "iteração" = rep(0,R), h = rep(0,R), h2 = rep(0,R)))

for(t in 1:ndraws){
  print(paste("Procedimento está na iteração", t))
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
  
  #GARCH  - ST
  garch.st <- garchFit(data=vy, cond.dist="sstd", trace=F)
  
  #GARCH  - T
  garch.t <- garchFit(data=vy, cond.dist="std", trace=F)
  
  #Gas Paramétrico - T-Student
  df=1
  gas.t <- optim.gas(vy, df)
  beta.t.til <-  gas.t$optimizer.new$par[3] - gas.t$optimizer.new$par[2]
  omega.t <- gas.t$optimizer.new$par[1]/(1 - gas.t$optimizer.new$par[3])

  paramEstimados_GARCH.st[,t] <- unlist(c(garch.st@fit$par,-garch.st@fit$llh))
  paramEstimados_GARCH.t[,t] <- unlist(c(garch.t@fit$par,-garch.t@fit$llh))
  paramEstimados_GAS[,t] <- unlist(c(omega, gas$optimizer.new$par[2], beta.til, gas$optimizer.new$par[4],-length(vy)*gas$optimizer.new$objective))
  paramEstimados_GAS.t[,t] <- unlist(c(omega.t, gas.t$optimizer.new$par[2], beta.t.til, gas.t$optimizer.new$par[4], gas.t$optimizer.new$par[5],-length(vy)*gas.t$optimizer.new$objective))
  paramEstimados_SPGAS[,t] <- unlist(c(omega.sp, sp_gas$optimizer.sp$par[2], beta.sp.til ,sp_gas$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas$optimizer.sp[2]), sp_gas$i, h, h2))
}





parametros.est <- rbind(paramEstimados_GARCH.st,paramEstimados_GARCH.t,paramEstimados_GAS,paramEstimados_GAS.t,paramEstimados_SPGAS)
parametros.est <- as.data.frame(parametros.est)
write.xlsx(parametros.est, file = "1_dgp_GARCH_st.xlsx", row.names = TRUE)

names <- c("mu.garch.t","omega.garch.t","alpha.garch.t","beta.garch.t","ddf.garch.t", "loglike.garch.t","mu.garch.st","omega.garch.st","alpha.garch.st", "beta.garch.st","ddf.garch.st","skew.garch.st","loglike.garch.st", 
           "omega.gas","alpha.gas","beta.gas","mu.gas", "loglike.gas","omega.gas.t","alpha.gas.t","beta.gas.t", "mu.gas.t","ddf.gas.t","loglike.gas.t",
           "omega.gas.st", "alpha.gas.st","beta.gas.st","mu.gas.st","loglike.gas.st", "it","h","h2")


