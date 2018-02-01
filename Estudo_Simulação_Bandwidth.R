####################
# DGP - GARCH - T ##
####################
rm(list=ls())

ndraws=100
horizon.simul=2000
alpha1=0.10
beta1 =0.85
mu=0.00
omega= 2
burnin = 101
model = garchSpec(model = list(omega = omega, alpha = alpha1, beta = beta1, skew = 2, shape = 5), cond.dist = "sstd")
R <- ndraws #Numero de repeticoes do estudo 

paramEstimados_GARCH.st <-  t(data.frame(mu.garch.st = rep(0,R) ,omega.garch.st = rep(0,R), alpha.garch.st = rep(0,R), beta.garch.st = rep(0,R), skew.garch.st = rep(0,R), ddf.garch.st = rep(0,R),loglike.garch.st = rep(0,R)))
paramEstimados_GAS <-  t(data.frame(omega.gas = rep(0,R) ,alpha.gas = rep(0,R), beta.gas = rep(0,R), mu.gas = rep(0,R), loglike.gas = rep(0,R)))
paramEstimados_SPGAS_priori <- t(data.frame(omega.gas.priori = rep(0,R) ,alpha.gas.priori = rep(0,R), beta.gas.priori = rep(0,R), mu.gas.priori = rep(0,R), loglike.gas.priori = rep(0,R), "iteração" = rep(0,R), h = rep(0,R), h2 = rep(0,R)))
paramEstimados_SPGAS_thumbrule <- t(data.frame(b0 = rep(0,R), b1 = rep(0,R), omega.gas.thumbrule = rep(0,R) ,alpha.gas.thumbrule = rep(0,R), beta.gas.thumbrule = rep(0,R), mu.gas.thumbrule = rep(0,R), loglike.gas.thumbrule = rep(0,R), "iteração" = rep(0,R)))
paramEstimados_SPGAS_loglike <- t(data.frame(b0_log = rep(0,R), b1_log = rep(0,R), omega.gas.loglike = rep(0,R) ,alpha.gas.loglike = rep(0,R), beta.gas.loglike = rep(0,R), mu.gas.loglike = rep(0,R), loglike.gas.loglike = rep(0,R), "iteração" = rep(0,R)))

for(t in 1:ndraws){
  print(paste("Estamos na iteração", t))
  vy <- as.numeric(garchSim(model, n = horizon.simul))
  
  #Gas Paramétrico - Normal
  df=0
  gas <- optim.gas(vy, df)
  beta.til <- gas$optimizer.new$par[3] - gas$optimizer.new$par[2]
  omega <- gas$optimizer.new$par[1]/(1 - gas$optimizer.new$par[3])
  
  #Obtenção dos resíduos padronizados: 
  list[vf.temp, sigma] = get.ft(vy, horizon = (length(vy)-1), gas)
  list[domega, vA, vB, dmu, ddf] = param(gas$optimizer.new)
  vp = c(domega, vA, vB, dmu)
  std.res <- (vy - dmu)/exp(vf.temp/2)
  
  #Gas Semi-Paramétrico Priori
  h = 0.5
  h2 = 0.5
  sp_gas_1 = sp.gas(vy, h, h2, horizon = length(vy), gas)
  beta.sp.til_1 <- sp_gas_1$optimizer.sp$par[3] - sp_gas_1$optimizer.sp$par[2]
  omega.sp_1 <- sp_gas_1$optimizer.sp$par[1]/(1 - sp_gas_1$optimizer.sp$par[3])
  
  #Gas Semi-Paramétrico Thumbrule
  band <- bw.nrd0(std.res)
  band1 <- h.bcv(std.res, deriv.order = 1)
  sp_gas_2 = sp.gas(vy, h = band, h2 = band1$h, horizon = length(vy), gas)
  beta.sp.til_2 <- sp_gas_2$optimizer.sp$par[3] - sp_gas_2$optimizer.sp$par[2]
  omega.sp_2 <- sp_gas_2$optimizer.sp$par[1]/(1 - sp_gas_2$optimizer.sp$par[3])
  
  #Gas Semi-Paramétrico Loglikehood
  vP = c(0.5,0.5) 
  optimizer.sp.band = nlminb(start = vP, objective = sp.gas.likelihood.bandwidth, 
                             gradient = NULL, hessian = NULL, 
                             vy = vy, std.res = std.res, gas,
                             control = list(trace = 1), 
                             lower = 0.01, upper = 0.9)
  
  h.livre1 = optimizer.sp.band$par[1]
  h.livre2 = optimizer.sp.band$par[2]
  sp_gas_3 = sp.gas(vy, h = h.livre1, h2 = h.livre2, horizon = length(vy), gas)
  beta.sp.til_3 <- sp_gas_3$optimizer.sp$par[3] - sp_gas_3$optimizer.sp$par[2]
  omega.sp_3 <- sp_gas_3$optimizer.sp$par[1]/(1 - sp_gas_3$optimizer.sp$par[3])
  
  #GARCH  - ST
  garch.st <- garchFit(data=vy, cond.dist="sstd", trace=F)

  paramEstimados_GARCH.st[,t] <- unlist(c(garch.st@fit$par,-garch.st@fit$llh))
  paramEstimados_GAS[,t] <- unlist(c(omega, gas$optimizer.new$par[2], beta.til, gas$optimizer.new$par[4],-length(vy)*gas$optimizer.new$objective))
  paramEstimados_SPGAS_priori[,t] <- unlist(c(omega.sp_1, sp_gas_1$optimizer.sp$par[2], beta.sp.til_1 ,sp_gas_1$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas_1$optimizer.sp[2]), sp_gas_1$i, h, h2))
  paramEstimados_SPGAS_thumbrule[,t] <- unlist(c(band, band1$h, omega.sp_2, sp_gas_2$optimizer.sp$par[2], beta.sp.til_2 ,sp_gas_2$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas_2$optimizer.sp[2]), sp_gas_2$i))
  paramEstimados_SPGAS_loglike[,t] <- unlist(c(h.livre1, h.livre2, omega.sp_3, sp_gas_3$optimizer.sp$par[2], beta.sp.til_3,sp_gas_3$optimizer.sp$par[4] ,-length(vy)*unlist(sp_gas_3$optimizer.sp[2]), sp_gas_3$i))
}


