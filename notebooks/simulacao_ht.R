randhelp <- function(
  horizon,
  h0 = 2e-4, 
  mu = 0, omega=2,
  alpha1 = 0.1,
  beta1  = 0.85
){
  ret <- zt <- et <- ht <- rep(0, horizon)
  ht[1] <- h0
  for(j in 1:horizon){
    zt[j]   <- rnorm(1,0,1)
    #zt[j]   <- rsged(1, mean = 0, sd = 1, nu = 2, xi = 1.5)
    et[j]   <- zt[j]*sqrt(ht[j])
    ret[j]  <- mu + et[j]
    if( j < horizon )
      ht[j+1] <- omega + alpha1*et[j]^2 + beta1*ht[j]
  }
  list(ret = ret, ht = ht)
}

ndraws = 10 

for(i in 1:ndraws){
burnin = 100
# Obtain simulated returns and volatility
ret <- randhelp(horizon = 3100)

# Setting in sample window 
vy.in.sample = ret$ret[burnin:1500]
vy.full = ret$ret[burnin:length(ret$ret)]
sigma_true = ret$ht[burnin:length(ret$ht)]

#Estimating GAS models: 
df = 0 
h = 0.5; h2 = 0.5
gas <- optim.gas(vy.in.sample, df)
sp_gas_n = sp.gas(vy.in.sample, h, h2, horizon = length(vy.in.sample), gas)
sp_gas_t = sp.gas_2(vy.in.sample, h, h2, horizon = length(vy.in.sample), gas)


# Get ht GAS
dmu = gas$optimizer.new$par[4]
list_gas_n = get.ft(vy.full, horizon = length(vy.full), gas)
vf.temp = list_gas_n$f_t; sigma_gas_n = list_gas_n$sigma[1:(length(list_gas_n$sigma)-1)]

# Get sigma GAS-N
vf.temp = vf.temp[1:(length(vf.temp)-1)]
std.res <- (vy.in.sample - dmu)/exp(vf.temp[1:length(vy.in.sample)]/2)
list_sp_n = get.ft.sp(vy.full, h = h, h2 = h2, horizon = length(vy.full), optimizer.sp = sp_gas_n$optimizer.sp, std.res = std.res, kernel = "gauss")
list_sp_t = get.ft.sp(vy.full, h = h, h2 = h2, horizon = length(vy.full), optimizer.sp = sp_gas_t$optimizer.sp, std.res = std.res, kernel = "t")
sigma_gas_sp_n = list_sp_n$sigma[1:(length(list_sp_n$sigma)-1)]
sigma_gas_sp_t = list_sp_t$sigma[1:(length(list_sp_t$sigma)-1)]

#Gas Paramétrico - T-Student
df=1
gas.t <- optim.gas(vy.in.sample, df)

# Get sigma GAS-t
dmu = gas.t$optimizer.new$par[4]
list_gas_t = get.ft(vy.full, horizon = length(vy.full), gas.t)
sigma_gas_t = list_gas_t$sigma[1:(length(list_gas_t$sigma)-1)]

#GARCH  - N
garch11.n.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "norm")

garch.fit_n = ugarchfit(spec = garch11.n.spec, data = vy.full, out.sample = 1500)
garch.fit_n@fit$sigma
garch.fit_n.forecast = ugarchforecast(garch.fit_n, n.roll=1500, n.ahead=1)

#GARCH  - T
garch11.t.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "std")

garch.fit_t = ugarchfit(spec = garch11.t.spec, data = vy.full, out.sample = 1500)
garch.fit_t@fit$sigma
garch.fit_t.forecast = ugarchforecast(garch.fit_t, n.roll=1500, n.ahead=1)

#GARCH  - ST
garch11.sstd.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "sstd")

garch.fit_sstd = ugarchfit(spec = garch11.sstd.spec, data = vy.full, out.sample = 1500)
garch.fit_sstd@fit$sigma
garch.fit_sstdforecast = ugarchforecast(garch.fit_sstd, n.roll=1500, n.ahead=1)

#GARCH  - SGED
garch11.sged.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "sged")

garch.fit_sged = ugarchfit(spec = garch11.sged.spec, data = vy.full, out.sample = 1500)
garch.fit_sged@fit$sigma
garch.fit_sged.forecast = ugarchforecast(garch.fit_sged, n.roll=1500, n.ahead=1)



in.sample = rbind(sqrt(sigma_true[1:1501]),sqrt(sigma_gas_n[1:1501]), sqrt(sigma_gas_t[1:1501]), sqrt(sigma_gas_sp_n[1:1501]), sqrt(sigma_gas_sp_t[1:1501]), 
                  garch.fit_n@fit$sigma, garch.fit_t@fit$sigma, garch.fit_sstd@fit$sigma, garch.fit_sged@fit$sigma)
                  
out.sample =  rbind(sqrt(sigma_true[1501:3001]),sqrt(sigma_gas_n[1501:3001]), sqrt(sigma_gas_t[1501:3001]), sqrt(sigma_gas_sp_n[1501:3001]), sqrt(sigma_gas_sp_t[1501:3001]),
                    as.numeric(garch.fit_n.forecast@forecast$sigmaFor), as.numeric(garch.fit_t.forecast@forecast$sigmaFor), as.numeric(garch.fit_sstdforecast@forecast$sigmaFor),
                    as.numeric(garch.fit_sged.forecast@forecast$sigmaFor))

if(i == 1){
  in.sample.1 = in.sample
  out.sample.1 = out.sample
} else if(i == 2){
  in.sample.array = abind(in.sample.1, in.sample, along = 3)
  out.sample.array = abind(out.sample.1, out.sample, along = 3)
} else{
  in.sample.array = abind(in.sample.array, in.sample, along = 3)
  out.sample.array = abind(out.sample.array, in.sample, along = 3)
}} 





