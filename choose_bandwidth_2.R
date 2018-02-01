sp.gas.likelihood.bandwidth = function(vp, vy, std.res, gas){
  vy <- as.numeric(vy)
  list[domega, vA, vB, dmu, ddf] <- param(gas$optimizer.new)
  cT = length(vy)
  dloglik = 0
  vf = rep(0, cT); vscore = rep(0, cT)
  q.store = as.matrix(rep(0,cT))
  q.line = 0
  ft = as.matrix(rep(0,cT))
  nabla = as.matrix(rep(0,cT)) 
  vf[1] = domega/(1-vB)  
  h = vp[1]
  h2 = vp[2]
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



