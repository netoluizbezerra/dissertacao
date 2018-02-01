choose.bandwidth = function(gas, horizon, vy, std.res, b0){
  vy <- as.numeric(vy)
  list[domega, vA, vB, dmu, ddf] <- param(gas$optimizer.new)
  cT = length(vy)
  dloglik = 0
  q.store = as.matrix(rep(0,cT))
  q.line = 0
  nabla = as.matrix(rep(0,cT))
  f_t = rep(0,(horizon+1)); vscore = rep(0,(horizon+1))
  f_t[1]  = domega/(1-vB)
  h.grid = cbind(seq(from = 0.2,to = 0.9,by = 0.05),seq(from = 0.2,to = 0.9,by = 0.05))
  likelihood = matrix(nrow = nrow(h.grid), ncol = nrow(h.grid), 0)
  for(i in 1:nrow(h.grid)){
    for(j in 1:nrow(h.grid)){
      h = h.grid[i,1]
      h2 = h.grid[j,2]
      for(t in 1:horizon){
        arg <- (vy[t] - dmu)/exp(f_t[t]/2)
        q.store[t] <- dens.np(x = std.res, x.i = arg, h = h, k = normal)
        q.line <- -ddens.np(x = std.res, x.i = arg, h = h2, k = d.normal)
        nabla[t] <- q.line/q.store[t]
        vscore[t] <-  -1/2 * (nabla[t]) * arg  - 1/2
        f_t[t+1] <- domega + 2*vscore[t]*vA + f_t[t]*vB
      }
      dloglik = 1/cT*(-1/2*(sum(f_t)) + sum(log(q.store)))
      if(is.finite(dloglik) == "FALSE" | is.complex(dloglik) == "TRUE"){
        dloglik <- -b0
      }
      likelihood[i,j] <- dloglik
    }
  }
  stack.grid <- expand.grid(seq(from = 0.2,to = 0.9,by = 0.05), seq(from = 0.2,to = 0.9,by = 0.05))
  like <- vec(likelihood)
  stack.like = data.frame(cbind(stack.grid, like))
  result <- stack.like %>% 
    filter(like == max(like))
  return(result)}


