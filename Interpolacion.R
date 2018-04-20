#vecindades <- function(i,j){}

#Full conditional para las v_i
v_i <- function(y_i, u_i, old_vi, lambda, c_i){
  cl <- makeCluster(detectCores())
  regiterDoParallel(cl)
  v_is <- foreach(i = 1:length(v_i), .combine = 'c') %dopar% {
  prop <- rnorm(1, old_vi[i], sigma_sqd)
  prob_num <- (-prop^2/(2*lambda))  - (c_i[i]*exp(u_i[i]+prop)) + (prop*y_i[i])
  prob_denom <-  (-old_vi[i]^2/(2*lambda))  - (c_i[i]*exp(u_i[i]+old_v_i[i])) + (old_vi[i]*y_i[i])
  prob <- exp(prob_num - prob-denom)
  prob <- min(1,prob)
  if(runif(1) < prob)
    prop
  else
    old_vi[i]
  }
  stopCluster(cl)
  return(v_is)
}

#Full condtional para las u_i
u_i <- function(y_i, v_i, old_ui, k, c_i, indices){
  cl <- makeCluster(detectCores())
  regiterDoParallel(cl)
  u_is <- foreach(i = 1:length(u_i), .combine = 'c') %dopar% {
    prop <- rnorm(1, old_ui[i], sigma_sqd)
    prob_num <- (-1/(2*k) * (prop-old_ui[indices[[i]] ])^2)  - 
      (c_i[i]*exp(v_i[i]+prop)) + (prop*y_i[i])
    
    prob_denom <- (-1/(2*k) * (old_ui[i]-old_ui[indices[[i]] ])^2)  - 
      (c_i[i]*exp(v_i[i]+old_u_i[i])) + (old_ui[i]*y_i[i])
    prob <- exp(prob_num - prob_denom)
    prob <- min(1,prob)
    if(runif(1) < prob)
      prop
    else
      old_ui[i]
  }
  stopCluster(cl)
  return(u_is)
}

#Full conditional para k
full_k <- function(epsilon, u_i, indices){
  alfa <- length(u_i)/2 - 1
  suma <- sum(sapply(1:length(u_i), function(i) (u_i[i] - u_i[indices[[i]] ])^2))
  beta <- (suma - epsilon)/2
  rinvgamma(1, alfa, scale=beta)
}

#Full conditional para lambda
full_lambda <- function(epsilon, v_i){
  alfa <- length(v_i)/2 - 1
  beta <- (sum(v_i^2)-epsilon)/2
  rinvgamma(1,alfa, scale=beta)
}

#Funcion principal desde la cual se corren las iteraciones
mcmc <- function(nsim, y_i, epsilon, sigma_sqd, c_i, init){
  library(foreach); library(parallel); library(invgamma)
  indices <- vecindades(y_i)
  v <- matrix(nrow=nrow(y_i), ncol=nsim)
  v[,1] <- init[['v']]
  
  u <- matrix(nrow=nrow(y_i), ncol=nsim)
  u[,1] <- init[['u']]
  
  lambda <- vector('numeric', nsim)
  lambda[1] <- init[['lambda']]
  
  k <- vector('numeric', nsim)
  k[1] <- init[['k']]
  
  for(i in 2:nsim){
    v <- v_i(y[i-1], u[[i-1]], v[[i-1]], lambda[i-1], c_i)
    lambda <- full_lambda(epsilon, v[[i]])
    u <- u_i(y[i-1], v[[i]], u[i-1]], k[i-1], c_i, indices)
    k <- full_k(epsilon, u[[i]], indices)
  }
  return(list(v,lambda,u,k))
}