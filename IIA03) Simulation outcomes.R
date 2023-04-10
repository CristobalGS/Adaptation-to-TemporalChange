# IIA 3) Function to plot mean population tolerance curve at equilibrium (from generation 401 onwards)
evo.out <- function(xvec, bvec, plot = TRUE, equilibrium = 401){
  
  # Initial parameters
  e <- 0
  e1 <- 0
  e2 <- 0
  b <- 1
  bvar <- 0.1
  w_b <- 1
  x <- 0
  xvar <- 1
  y <- 0
  yvar <- 0.1
  w_t <- 1
  wmax <- 1
  A <- 0
  B <- 1
  theta <- A + B*e
  tvar <- xvar + bvar*e^2 + yvar
  h_b <- 0.5
  
  start <- equilibrium # 801 # 401
  end <- length(xvec) #same as length(bvec)
  x <- mean(xvec[start:end])
  b <- mean(bvec[start:end])
  seq.ini <- -10 #-10 / -5
  seq.end <- 25 # 25 / 15
  environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
  trait <- x + b*environment + y
  theta <- environment
  fitness <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
  if(plot == TRUE){
    plot(environment, fitness, type = "l", xlab = "Environment", ylab = "Fitness") #tolerance curve
    abline(v = environment[which.max(fitness)], col = "blue") 
  }
  print(c(paste0("Mean breeding value = ", x), paste0("Mean plasticity = ", b)))
  return(list(x, b, cbind(env = environment, fitness = fitness, trait = trait)))
}
