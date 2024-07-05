# Function to plot tolerance curve at equilibrium with variation in evolved plasticity (breadth)

plot.tolcurve <- function(xvec, bvec, reacNorm = FALSE, out = TRUE, env.i = -10, env.f = 25, equilibrium = 401, plot = TRUE){
  
  # Parameters (for plotting)
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
  
  # Tolerance curve at equilibrium (from generation 401 onwards)
  start <- equilibrium # 801 for trend, 401 for all others
  end <- length(xvec) #same as length(bvec)
  x <- mean(xvec[start:end])
  b <- mean(bvec[start:end])
  x
  b
  seq.ini <- env.i #-10 / -5
  seq.end <- env.f # 25 / 15
  environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
  trait <- x + b*environment + y
  theta <- environment
  fitness <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))

  # Max plasticity
  bmax <- max(bvec[start:end])
  x.bmax <- xvec[start:end][which.max(bvec[start:end])]
  environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
  trait <- x.bmax + bmax*environment + y
  theta <- environment
  fitness.bmax <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
  
  # Min plasticity
  bmin <- min(bvec[start:end])
  x.bmin <- xvec[start:end][which.min(bvec[start:end])]
  environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
  trait <- x.bmin + bmin*environment + y
  theta <- environment
  fitness.bmin <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
  
  # Plot
  height <- max(fitness)
  center <- environment[which.max(fitness)]
  center.bmax <- environment[which.max(fitness.bmax)]
  center.bmin <- environment[which.max(fitness.bmin)]
  width.bmax <- abs(environment[which.min(abs(fitness.bmax-(height/2)))] - center.bmax)/1.1775
  width.bmin <- abs(environment[which.min(abs(fitness.bmin-(height/2)))] - center.bmin)/1.1775
  f.bmax <- height*exp(-((environment-center)^2)/(2*width.bmax^2))
  f.bmin <- height*exp(-((environment-center)^2)/(2*width.bmin^2))
  if(plot == TRUE){
    par(mar=c(4.3,4.3,0.5,0.5))
    plot(environment, fitness, type = "l", xlab = "Environment", ylab = "Fitness", lwd = 4) 
    segments(x0 = environment[which.max(fitness)], y0 = -1, x1 = environment[which.max(fitness)], y1 = max(fitness), lty = 2)
    polygon(c(environment, rev(environment)), c(f.bmin, rev(f.bmax)), col = grey(0.7, alpha = 0.3), border = NA)
    
  }

  # Also plot reaction norm with variation in evolved plasticity if so desired
  if(reacNorm == TRUE){
    par(mar=c(4.3,4.3,0.5,0.5))
    trait <- x + b*environment + y
    trait.bmax <- x + bmax*environment + y
    trait.bmin <- x + bmin*environment + y
    plot(environment, trait, type = "l", xlab = "Environment", ylab = "Trait", ylim = c(0, 10), lwd = 4)
    polygon(c(environment, rev(environment)), c(trait.bmin, rev(trait.bmax)), col = grey(0.7, alpha = 0.3), border = NA)
  }
  
  # Return outcomes for further plotting if needed
  if(out == TRUE){
    return(cbind(fitness, f.bmax, f.bmin))
  }
}

