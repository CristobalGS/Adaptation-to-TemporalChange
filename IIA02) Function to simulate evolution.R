# Function to simulate evolution under each environment

#----------------------------------------  Equations --------------------------------------------------------------

# Trait value (phenotype) under an arbitrary reference environment (e = 0)
# t = x + b*e + y
# t: trait
# x: additive genetic breeding value
# b: plasticity (environmental effect, could be adaptive, neutral, or maladaptive)
# e: environment
# y: residual variation
# x, b, and y are all independent normal random variables

# Fitness is a Gaussian function of trait expression
# w = wmax*exp(-((theta - t)^2)/(2*w_t^2) - (b^2)/(2*w_b^2))
# wmax: maximum individual fitness
# theta: optimal trait expression in the current generation
# 1/w_t: strength of stabilizing selection on the trait (wt is the width of fitness function)
# 1/w_b: cost of plasticity modeled as strength of stabilizing selection acting against the reaction norm slope, favoring no plasticity

# Change in b and var(b) after one generation of selection
# b.new <- (b*w_b^2)/(w_b^2 + bvar)
# bvar.new <- (bvar*w_b^2)/(w_b^2 + bvar)

# Change in t and var(t) after one generation of selection
# t.new <- x + b.new*e
# tvar.new <- xvar + bvar.new*e^2 + yvar
#

# Mean lifetime fitness after one generaiton of selection
# w <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar.new)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - t.new)^2)/(2*w_t^2 + 2*tvar.new))
#

# The optimum phenotype (theta) is a linear function of the environment (Lande 2014)
# theta <- A + B*e
# A: optimum phenotype in the average o reference environment (e = 0) => A = 0
# B: slope of the optimum phenotype as a function of environment => B = 1
# -----------------------------------------------------------------------------------------------------

# Function to simulate evolution
sim.evo <- function(gt = 11, env, time = 1:10000, plot = TRUE){
  
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
  
  # Fitness derivatives
  w <- expression(log(wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + (xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar))))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - (x + ((b*w_b^2)/(w_b^2 + bvar))*e1))^2)/(2*w_t^2 + 2*(xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar)))))
  wx <- D(w, "x")
  wb <- D(w, "b")
  
  # Selection events based on generation time
  timesteps <- length(time)
  selevents <- trunc(timesteps/gt)
  
  # Simulation
  xvec <- c() 
  bvec <- c()
  for(i in 1:selevents){
    e1 <- env[gt*i - (gt - 1)] #environment of development (determines trait expression)
    e2 <- env[gt*i] #environment where selection occurs (determines theta)
    
    theta <- A + B*e2 #theta = 0 + 1*e = e2
    
    delta.x <- eval(wx)*xvar
    delta.b <- eval(wb)*bvar*h_b
    
    x <- x + delta.x
    b <- b + delta.b
    xvec <- c(xvec, x)
    bvec <- c(bvec, b)
  }
  
  if(plot == TRUE){
    plot(time[1:length(xvec)], xvec, type = "l", xlab = "Generations", ylab = expression(bar(x)))
    plot(time[1:length(bvec)], bvec, type = "l", xlab = "Generations", ylab = expression(bar(b))) 
  }
  return(cbind(x = xvec, b = bvec))
}
