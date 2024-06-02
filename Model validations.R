# Reality checks

time <- 1:10000

# Function to set up environment
setup.env <- function(e0, nt, nf, nn, P, time, rho){
  set.seed(2)
  env <- c()
  gamma <- rnorm(length(time), mean = 0, sd = 1)
  noise <- c(gamma[1])
  for(i in time){
    noisei <- rho*noise[i] + gamma[i]*sqrt(1-rho^2)
    noise <- c(noise, noisei)
  }
  noise <- noise[time]
  for (t in time){
    e_t <- e0 + nt*t + nf*sin(2*pi*t/P) + nn*noise[t]
    env <- c(env, e_t)
  }
  plot(time, env, type = "l", xlab = "Time", ylab = "Environment")
  return(env)
}

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

# Function to run simulation
run.simulation <- function(gt, env){
  
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
  
  #xvec <- rep(xvec, each = gt)
  #bvec <- rep(bvec, each = gt)
  
  plot(time[1:length(xvec)], xvec, type = "l", xlab = "Generations", ylab = expression(bar(x)))
  plot(time[1:length(bvec)], bvec, type = "l", xlab = "Generations", ylab = expression(bar(b)))
  return(cbind(x = xvec, b = bvec))
}

#Function to plot mean population tolerance curve at equilibrium (from generation 401 onwards)
tol.curve <- function(xvec, bvec){
  start <- 800 # 401
  end <- length(xvec) #same as length(bvec)
  x <- mean(xvec[start:end])
  b <- mean(bvec[start:end])
  x
  b
  seq.ini <- -10 #-10 / -5
  seq.end <- 25 # 25 / 15
  environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
  trait <- x + b*environment + y
  theta <- environment
  fitness <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
  plot(environment, fitness, type = "l", xlab = "Environment", ylab = "Fitness") #tolerance curve
  #abline(v = 5, col = "blue")
  abline(v = environment[which.max(fitness)], col = "blue")
}

#####################################################################################

# Lande and Shannon 1996 #
# Initial parameters
e <- 0
e1 <- 0
e2 <- 0
b <- 0
bvar <- 0
w_b <- 1
x <- 0
xvar <- 0.001
y <- 0
yvar <- 0.1
w_t <- 1
wmax <- 1
A <- 0
B <- 1
theta <- A + B*e
tvar <- xvar + bvar*e^2 + yvar
h_b <- 0.5


# 1) Deterministic change

# 1a) Constant environment: mean trait evolves to optimum
env <- setup.env(e0 = 10, nt = 0, nf = 0, nn = 0, P = 10, time = time, rho = 0)
xb <- run.simulation(gt = 1, env = env)
plot(time, env, type = "l", ylim = c(0, 12))
lines(time, xb[ , 1], col = "blue")
# Check!

# 1b) Directional change of environment: mean trait tracks optimum with constant lag at equilibrium
# lag = (nt*w_t)/xvar
env <- setup.env(e0 = 0, nt = 0.001, nf = 0, nn = 0, P = 10, time = time, rho = 0)
xb <- run.simulation(gt = 1, env = env)
(0.001*w_t)/xvar #expected lag at equilibrium: 1 
plot(time, env, type = "l")
lines(time, xb[ , 1], col = "blue")
mean(env[4000:10000] - xb[ , 1][4000:10000]) #observed lag at equilibrium: 1.094697 
# Check!

# 1c) Cyclic environment: 
env <- setup.env(e0 = 10, nt = 0, nf = 1, nn = 0, P = 1000, time = time, rho = 0)
xb <- run.simulation(gt = 1, env = env)
xvec <- xb[ , 1] 
plot(time, env, type = "l", ylim = c(0, 12))
lines(time, xb[ , 1], col = "blue")
sum(((xvec[6000:6001] - env[6000:6001])^2)/(2*w_t))/(2*pi/1000) #obs load
((1000^2)/(4*w_t))/((1000^2) + (xvar^2)/(w_t^2)) #exp load
# ?

# 2) Stochastic fluctuations

# 2a) Random (uncorrelated) change
env <- setup.env(e0 = 10, nt = 0, nf = 0, nn = 1, P = 1000, time = time, rho = 0)
xb <- run.simulation(gt = 1, env = env)
plot(time, env, type = "l", ylim = c(0, 14))
lines(time, xb[ , 1], col = "blue")
# ?

# 2a) Random (positively autocorrelated) change
env <- setup.env(e0 = 10, nt = 0, nf = 0, nn = 1, P = 1000, time = time, rho = 0.9)
xb <- run.simulation(gt = 1, env = env)
plot(time, env, type = "l", ylim = c(0, 14))
lines(time, xb[ , 1], col = "blue")
# ?


# Lande 2009 #

# Initial parameters
time <- 1:50000
e <- 0
e1 <- 0
e2 <- 0
b <- 0
bvar <- 0.1
w_b <- 10
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

env <- setup.env(e0 = 10, nt = 0, nf = 0, nn = 0.5, P = 5, time = time, rho = 0.9)
xb <- run.simulation(gt = 2, env = env)
bvec <- xb[ , 2]
mean(bvec[20000:25000])

## Same dynamics as in Lande 2009: rapid adaptation to sudden change in the environment by evolution of plasticity first, 
# which then slowly goes down to final equilibrium value (>0), with mean breeding value increasing as plasticity decreases.
# Optimal plasticity at equilibrium increases with correlation between environment of development and selection. 




