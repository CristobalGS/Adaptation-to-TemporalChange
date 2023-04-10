# IIA 1) Function to simulate environments

# Environment
# e = e0 + nt*t + nc*sin(2*pi*t/P) + nn*y
# e0: initial mean environment (intercept)
# nt: scaling constant reflecting relative importance of temporal trend (slope)
# nc: scaling constant reflecting relative importance of cyclic change
# nn: scaling constant reflecting relative importance of noise
# P: period of fluctuation (2*pi/P = constant)
# t: time
# gamma: random noise

sim.env <- function(e0 = 5, nt = 0, nc = 1, nn = 1, P = 10, time = 1:10000, rho = 0, plot = TRUE){
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
    e_t <- e0 + nt*t + nc*sin(2*pi*t/P) + nn*noise[t]
    env <- c(env, e_t)
  }
  if(plot == TRUE){
    plot(time, env, type = "l", xlab = "Time", ylab = "Environment") 
  }
  return(env)
}
