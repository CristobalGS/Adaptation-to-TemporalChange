# Simulate evolution under each environmental change scenario and plot outcomes

# Data frame to store simulated tolerance curves (to plot baseline and change together)
curves <- data.frame(Environment = seq(from = -10, to = 25, length.out = 1000))

# Get baseline environment and tolerance curve for plotting
env.base <- sim.env(plot = FALSE)
xb <- sim.evo(env = env.base, plot = FALSE)
out <- evo.out(xvec = xb[ , 1], bvec = xb[ , 2], plot = FALSE)
base <- plot.tolcurve(xvec = xb[ , 1], bvec = xb[ , 2], reacNorm = FALSE, out = TRUE, env.i = -10, env.f = 25, equilibrium = 401)
base <- as.data.frame(base)
curves$base <- base$fitness
curves$base.bmax <- base$f.bmax
curves$base.bmin <- base$f.bmin

# Wrapper function to plot baseline versus alternative environments and evolved tolerance curves
plot.conrast <- function(env, interval = 1:200, equilibrium = 401){
  xb <- sim.evo(env = env, plot = FALSE)
  out <- evo.out(xvec = xb[ , 1], bvec = xb[ , 2], plot = FALSE, equilibrium = equilibrium)
  
  # Plot environments
  plot(time[interval], env.base[interval], ylim = c(-1, 20), type = "l", xlab = "Time", ylab = "Environment", col = "black")
  lines(time[interval], env[interval], col = "red")
  
  # Plot tolerance curves
  change <- plot.tolcurve(xvec = xb[ , 1], bvec = xb[ , 2], plot = FALSE, equilibrium = equilibrium)
  change <- as.data.frame(change)
  curves$change <- change$fitness
  curves$change.bmax <- change$f.bmax
  curves$change.bmin <- change$f.bmin
  
  environment <- seq(from = -10, to = 25, length.out = 1000)
  plot(environment, curves$base, type = "l", xlab = "Environment", ylab = "Fitness", lwd = 4, ylim = c(0, 0.7))
  segments(x0 = environment[which.max(curves$base)], y0 = -1, x1 = environment[which.max(curves$base)], y1 = max(curves$base), lty = 2)
  polygon(c(environment, rev(environment)), c(curves$base.bmin, rev(curves$base.bmax)), col = grey(0.7, alpha = 0.3), border = NA)
  
  lines(environment, curves$change, lwd = 4, col = "red")
  segments(x0 = environment[which.max(curves$change)], y0 = -1, x1 = environment[which.max(curves$change)], y1 = max(curves$change), lty = 2, col = "red") #col = "red"
  polygon(c(environment, rev(environment)), c(curves$change.bmax, rev(curves$change.bmin)), col = rgb(1, 0, 0, alpha = 0.2), border = NA) #col = rgb(1, 0, 0, alpha = 0.2)
}

# A) Change in mean
plot.conrast(env = sim.env(plot = TRUE, e0 = 10))

# B) Change in trend
plot.conrast(env = sim.env(plot = TRUE, nt = 0.001), interval = 1:10000, equilibrium = 800)

# C) Cyclic change
plot.conrast(env = sim.env(plot = TRUE, nc = 4))

# D) Change in noise magnitude
plot.conrast(env = sim.env(plot = TRUE, nn = 2))

# E) Change in noise colour
plot.conrast(env = sim.env(plot = TRUE, rho = 0.65))

# F) Change in mean + trend
plot.conrast(env = sim.env(plot = TRUE, e0 = 10, nt = 0.001), interval = 1:10000, equilibrium = 800)

# G) Change in mean + cyclic
plot.conrast(env = sim.env(plot = TRUE, e0 = 10, nc = 4))

# H) Change in mean + noise magnitude
plot.conrast(env = sim.env(plot = TRUE, e0 = 10, nn = 2))

# I) Change in mean + noise colour
plot.conrast(env = sim.env(plot = TRUE, e0 = 10, rho = 0.65))

# J) Change in trend + cyclic
plot.conrast(env = sim.env(plot = TRUE, nt = 0.001, nc = 4), interval = 1:10000, equilibrium = 800)

# K) Change in trend + noise magnitude
plot.conrast(env = sim.env(plot = TRUE, nt = 0.001, nn = 2), interval = 1:10000, equilibrium = 800)

# L) Change in trend + noise colour
plot.conrast(env = sim.env(plot = TRUE, nt = 0.001, rho = 0.65), interval = 1:10000, equilibrium = 800)

# M) Cyclic change + noise magnitude
plot.conrast(env = sim.env(plot = TRUE, nc = 4, nn = 2))

# N) Cyclic change + noise colour
plot.conrast(env = sim.env(plot = TRUE, nc = 4, rho = 0.65))

# O) Change in noise magnitude + colour
plot.conrast(env = sim.env(plot = TRUE, nn = 2, rho = 0.65))
