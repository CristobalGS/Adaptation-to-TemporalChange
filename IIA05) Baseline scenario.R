# Simulate evolution under baseline scenario and plot outcomes
env <- sim.env() # defaults to baseline
xb <- sim.evo(env = env) # Also plots mean breeding values and plasticities
out <- evo.out(xvec = xb[ , 1], bvec = xb[ , 2]) #"Mean breeding value = 4.52156105088215" "Mean plasticity = 0.0720575447107392"

# Plot environment
time <- 1:10000
interval <- 1:500
plot(time[interval], env[interval], type = "l", xlab = "Time-steps", ylab = "Environment")

# Plot reaction norm with variation in evolved plasticity (slope)
par(mar=c(4.3,4.3,0.5,0.5))
trait <- x + b*environment + y
trait.bmax <- x + bmax*environment + y
trait.bmin <- x + bmin*environment + y
plot(environment, trait, type = "l", xlab = "Environment", ylab = "Trait", ylim = c(0, 10), lwd = 4)
polygon(c(environment, rev(environment)), c(trait.bmin, rev(trait.bmax)), col = grey(0.7, alpha = 0.3), border = NA)

# Plot tolerance curve and reaction norm with variation in evolved plasticity (breadth/slope respectively)
plot.tolcurve(xvec = xb[ , 1], bvec = xb[ , 2], reacNorm = TRUE, out = FALSE, env.i = -5, env.f = 15, equilibrium = 401)
