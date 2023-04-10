# IIB) Rates of evolutionary adaptation

# 1) Change in mean
env1 <- sim.env(plot = FALSE, e0 = 10)
xb1 <- sim.evo(env = env1, plot = FALSE)
xvec.1 <- xb1[ , 1]; bvec.1 <- xb1[ , 2]

# 2) Change in mean + cycle
env2 <- sim.env(plot = FALSE, e0 = 10, nc = 4)
xb2 <- sim.evo(env = env2, plot = FALSE)
xvec.2 <- xb2[ , 1]; bvec.2 <- xb2[ , 2]

# 3) Change in mean + noise magnitude
env3 <- sim.env(plot = FALSE, e0 = 10, nn = 2)
xb3 <- sim.evo(env = env3, plot = FALSE)
xvec.3 <- xb3[ , 1]; bvec.3 <- xb3[ , 2]

# 4) Change in mean + noise colour
env4 <- sim.env(plot = FALSE, e0 = 10, rho = 0.65)
xb4 <- sim.evo(env = env4, plot = FALSE)
xvec.4 <- xb4[ , 1]; bvec.4 <- xb4[ , 2]

# Set up colours for plotting
col.black <- rgb(0, 0, 0, max = 255, alpha = 125)
col.blue <- rgb(0, 0, 255, max = 255, alpha = 125)
col.green <- rgb(0, 255, 0, max = 255, alpha = 125)
col.red <- rgb(255, 0, 0, max = 255, alpha = 125)

####### Evolutionary trajectories and rates of adaptation of mean breeding values ####### 
interval <- 1:600
time <- interval
fit1.1 <- lm(xvec.1[interval] ~ 1 + time + I(time^2) + I(time^3) + I(time^4))
plot(time, xvec.1[interval], col = "blue", type = "l")
lines(time, predict(fit1.1))
c1.1 <- unname(fit1.1$coefficients[1])
c2.1 <- unname(fit1.1$coefficients[2])
c3.1 <- unname(fit1.1$coefficients[3])
c4.1 <- unname(fit1.1$coefficients[4])
c5.1 <- unname(fit1.1$coefficients[5])
x1 <- expression(c1.1 + c2.1*time + c3.1*time^2 + c4.1*time^3 + c5.1*time^4)
x1.t.1 <- D(x1, "time")
plot(time, eval(x1.t.1), type = "l")

fit2.1 <- lm(xvec.2[interval] ~ 1 + time + I(time^2) + I(time^3)+ I(time^4)) # + I(time^5)
plot(time, xvec.2[interval], col = "red", type = "l")
lines(time, predict(fit2.1))
c1.2 <- unname(fit2.1$coefficients[1])
c2.2 <- unname(fit2.1$coefficients[2])
c3.2 <- unname(fit2.1$coefficients[3])
c4.2 <- unname(fit2.1$coefficients[4])
c5.2 <- unname(fit2.1$coefficients[5])
# c6.2 <- unname(fit2.1$coefficients[6])
x2 <- expression(c1.2 + c2.2*time + c3.2*time^2 + c4.2*time^3 + c5.2*time^4) #+ c6.2*time^5
x2.t.1 <- D(x2, "time")
plot(time, eval(x2.t.1), type = "l")

fit3.1 <- lm(xvec.3[interval] ~ 1 + time + I(time^2) + I(time^3) + I(time^4))
plot(time, xvec.3[interval], col = "black", type = "l")
lines(time, predict(fit3.1))
c1.3 <- unname(fit3.1$coefficients[1])
c2.3 <- unname(fit3.1$coefficients[2])
c3.3 <- unname(fit3.1$coefficients[3])
c4.3 <- unname(fit3.1$coefficients[4])
c5.3 <- unname(fit3.1$coefficients[5])
x3 <- expression(c1.3 + c2.3*time + c3.3*time^2 + c4.3*time^3 + c5.3*time^4)
x3.t.1 <- D(x3, "time")
plot(time, eval(x3.t.1), type = "l")

fit4.1 <- lm(xvec.4[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, xvec.4[interval], col = "green", type = "l")
lines(time, predict(fit4.1))
c1.4 <- unname(fit4.1$coefficients[1])
c2.4 <- unname(fit4.1$coefficients[2])
c3.4 <- unname(fit4.1$coefficients[3])
c4.4 <- unname(fit4.1$coefficients[4])
c5.4 <- unname(fit4.1$coefficients[5])
x4 <- expression(c1.4 + c2.4*time + c3.4*time^2 + c4.4*time^3 + c5.4*time^4)
x4.t.1 <- D(x4, "time")
plot(time, eval(x4.t.1), type = "l")

# Plot all evolutionary trajectories together
interval <- 1:400
time <- interval
par(mar = c(4, 4, 0.5, 0.5), mgp=c(2.5, 1, 0))

# Change 1
plot(time, xvec.1[interval], type = "l", xlab = "Generations", ylab = expression(Mean~breeding~value~(italic(bar(x)))), ylim = c(0, 10), cex.lab = 1.5, col = col.blue) #ylim = c(0, 10) / ylim = c(0, 20) #col = col.black
lines(time, predict(fit1.1)[interval], lwd = 2, col = "blue")

# Change 2
lines(time, xvec.2[interval], col = col.red)
lines(time, predict(fit2.1)[interval], lwd = 2, col = "red")

# Change 3
lines(time, xvec.3[interval], col = col.black)
lines(time, predict(fit3.1)[interval], lwd = 2, col = "black")

# Change 4
lines(time, xvec.4[interval], col = col.green)
lines(time, predict(fit4.1)[interval], lwd = 2, col = "green")

legend("bottomright", c("Mean", "Mean + Cycle", "Mean + Noise", "Mean + Noise colour"),
       col = c("blue", "red", "black", "green"), bty = "n", lwd = 2)

# Plot rates of adaptation
par(mar = c(4, 6, 0.5, 0.5), mgp=c(2.5, 1, 0))
plot(time, eval(x1.t.1), type = "l", xlab = "Generations", ylab = expression(atop("Rate of adaptation", (Delta*bar(x)/generation))), lwd = 2, ylim = c(0, 0.08), cex.lab = 1.5, col = "blue")
lines(time, eval(x2.t.1), col = "red", lwd = 2)
lines(time, eval(x3.t.1), col = "black", lwd = 2)
lines(time, eval(x4.t.1), col = "green", lwd = 2)
abline(h = 0, lty = 2)
legend("topright", c("Mean", "Mean + Cycle", "Mean + Noise", "Mean + Noise colour"),
       col = c("blue", "red", "black", "green"), bty = "n", lwd = 2)


####### Evolutionary trajectories and rates of adaptation of mean plasticities ####### 
interval <- 1:600
time <- interval
fit1.2 <- lm(bvec.1[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, bvec.1[interval], col = "blue", type = "l")
lines(time, predict(fit1.2))
c1.1 <- unname(fit1.2$coefficients[1])
c2.1 <- unname(fit1.2$coefficients[2])
c3.1 <- unname(fit1.2$coefficients[3])
c4.1 <- unname(fit1.2$coefficients[4])
c5.1 <- unname(fit1.2$coefficients[5])
x1 <- expression(c1.1 + c2.1*time + c3.1*time^2 + c4.1*time^3 + c5.1*time^4)
x1.t.2 <- D(x1, "time")
plot(time, eval(x1.t.2), type = "l")

fit2.2 <- lm(bvec.2[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, bvec.2[interval], col = "red", type = "l")
lines(time, predict(fit2.2))
c1.2 <- unname(fit2.2$coefficients[1])
c2.2 <- unname(fit2.2$coefficients[2])
c3.2 <- unname(fit2.2$coefficients[3])
c4.2 <- unname(fit2.2$coefficients[4])
c5.2 <- unname(fit2.2$coefficients[5])
x2 <- expression(c1.2 + c2.2*time + c3.2*time^2 + c4.2*time^3 + c5.2*time^4)
x2.t.2 <- D(x2, "time")
plot(time, eval(x2.t.2), type = "l")

fit3.2 <- lm(bvec.3[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, bvec.3[interval], col = "black", type = "l")
lines(time, predict(fit3.2))
c1.3 <- unname(fit3.2$coefficients[1])
c2.3 <- unname(fit3.2$coefficients[2])
c3.3 <- unname(fit3.2$coefficients[3])
c4.3 <- unname(fit3.2$coefficients[4])
c5.3 <- unname(fit3.2$coefficients[5])
x3 <- expression(c1.3 + c2.3*time + c3.3*time^2 + c4.3*time^3 + c5.3*time^4)
x3.t.2 <- D(x3, "time")
plot(time, eval(x3.t.2), type = "l")

fit4.2 <- lm(bvec.4[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, bvec.4[interval], col = "green", type = "l")
lines(time, predict(fit4.2))
c1.4 <- unname(fit4.2$coefficients[1])
c2.4 <- unname(fit4.2$coefficients[2])
c3.4 <- unname(fit4.2$coefficients[3])
c4.4 <- unname(fit4.2$coefficients[4])
c5.4 <- unname(fit4.2$coefficients[5])
x4 <- expression(c1.4 + c2.4*time + c3.4*time^2 + c4.4*time^3 + c5.4*time^4)
x4.t.2 <- D(x4, "time")
plot(time, eval(x4.t.2), type = "l")

# Plot all evolutionary trajectories together
interval <- 1:400
time <- interval
par(mar = c(4, 4.1, 0.5, 0.5), mgp=c(2.5, 1, 0))

# Mean
plot(time, bvec.1[interval], type = "l", xlab = "Generations", ylab = expression(Mean~plasticity~(italic(bar(b)))), ylim = c(0, 1), cex.lab = 1.5, col = col.blue)
lines(time, predict(fit1.2)[interval], lwd = 2, col = "blue")

# Mean + cycle
lines(time, bvec.2[interval], col = col.red)
lines(time, predict(fit2.2)[interval], lwd = 2, col = "red")

# Mean + Noise
lines(time, bvec.3[interval], col = col.black)
lines(time, predict(fit3.2)[interval], lwd = 2, col = "black")

# Mean + Noise colour
lines(time, bvec.4[interval], col = col.green)
lines(time, predict(fit4.2)[interval], lwd = 2, col = "green")

legend("topright", c("Mean", "Mean + Cycle", "Mean + Noise", "Mean + Noise colour"),
       col = c("blue", "red", "black", "green"), bty = "n", lwd = 2)

# Plot rates of adaptation
par(mar = c(4, 6, 0.5, 0.5), mgp=c(2.5, 1, 0))
plot(time, eval(x1.t.2), type = "l", xlab = "Generations", ylab = expression(atop("Rate of adaptation", (Delta*bar(b)/generation))), lwd = 2, cex.lab = 1.5, ylim = c(-0.009, 0.0001), col = "blue")
lines(time, eval(x2.t.2), col = "red", lwd = 2)
lines(time, eval(x3.t.2), col = "black", lwd = 2)
lines(time, eval(x4.t.2), col = "green", lwd = 2)
abline(h = 0, lty = 2)
legend("bottomright", c("Mean", "Mean + Cycle", "Mean + Noise", "Mean + Noise colour"),
       col = c("blue", "red", "black", "green"), bty = "n", lwd = 2)





