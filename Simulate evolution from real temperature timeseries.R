# Extracting temperature timeseries for different locations in the south east australian sea
# and predicting evolved tolerance curves for each location

# Locations: Southwest Tasmania (-43.43251, 145.76825), Port Phillip Bay in Victoria (-38.08977, 144.87799), Southern New South Wales (-36.25265, 151.6857)

# 1) Stack all rasters / run this in the cluster
library(raster)
stack <- stack()
for (year in 2010:2019){
  if(year == 2010){
    ini <- 7
    end <- 12
  } else if (year == 2019) {
    ini <- 1
    end <- 8
  } else {
    ini <- 1
    end <- 12
  }
  for (j in ini:end){
    month <- j
    month <- formatC(month, width = 2, format = "d", flag = "0")
    days <- formatC(1:31, width = 2, format = "d", flag = "0")
    ncfolder <- paste("~/Small ncFiles/", year, "/", year, "-", month, sep = "")
    files <- list.files(ncfolder, pattern = "*.nc$", full.names = T)
    map <- stack(files[1:length(files)])
    #val <- values(map) #This lines check if there are any values less than 0 (an error in the data that sometimes occur)
    #val[val < 0] <- NA #and replaces them with NAs, before adding the layers to the main stack
    #values(map) <- val
    #rm(val)
    stack <- addLayer(stack, map)
    rm(map)
    print(paste("year", year, "month", month, "done"))
  }
}


# 2) Extract SST for the three locations and organise in data frame / run this in the cluster
detach("package:tidyr", unload=TRUE)
sst <- extract(stack, cbind(c(145.76825, 144.87799, 151.6857), c(-43.43251, -38.08977, -36.25265)))
sst <- as.data.frame(sst)
Dates <- as.character(seq(as.Date("2010-07-01"), as.Date("2019-08-25"), by="days"))
names(sst) <- Dates
library(tidyr)
sst <- gather(data = sst, key = Date, value = SST, Dates)
sst$Location <- rep(c("Southwest Tasmania", "Port Phillip Bay", "Southern NSW"), nrow(sst)/3)
sst <- sst[ , c(3, 1, 2)]
#write.csv(sst, "~/case_study_sst.csv", row.names = FALSE)
sst <- read.csv("~/case_study_sst.csv")
sst$Location <- as.character(sst$Location)
sst$Date <- as.character(sst$Date)

# * There seems to be a problem witht the data between ~2018-05-01 and 2018-08-29, so I'll filter out the data from 2018-05-01 onwards
sst <- sst[1:unname(which(sst == "2018-05-01", arr.ind=TRUE)[1, 1])-1, ]
Dates <- as.character(seq(as.Date("2010-07-01"), as.Date("2018-04-30"), by="days"))

# 3) Extract time series and fill NAs
library(imputeTS)
tas <- sst[sst$Location == "Southwest Tasmania" , "SST"]
statsNA(tas)
tas <- na_interpolation(tas, option = "linear") 
statsNA(tas)

vic <- sst[sst$Location == "Port Phillip Bay" , "SST"]
statsNA(vic)
vic <- na_interpolation(vic, option = "linear") 
statsNA(vic)

nsw <- sst[sst$Location == "Southern NSW" , "SST"]
statsNA(nsw)
nsw <- na_interpolation(nsw, option = "linear") 
statsNA(nsw)


# 4) Plot sst data
#par(mar=c(4.5, 4.1, 2.1, 10.5), xpd = FALSE)
pdf("~/Case study/environments.pdf", width = 6.5, height = 4.3)
par(mar = c(4, 4, 0.5, 0.5), mgp=c(2.5, 1, 0))
#par(mar=c(4.5,4.3,0.5,0.5))
plot(as.Date(Dates), tas, type = "l", col = "blue", xlab = "Date", ylab = "Temperature (?C)", lwd = 1, ylim = c(min(sst$SST, na.rm = T), max(sst$SST, na.rm = T)), cex.lab = 1.5)
lines(as.Date(Dates), vic, col = "black", lwd = 1)
lines(as.Date(Dates), nsw, col = "red", lwd = 1)
#axis(2, at = mean(tas), labels = F, col = "blue", lty = 1, tcl=0.5)
#axis(2, at = mean(vic), labels = F, col = "black", lty = 1, tcl=0.5)
#axis(2, at = mean(nsw), labels = F, col = "red", lty = 1, tcl=0.5)
abline(h = mean(tas, na.rm = T), lty = 2, col = "blue")
abline(h = mean(vic, na.rm = T), lty = 2, col = "black")
abline(h = mean(nsw, na.rm = T), lty = 2, col = "red")
dev.off()

#par(xpd = TRUE)
legend("topright", c("NSW", "VIC", "TAS"), col = c("red", "black", "blue"), lwd = 4, bty = "n", cex = 0.8) #, inset=c(-0.53, 0)
dev.off()

# 5) Simulate evolved tolerance curve for each time series
location <- "Southwest Tasmania"
location <- "Port Phillip Bay"
location <- "Southern NSW"
env <- rep(tas, 20) # Repeat time series 20 times to allow evolutionary equilibrium
env <- rep(vic, 20) # Run this for each location separately
env <- rep(nsw, 20)

# Initial parameters
time <- 1:length(env)
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

# Derivatives:
#This funciton defines mean lifetime fitness
w <- expression(log(wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + (xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar))))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - (x + ((b*w_b^2)/(w_b^2 + bvar))*e1))^2)/(2*w_t^2 + 2*(xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar)))))
#dlog(w)/dx
wx <- D(w, "x")
#dlog(w)/db
wb <- D(w, "b")

# Specify generation time in days
gt <- 30
timesteps <- length(time)
selevents <- trunc(timesteps/gt)

# Run simulation
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

# Store evolution of mean breeding values and plasticities for each location run
xvec.1 <- xvec # tas
bvec.1 <- bvec

xvec.2 <- xvec # vic
bvec.2 <- bvec

xvec.3 <- xvec # nsw
bvec.3 <- bvec

# Plot results
par(mar=c(4.3,4.3,0.5,0.5))
plot(1:length(xvec), xvec, type = "l", xlab = "Generations", ylab = expression(bar(x)))
plot(1:length(bvec), bvec, type = "l", xlab = "Generations", ylab = expression(bar(b)))

# nsw (the one that takes longer) takes about 1,500 generations to reach equilibrium
start <- 1501
end <- length(xvec)
x <- mean(xvec[start:end])
b <- mean(bvec[start:end])
x
b
min(xvec[start:end]); max(xvec[start:end])
min(bvec[start:end]); max(bvec[start:end])
seq.ini <- 0
seq.end <- 40
environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
trait <- x + b*environment + y
theta <- environment
fitness <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
plot(environment, fitness, type = "l", xlab = "Environment", ylab = "Fitness") #tolerance curve
#abline(v = 5, col = "blue")
abline(v = environment[which.max(fitness)], col = "blue")

max(xvec[start:end])
min(xvec[start:end])

#Tolerance curve at bmax
bmax <- max(bvec[start:end])
x.bmax <- xvec[start:end][which.max(bvec[start:end])]
environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
trait <- x.bmax + bmax*environment + y
theta <- environment
fitness.bmax <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
plot(environment, fitness.bmax, type = "l", xlab = "Environment", ylab = "Fitness")

#Tolerance curve at bmin
bmin <- min(bvec[start:end])
x.bmin <- xvec[start:end][which.min(bvec[start:end])]
environment <- seq(from = seq.ini, to = seq.end, length.out = 1000)
trait <- x.bmin + bmin*environment + y
theta <- environment
fitness.bmin <- wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + tvar)))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - trait)^2)/(2*w_t^2 + 2*tvar))
plot(environment, fitness.bmin, type = "l", xlab = "Environment", ylab = "Fitness")

# Mean curve with shaded width variation (bmax - bmin)
height <- max(fitness)
center <- environment[which.max(fitness)]
center.bmax <- environment[which.max(fitness.bmax)]
center.bmin <- environment[which.max(fitness.bmin)]
width.bmax <- abs(environment[which.min(abs(fitness.bmax-(height/2)))] - center.bmax)/1.1775
width.bmin <- abs(environment[which.min(abs(fitness.bmin-(height/2)))] - center.bmin)/1.1775
f.bmax <- height*exp(-((environment-center)^2)/(2*width.bmax^2))
f.bmin <- height*exp(-((environment-center)^2)/(2*width.bmin^2))
par(mar=c(4.3,4.3,0.5,0.5))
plot(environment, fitness, type = "l", xlab = "Environment", ylab = "Fitness", lwd = 4) 
segments(x0 = environment[which.max(fitness)], y0 = -1, x1 = environment[which.max(fitness)], y1 = max(fitness), lty = 2)
#lines(environment, f.bmax, lty = 2)
#lines(environment, f.bmin, lty = 2)
polygon(c(environment, rev(environment)), c(f.bmin, rev(f.bmax)), col = grey(0.7, alpha = 0.3), border = NA)

# Plot curves in same graph
# Data frame to store simulated tolerance curves
curves <- data.frame(Environment = environment)

curves$tas <- fitness
curves$tas.bmax <- f.bmax
curves$tas.bmin <- f.bmin

curves$vic <- fitness
curves$vic.bmax <- f.bmax
curves$vic.bmin <- f.bmin

curves$nsw <- fitness
curves$nsw.bmax <- f.bmax
curves$nsw.bmin <- f.bmin

pdf("~/Case study/tol curves.pdf", width = 6.5, height = 4.3)
par(mar = c(4, 4, 0.5, 0.5), mgp=c(2.5, 1, 0))
plot(environment, curves$tas, type = "l", xlab = "Temperature (?C)", ylab = "Fitness", lwd = 4, ylim = c(0, 0.7), col = "black", xaxt="n", xlim = c(5, 30), cex.lab = 1.5)
axis(1, at = seq(0, 40, by = 5))
segments(x0 = environment[which.max(curves$tas)], y0 = -1, x1 = environment[which.max(curves$tas)], y1 = max(curves$tas), lty = 2, col = "black")
polygon(c(environment, rev(environment)), c(curves$tas.bmin, rev(curves$tas.bmax)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)

lines(environment, curves$vic, lwd = 4, col = "red")
segments(x0 = environment[which.max(curves$vic)], y0 = -1, x1 = environment[which.max(curves$vic)], y1 = max(fitness), lty = 2, col = "red")
polygon(c(environment, rev(environment)), c(curves$vic.bmax, rev(curves$vic.bmin)), col = rgb(0, 0, 0, alpha = 0.2), border = NA)

lines(environment, curves$nsw, lwd = 4, col = "blue")
segments(x0 = environment[which.max(curves$nsw)], y0 = -1, x1 = environment[which.max(curves$nsw)], y1 = max(fitness), lty = 2, col = "blue")
polygon(c(environment, rev(environment)), c(curves$nsw.bmax, rev(curves$nsw.bmin)), col = rgb(1, 0, 0, alpha = 0.2), border = NA)
#legend("topright", legend = c("NSW", "VIC", "TAS"), lwd = 4, col = c("blue", "red", "black"), bty = "n", cex = 0.8)
dev.off()

# Plot evolutionary trajectories of mean breeding value and plasticities

# Fit curves to adaptation of mean breeding value and mean plasticity
interval <- 1:1600 
time <- interval

col.black <- rgb(0, 0, 0, max = 255, alpha = 125)
col.blue <- rgb(0, 0, 255, max = 255, alpha = 125)
col.red <- rgb(255, 0, 0, max = 255, alpha = 125)

# 1) Mean breeding value
fit1.1 <- lm(xvec.1[interval] ~ 1 + time + I(time^2) + I(time^3) + I(time^4))
plot(time, xvec.1[interval], col = "black", type = "l")
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
plot(time, xvec.3[interval], col = "blue", type = "l")
lines(time, predict(fit3.1))
c1.3 <- unname(fit3.1$coefficients[1])
c2.3 <- unname(fit3.1$coefficients[2])
c3.3 <- unname(fit3.1$coefficients[3])
c4.3 <- unname(fit3.1$coefficients[4])
c5.3 <- unname(fit3.1$coefficients[5])
x3 <- expression(c1.3 + c2.3*time + c3.3*time^2 + c4.3*time^3 + c5.3*time^4)
x3.t.1 <- D(x3, "time")
plot(time, eval(x3.t.1), type = "l")

# Plot adaptation of mean breeding value with fitted lines
pdf("~/case study_breeding value.pdf", width = 7, height = 4.3)
interval <- 1:1500
time <- interval
par(mar = c(4, 4, 0.5, 0.5), mgp=c(2.5, 1, 0))

# tas
plot(time, xvec.1[interval], type = "l", xlab = "Generations", ylab = expression(Mean~breeding~value~(italic(bar(x)))), ylim = c(0, 10), cex.lab = 1.5, col = col.black) #ylim = c(0, 10) / ylim = c(0, 20) #col = col.black
lines(time, predict(fit1.1)[interval], lwd = 2, col = "black") #col = "black"

# vic
lines(time, xvec.2[interval], col = col.red) #col = col.blue
lines(time, predict(fit2.1)[interval], lwd = 2, col = "red") #col = "blue"

# nsw
lines(time, xvec.3[interval], col = col.blue) #col = col.green
lines(time, predict(fit3.1)[interval], lwd = 2, col = "blue") #col = "green"

legend("bottomright", c("NSW", "VIC", "TAS"), col = c("blue", "red", "black"), bty = "n", lwd = 2)
dev.off()

# 2) Mean plasticity
interval <- 1:1600
time <- interval
fit1.2 <- lm(bvec.1[interval] ~ time + I(time^2) + I(time^3) + I(time^4))
plot(time, bvec.1[interval], col = "black", type = "l")
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
plot(time, bvec.3[interval], col = "blue", type = "l")
lines(time, predict(fit3.2))
c1.3 <- unname(fit3.2$coefficients[1])
c2.3 <- unname(fit3.2$coefficients[2])
c3.3 <- unname(fit3.2$coefficients[3])
c4.3 <- unname(fit3.2$coefficients[4])
c5.3 <- unname(fit3.2$coefficients[5])
x3 <- expression(c1.3 + c2.3*time + c3.3*time^2 + c4.3*time^3 + c5.3*time^4)
x3.t.2 <- D(x3, "time")
plot(time, eval(x3.t.2), type = "l")

# Plot adaptation of mean plasticity with fitted lines
pdf("~/case study_plasticity.pdf", width = 7, height = 4.3)
interval <- 1:1500
time <- interval
par(mar = c(4, 4.1, 0.5, 0.5), mgp=c(2.5, 1, 0))

# tas
plot(time, bvec.1[interval], type = "l", xlab = "Generations", ylab = expression(Mean~plasticity~(italic(bar(b)))), ylim = c(0, 1), cex.lab = 1.5, col = col.black) #col = col.black, ylim = c(-0.1, 1) / ylim = c(-0.2, 1.1)
lines(time, predict(fit1.2)[interval], lwd = 2, col = "black") #col = "black"

# vic
lines(time, bvec.2[interval], col = col.red) #col = col.blue
lines(time, predict(fit2.2)[interval], lwd = 2, col = "red") #col = "blue"

# nsw
lines(time, bvec.3[interval], col = col.blue) #col = col.green
lines(time, predict(fit3.2)[interval], lwd = 2, col = "blue") #col = "green"

legend("topright", c("NSW", "VIC", "TAS"), col = c("blue", "red", "black"), bty = "n", lwd = 2)
dev.off()
