# IIC) Evolution of plasticity under different forms of environmental predictability

# A) Cyclic change dominates
predANDplas <- data.frame(gt = NA, rho = NA, cor = NA, b = NA, ac = NA)[numeric(0), ]
for(gt_i in c(1, 2, 11)){
  rhos <- c() # autocorrelation of noise
  cors <- c() # correlation between environments of development and selection
  bs <- c() # evolved mean plasticities
  acs <- c() # environmental time series autocorrelation with all components put together
  for(i in 0:9){
    
    rho <- i/10 
    env <- sim.env(nc = 4, rho = rho, plot = FALSE)
    
    rhos <- c(rhos, rho)
    acs <- c(acs, pacf(env, plot = FALSE)$acf[1])
    
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
    #
    
    # Derivatives:
    #This funciton defines mean lifetime fitness
    w <- expression(log(wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + (xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar))))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - (x + ((b*w_b^2)/(w_b^2 + bvar))*e1))^2)/(2*w_t^2 + 2*(xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar)))))
    #dlog(w)/dx
    wx <- D(w, "x")
    #dlog(w)/db
    wb <- D(w, "b")
    
    gt <- gt_i
    timesteps <- length(time)
    selevents <- trunc(timesteps/gt)
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    
    # ---------------------------------------------- Simulation -----------------------------------------------------------------------
    xvec <- c() 
    bvec <- c()
    env.dev <- c()
    env.sel <- c()
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
      
      env.dev <- c(env.dev, e1)
      env.sel <- c(env.sel, e2)
    }
    
    cors <- c(cors, cor(env.dev, env.sel))
    
    xvec <- rep(xvec, each = gt)
    bvec <- rep(bvec, each = gt)
    
    start <- 401
    end <- length(xvec) #same as length(bvec)
    x <- mean(xvec[start:end])
    b <- mean(bvec[start:end])
    
    bs <- c(bs, b)
    
  }
  predANDplas <- rbind(predANDplas, data.frame(gt = rep(gt_i, 10), rho = rhos, cor = cors, b = bs, ac = acs))
}

# Save results for when cyclic change dominates
predANDplas_cyclic <- predANDplas
predANDplas_cyclic$dominates <- rep("cyclic", nrow(predANDplas_cyclic))


# B) Noise dominates
predANDplas <- data.frame(gt = NA, rho = NA, cor = NA, b = NA, ac = NA)[numeric(0), ]
for(gt_i in c(1, 2, 11)){
  rhos <- c() # autocorrelation of noise
  cors <- c() # correlation between environments of development and selection
  bs <- c() # evolved mean plasticities
  acs <- c() # environmental time series autocorrelation with all components put together
  for(i in 0:9){
    
    rho <- i/10 
    env <- sim.env(nn = 2, rho = rho, plot = FALSE)
    
    rhos <- c(rhos, rho)
    acs <- c(acs, pacf(env, plot = FALSE)$acf[1])
    
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
    #
    
    # Derivatives:
    #This funciton defines mean lifetime fitness
    w <- expression(log(wmax*sqrt((w_b^2*w_t^2)/((w_b^2 + bvar)*(w_t^2 + (xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar))))*exp(-(b^2)/(2*w_b^2 + 2*bvar) - ((theta - (x + ((b*w_b^2)/(w_b^2 + bvar))*e1))^2)/(2*w_t^2 + 2*(xvar + ((bvar*w_b^2)/(w_b^2 + bvar))*e1^2 + yvar)))))
    #dlog(w)/dx
    wx <- D(w, "x")
    #dlog(w)/db
    wb <- D(w, "b")
    
    gt <- gt_i
    timesteps <- length(time)
    selevents <- trunc(timesteps/gt)
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    
    # ---------------------------------------------- Simulation -----------------------------------------------------------------------
    xvec <- c() 
    bvec <- c()
    env.dev <- c()
    env.sel <- c()
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
      
      env.dev <- c(env.dev, e1)
      env.sel <- c(env.sel, e2)
    }
    
    cors <- c(cors, cor(env.dev, env.sel))
    
    xvec <- rep(xvec, each = gt)
    bvec <- rep(bvec, each = gt)
    
    start <- 401
    end <- length(xvec) #same as length(bvec)
    x <- mean(xvec[start:end])
    b <- mean(bvec[start:end])
    
    bs <- c(bs, b)
    
  }
  predANDplas <- rbind(predANDplas, data.frame(gt = rep(gt_i, 10), rho = rhos, cor = cors, b = bs, ac = acs))
}

# Save results for when noise dominates
predANDplas_noise <- predANDplas
predANDplas_noise$dominates <- rep("noise", nrow(predANDplas_noise))

# Combine and reorder
predANDplas_full <- rbind(predANDplas_cyclic, predANDplas_noise)
predANDplas_full$gt <- as.factor(predANDplas_full$gt)
predANDplas_full <- predANDplas_full[predANDplas_full$gt %in% c(1, 2, 11), ]

# Plot environments
env_cyclic_white <- sim.env(nc = 4)
env_cyclic_red <- sim.env(nc = 4, rho = 0.9)
env_noise_white <- sim.env(nn = 2)
env_noise_red <- sim.env(nn = 2, rho = 0.9)

interval <- 1:200
par(mar=c(4.3,4.3,0.5,0.5))
plot(time[interval], env_cyclic_white[interval], type = "l", xlab = "Time-steps", ylab = "Environment")
plot(time[interval], env_cyclic_red[interval], type = "l", xlab = "Time-steps", ylab = "Environment")
plot(time[interval], env_noise_white[interval], type = "l", xlab = "Time-steps", ylab = "Environment")
plot(time[interval], env_noise_red[interval], type = "l", xlab = "Time-steps", ylab = "Environment")

# Plot evolved plasticities
library(viridis)
library(ggplot2)
library(lemon)

env.labs <- c("A) Cyclic change dominates", "B) Noise dominates")
names(env.labs) <- c("cyclic", "noise")

gg <- ggplot(predANDplas_full, aes(x = rho, y = b, group = gt)) + #, colour = gt
  geom_point(aes(fill = cor), size = 4, shape = 21, alpha = 0.8) +
  xlim(0, 1) +
  ylim(-0.5, 1) +
  scale_fill_viridis(option = "D", name = "Correlation between\nenvironments of\ndevelopment\nand selection") +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, lwd = 0.3, lty = 2, color = "black") + #, color = "black"
  labs(y = "Mean plasticity (b)", x = "Autocorrelation or colour of noise (rho)") +
  facet_rep_grid(. ~ dominates, labeller = labeller(dominates = env.labs), repeat.tick.labels = FALSE) +
  theme_classic() +
  theme(panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, vjust = 0, hjust = 0), #facet titles
        legend.box.margin=margin(0, 0, 0, 20),
        axis.title = element_text(size = 14), # axis titles
        axis.text = element_text(size = 12),
        legend.text.align = 0,
        legend.text = element_text(size = 12))
gg


# Autocorrelation of the entire time series versus autocorrelation of noise

# 1) When cyclic change dominates
cyclic <- predANDplas_full[predANDplas_full$dominates == "cyclic", ]
plot(cyclic$rho, cyclic$ac, ylim = c(0, 1), xlim = c(0, 1))
range(cyclic$ac) # autocorrelation of the entire time series goes from 0.72 to 0.82
range(cyclic$rho) # noise autocorrelation goes from 0 to 0.9 (as defined by me)
cor(cyclic$rho, cyclic$ac) # Their correlation is 0.9999
summary(lm(cyclic$ac ~ cyclic$rho)) # Their slope is 0.11 and R2 is 1 (as should be)

# 2) When noise dominates
noise <- predANDplas_full[predANDplas_full$dominates == "noise", ]
plot(noise$rho, noise$ac, ylim = c(0, 1), xlim = c(0, 1))
range(noise$ac) # autocorrelation of the entire time series goes from 0.1 to 0.9
range(noise$rho) # noise autocorrelation goes from 0 to 0.9 (as defined by me)
cor(noise$rho, noise$ac) # Their correlation is 0.9999
summary(lm(noise$ac ~ noise$rho)) # Their slope is 0.9 and R2 is 0.9999 (as should be)
