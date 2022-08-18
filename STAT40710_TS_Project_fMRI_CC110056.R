library(TSA)
library(fGarch)
library(forecast)
library(tseries)
library(reshape2)
library(ggplot2)
library(rugarch)

acf_identifier <- function(x, lags){
  # ACF
  ro <- acf(x, lag.max = lags, plot = FALSE)
  # PACF
  phi <- pacf(x, lag.max = lags, plot = FALSE)
  
  return(cbind(acf = ro$acf, pacf = phi$acf))
}

order_identifier <- function(acf_pacf, lags){
  p_hat <- 1
  q_hat <- 1
  
  for(l in (2:lags)){
    if(abs(acf_pacf[,'acf'][l]) > ci[2]){
      q_hat = q_hat + 1
    }
    else{
      break
    }
  }
  
  for(l in (2:lags)){
    if(abs(acf_pacf[,'pacf'][l]) > ci[2]){
      p_hat = p_hat + 1
    }
    else{
      break
    }
  }  
  return(list(p_hat = p_hat, q_hat = q_hat))
}

# Create a function to estimate AR orders and get residuals
ar_estimator <- function(series, K, N, lags){
  
  acfs <- matrix(NA, ncol = 2, nrow = lags, dimnames = list(NULL,c('acf','pacf')))
  orders <- matrix(NA, ncol = 2, nrow = K, dimnames = list(NULL,c('p_hat','q_hat')))
  residuals_matrix <- matrix(NA, ncol = N, nrow = K)
  phis <- matrix(NA, ncol = 2, nrow = K)
  phis_se <- matrix(NA, ncol = 2, nrow = K)  
  results <- list(NA)
  
  for(index in (1:K)){
    tar_ts <- series[index,]
    
    # ACF
    ro <- acf(tar_ts, lag.max = lags, plot = FALSE)
    # PACF
    phi <- pacf(tar_ts, lag.max = lags, plot = FALSE)
    
    acfs[,'acf'] <- as.vector(ro$acf)
    acfs[,'pacf'] <- as.vector(phi$acf)
    arma_orders <- order_identifier(acfs, lags)  
    orders[index,] <- c(arma_orders$p_hat, arma_orders$q_hat)
    
    # assumption: series has AR component
    p_hat <- arma_orders$p_hat
    fit_ts = arima(tar_ts, order = c(p_hat,0,0), method = c("ML"))
    x = fit_ts$residuals  
    residuals_matrix[index,] <- x
    
    if (p_hat > 2){
      p_hat <- 2
    }
    phis[index,1:p_hat] <- fit_ts$coef[1:p_hat]
    phis_se[index,1:p_hat] <- sqrt(diag(fit_ts$var.coef))[1:p_hat]
  }
  results <- list(orders = orders, 
                  residuals = residuals_matrix, 
                  phi = phis, se = phis_se,
                  fit_ts = fit_ts)
  return(results) 
}

#################### CC110056 WIP BEGIN ####################
# load fMRI data
# CC110056_schaefer_400ROI_rest.txt - showed high movement
CC110056 <- read.delim("CC110056_schaefer_400ROI_rest.txt")
CC110056 <- CC110056[,2:400]
# CC110045_schaefer_400ROI_rest.txt - moved very little
# plot fMRI data
matplot(CC110056, type="l", lty = "solid")

X <- data.matrix(CC110056)
X <- t(X)

lags = 15
K <- nrow(X)
N <- ncol(X)
ci = c(-1.96/sqrt(N), 1.96/sqrt(N))

arma.est1 <- ar_estimator(X, K, N, lags)

# estiamted orders
arma.est1$orders[,1]

# all estimated phi values
phi.est1 <- arma.est1$phi
# phi 1 standardised error
phi.est1.se <- arma.est1$se

phi.est1[is.na(phi.est1)] <- 0
phi.est1.se[is.na(phi.est1.se)] <- 0


# w0 = 1/phi
w0 <- 1/rowSums(phi.est1)
# w = w0/sum(w0)
w <- w0/(sum(w0))

# weights
#w <- c(w)
# all residuals
est_residual <- arma.est1$residuals
# calculate estimated eta with weights
eta_hat <- w%*%est_residual
matplot(c(eta_hat), type="l")

x <- c(eta_hat)

# examine cross-correlation
res.mat <- est_residual%*%t(est_residual)
res.cor <- cor(res.mat) 
res.cor_lower <- lower.tri(res.cor, diag = FALSE)
cross_corr <- c(res.cor[res.cor_lower])
hist(cross_corr, # histogram
     col="#c2e7cd", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "Cross-correlation",
     main = "Cross-correlation of squared residuals")
lines(density(cross_corr), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")

# plot distribution of cross-correlation values
residual_vec <- as.numeric(rowSums(est_residual))
ccf(residual_vec, residual_vec, 
    ylab = "cross-correlation", main = "")

# ACF and PACF
arch_identify <- acf_identifier(x, lags)
garch_orders_hat <- order_identifier(arch_identify, lags)
garch_orders_hat

# McLeod.Li test for conditional heteroscedascity 
x.arma <- arima(x, order=c(1,0,1))
x.res <- residuals(x.arma)
McLeod.Li.test(y=x.res)

#Fit GARCH(1,1)
fitted.res1<- garchFit(~ garch(1,1), 
                       data=x, 
                       cond.dist="norm", 
                       include.mean = FALSE, 
                       trace = FALSE)
llh <- -unname(fitted.res1@fit$llh)
AIC_1 <- (-2*(llh))/N + 2*(length(fitted.res1@fit$par))/N
AIC_1

#Fit GARCH(2,2)
fitted.res2 <- garchFit(~ garch(2,2), 
                        data=x, 
                        cond.dist="norm", 
                        include.mean = FALSE, 
                        trace = FALSE)
llh <- -unname(fitted.res2@fit$llh)
AIC_2 <- (-2*(llh))/N + 2*(length(fitted.res2@fit$par))/N
AIC_2

if(AIC_1 < AIC_2){
  fitted.residual <- fitted.res1 
} else {
  fitted.residual <- fitted.res2
}

# Check Standardised Residuals of GARCH (epsilon_t)
garch_res = eta_hat / fitted.residual@sigma.t

garch_res <- as.matrix(garch_res)
qqnorm(garch_res, main='GARCH Standard Residuals')
qqline(colMeans(garch_res))

library(car)
qqPlot(garch_res, main='GARCH Standard Residuals',
       ylab='Sample Quantiles',
       xlab='Theoretical Quantiles')

# Check statistics of Standardised Residuals 
# Ljung-Box Test on R and R^2
summary(fitted.residual)

# Extract estimator sigma
sigma.est <- fitted.residual@sigma.t

# remove eta_hat from series
X_hat <- sweep(X, 2, c(eta_hat)) 

phi.est2 <- matrix(NA, nrow = K, ncol = 2)
phi.est2.se <- matrix(NA, nrow = K, ncol = 2)

for(j in (1:K)){
  u <- arma.est1$orders[j,1]
  arma_fit <- arima(X_hat[j,], 
                  order = c(u,0,0), 
                  method = "ML")
  # extract the estimated coefficients
  if (u > 2){
    u = 2
  } else {
    phi.est2[j,2] <- 0.0
    phi.est2.se[j,2] <- 0.0
  }
  phi.est2[j,1:u] <- arma_fit$coef[1:u]
  phi.est2.se[j,1:u] <- sqrt(diag(arma_fit$var.coef))[1:u]
}

head(phi.est2)
head(phi.est2.se)


### Modelling each series one by one ###
phi1.old <- matrix(NA, nrow = K, ncol = 2)
phi1.old.se <- matrix(NA, nrow = K, ncol = 2)
sigma.true <- matrix(NA, nrow = K, ncol = N)

# method 1
# for(i in (1:K)){
#   y <- X[i,]
#   u <- arma.est1$orders[i,1]
# 
#   if (u > 2){
#     u = 2
#     fit <- garchFit(~ arma(2,0) + garch(1,1),
#                     data = y,
#                     include.mean=FALSE,
#                     trace = FALSE)
#   } else {
#     fit <- garchFit(~ arma(1,0) + garch(1,1),
#                     data = y,
#                     include.mean=FALSE,
#                     trace = FALSE)
# 
#     phi.est2[j,2] <- 0.0
#     phi.est2.se[j,2] <- 0.0
#   }
# 
#   phi1.old[i,1:u] <- fit@fit$coef[1:u]
#   phi1.old.se[i,1:u] <- fit@fit$se.coef[1:u]
#   sigma.true[i,] <- fit@sigma.t
# }

# method 2
AIC <- matrix(NA, nrow = K, ncol = 1)
for(i in (1:K)){
  y <- X[i,]
  u <- arma.est1$orders[i,1]
  
  if (u < 2){
    phi1.old[i,2] <- 0.0
    phi1.old.se[i,2] <- 0.0       
  } else {
    u = 2
  }
  
  ar_fit <- arima(y, order = c(u,0,0), method = "ML")
  phi1.old[i,1:u] <- ar_fit$coef[1:u]
  phi1.old.se[i,1:u] <- sqrt(diag(ar_fit$var.coef))[1:u]
  eta = ar_fit$residuals  
  
  garch_fit <- garchFit(~ garch(1,1), 
                  data = eta, 
                  include.mean=FALSE,
                  trace = FALSE)
  
  llh <- -unname(garch_fit@fit$llh)
  AIC[i,] <- (-2*(llh))/N + 2*(length(garch_fit@fit$par))/N
  sigma.true[i,] <- garch_fit@sigma.t
}

# How is the parameter sigma estimated ?
matplot(t(sigma.true), type="l", main = "CC110056 Sigma",
        xlab="Time Points", 
        ylab="Sigma", 
        lty = "solid",
        col = alpha("#f4766d", 0.4))
matlines(sigma.est, type="l", lwd = 2,
        col = alpha("#47bfc4", 0.8))


# How is the parameter phi estimated ?
# average two phi estimator
# check standard error of phi
########### phi 1 ###########
phi1.se <- cbind(phi.est1.se[,1], phi1.old.se[,1], phi.est2.se[,1])

########### phi 2 ###########
phi2.se <- cbind(phi.est1.se[,2], phi1.old.se[,2], phi.est2.se[,2])

#### f4766d red
#### 47bfc4 green
par(mfrow=c(1,2))
matplot(phi1.se, main = "CC110056 Phi 1",
        xlab="ROI", 
        ylab="Standard Error", 
        type = "p",
        pch = 1,
        lwd = 0.8,
        col = c("#f4766d","#47bfc4"))

matplot(phi2.se, main = "CC110056 Phi 2",
        xlab="ROI", 
        ylab="Standard Error", 
        type = "p",
        pch = 1,
        lwd = 0.8,
        col = c("#f4766d","#47bfc4"))
#############################

phi.est <- (phi.est1 + phi.est2)/2

### standard error ### 
phi.est.se <- 0.5*sqrt((phi.est1.se^2 + phi.est2.se^2))

### phi 1 ### 
phi1.comp <- melt(cbind(phi.est.se[,1], phi1.old.se[,1]))
phi1.comp$Var2 <- as.factor(phi1.comp$Var2)
levels(phi1.comp$Var2) <- c("phi1.new.est", "phi1.old.est")
phi1.comp$Var1 <- as.factor(phi1.comp$Var1)
# compute Mean Squared Error
# phi1.mse <- (1/(2*K))*sum((phi1.old[,1] - phi.est[,1])^2)
# phi1.mse

### phi 2 ### 
phi2.comp <- melt(cbind(phi.est.se[,2], phi1.old.se[,2]))
phi2.comp$Var2 <- as.factor(phi2.comp$Var2)
levels(phi2.comp$Var2) <- c("phi2.new.est", "phi2.old.est")
phi2.comp$Var1 <- as.factor(phi2.comp$Var1)
# compute Mean Squared Error
# phi2.mse <- (1/(2*K))*sum((phi1.old[,2] - phi.est[,2])^2)
# phi2.mse

ggplot(data=phi1.comp, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "Old Estimated Phi 1 (Green) vs New Estimated Phi 1 (Red)") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.2) +
  scale_x_discrete(breaks = c("1","50","100","150","200","250","300","350","400")) +
  theme(legend.position="none")


ggplot(data=phi2.comp, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "Old Estimated Phi 2 (Green) vs New Estimated Phi 2 (Red)") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.2) +
  scale_x_discrete(breaks = c("1","50","100","150","200","250","300","350","400")) +
  theme(legend.position="none")


#################### CC110056 WIP END ####################

# plot(fitted.residual@residuals/fitted.residual@sigma.t,
#      type='h',
#      ylab='standard residuals')
# 
# matplot(fitted.residual@residuals, type="l")
# matplot(c(garch_res), type="l")
# 
# matplot(fitted.residual@residuals, type="l", col="#f4766d")
# matlines(sigma.est, type="l")
# 
# 
# matplot(sigma.est^2, type="l")
# matlines(c(eta_hat), type="l")

matplot(fitted.residual@residuals, type="l", 
        xlab="Time Points",
        ylab="Returns (eta_hat)", 
        lty = "solid",
        lwd = 2)

# Goodness of fit for the variance prediction
e <- fitted.residual@residuals
d <- e^2 - fitted.residual@sigma.t^2
mean(d^2)


#################### #################### ####################
garch_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
                         mean.model=list(armaOrder=c(0,0), include.mean = FALSE),  
                         distribution.model="norm")
ugfit <- ugarchfit(spec = garch_spec, data = c(eta_hat))
ugfit
plot(ugfit)
plot(ugfit, which = 'all')



