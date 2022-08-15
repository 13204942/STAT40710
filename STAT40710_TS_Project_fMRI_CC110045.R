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
  phis <- matrix(NA, ncol = 3, nrow = K)
  phis_se <- matrix(NA, ncol = 3, nrow = K)
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
    
    if (p_hat > 3){
      p_hat <- 3
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
# CC110045_schaefer_400ROI_rest.txt - moved very little
CC110045 <- read.delim("CC110045_schaefer_400ROI_rest.txt")
CC110045 <- CC110045[,2:400]
# plot fMRI data
matplot(CC110045, type="l", lty = "solid")

X <- data.matrix(CC110045)
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
# McLeod.Li test for conditional heteroscedascity 
McLeod.Li.test(y=x)

# ACF and PACF
arch_identify <- acf_identifier(x, lags)
garch_orders_hat <- order_identifier(arch_identify, lags)
garch_orders_hat

#Fit GARCH(1,1)
fitted.res1<- garchFit(~ garch(1,1), 
                       data=x, 
                       cond.dist="norm", 
                       include.mean = FALSE, 
                       trace = FALSE)
llh <- -unname(fitted.res1@fit$llh)
AIC_1 <- (-2*(llh))/N + 2*(length(fitted.res1@fit$par))/N
AIC_1

#Fit GARCH(2,1)
fitted.res2 <- garchFit(~ garch(2,1), 
                        data=x, 
                        cond.dist="norm", 
                        include.mean = FALSE, 
                        trace = FALSE)
llh <- -unname(fitted.res2@fit$llh)
AIC_2 <- (-2*(llh))/N + 2*(length(fitted.res2@fit$par))/N
AIC_2

#Fit GARCH(2,2)
fitted.res3 <- garchFit(~ garch(2,2), 
                        data=x, 
                        cond.dist="norm", 
                        include.mean = FALSE, 
                        trace = FALSE)
llh <- -unname(fitted.res3@fit$llh)
AIC_3 <- (-2*(llh))/N + 2*(length(fitted.res3@fit$par))/N
AIC_3


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

phi.est2 <- matrix(NA, nrow = K, ncol = 3)
phi.est2.se <- matrix(NA, nrow = K, ncol = 3)

for(j in (1:K)){
  u <- arma.est1$orders[j,1]
  arma_fit <- arima(X_hat[j,], 
                    order = c(u,0,0), 
                    method = "ML")
  # extract the estimated coefficients
  if (u > 3){
    u = 3
  } else if (u == 2) {
    phi.est2[j,3] <- 0.0
    phi.est2.se[j,3] <- 0.0
  }  else {
    phi.est2[j,2:3] <- c(0.0, 0.0)
    phi.est2.se[j,2:3] <- c(0.0, 0.0)
  }  
  # extract the estimated coefficients
  phi.est2[j,1:u] <- arma_fit$coef[1:u]
  phi.est2.se[j,1:u] <- sqrt(diag(arma_fit$var.coef))[1:u]
}

head(phi.est2)


### Modelling each series one by one ###
phi.true <- matrix(NA, nrow = K, ncol = 3)
phi.true.se <- matrix(NA, nrow = K, ncol = 3)
sigma.true <- matrix(NA, nrow = K, ncol = N)

# method 1
# for(i in (1:K)){
#   y <- X[i,]
#   u <- arma.est1$orders[i,1]
#   
#   if (u > 3){
#     u = 3
#     fit <- garchFit(~ arma(3,0) + garch(2,1), 
#                     data = y, 
#                     include.mean=FALSE,
#                     trace = FALSE) 
#   } else if (u == 1) {
#     fit <- garchFit(~ arma(1,0) + garch(2,1),
#                     data = y,
#                     include.mean=FALSE,
#                     trace = FALSE)
#     phi.est2[j,2:3] <- 0.0
#     phi.est2.se[j,2:3] <- 0.0
#   } else {
#     fit <- garchFit(~ arma(2,0) + garch(2,1),
#                     data = y,
#                     include.mean=FALSE,
#                     trace = FALSE)
#     phi.est2[j,3] <- 0.0
#     phi.est2.se[j,3] <- 0.0     
#   }
#     phi.true[i,1:u] <- fit@fit$coef[1:u]
#     phi.true.se[i,1:u] <- fit@fit$se.coef[1:u]
#     sigma.true[i,] <- fit@sigma.t
# }


# method 2
AIC <- matrix(NA, nrow = K, ncol = 1)
for(i in (1:K)){
  y <- X[i,]
  u <- arma.est1$orders[i,1]

  if (u == 1){
    phi.true[i,2:3] <- 0.0
    phi.true.se[i,2:3] <- 0.0       
  } else if (u == 2) {
    phi.true[i,3] <- 0.0
    phi.true.se[i,3] <- 0.0   
  }
  else {
    u = 3
  }
  
  ar_fit = arima(y, order = c(u,0,0), method = c("ML"))
  phi.true[i,1:u] <- ar_fit$coef[1:u]
  phi.true.se[i,1:u] <- sqrt(diag(ar_fit$var.coef))[1:u]  
  eta = ar_fit$residuals  
  
  garch_fit <- garchFit(~ garch(2,1), 
                  data = eta, 
                  include.mean=FALSE,
                  trace = FALSE)
  
  llh <- -unname(garch_fit@fit$llh)
  AIC[i,] <- (-2*(llh))/N + 2*(length(garch_fit@fit$par))/N
  sigma.true[i,] <- garch_fit@sigma.t
}


# How is the parameter sigma estimated ?
matplot(t(sigma.true), type="l", main = "CC110045 Sigma",
        xlab="ROI", 
        ylab="Sigma", 
        col = alpha("#47bfc4", 0.4)) 
matlines(sigma.est, type="l", lwd = 2,
        col = alpha("#f4766d", 0.8)) 

# How is the parameter phi estimated ?
# average two phi estimator
# check standard error of phi
########### phi 1 ###########
phi1.se <- cbind(phi.est1.se[,1], phi.true.se[,1], phi.est2.se[,1])

########### phi 2 ###########
phi2.se <- cbind(phi.est1.se[,2], phi.est2.se[,2], phi.true.se[,2])

########### phi 3 ###########
phi3.se <- cbind(phi.est1.se[,3], phi.true.se[,3], phi.est2.se[,3])

par(mfrow=c(1,3))
matplot(phi1.se, main = "CC110045 Phi 1",
        xlab="ROI", 
        ylab="Standard Error", 
        type = "p",
        pch = 1,
        lwd = 0.8,
        col = c("#f4766d","#47bfc4"))

matplot(phi2.se, main = "CC110045 Phi 2",
        xlab="ROI", 
        ylab="Standard Error", 
        type = "p",
        pch = 1,
        lwd = 0.8,
        col = c("#f4766d","#47bfc4"))

matplot(phi3.se, main = "CC110045 Phi 3",
        xlab="ROI", 
        ylab="Standard Error", 
        type = "p",
        pch = 1,
        lwd = 0.8,
        col = c("#f4766d","#47bfc4"))
#############################

phi.est <- (phi.est1 + phi.est2)/2

### phi 1 ### 
phi1.comp <- melt(cbind(phi.est[,1], phi.true[,1]))
phi1.comp$Var2 <- as.factor(phi1.comp$Var2)
levels(phi1.comp$Var2) <- c("phi1.est", "phi1.true")
phi1.comp$Var1 <- as.factor(phi1.comp$Var1)
# compute Mean Squared Error
phi1.mse <- (1/(2*K))*sum((phi.true[,1] - phi.est[,1])^2)
phi1.mse

### phi 2 ### 
phi2.comp <- melt(cbind(phi.est[,2], phi.true[,2]))
phi2.comp$Var2 <- as.factor(phi2.comp$Var2)
levels(phi2.comp$Var2) <- c("phi2.est", "phi2.true")
phi2.comp$Var1 <- as.factor(phi2.comp$Var1)
# compute Mean Squared Error
phi2.mse <- (1/(2*K))*sum((phi.true[,2] - phi.est[,2])^2)
phi2.mse

### phi 3 ### 
phi3.comp <- melt(cbind(phi.est[,3], phi.true[,3]))
phi3.comp$Var2 <- as.factor(phi3.comp$Var2)
levels(phi3.comp$Var2) <- c("phi3.est", "phi3.true")
phi3.comp$Var1 <- as.factor(phi3.comp$Var1)
# compute Mean Squared Error
phi3.mse <- (1/(2*K))*sum((phi.true[,3] - phi.est[,3])^2)
phi3.mse

# Plotting
par(mfrow=c(1,3))
ggplot(data=phi1.comp, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "True Phi 1 (Green) vs Estimated Phi 1 (Red)") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.2) +
  scale_x_discrete(breaks = c("1","50","100","150","200","250","300","350","400")) +
  theme(legend.position="none")

ggplot(data=phi2.comp, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "True Phi 2 (Green) vs Estimated Phi 2 (Red)") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.2) +
  scale_x_discrete(breaks = c("1","50","100","150","200","250","300","350","400")) +
  theme(legend.position="none")

ggplot(data=phi3.comp, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "True Phi 3 (Green) vs Estimated Phi 3 (Red)") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.3) +
  scale_x_discrete(breaks = c("1","50","100","150","200","250","300","350","400")) +
  theme(legend.position="none")
#################### CC110056 WIP END ####################

# e <- fit@residuals
# d <- e^2 - fit@sigma.t^2
# mean(d^2)

# test <- as.data.frame(fit@residuals)
# ind <- seq(1, 400, by=20)
# ggplot(data=test, aes(x=ind, y=fit@residuals, group=1)) +
#   geom_line()

# matplot(e, type="l", main = "CC110045 Residuals",
#         xlab="", 
#         ylab="Residuals", 
#         lty = "solid",
#         lwd = 1)


matplot(fitted.residual@residuals, type="l", 
        xlab="Time Points",
        ylab="Returns (eta_hat)", 
        lty = "solid",
        lwd = 2)

# matplot(sqrt(e), type="l")
# matplot(fitted.residual@h.t, type="l")
# matplot(sigma.est^2, type="l")

# Goodness of fit for the variance prediction
e <- fitted.residual@residuals
d <- e^2 - fitted.residual@sigma.t^2
mean(d^2)


#################### #################### ####################
garch_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(2,1)), 
                   mean.model=list(armaOrder=c(0,0), include.mean = FALSE),  
                   distribution.model="norm")
ugfit <- ugarchfit(spec = garch_spec, data = c(eta_hat))
ugfit
plot(ugfit)
plot(ugfit, which = 'all')

