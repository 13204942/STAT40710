library(TSA)
library(fGarch)
library(forecast)
library(tseries)
library(reshape2)
library(ggplot2)


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

sim_series_dif_phi <- function(u, K, n, lags, ci, eta){
  multi_series = matrix(NA, K, n)
  series_param = matrix(NA, nrow = K, ncol = 4)  
  # multiple time series share the same eta_t part
  for(k in 1:K){
    # generate different phi each time 
    # phi needs to be close to 1.0
    phi = runif(1, min = 0.7, max = 0.9)
    
    skip_to_next <- FALSE
    
    series_param[k,1] = u
    
    # the values of AR coefficients (ARIMA)
    series_param[k,2] = phi
    
    possibleError <- tryCatch({
      # use shared eta to simulate a AR(u) model
      multi_series[k,] = arima.sim(list(ar = phi), n, innov = eta$garch)
      
      # ACF
      series_ro <- pacf(multi_series[k,], lag.max = lags, plot = FALSE)
      series_param[k,3] <- sum(abs(series_ro$acf) > ci[2])
      series_param[k,4] <- 1    
    }, 
    error = function(e) {
      e
      print(paste("Skip the non-stationary series in loop ",k,sep = ""))
    })
    
    if(inherits(possibleError, "error")) next
  }
  results <- list(series = multi_series, params = series_param)
  return(results)  
}

sim_series_fix_phi <- function(phi, u, K, n, lags, ci, eta){
  multi_series = matrix(NA, K, n)
  series_param = matrix(NA, nrow = K, ncol = 4)  
  # multiple time series share the same eta_t part
  for(k in 1:K){
    skip_to_next <- FALSE
    
    series_param[k,1] = u
    
    # the values of AR coefficients (ARIMA)
    series_param[k,2] = phi
    
    possibleError <- tryCatch({
      # use shared eta to simulate a AR(u) model
      multi_series[k,] = arima.sim(list(ar = phi), n, innov = eta$garch)
      
      # ACF
      series_ro <- pacf(multi_series[k,], lag.max = lags, plot = FALSE)
      series_param[k,3] <- sum(abs(series_ro$acf) > ci[2])
      series_param[k,4] <- 1    
    }, 
    error = function(e) {
      e
      print(paste("Skip the non-stationary series in loop ",k,sep = ""))
    })
    
    if(inherits(possibleError, "error")) next
  }
  results <- list(series = multi_series, params = series_param)
  return(results)  
}


u = 1 # AR order
v = 0 # MA order
p = 1 # GARCH order
q = 1 # GARCH order

mu = 0 # (the linear part) the mean value, by default NULL
omega = 0.01 # the constant coefficient of the variance equation (default 1e-6)
alpha = 0.2
beta = 0.5
lags = 20
phi = 0.05

# length of time series
n = 300

# n = 300
K = 400  # simulate 400 time series
acf_matrix <- matrix(NA, ncol = 2, nrow = lags, dimnames = list(NULL,c('acf','pacf')))
# confidence interval
ci = c(-1.96/sqrt(n), 1.96/sqrt(n))

AIC_matrix <- matrix(NA, ncol = 2, nrow = 2, dimnames = list(NULL,c('garch_11','garch_22')))


# GARCH(1,1)
spec = garchSpec(model = list(mu = mu, omega = omega, alpha = alpha, beta = beta), cond.dist = "norm")
eta = garchSim(spec, n)
fitted.true <- garchFit(~ garch(1,1), 
                        data=eta,
                        cond.dist="norm",
                        include.mean = FALSE,
                        trace = FALSE)

d_phi.res <- sim_series_dif_phi(u, K, n, lags, ci, eta)
d_phi.multi_series <- d_phi.res$series
d_phi.series_param <- d_phi.res$params
colnames(d_phi.series_param) <- c("u", "phi1", "autocorrelation", "accepted")

f_phi.res <- sim_series_fix_phi(phi, u, K, n, lags, ci, eta)
f_phi.multi_series <- f_phi.res$series
f_phi.series_param <- f_phi.res$params
colnames(f_phi.series_param) <- c("u", "phi1", "autocorrelation", "accepted")

################  working on different phi data begin ################
#Identify $u$ for each series separately.  
target_series <- d_phi.multi_series
N <- nrow(target_series)
lags <- 15
acf_matrix <- matrix(NA, ncol = 2, nrow = lags, dimnames = list(NULL,c('acf','pacf')))
estimate_orders <- matrix(NA, ncol = 2, nrow = N, dimnames = list(NULL,c('p_hat','q_hat')))
residuals_matrix <- matrix(NA, ncol = n, nrow = N)

for(index in (1:N)){
  tar_ts <- target_series[index,]
  
  # ACF
  ro <- acf(tar_ts, lag.max = lags, plot = FALSE)
  # PACF
  phi <- pacf(tar_ts, lag.max = lags, plot = FALSE)
  
  acf_matrix[,'acf'] <- as.vector(ro$acf)
  acf_matrix[,'pacf'] <- as.vector(phi$acf)
  arma_orders <- order_identifier(acf_matrix, lags)  
  estimate_orders[index,] <- c(arma_orders$p_hat, arma_orders$q_hat)
  
  # assumption: series has AR component
  fit_ts = arima(tar_ts, order = c(arma_orders$p_hat,0,0))
  x = residuals(fit_ts)  
  residuals_matrix[index,] <- x
}

d_phi.eta_hat <- colMeans(residuals_matrix)
x <- d_phi.eta_hat
# ACF and PACF
arch_identify <- acf_identifier(x, lags)
garch_orders_hat <- order_identifier(arch_identify, lags)
garch_orders_hat$p_hat
garch_orders_hat$q_hat

fit_res_square <- garchFit(~ garch(2,2), 
                           data=x, 
                           cond.dist="norm", 
                           include.mean = FALSE, 
                           trace = FALSE)
llh <- -unname(fit_res_square@fit$llh)
AIC_matrix[1,2] <- (-2*(llh))/n + 2*(length(fit_res_square@fit$par))/n

fit_res_square <- garchFit(~ garch(1,1), 
                           data=x, 
                           cond.dist="norm", 
                           include.mean = FALSE, 
                           trace = FALSE)
llh <- -unname(fit_res_square@fit$llh)
AIC_matrix[1,1] <- (-2*(llh))/n + 2*(length(fit_res_square@fit$par))/n

d_phi.fitted <- fit_res_square

garch_res <- d_phi.eta_hat / fit_res_square@sigma.t
garch_res <- as.matrix(garch_res)
qqnorm(garch_res, main='Standard Residuals with Dynamic phi')
qqline(garch_res)

d_phi.est <- matrix(NA, nrow = nrow(target_series), ncol = 2)
d_phi.est[,1] <- estimate_orders[,1]
# focus on the selected simulation data
# remove shared GARCH from each series
estimated_ar <- target_series - d_phi.eta_hat

for(j in (1:nrow(estimated_ar))){
  ar_fit <- arima(estimated_ar[j,], 
                  order = c(estimate_orders[j,1],0,0), 
                  method = "ML")
  # extract the estimated coefficients
  est_coef <- ar_fit$coef[1:estimate_orders[j,1]]
  d_phi.est[j,2:(length(est_coef)+1)] <- est_coef
}
colnames(d_phi.est) <- c("est_orders","est_phi_1")

d_phi.true <- d_phi.series_param[,1:2]
colnames(d_phi.true) <- c("true_orders","true_phi_1")
################  working on different phi data end ################

################  working on fixed phi data begin ################
target_series <- f_phi.multi_series
N <- nrow(target_series)
lags <- 15
acf_matrix <- matrix(NA, ncol = 2, nrow = lags, dimnames = list(NULL,c('acf','pacf')))
estimate_orders <- matrix(NA, ncol = 2, nrow = N, dimnames = list(NULL,c('p_hat','q_hat')))
residuals_matrix <- matrix(NA, ncol = n, nrow = N)

for(index in (1:N)){
  tar_ts <- target_series[index,]
  
  # ACF
  ro <- acf(tar_ts, lag.max = lags, plot = FALSE)
  # PACF
  phi <- pacf(tar_ts, lag.max = lags, plot = FALSE)
  
  acf_matrix[,'acf'] <- as.vector(ro$acf)
  acf_matrix[,'pacf'] <- as.vector(phi$acf)
  arma_orders <- order_identifier(acf_matrix, lags)  
  estimate_orders[index,] <- c(arma_orders$p_hat, arma_orders$q_hat)
  
  # assumption: series has AR component
  fit_ts = arima(tar_ts, order = c(arma_orders$p_hat,0,0))
  x = residuals(fit_ts)  
  residuals_matrix[index,] <- x
}

f_phi.eta_hat <- colMeans(residuals_matrix)
x <- f_phi.eta_hat
# ACF and PACF
arch_identify <- acf_identifier(x, lags)
garch_orders_hat <- order_identifier(arch_identify, lags)
garch_orders_hat$p_hat
garch_orders_hat$q_hat

fit_res_square <- garchFit(~ garch(2,2), 
                           data=x, 
                           cond.dist="norm", 
                           include.mean = FALSE, 
                           trace = FALSE)
llh <- -unname(fit_res_square@fit$llh)
AIC_matrix[2,2] <- (-2*(llh))/n + 2*(length(fit_res_square@fit$par))/n

fit_res_square <- garchFit(~ garch(1,1), 
                           data=x, 
                           cond.dist="norm", 
                           include.mean = FALSE, 
                           trace = FALSE)
llh <- -unname(fit_res_square@fit$llh)
AIC_matrix[2,1] <- (-2*(llh))/n + 2*(length(fit_res_square@fit$par))/n

f_phi.fitted <- fit_res_square

garch_res <- f_phi.eta_hat / fit_res_square@sigma.t
garch_res <- as.matrix(garch_res)
qqnorm(garch_res, main='Standard Residuals with Fixed phi')
qqline(garch_res)

f_phi.est <- matrix(NA, nrow = nrow(target_series), ncol = 2)
f_phi.est[,1] <- estimate_orders[,1]
# focus on the selected simulation data
# remove shared GARCH from each series
estimated_ar <- target_series - f_phi.eta_hat

for(j in (1:nrow(estimated_ar))){
  ar_fit <- arima(estimated_ar[j,], 
                  order = c(estimate_orders[j,1],0,0), 
                  method = "ML")
  # extract the estimated coefficients
  est_coef <- ar_fit$coef[1:estimate_orders[j,1]]
  f_phi.est[j,2:(length(est_coef)+1)] <- est_coef
}
colnames(f_phi.est) <- c("est_orders","est_phi_1")

f_phi.true <- f_phi.series_param[,1:2]
colnames(f_phi.true) <- c("true_orders","true_phi_1")
################  working on fixed phi data end ################

################ Compare sigma
sig.true <- fitted.true@h.t     #- this is true sigma
sig.fixed.est <- f_phi.fitted@h.t    #- estimated sigma fixed phi
sig.diff.est <- d_phi.fitted@h.t    #- estimated sigma diff phi

sig.comparison <- melt(cbind(sig_400.fixed.est, sig_400.diff.est, sig.true))
sig.comparison$Var2 <- as.factor(sig.comparison$Var2)
ggplot(data=sig.comparison, aes(x=Var1, y=value, colour=Var2)) +
  labs(title = "True Sigma (Blue) vs Estimated Sigma") +
  xlab("") +
  ylab("Sigma") + 
  labs(color="Sigma") +  
  geom_line() 


################ Compare phi
phi.true <- rbind(f_phi.true, d_phi.true)
phi.est <- rbind(f_phi.est, d_phi.est)

phi.comparison <- melt(cbind(phi.est[,2], phi.true[,2]))
phi.comparison$Var2 <- as.factor(phi.comparison$Var2)
levels(phi.comparison$Var2) <- c("phi.est", "phi.true")
phi.comparison <- data.frame(phi.comparison)
phi.comparison$Var1 <- as.factor(phi.comparison$Var1)

ggplot(data=phi.comparison, aes(x=Var1, y=value), position=position_dodge(0.5)) +
  labs(title = "True Phi (Blue) vs Estimated Phi") +
  xlab("") +
  ylab("Value") + 
  labs(color="Phi") +  
  geom_point(aes(colour=Var2)) +
  geom_line(arrow = arrow(length=unit(0.30, "cm")), size = 0.2) +
  scale_x_discrete(breaks = c("1","40","80","120","160","200","240","280","320","360","400")) +
  theme(legend.position="none")







