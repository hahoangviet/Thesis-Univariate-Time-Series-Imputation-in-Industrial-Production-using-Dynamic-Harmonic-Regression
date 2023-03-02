library(forecast)
library(imputeTS)
library(fpp3)

dat <- co2
n <- length(dat)
dominant_freq <- findfrequency(dat)
max_freq <- frequency(dat)

# visualize decomposition
autoplot(dat)
autoplot(decompose(dat)) # seasonal MA decomposition
autoplot(mstl(dat)) # MSTL decomposition shows seasonal 48 and seasonal 336

# dataset investigation
acf(dat) # acf shows cycle in the data

# missing data simulation
se <- 111
set.seed(se)
misrate <- 0.1 # missing rate
# MAR simulation
# a gap at random index
idx1 <- floor(runif(1, 0, n))
idx2 <- idx1 + floor(misrate*n)
if(idx2 > n) {idx2 <- idx1 - floor(misrate*n)}
mar <- c(1:n)
dat_mar <- dat
if(idx1 < idx2) {dat_mar[idx1:idx2] <- NA} else {dat_mar[idx2:idx1] <- NA}
ggplot_na_distribution(dat_mar)
statsNA(dat_mar)
which(is.na(dat_mar))

# get all missing data index
id.na <- which(is.na(dat_mar))
id.notna <- which(!is.na(dat_mar))
n_gap <- length(id.na)

# visualization MAR
viz_imp_1 <-function(dat_mis, imp) {
  plot(c(1:n_gap), imp, col='red') # imputed values
  points(c(1:n_gap), dat[id.na], col='blue') # true values
  lines(c(1:n_gap), imp, col='red')
  lines(c(1:n_gap), dat[id.na], col='blue')
  mean_imp <- rep(mean(dat_mis, na.rm=TRUE), n_gap)
  legend('topleft', legend=c('true values', 'imputed values'), 
         col=c('blue', 'red'), pch=c(1, 1))
  lines(c(1:n_gap), mean_imp, col='green')
}
viz_imp_2 <- function(dat_mis, imp, s, e) {
  dat_mis_imp <- dat_mis
  dat_mis_imp[id.na] <- imp
  ggplot_na_imputations(dat_mis[s:e], dat_mis_imp[s:e], dat[s:e])
}

# preapare data before and after gap
if(idx2 > idx1) {
  s <- idx2+1
  e <- n
  dat_af_gap <- dat[s:e]
  s <- 0
  e <- idx1-1
  dat_be_gap <- dat[s:e]
  start <- idx1
  end <- idx2
} else {
  s <- idx1+1
  e <- n  
  dat_af_gap <- dat[s:e]
  s <- 0
  e <- idx2-1
  dat_be_gap <- dat[s:e]
  start <- idx2
  end <- idx1
}

# dhr - step-by-step
x <- dat_mar
autoplot(x)

# step 1
# linear regression with Fourier terms, exclude missing data
ft <- fourier(x, K=5)
tslm <- tslm(x ~ ft)
summary(tslm)
# there is still information in remainder remainder component (residuals is not WN)
plot(tslm$residuals)
# seasonality component
sea_comp <- tslm$fitted.values
plot(sea_comp)
checkresiduals(tslm$residuals)
# impute the seasonal component
new <- data.frame(ft)
sea_imp <- predict.lm(tslm, new)
# seasonality imputation
plot(c(dat_mar))
lines(c(dat_mar), col='green')
lines(sea_imp, col='red')

# step 2
desea_res <- tslm$residuals - mstl(tslm$residuals)[,3]
x <- desea_res
# fit seasonally adjusted remainder to ARIMA
res_arima <- auto.arima(x, seasonal=F, stepwise=T, trace=T)
summary(res_arima)
# residuals of remainder component
plot(res_arima$residuals)
# fitted values by arima, capture autocorrelations
adj_dynamics_comp <- res_arima$fitted
plot(tslm$residuals)
lines(adj_dynamics_comp, col='red')
# residuals of residuals
plot(res_arima$residuals, col='red')
# there are still information in residuals of residuals
checkresiduals(res_arima$residuals)

# step 3
# impute the dynamic component
adj_dynamics_imp <- adj_dynamics_comp
# use Kalman filter and Smoother
adj_dynamics_imp <- na_kalman(adj_dynamics_comp, res_arima$model)
plot(adj_dynamics_imp)
# get the imputed dynamics
dynamics_imp <- tslm$residuals
dynamics_imp[id.na] <- adj_dynamics_imp[id.na] + mstl(tslm$residuals)[,3][id.na]
# dynamics component, red is the imputation for dynamics component
plot(dynamics_imp, col='red')
lines(tslm$residuals)
plot(dynamics_imp[id.na], col='red')
lines(dynamics_imp[id.na])

# step 4
# dhr imputation, red is for imputation
dhr_imp <- dat_mar
dhr_imp[id.na] <- sea_imp[id.na]+dynamics_imp[id.na]
plot(dhr_imp, col='red')
lines(dat, col='blue')

# seasonal component
viz_imp_2(dat_mar, sea_imp[id.na], start-50, end+50)
# trend component
viz_imp_2(dat_mar, dynamics_imp[id.na], start-50, end+50)
# visualization for imputation
viz_imp_1(dat_mar, dhr_imp[id.na])
viz_imp_2(dat_mar, dhr_imp[id.na], start-50, end+50)

plot(dhr_imp, col='red')
lines(dat_mar, col='blue')

# evaluation
forecast::accuracy(dhr_imp, dat)