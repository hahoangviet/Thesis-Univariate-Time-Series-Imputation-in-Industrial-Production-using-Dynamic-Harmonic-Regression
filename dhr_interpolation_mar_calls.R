library(forecast)
library(imputeTS)
library(fpp3)
library(fpp2)
library(zoo)

dat <- calls
n <- length(dat)
dominant_freq <- findfrequency(dat)
max_freq <- frequency(dat)

# visualize decomposition
autoplot(dat)
autoplot(decompose(dat)) # seasonal MA decomposition
autoplot(mstl(dat)) # MSTL decomposition shows seasonal 169 and seasonal 845

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

# index for visualization
if(idx2 > idx1) {
  start <- idx1
  end <- idx2
} else {
  start <- idx2
  end <- idx1
}

# dhr - step-by-step
x <- dat_mar
autoplot(x)

# step 1
# linear regression with Fourier terms, exclude missing data
ft <- fourier(x, K=c(10,20))
tslm <- tslm(x ~ ft)
summary(tslm)
# there is still information in remainder remainder component (residuals is not WN)
plot(tslm$residuals)
checkresiduals(tslm$residuals)
# seasonality component
sea_comp <- tslm$fitted.values
# impute the seasonal component
new <- data.frame(ft)
sea_imp <- predict.lm(tslm, new)
# after imputation
plot(c(dat_mar), col='blue')
lines(sea_imp)
# seasonality imputation
plot(sea_imp[id.na])
lines(sea_imp[id.na])

# step 2
tslm_res <- tslm$residuals %>% msts(seasonal.periods=c(169, 845))
# mstl(tslm_res)[,3] (seasonal169), mstl(tslm_res)[,4] (seasonal845)
plot(mstl(tslm_res))
desea_res <- tslm_res - mstl(tslm_res)[,3] - mstl(tslm_res)[,4]
x <- desea_res
# fit seasonally adjusted remainder to ARIMA
res_arima <- auto.arima(x, seasonal=F, stepwise=T, trace=T, lambda='auto')
summary(res_arima)
# residuals of remainder component
plot(res_arima$residuals)
# there are still information in residuals of residuals
checkresiduals(res_arima)
# fitted values by arima, capture autocorrelations
adj_dynamics_comp <- res_arima$fitted
plot(tslm$residuals)
lines(adj_dynamics_comp, col='red')

# step 3
# impute the dynamics component
# use Kalman filter and Smoother
adj_dynamics_imp <- na_kalman(adj_dynamics_comp, res_arima$model)
adj_dynamics_imp <- msts(adj_dynamics_imp, c(169, 845))
autoplot(adj_dynamics_imp)

# get the imputed dynamics
dynamics_imp <- tslm_res
dynamics_imp[id.na] <- adj_dynamics_imp[id.na] +
  mstl(tslm_res)[,3][id.na] + mstl(tslm_res)[,4][id.na]

# dynamics component, red is the imputation for dynamics component
plot(dynamics_imp, col='red')
lines(tslm_res)
plot(dynamics_imp[id.na], col='red')
lines(dynamics_imp[id.na])

plot(dynamics_imp[id.na][1:100])
lines(dynamics_imp[id.na][1:100])

# dhr imputation, red is for imputation
dhr_imp <- dat_mar
# step 4: add imputated seasonal to imputated dynamics
dhr_imp[id.na] <- sea_imp[id.na]+dynamics_imp[id.na]
plot(dhr_imp, col='red')
lines(dat_mar)

plot(dat[id.na][1:500])
lines(dhr_imp[id.na][1:500], col='red')

# dynamics component
viz_imp_2(dat_mar, dynamics_imp[id.na], start-50, end+50)
# seasonal component
viz_imp_2(dat_mar, sea_imp[id.na], start-50, end+50)
# visualization imputation (dynamics component + seasonal component)
viz_imp_1(dat_mar, dhr_imp[id.na])
viz_imp_2(dat_mar, dhr_imp[id.na], start-50, end+50)

# evaluation
forecast::accuracy(dhr_imp, dat)