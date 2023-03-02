library(forecast)
library(imputeTS)
library(fpp3)

dat <- taylor
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
misrate <- 0.05 # missing rate
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
  e <- idx2+max_freq*3
  dat_af_gap <- dat[s:e]
  s <- idx1-max_freq*3
  if(idx1 < max_freq*3) {
    s = 0
    print("not enough data to forecast")
  }
  e <- idx1-1
  dat_be_gap <- dat[s:e]
  start <- idx1
  end <- idx2
} else {
  s <- idx1+1
  e <- idx1+max_freq*3
  dat_af_gap <- dat[s:e]
  s <- idx2-max_freq*3
  if(idx2 < max_freq*3) {
    s = 0
    print("not enough data to forecast")
  }
  e <- idx2-1
  dat_be_gap <- dat[s:e]
  start <- idx2
  end <- idx1
}

# sarima - imputation
x <- dat_be_gap %>% msts(seasonal.periods=c(48, 336))
sarima <- auto.arima(x, seasonal=TRUE, stepwise=TRUE, trace=TRUE)
fc <- forecast::forecast(x, h=n_gap)
imp_sarima <- dat_mar
imp_sarima[id.na] <- fc$mean
autoplot(imp_sarima)
checkresiduals(sarima)
summary(sarima)
# visualize sarima
viz_imp_1(dat_mar, fc$mean)
viz_imp_2(dat_mar, fc$mean, start-50, end+50)

# dhr - forecast
x <- dat_be_gap %>% msts(seasonal.periods=c(48, 336))
dhr <- auto.arima(x, xreg=cbind(fourier(x, K=c(24, 56))), seasonal=F, trace=TRUE)
fc <- forecast::forecast(dhr, xreg=cbind(fourier(x, K=c(24, 56))))
imp_dhr <- dat_mar
imp_dhr[id.na] <- fc$mean[1:n_gap]
autoplot(imp_dhr)
# check residuals & visualization
checkresiduals(dhr)
summary(dhr)
viz_imp_1(dat_mar, fc$mean[1:n_gap])
viz_imp_2(dat_mar, fc$mean[1:n_gap], start-50, end+50)

# STL+arima
x <- dat_be_gap %>% msts(seasonal.periods=c(48, 336))
stlm <- stlm(x, s.window='periodic', method='arima', trace=T)
fc <- forecast::forecast(stlm, h=n_gap)
imp_stlm <- dat_mar
imp_stlm[id.na] <- fc$mean
autoplot(imp_stlm)
checkresiduals(stlm)
viz_imp_1(dat_mar, fc$mean)
viz_imp_2(dat_mar, fc$mean, start-50, end+50)

# regression with dummy seasonal trend
x <- dat_be_gap %>% msts(seasonal.periods=c(48, 336))
tslm <- tslm(x ~ trend + season)
fc <- forecast::forecast(tslm, h=n_gap)
imp_sreg <- dat_mar
imp_sreg[id.na] <- fc$mean
autoplot(imp_sreg)
checkresiduals(tslm)
summary(tslm)
viz_imp_1(dat_mar, fc$mean)
viz_imp_2(dat_mar, fc$mean, start-50, end+50)

# regression with deterministic cosine trend
tslm2 <- tslm(x ~ trend + I(sin(2*pi*trend/336)) + I(cos(2*pi*trend/336))
              + I(sin(2*pi*trend/48)) + I(cos(2*pi*trend/48))
              + I(sin(2*pi*trend/24)) + I(cos(2*pi*trend/24))
              + I(sin(2*pi*trend/12)) + I(cos(2*pi*trend/12)))
fc <- forecast::forecast(tslm2, h=n_gap)
imp_sreg2 <- dat_mar
imp_sreg2[id.na] <- fc$mean
autoplot(imp_sreg2)
checkresiduals(tslm2)
summary(tslm2)
viz_imp_1(dat_mar, fc$mean)
viz_imp_2(dat_mar, fc$mean, start-50, end+50)

# regression with fourier terms
x <- dat_be_gap %>% msts(seasonal.periods=c(48, 336))
tslm3 <- tslm(x ~ trend + fourier(x, K=c(24,56)))
fc <- forecast::forecast(tslm3, h=336*3)
imp_sreg3 <- dat_mar
imp_sreg3[id.na] <- fc$mean[1:n_gap]
autoplot(imp_sreg3)
checkresiduals(tslm3)
summary(tslm3)
viz_imp_1(dat_mar, fc$mean[1:n_gap])
viz_imp_2(dat_mar, fc$mean[1:n_gap], start-50, end+50)

# evaluation
forecast::accuracy(imp_dhr, dat)
forecast::accuracy(imp_stlm, dat)
forecast::accuracy(imp_sarima, dat)
forecast::accuracy(imp_sreg, dat)
forecast::accuracy(imp_sreg2, dat)
forecast::accuracy(imp_sreg3, dat)