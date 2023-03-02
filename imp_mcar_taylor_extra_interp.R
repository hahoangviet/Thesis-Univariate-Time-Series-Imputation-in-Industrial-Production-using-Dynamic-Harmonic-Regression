library(forecast)
library(imputeTS)
library(fpp3)
library(imputeTestbench)

dat <- taylor
n <- length(dat)

# visualize decomposition
autoplot(dat)
autoplot(decompose(dat)) # seasonal MA decomposition
autoplot(mstl(dat)) # MSTL decomposition shows seasonal 48 and seasonal 336

# dataset investigation
acf(dat) # acf shows cycle in the data

# missing data simulation
set.seed(100)
misrate <- 0.1 # 10% missing rate
# MCAR simulation
mcar = runif(n, min=0, max=1)
dat_mcar <- ifelse(mcar<misrate, NA, dat)
# distribution of missing values
ggplot_na_distribution(dat_mcar[1:336])
statsNA(dat_mcar)
which(is.na(dat_mcar))

# get all missing data index
x <- dat_mcar
id.na <- which(is.na(dat_mcar))
n_gap <- length(id.na)


# visualization MCAR
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
viz_imp_2 <- function(dat_mis, imp) {
  dat_mis_imp <- dat_mis
  dat_mis_imp[id.na] <- imp
  ggplot_na_imputations(dat_mis, dat_mis_imp, dat)
  ggplot_na_imputations(dat_mis[1200:1800], dat_mis_imp[1200:1800], dat[1200:1800])
}

# common simple imputation methods
a <- impute_errors(taylor)
plot_errors(a, plotType = 'line')
plot_impute(dataIn = taylor[1:336], showmiss = T)

# mean imputation
x <- dat_mcar
imp_mean <- na_mean(x)[id.na]
viz_imp_1(dat_mcar, imp_mean)
viz_imp_2(dat_mcar, imp_mean)

# na.locf
x <- dat_mcar
imp_locf <- na_locf(x)[id.na]
viz_imp_1(dat_mcar, imp_locf)
viz_imp_2(dat_mcar, imp_locf)


# na.interpolation
x <- dat_mcar
imp_interpolation <- na.interpolation(x)[id.na]
viz_imp_1(dat_mcar, imp_interpolation)
viz_imp_2(dat_mcar, imp_interpolation)

# na_interp
x <- dat_mcar
imp_interp <- na.interp(x)[id.na]
viz_imp_1(dat_mar, imp_interp)
viz_imp_2(dat_mcar, imp_interp)


# arima
x <- dat_mcar
arima <- auto.arima(x, trace=TRUE)
imp_arima_km <- na_kalman(x, arima$model)[id.na]
checkresiduals(imp_arima_km)
summary(arima)
viz_imp_1(dat_mcar, imp_arima_km)
viz_imp_2(dat_mcar, imp_arima_km)

# dhr
x %>% msts(seasonal.periods=c(48, 336)) -> x
dhr <- auto.arima(x, xreg=cbind(fourier(x, K=c(10, 15))), seasonal=FALSE, trace=TRUE)
imp_dhr_km <- na_kalman(x, dhr$model)[id.na]
checkresiduals(imp_dhr_km)
summary(dhr)
viz_imp_1(dat_mcar, imp_dhr_km)
viz_imp_2(dat_mcar, imp_dhr_km)

# evaluation
forecast::accuracy(imp_mean, dat[id.na])
forecast::accuracy(imp_locf, dat[id.na])
forecast::accuracy(imp_interpolation, dat[id.na])
forecast::accuracy(imp_interp, dat[id.na])
forecast::accuracy(imp_arima_km, dat[id.na])
forecast::accuracy(imp_dhr_km, dat[id.na])