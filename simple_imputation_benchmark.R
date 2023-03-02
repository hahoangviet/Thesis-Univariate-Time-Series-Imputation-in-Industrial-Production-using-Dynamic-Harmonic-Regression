library(TSA)
library(imputeTestbench)
set.seed(100)
a <- impute_errors(taylor)
plot_errors(a, plotType = 'line')
plot_impute(dataIn = taylor[1:336], showmiss = T)
mcar <- impute_errors(taylor, smps='mcar', 
             methods=c('na.ma','na.interp','na.interpolation', 'na.locf', 'na.mean'),
             errorParameter = 'rmse',
             missPercentFrom = 10, missPercentTo = 90)
plot_errors(mcar, plotType = 'line')