rm(list = ls())
setwd("C:/Users/Neel/Desktop/2023-2024/Spring 2024/Global Economics/R Functions")
pwd = getwd() # getting the path for present working directory (pwd)

##########################################################
# Adding required packages and functions
##########################################################
source("model_selection_function.R") # function for model selection
source("ols_function.R") # function for OLS estimation
source("t_test_function.R") # fucntion for t test
source("expanding_window_forecast_function.R") # function for forecasting using expanding windows
source("rolling_window_forecast_function.R") # function for forecasting using rolling windows
library("quantmod") # add quantmod to the list of Packages
install.packages("fBasics")
library("fBasics") # add fBasics to the list of Packages

#########################################################################################################
# Fetch Data for Housing Price index in the US
#########################################################################################################

getSymbols(Symbols ="CSUSHPISA",src = "FRED",warnings = FALSE) # Download Quarterly data for S&P/Case-Shiller U.S. National Home Price Index

Y <- as.matrix(CSUSHPISA[,1])
DY <- diff(Y)/Y[1:(dim(Y)[1]-1),]
head(DY)
tail(DY)

keep_data <- seq(from = as.Date("1990-01-01"), to = as.Date("2023-11-1"), by = "month")
DY_new = as.matrix(DY[as.Date(rownames(DY)) %in% keep_data,])
colnames(DY_new) = "House Price Index Growth"
n_obs = dim(DY_new)[1]
DY_new_date = as.Date(row.names(DY_new))

plot(x = DY_new_date, y = DY_new,xlab='time',ylab='House Price Index Growth',type='l',col="black")

basicStats(DY_new)

acf(DY_new,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data

Box.test(DY_new, lag = round(n_obs^(1/3)), type = "Ljung-Box") # applying Ljung and Box (1978) joint test of auto correlations

pacf(DY_new,lag=round(n_obs^(1/3)),main="House Price Index Growth") # command to obtain sample PACF of the data

#########################################################################################################
# select the number of lags and model checking
#########################################################################################################

results <- model_selection(round(n_obs^(1/3)),DY_new)
aic_values = results$AIC
bic_values = results$BIC
num_lags_aic = results$op_lag_AIC  
num_lags_bic = results$op_lag_BIC
num_lags_aic
num_lags_bic

num_lags = num_lags_aic
lags_DY_new = matrix(NA,nrow = n_obs, ncol = num_lags)
for (i in 1:num_lags) {
  lags_DY_new[(i+1):n_obs,i] = as.matrix(DY_new[1:(n_obs-i),1])
}
intercept = matrix(1,n_obs)
X = cbind(intercept,lags_DY_new)
y = DY_new
reg_result = ols(X[(num_lags+1):n_obs,],as.matrix(y[(num_lags+1):n_obs,1]))
beta_hat = reg_result$beta_hat

ar_coeff <- as.numeric(beta_hat[2:(num_lags+1)])
ma_coeff <- 0
ACF = acf(DY_new,lag=round(n_obs^(1/3)),plot = FALSE) # command to obtain sample ACF of the data
TACF <- ARMAacf(ar_coeff, ma_coeff, lag.max = round(n_obs^(1/3))) # command to obtain theorical ACF
plot(c(0:round(n_obs^(1/3))),ACF$acf,type='l',xlab='Lag',ylab='ACF',ylim=c(-1,1))
lines(0:round(n_obs^(1/3)),TACF,lty=2)
grid(nx = 4, ny = 4)

residuals = reg_result$u_hat # get the AR model residuals
acf(residuals,lag=round(n_obs^(1/3)),main = "residuals of House Price Index Growth") # command to obtain sample ACF of the data

Box.test(residuals, lag = round(n_obs^(1/3)), type = "Ljung-Box") # applying Ljung and Box (1978) joint test of auto correlations

#########################################################################################################
# Forecasting house price index using expanding windows
#########################################################################################################

lag_choice = NA
init_win_len = 120 # the first 10 years
num_step_ahead = 24 # 1 to 24 steps ahead forecastes 
prediction_results = expanding_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'aic')
yhat_f_aic <- prediction_results$forecast

y_f_aic <- prediction_results$actual_value
plot(x = DY_new_date[121:n_obs], y = y_f_aic,xlab='time',ylab='House Price Index Growth',type='l',col="black")
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,24],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_aic) - yhat_f_aic
rmsfe_ar_aic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ar_aic

lag_choice = NA
init_win_len = 120 # the first 10 years
num_step_ahead = 24 # 1 to 24 steps ahead forecastes 
prediction_results = expanding_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'bic')
yhat_f_bic <- prediction_results$forecast

y_f_bic <- prediction_results$actual_value
plot(x = DY_new_date[121:n_obs], y = y_f_bic,xlab='time',ylab='House Price Index Growth',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,24],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_bic
rmsfe_ar_bic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ar_bic

yhat_f_ave = (yhat_f_aic + yhat_f_bic)/2
forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_ave
rmsfe_ave = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ave

rmsfe_all_expanding = rbind(rmsfe_ar_aic,rmsfe_ar_bic,rmsfe_ave)
rmsfe_all_expanding

#########################################################################################################
# Forecasting house price index using rolling windows
#########################################################################################################

lag_choice = NA
init_win_len = 120 # the first 10 years
num_step_ahead = 24 # 1 to 24 steps ahead forecastes 
prediction_results = rolling_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'aic')
yhat_f_aic <- prediction_results$forecast

y_f_aic <- prediction_results$actual_value
plot(x = DY_new_date[121:n_obs], y = y_f_aic,xlab='time',ylab='House Price Index Growth',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,24],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_aic) - yhat_f_aic
rmsfe_ar_aic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ar_aic

lag_choice = NA
init_win_len = 120 # the first 10 years
num_step_ahead = 24 # 1 to 24 steps ahead forecastes 
prediction_results = rolling_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'bic')
yhat_f_bic <- prediction_results$forecast

y_f_bic <- prediction_results$actual_value
plot(x = DY_new_date[121:n_obs], y = y_f_bic,xlab='time',ylab='House Price Index Growth',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,24],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_bic
rmsfe_ar_bic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ar_bic

yhat_f_ave = (yhat_f_aic + yhat_f_bic)/2
forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_ave
rmsfe_ave = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))
rmsfe_ave

rmsfe_all_rolling = rbind(rmsfe_ar_aic,rmsfe_ar_bic,rmsfe_ave)
rmsfe_all_rolling


rmsfe_all = rbind(rmsfe_all_expanding,rmsfe_all_rolling)
rmsfe_all

