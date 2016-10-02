data("AirPassengers")
class(AirPassengers)#shows that data is in a time series format
start(AirPassengers)
end(AirPassengers)

help("AirPassengers")
summary(AirPassengers)

#plot
plot(AirPassengers)
help(abline)
abline(reg = lm(AirPassengers~time(AirPassengers)))

#sanity check
head(AirPassengers)
tail(AirPassengers)
head(time(AirPassengers))
tail(time(AirPassengers))

cycle(AirPassengers)
#aggregate the cycles and display a year on year trend
plot(aggregate(AirPassengers, FUN = mean))
#boxplot across months
boxplot(AirPassengers~cycle(AirPassengers))

#year on year trend shows #passengers increases with time
#variance and mean is higher in summer months
#variance of the mean value of each month is small
#time series appears to be nonstationary

#We don't have constant variance, so let's do a log transform:
plot(log(AirPassengers))


#######test if stationary########
#install.packages('tseries')
library(tseries)
#augmented Dickey-Fuller test statistic is always negative.
#The more negative, the stronger the rejection of the null hypothesis that there is a unit root in a ts
#help(adf.test)
adf.test(diff(log(AirPassengers)), alternative = 'stationary', k=0)
#since p=.01 < .05, we reject the null hypothesis, so our ts is stationary.

#######find parameters for ARIMA model#######
#d=1 since we need 1 difference to make the series stationary
acf(log(AirPassengers))
#decay is slow which shows that population is not stationary
par(mfcol = c(1,2))
acf(diff(log(AirPassengers)))
pacf(diff(log(AirPassengers)))

#?somehow p,d,q = 0,1,1

#######fit an ARIMA model########
fit =arima(log(AirPassengers), order = c(0,1,1), seasonal = list(order = c(0,1,1), period = 12))
summary(fit)
#predict the future in 10 years
pred = predict(fit, n.ahead = 10*12)
par(mfcol = c(1,1))
ts.plot(AirPassengers, 2.718^pred$pred, log = 'y', lty = c(1,3))
