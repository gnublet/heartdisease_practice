sample(10)#w/o replacement
sample(10, replace = T)
#boostrap sample from same seq w/ probabilities that favor #s 1-2
prob1 = c(rep(.2, 2), rep(.075,8))
sample(1:10, replace=T, prob=prob1)

#matrix
y1 = matrix(round(rnorm(25,5)), nrow = 5)
y1

x1 = y1[sample(25,5)]
x1

y2 = matrix(round(rnorm(40,5)), ncol = 5)
y2

x2 = y2[sample(8,3), ]#3 rows
x2

##########
data = round(rnorm(100,5,3))
data[1:10]

#obtain 20 bootstrap samples
resamples = lapply(1:20, function(i)
  sample(data, replace = T))
#resample from data 20 times.

#display first boostrap sample
resamples[1]

#find median of each bootstrap sample
r.median = sapply(resamples, median)
r.median

r.mean=sapply(resamples,mean);hist(r.mean)

sqrt(var(r.median))
hist(r.median)

#all in 1
b.median = function(data,num){
  resamples = lapply(1:num, function(i)sample(data,replace=T))#get num resamples
  r.median = sapply(resamples, median)#find median of each resample
  r.mean = sapply(resamples, mean)
  std.err=sqrt(var(r.median))#std error
  list(std.err = std.err, resamples = resamples, medians = r.median, means = r.mean)
}

data1 = round(rnorm(500,5,3))#generate 100 samples, drawn from norm dist mean 5, sd 3
b1 = b.median(data1, 50)
b1$resamples[1]
b1$std.err
hist(b1$medians, 30)
hist(b1$means, 30) #recall Central Limit Theorem


b.median(rnorm(100,5,2), 50)$std.err

#built in boostrapping functions
library(boot)
data(bigcity)#random sampling of population (in 1000s of 49 cities in 1929 (u) and 1930 (x))
#49 cities are a random sample taken from 196 largest cities in the US.

#define ratio function
ratio = function(d,w) sum(d$x*w)/sum(d$u*w)

#using the built in boot function
boot(city, ratio, R=999, stype='w')
#data,statistic, # of bootstrap replicates, second argument of the statisticw = weights

samplemean = function(x,d){
  return(mean(x[d]))
}
b = boot(data1, statistic = samplemean, R = 50)
str(b)
print(b$t[,1])
print(sd(b$t[,1]))

b2 = boot(data1, statistic = samplemean, R = 5)
print(b2$t[,1])
print(sd(b2$t[,1]))#std error goes up with smaller number of resamples which is one of the reasons we wish to use bootstrapping.
