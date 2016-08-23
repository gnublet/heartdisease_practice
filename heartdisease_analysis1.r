#####want to predict heart disease or not#####

#####load xgboost package#####
library(xgboost)

#####load heart dataset#####
data = read.table('heart.dat', header = F)
datared = data[,c(1,2,4,8,12,14)]
names(datared) = c('age','sex', 'rbp', 'maxhr', 'numvessels', 'hd')
head(datared)
train.x = datared[1:200,1:5]
train.y = datared[1:200, 6]-1
test.x = datared[201:270,1:5]
test.y = datared[201:270, 6]-1

#data(agaricus.train, package='xgboost')
#data(agaricus.test, package='xgboost')
#train <- agaricus.train
#test <- agaricus.test

#caret package may help spliting data in future.

str(train.x)
#dim(train$data)
#dim(test$data)

#####Explore the data######
head(train.x)
pairs(train.x)

library(GGally)
ggpairs(data)

#let's look at a reduced data set of only age, sex, blood pressure, max heart rate, # of major vessels, heartdisease
datared = data[,c(1,2,4,8,12,14)]
names(datared) = c('age','sex', 'rbp', 'maxhr', 'numvessels', 'hd')
head(datared)
g <- ggplot(datared, aes(age, hd, alpha = .25))
g + geom_point(aes(colour = factor(sex), shape = factor(sex)) )

ggpairs(datared)
#####Training######

dtrain = xgb.DMatrix(data = as.matrix(train.x), label = as.matrix(train.y))
bst = xgboost(data = dtrain, max.depth = 2, eta = 1, nthread = 2, nround = 2, objective = 'binary:logistic', verbose = 2)

#####Prediction#####
#now that we have trained models, we can classify new data
pred = predict(bst, as.matrix(test.x))

print(length(pred))
#dim(test$data)
print(head(pred))
#xgboost only does regression. Need to map this to {0,1} for classification

prediction = as.numeric(pred>.5)
print(head(prediction))

mean(prediction)
#####Mesuring model performance#####
#here we use average error as the metric
err = mean(as.numeric(pred>.5) != test.y)#need to examine splitting precedure 27% test error is high. (might be due to small dataset)
print(paste('test-error', err))

mean(prediction != test.y)
#very low error of 2 percent
#high error of 27 percent

#####advanced features#####
dtrain = xgb.DMatrix(data = as.matrix(train.x), label = as.matrix(train.y))
dtest = xgb.DMatrix(data = as.matrix(test.x), label = (test.y))
######Measuring learning progress#####
#

watchlist = list(train = dtrain, test = dtest)
bst = xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 2, nround = 2, watchlist = watchlist, eval.metric = 'error', eval.metric = 'logloss',objective = 'binary:logistic', verbose = 2)
#if train and test error don't match, should look for problems in how dataset was split (caret)

#####linear boosting#####
bst = xgb.train(data = dtrain, booster = 'gblinear', max.depth = 2, nthread = 2, nround = 2, watchlist = watchlist, eval.metric = 'error', eval.metric = 'logloss', objective = 'binary:logistic', verbose = 2)#no eta param
#does better than decision trees since there's nothing better than a linear algorithm to catch a linear link. Decision trees are better to catch nonlinear link between predictors and outcome
#linear booster does slightly better

#####manipulating xgb.Dmatrix
xgb.DMatrix.save(dtrain, 'dtrain.buffer')
dtrain2 = xgb.DMatrix('dtrain.buffer')
bst <- xgb.train(data=dtrain2, max.depth=2, eta=1, nthread = 2, nround=2, watchlist=watchlist, objective = "binary:logistic")

#information extraction
label = getinfo(dtest, 'label')#getinfo
pred = predict(bst, dtest)
err = as.numeric(sum(as.integer(pred>.5)!=label)/length(label))
print(paste('test-error=', err))


#FEATURE IMPORTANCE
importance_matrix = xgb.importance(names(datared), model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)
#number of vessels is the most important ~.5


xgb.dump(bst, with.stats = T)
xgb.plot.tree(model = bst)


#SAVE/LOAD models
# save model to binary local file
xgb.save(bst, "xgboost.model")

# load binary model to R
bst2 <- xgb.load("xgboost.model")
pred2 <- predict(bst2, test$data)

# And now the test
print(paste("sum(abs(pred2-pred))=", sum(abs(pred2-pred))))



# save model to R's raw vector
rawVec <- xgb.save.raw(bst)

# print class
print(class(rawVec))

# load binary model to R
bst3 <- xgb.load(rawVec)
pred3 <- predict(bst3, test$data)

# pred2 should be identical to pred
print(paste("sum(abs(pred3-pred))=", sum(abs(pred2-pred))))





#####simple linear regression
mod = lm(hd~., datared)
summary(mod)

#linear model shows that age is not significant, but the other predictors are, particularly number of vessels, mas heart rate, and sex.
plot(mod)

qqnorm(datared$hd)#definitely not normally distributed (since logistic)
