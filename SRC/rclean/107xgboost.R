#######################################
# Project: DuchenneConnect2016
# Purpose: XGBoost analysis of DuchenneConnect Data from October 2016 & October 2017
# Module: DuchenneConnect2016/SRC/rclean/104do.R
# Author: Richard T Wang
#
#
# Try to find best predictors of age at loss of ambulation LOA
#
# XGBoost only works on NUMERIC vectors! Categorical vectors need to be one-hot encoded
#
# https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/
# Up to this point, we dealt with basic data cleaning and data inconsistencies. To use xgboost package, keep these things in mind:
#   Convert the categorical variables into numeric using one hot encoding
# For classification, if the dependent variable belongs to class factor, convert it to numeric

library(xgboost)
library(caret)
library(Matrix)


# what's our missingness
missingness = sapply(all, function(x) sum(is.na(x))/length(x))
barplot(missingness)


# how many columns are factors?
lapply(all, class)
sapply(all, function(x) {
  if (class(x) == "factor") {
    return("factor")
  } else {
    return("notfactor")
  }
}) %>% table()


# ----- what features predict age LOA? -----
# use time_to_wheelchair and walking_01 to filter data
all.anymutation %>% xtabs(~ walking_01 + time_to_wheelchair , data =.)

# ----- data -----
# input data is all NOT walking with time_to_wheelchair
dat = all.anymutation %>% filter(walking_01==0, !is.na(time_to_wheelchair))

# ---- data cleanup -----
dat[is.na(dat)] <- "Missing"

# NOT FIXED attempt to find correlated variables -begin
cols.num = sapply(dat, is.numeric)  #logical vector
tmp = dat[, cols.num, with=FALSE] #with=FALSE allows logical vector to select
cols.na = sapply(tmp, is.na)
apply(cols.na, 1, sum)
tmp = tmp[, -cols.na, with=FALSE]
cor.mat = cor(x=tmp, use="pairwise.complete.obs")
findCorrelation(tmp, cutoff=0.9)

cor.mat[upper.tri(cor.mat)] <- 0
diag(cor.mat) <- 0
data.new <- tmp[,!apply(cor.mat,2,function(x) any(x > 0.99))]
head(data.new)
# attempt to find correlated variables -end


# how to count %of factors per column?
# we need 2 or more factors per column, otherwise remove the column
zzz = lapply(dat, function(x){  #not correct
  prop.table(table(x))
})


# Each columns must have 2 or more factors. But some columns are sparesely 
# populated, so after splitting into TRAIN and TEST sets, the TEST set
# might have only 1 factor for some columns.
#
# 1. Remove sparsely populated columns? This does not work well if 1 factor is very rare/
# what columns have at least 1 factor
gt1 = lapply(dat, function(x) { table(x) %>% length }) %>% as.data.frame  
col.del = which(gt1 == 1)
# remove cols
dat[,col.del] <- NULL

# ---- split data into TRAIN and TEST -----
set.seed(107)
# from caret, split 75% of data for training
inTrain <- createDataPartition(y = dat$time_to_wheelchair,
                               p = 0.75,
                               list=FALSE
)
str(inTrain)
#train = dat[inTrain, -c("time_to_wheelchair")]
train = dat[inTrain, ]
label = dat[inTrain, c("time_to_wheelchair")] %>% unlist(use.names=F)
#test = dat[-inTrain, -c("time_to_wheelchair")]
test = dat[-inTrain, ]
test.label = dat[-inTrain, c("time_to_wheelchair")] %>% unlist(use.names=F)

# Remove columns that have only 1 factor in TRAIN and TEST
# combine indicies and remvoe from both TRAIN and TEST so they have same number
# of columns.

gt1 = lapply(train, function(x) { table(x) %>% length }) %>% as.data.frame  
col.del = which(gt1 == 1)
gt2 = lapply(test, function(x) { table(x) %>% length }) %>% as.data.frame  
col.del2 = which(gt2 == 1)

col.del3 = c(col.del, col.del2) %>% sort %>% unique

train[,col.del3] <- NULL  # remove cols
test[,col.del3] <- NULL  # remove cols

#### Other columns to remove
# age at LOA is same as time to wheelchair, so remove!!
train[,'MusQ3.2'] <- NULL
test[,'MusQ3.2'] <- NULL

# Patient.ID is useless
train[,'Patient.ID'] <- NULL
test[,'Patient.ID'] <- NULL

# Remove Dates
train[,grep("Time",names(train))] <- NULL
test[,grep("Time",names(test))] <- NULL

# one hot encode
sparse_mat_train = sparse.model.matrix(time_to_wheelchair ~ .-1, data = train)
sparse_mat_test  = sparse.model.matrix(time_to_wheelchair ~ .-1, data = test)

dtrain <- xgb.DMatrix(data = sparse_mat_train, label = label) 
dtest <- xgb.DMatrix(data = sparse_mat_test, label=test.label)

# gbtree
params <- list(booster = "gbtree", 
               objective = "reg:linear", 
               eta=1.0, 
               gamma=1.0, 
               max_depth=6, 
               min_child_weight=1, 
               subsample=0.6, 
               colsample_bytree=0.6)
# Cross Vaidation
xgbcv <- xgb.cv( params = params, 
                 data = dtrain, 
                 nrounds = 100, 
                 nfold = 5, 
                 showsd = T, 
                 stratified = T, 
                 print.every.n = 10, 
                 early.stop.round = 20, 
                 maximize = F)


xgb1 = xgb.train (params = params, 
                  data = dtrain, 
                  nrounds = 100, 
                  watchlist = list(val=dtest,train=dtrain), 
                  print_every_n = 10, 
                  early_stop_round = 10, 
                  maximize = F , 
                  eval_metric = "rmse")



xgbpred <- predict (xgb1,dtest)
xgbpred

plot(xgbpred, test.label)
#confusion matrix
# library(caret)
# confusionMatrix (xgbpred, test.label)

#view variable importance plot
mat <- xgb.importance (feature_names = colnames(sparse_mat_train), model = xgb1)
print(mat)
png("RESULTS/xgboost/importance.png")
xgb.plot.importance (importance_matrix = mat[1:20]) 
dev.off()

# view tree
xgb.plot.tree(feature_names = colnames(sparse_mat_train), model = xgb1)
xgb.plot.deepness(model = xgb1)


#
# ------- which columns to use ---------

# ----- grid search param optimzation -----
library(dummy)
library(mlr)
fact_col = colnames(train)[sapply(train, is.character)]
for(i in fact_col) set(train,j=i,value = factor(train[[i]]))
for (i in fact_col) set(test,j=i,value = factor(test[[i]]))
#create tasks
traintask <- makeRegrTask (data = train, target = "time_to_wheelchair")
testtask <- makeRegrTask (data = test, target = "time_to_wheelchair")
#do one hot encoding
#traintask <- createDummyFeatures (obj = traintask,target = "time_to_wheelchair")
traintask <- createDummyFeatures (obj = traintask) 
#testtask <- createDummyFeatures (obj = testtask,target = "time_to_wheelchair")
testtask <- createDummyFeatures (obj = testtask)
#create learner
lrn <- makeLearner("regr.xgboost",predict.type = "response", fix.factors.prediction=TRUE)
lrn$par.vals <- list( objective="reg:linear", eval_metric="error", nrounds=100L, eta=0.1)

#set parameter space
params <- makeParamSet( makeDiscreteParam("booster",values = c("gbtree","gblinear")), 
                        makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                        makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                        makeNumericParam("subsample",lower = 0.5,upper = 1), 
                        makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))

#set resampling strategy
rdesc <- makeResampleDesc("CV",stratify = F,iters=5L)
#search strategy
ctrl <- makeTuneControlRandom(maxit = 10L)
library(parallel)
library(parallelMap) 
parallelStartSocket(cpus = detectCores())
#parameter tuning
mytune <- mlr::tuneParams(learner = lrn, task = traintask, resampling = rdesc, measures = acc, 
                     par.set = params, control = ctrl, show.info = T)
