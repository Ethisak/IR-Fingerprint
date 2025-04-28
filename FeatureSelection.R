#Libraries 
library(randomForest)
library(shapper)
library(readxl)
library(Metrics)
library(GA)
library(rpart)
library(caret)
library(dplyr)
library(stringr)
library(e1071)

#Random Forest - feature selection
dt <- read.csv("FPteor.csv")
set.seed(1)
col_names=colnames(dt[,2:102])
cv <- sample(1:10,nrow(dt),replace = T, prob = rep(0.1,10))

#Objective function
tree <- function(string){
  inc <- which(string == 1)
  if(length(inc)<20)return(-Inf)
  val <- c() 
  for (ii in 1:10) {
    test   <- dt[which(cv == ii), ]
    train  <- dt[which(cv != ii), ]
    xx <- train[2:102]
    yy <- train$LogP  
    X <- xx[,inc]
    
    mod <- randomForest(yy~.,data = X, ntree = 1000, mtry = 20, sampsize = nrow(X), replace = T, nodesize = 2)
    val <- c(val, -rmse(test$LogP, predict(mod,test)))
  }
  mean(val)
}

#Genetic algorithm
ga_GA_1 = ga(fitness = function(string) tree(string = string), 
             type = "binary", # optimization data type
             population = gabin_Population,
             crossover=gabin_uCrossover,  # cross-over method
             selection = gabin_lrSelection,
             elitism = 10, # number of best ind. to pass to next iteration
             pcrossover = 0.69,
             pmutation = function(...) ga_pmutation(...,p0 = 0.075, p = 0.5), # mutation rate prob
             popSize = 100, # the number of indivduals/solutions
             nBits = 101, # total number of variables
             names = col_names, # variable name
             run = 20, # max iter without improvement (stopping criteria)
             maxiter = 1000, # total runs or generations
             keepBest = TRUE, # keep the best solution at the end
             parallel = T, # allow parallel procesing
             monitor = T 
             
)

#Best solution
bestsol <- as.vector(which(ga_GA_1@solution[1,]==1))


#Random Forest - Testing solution on 100 different seeds
set.seed(1)
sid <- sample(1:10000000,100,replace = T)
val <- c()

for(i in sid){
  set.seed(i)
  cv <- sample(1:10,nrow(dt),replace = T, prob = rep(0.1,10))
  wynik <- c()
  for (ii in 1:10) {
    test   <- dt[which(cv == ii), ]
    train  <- dt[which(cv != ii), ]
    yy <- train$LogP
    xx <- train[,2:102]
    X <- xx[,bestsol]
    
    mod <- lm(yy~.,data = X)
    pr <- predict(mod,test)
    wynik <-  c(wynik, rmse(pr,test$LogP))
  }
  val <- c(val, mean(wynik))
}





#Support Vector Regression - feature selection
dt <- read.csv("FPteor.csv")
dt <- dt[,-c(2,79,which(colSums(dt) == 0))]
col_names=colnames(dt[,2:88])
set.seed(1)
cv <- sample(1:10,nrow(dt),replace = T, prob = rep(0.1,10))

#Objective function
tree <- function(string){
  inc <- which(string == 1)
  if(length(inc)<20)return(-Inf)
  val <- c() 
  for (ii in 1:10) {
    test   <- dt[which(cv == ii), ]
    train  <- dt[which(cv != ii), ]
    xx <- train[2:88]
    yy <- train$LogP  
    xx <- xx[,inc]
    
    if(length(which(colSums(xx) == 0)) > 0 | length(which(colSums(xx) == nrow(xx))) > 0){
      xyz <- c(which(colSums(xx) == 0),which(colSums(xx) == nrow(xx)))
      X <- xx[,-xyz]
    }else{
      X <- xx
    }
    
    mod <- svm(yy~., data = X, epsilon = 0.39, cost = 3, kernel = "radial")
    val <- c(val, -rmse(test$LogP, predict(mod,test)))
  }
  mean(val)
}


#Genetic algorithm
ga_GA_1 = ga(fitness = function(string) tree(string = string), 
             type = "binary", # optimization data type
             population = gabin_Population,
             crossover=gabin_uCrossover,  # cross-over method
             selection = gabin_lrSelection,
             elitism = 10, # number of best ind. to pass to next iteration
             pcrossover = 0.69,
             pmutation = function(...) ga_pmutation(...,p0 = 0.075, p = 0.5), # mutation rate prob
             popSize = 100, # the number of indivduals/solutions
             nBits = 87, # total number of variables
             names = col_names, # variable name
             run = 20, # max iter without improvement (stopping criteria)
             maxiter = 1000, # total runs or generations
             keepBest = TRUE, # keep the best solution at the end
             parallel = T, # allow parallel procesing
             monitor = T 
)

#Best solution
bestsol <- as.vector(which(ga_GA_1@solution[1,]==1))


#Testing solution on 100 different seeds
set.seed(1)
sid <- sample(1:10000000,10,replace = T)
val <- c()

for(i in sid){
  set.seed(i)
  cv <- sample(1:10,nrow(dt2),replace = T, prob = rep(0.1,10))
  wynik <- c()
  for (ii in 1:10) {
    test   <- dt[which(cv == ii), ]
    train  <- dt[which(cv != ii), ]
    yy <- train$LogP
    xx <- train[,2:88]
    xx <- xx[,bestsol]
    
    if(length(which(colSums(xx) == 0)) > 0 | length(which(colSums(xx) == nrow(xx))) > 0){
      xyz <- c(which(colSums(xx) == 0),which(colSums(xx) == nrow(xx)))
      X <- xx[,-xyz]
    }else{
      X <- xx
    }
    
    mod <- svm(yy~., data = X, epsilon = 0.39, cost = 3, kernel = "radial")
    pr <- predict(mod,test)
    wynik <-  c(wynik, rmse(pr,test$LogP))
  }
  val <- c(val, mean(wynik))
}
