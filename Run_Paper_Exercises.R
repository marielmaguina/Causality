install.packages("timeSeries")
install.packages("NlinTS")
install.packages('RTransferEntropy')

library(timeSeries)
library(NlinTS)
library(RTransferEntropy)

# 1. The Granger Causality Test ------------------------------------------------------------------------------------------------

data=LPP2005REC

# Construct the causality model from the second column to the first one, 
# with a lag equal to 2,and without taking into account stationarity

model=causality.test(data[,1],data[,2],2,FALSE)
#This tests whether column 2 causes column 1 

#Compute the causality index,the F-test,and the p-value of the test

model$summary()
model$gci 
model$Ftest
model$pvalue

data[,1]
data[,2]
View(data)

# 2.. The VARNN model ------------------------------------------------------------------------------------------------------------

#The lag parameter 
lag=1 
#The training set 
train_data=data[1:(nrow(data)-1),] 
#Build and train the model 
model=varmlp(train_data,1,c(10,5),100)


#Predict the last row of the data 
predictions=model$forecast(train_data) 
#Show the predictions 
print(predictions[nrow(predictions),]) 

#Update the model (two observations are required at least since lag=1) 
#model$train(data[nrow(data)-lag:nrow(data)])
#model$train (data[nrow (data) - lag: nrow (data)])

# 3. The non-linear Granger causality test ------------------------------------------------------------------------------------------------------------

#Build and train the model 
model=nlin_causality.test(data[,1],data[,2],2,c(2),c(4))

#Compute the causality index, the F-test, and the p-value of the test 
model$summary() 
model$gci 
model$Ftest 
model$pvalue

#4. The discrete entropy ------------------------------------------------------------------------------------------------------------

#The entropy of an integer vector 
print(entropy_disc(c(3,2,4,4,3)))

#5. The continuous entropy ------------------------------------------------------------------------------------------------------------

#The entropy of the first column with k=3 
print(entropy_cont(data[,1],3))

#6. The discrete mutual information ------------------------------------------------------------------------------------------------------------

#Construct an integer data frame with 2columns 
df=data.frame(c(3,2,4,4,3),c(1,4,4,3,3)) 
#The mutual information between columns of df 
mi=mi_disc(df) 
print(mi)

#7. The continuous estimation of the mutual information ------------------------------------------------------------------------------------------------------------
#The mutual information between the two first columns of the data with k=3 
print(mi_cont(data[,1],data[,2],3))

#8. The discrete Transfer entropy  ------------------------------------------------------------------------------------------------------------
#The transfer entropy between two integer vectors with lag=1 
te=te_disc(c(3,2,4,4,3),c(1,4,4,3,3),1,1) 
print(te)

#9. The continuous estimation of the Transfer entropy ------------------------------------------------------------------------------------------------------------
#The transfer entropy between two columns with lag=1 and k=3 
te=te_cont(data[,1],data[,2],1,1,3) 
print(te)

# What kind of data fits the functions better? There's cont. and discrete.
# Q: But I am not sure what makes certain data discrete.

# When transfer entropy will perform better than Granger Causality?
# Transfer entropy is considered as a non-linear alternative for the Granger causality, since it does not model the 
# relationships between variables using a statistical model,instead,it is based on information theory.
# But this paper extends the Granger Causality test to non-linear cases. 

# How to simulate data think about it
# Come up with a function of t, add noise

#Email Laury, submit report by Dec 23rd


# Linear Data & Granger Causality ----------------------------------------------------------------------------------------------------------------------------------------------

set.seed(123) # Set seed
n <- 100 # Select number of observations

# Construct two time-series that are independent of each other

x <- numeric(n)
y <- numeric(n)

for (i in 2:n) { 
  y[i] <- 0.5 * y[i - 1] + rnorm(2)
  x[i] <- 0.7 * x[i - 1] + rnorm(1)
}

# Construct a time-series that makes x and y dependent on each other

z <- numeric(n)

for (i in 2:n) {
  z[i] <- y[i] + x[i-1]
}

z

# Perform a Granger Causality test 

model=causality.test(z,x,1,FALSE)
model$summary() # The p-value is very small so reject the null => x does cause z
model$pvalue

# Non-Linear Data & Granger Causality ----------------------------------------------------------------------------------------------------------------------------------------------

# Construct two time-series that are independent of each other

x_nl <- numeric(n)
y_nl <- numeric(n)

for (i in 2:n) { 
  y_nl[i] <- 0.5 * 1/exp (y_nl[i - 1]) + rnorm(2)
  x_nl[i] <- 0.7 * 1/exp (x_nl[i - 1]) + rnorm(1)
}

# Construct a time-series that makes x and y dependent on each other

z_nl <- numeric(n)

for (i in 2:n) {
  z_nl[i] <- (1/exp(y_nl[i])) * (1/ exp(x_nl[i-1]))
}

# Perform a Granger Causality test 

model=causality.test(z_nl,x_nl,1,FALSE)
model$summary() # The p-value is very small so reject the null => x does cause z

#The data is not linear but the test rejected the null. Somehow.
#Use sine or quadratic that is centered around the x, use a multiplicative error instead

# Linear Data & Transfer Entropy ----------------------------------------------------------------------------------------------------------------------------------------------

te_cont (z, x, 1, 1)
transfer_entropy(x,z)

# Non-Linear Data & Transfer Entropy ----------------------------------------------------------------------------------------------------------------------------------------------

te_cont (z_nl, x_nl, 1, 1)
#This is odd, the entropy for the non-linear data is lower. But not by a lot. 
transfer_entropy(x_nl,z_nl)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


