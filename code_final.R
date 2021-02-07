### ST326 Financial Statistics Project ###

### Q1 ####

library(tidyquant)
library(quantmod)
library(bbmle) # Library for MLE

tickers <- c('MSFT','JPM','PG','AMZN','CVX','XOM','MS','JNJ','BAC','INTC','^GSPC') 

getSymbols(tickers, from="2010-11-27",to='2020-11-27') 

MSFTdata <- MSFT[,4] # Getting the closing prices for MSFT
JPMdata <- JPM[,4]
PGdata <- PG[,4]
AMZNdata <- AMZN[,4]
CVXdata <- CVX[,4]
XOMdata <- XOM[,4]
MSdata <- MS[,4]
JNJdata <- JNJ[,4]
BACdata <- BAC[,4]
INTCdata <- INTC[,4]
SPdata <- GSPC[,4]


df.prices <- data.frame(MSFTdata , JPMdata , PGdata , AMZNdata , CVXdata , XOMdata , MSdata , JNJdata ,BACdata , INTCdata ,SPdata) # Create a DF of prices
df.prices <- setNames(df.prices,c('MSFT','JPM','PG','AMZN','CVX','XOM','MS','JNJ','BAC','INTC','S&P500')) # Change DF column names
head(df.prices)

# Let's check if there are any missing data in our prices

total.na <- 0 # Setting a counter for missing values

for (stock in df.prices)
{
  total.na <- total.na + sum(is.na(stock)) # Note that FALSE is equivalent to 0 in R
}

sprintf('There are %s missing prices in total',total.na)


# Now that we confirmed that there are no missing values, we convert prices in log-prices

df.log.price <- log(df.prices)
head(df.log.price)

# Plotting log prices
ts.plot(df.log.price, col = 1:11,ylim=c(min(df.log.price),max(df.log.price)+0.3), ylab='Log price')
legend('top', c('MSFT','JPM','PG','AMZN','CVX','XOM','MS','JNJ','BAC','INTL','S&P500'), lty=c(rep(1,8),rep(3,3)), col=1:11 , cex=0.5 , horiz = T)

returns <- matrix(0 , nrow=dim(df.log.price)[1] - 1 , ncol=11 ) # Return matrix of dimensions (n_days - 1) by (n_stocks) since day 1 return is NaN

# Creating a returns matrix using the diff() function

for (i in seq(1,11)){
  returns[,i] <- diff(df.log.price[,i])
}


### Q2 ###


# Split data in training, validation and test sets

end.training <- dim(returns)[1] / 2
end.valid <- end.training + dim(returns)[1] / 4
end.test <- end.valid + dim(returns)[1] / 4

data.train <- returns[1 : end.training , ]
data.valid <- returns[(end.training+1) : end.valid , ]
data.test <- returns[(end.valid+1) : end.test, ]



## This function calculates lambda values for each stock in the TRAINING set and returns a Data Frame containing the lambda values

get.train.lambdas <- function(){
  
  lambdas <- t( c(rep(0,11)) ) # Create a row vector to store lambdas
  
  stocks <- seq(1,11) # Iterable to be used in the for loop 
  
  # This loop is for repeating the procedure for each stock
  for (s in stocks){
    
    # This function calculates the negative of the log-likelihood for each stock because mle2() works by minimising the function instead of maximising
    
    neg.LL <- function(lambda) {
      
      x <- data.train[,s] # Stock 's' returns in the training set
      n <- length(x) 
      sigma2 <- x^2
      
      for (t in 2:n){
        
        sigma2[t] <- (1-lambda)*x[t-1]^2 + lambda*sigma2[t-1]
      }
      
      - sum( log( 1 / (2*pi*sigma2)^0.5 ) - (x^2 / (2*sigma2)) ) # Return the negative of the log-likelihood function 
    }
    
    lambdas[s] <- coef( mle2(neg.LL, start = c(lambda=0.9)) )[1] # Add the lambda coefficient for stock 's' to the vector of lambdas
  }
  # Create a data frame consisting of the lambda values
  
  df.lambda <- data.frame(lambdas)
  df.lambda <- setNames(df.lambda , c('MSFT','JPM','PG','AMZN','CVX','XOM','MS','JNJ','BAC','INTL','S&P500'))
  df.lambda
}


train.lambdas <- get.train.lambdas()


### Function to check that results found above are indeed optimal ###

lambda.check <- function(stock_n = 1){ # Change for different stock
  
  param.space <- seq(0.01,0.99,by = 0.01)
  results <- matrix(0,1,99)
  x <- data.train[,stock_n] 
  sigma2 = x^2
  n = length(x)
  i <- 1
  
  for (lambda in param.space){
    
    s = 0
    
    for (t in 2:n){
      
      sigma2[t] = (1-lambda)*x[t-1]^2 + lambda*sigma2[t-1]
      
      s = s + log( 1 / (2*pi*sigma2[t])^0.5 ) - ( x[t]^2 / (2*sigma2[t]) )  # This is the LOG-Likelihood
    }
    
    results[1,i] = s
    i <- i + 1
  }
  
  plot(param.space,results)
  results
}

lambdas.check <- lambda.check(11)

##################################


vol.exp.sm <- function(x, lambda) {
  
  # Exponential smoothing of x^2 with parameter lambda
  
  sigma2 <- x^2 # Used to create of a vector of conditional variances (and also initiate first entry by setting it to x^2)
  n <- length(x)
  
  for (i in 2:n)
    sigma2[i] <- sigma2[i-1] * lambda + x[i-1]^2 * (1-lambda)
  
  sigma <- sqrt(sigma2)
  
  resid <- x/sigma
  resid[is.na(resid)] <- 0 # Set 'na' values to 0
  sq.resid <- resid^2
  
  
  list(sigma2=sigma2, sigma=sigma, resid = resid, sq.resid = sq.resid)
}


vol.train <- matrix(0 , dim(data.train)[1] , dim(data.train)[2] ) # A matrix to allocate the volatilities of the training set using the exponential smoothing function
vol.valid <- matrix(0 , dim(data.valid)[1] , dim(data.valid)[2] )
vol.test <- matrix(0 , dim(data.test)[1] , dim(data.test)[2] )



for (i in 1:11){
  vol.train[,i] <- vol.exp.sm( data.train[,i] , train.lambdas[1,i] )$sigma # Note also that by construction, conditional vol. at days 1 and 2 is the same for each asset
  vol.valid[,i] <- vol.exp.sm( data.valid[,i] , train.lambdas[1,i] )$sigma 
  vol.test[,i] <- vol.exp.sm( data.test[,i] , train.lambdas[1,i] )$sigma 
}



### Q4 ###

# Normalising our returns series using the volatilities obtained from exponential smoothing

data.train <- data.train / vol.train 
data.valid <- data.valid / vol.valid
data.test <- data.test / vol.test

x.train <- data.train[ , 1:10]
x.valid <- data.valid[ , 1:10]
x.test <- data.test[ , 1:10]


y.train <- data.train[ , 11]
y.valid <- data.valid[ , 11]
y.test <- data.test[ , 11]



# This function gives the OLS solution for a given matrix X and vector Y 
# This function also takes care of the particular case where 'X' is actually a vector (can occur when using only 1 PC as a regressor, in Q5)

beta.ols <- function(x, y){
  
  d <- dim(x)
  
  if (length(d)){ # Executes if d is not 'NULL', i.e has an integer length
    
    new.x <- matrix( c( rep(1, d[1]), x) , d[1], (d[2]+1) ) # Adding intercept term
    
  } else {
    new.x <- matrix( c( rep(1, length(x)), x) , length(x), 2) # length(x) by 2 matrix 
  }
  
  gram = t(new.x) %*% new.x # Gram matrix 
  
  if (length(d)){
    beta <- solve(gram) %*% t(new.x) %*% matrix(y, d[1], 1) # OLS solution for beta
    
  } else {
    beta <- solve(gram) %*% t(new.x) %*% matrix(y, length(x), 1)
  }
  
  beta
  
}


# This function gets the sharpe ratio for a single window length size

get.sharpe <- function(x, y, win = 150 , t0=250 , q=0){
  
  # 'win' is the rolling window size to be used (Note that win must be less than t0)
  # 't0' is the starting point in the training set (same as how it is defined in the notes)
  # 'q' is the number of S&P500 lags to include
  
  if (win>t0){
    print('Window size must be smaller than t0!!')
  }
  
  n <- length(y)
  
  y.pred <- rep(0, n-t0) # Vector of predicted S&P returns 
  y.true <- rep(0, n-t0) # Vector of actual S&P returns
  ret <- rep(0, n-t0) # Vector of returns
  position <- rep(1 , n-t0) # vector of signs to calculate our actual return 
  
  
  for (i in t0:n-1) { # (n-t0) iterations    
    
    target <- y[(i-win+1): i] # Gives a vector of length win
    
    
    if (q==1){
      predictors <- cbind( x[(i-win):(i-1) , ] , y[ (i-win):(i-1) ] ) # This part includes the S&P500 lag
      
    } else {
      
      predictors <- x[(i-win):(i-1) , ] 
    }
    
    beta <- beta.ols(predictors,target) # Get beta OLS solution of regressing Y(t) on X(t-1)
    
    if (q==1){
      pred <- sum( c(1,x[i,],y[i]) * beta ) # This is Y(t+1) obtained by multiplying X(t) and â OLS 
      y.pred[i-t0+1] <- pred
      
    } else { # This part executes when q=0
      
      pred <- sum( c(1, x[i,]) * beta )
      y.pred[i-t0+1] <- pred
      
    }
    
    y.true[i-t0+1] <- y[i+1]
    
    if ( pred < 0){ # If predicted Y(t+1) is negative, add a position of -1 to the position vector
      position[i-t0+1] <- -1
    }
    
  }
  
  ret <- position * y.true # This is element-wise multiplication of the two vectors
  
  sharpe.ratio <- sqrt(250) * ( mean(ret) / sqrt(var(ret)) ) 
  

  sharpe.ratio
}



# Function to compute the sharpe ratio for each window size on a given set
# This function will be used to find the optimal (aka tune) window size

tune.win <- function(x , y, window = seq(50,350,by=20) , q=0 ){
  
  # 'window' is a vector of different window sizes to try
  # 'q' is the number of S&P500 lags to use 
  
  
  sharpes <- c(rep(0,length(window)))
  i <- 1 # Set a counter
  
  for (d in window){
    sharpes[i] <- get.sharpe(x , y, win=d, t0=d+1 , q=q)
    i <- i + 1
  }
  
  s.mean <- mean(sharpes)
  
  # Presenting the results graphically
  plot(window , sharpes, pch=4, main='Sharpe ratio for different window lengths', xlab='Window length' , ylab='Sharpe ratio' )
  points(mean(window), s.mean, col = 4, pch = 19) # Plotting the mean sharpe ratio
  abline(lm(sharpes ~ window),col='red',lwd=1)
  
  list(sharpe = sharpes , mean.sharpe = s.mean)
}



### Finally, testing the functions above on the validation and test sets using different window sizes


# q=0 results

tune.win(x.train, y.train, seq(50,500,by=20) , q=0)

tune.win(x.valid, y.valid, seq(50,500,by=20) , q=0)

tune.win(x.test,y.test, seq(50,500,by=20) , q=0)

# q=1 results

tune.win(x.train, y.train, seq(50,500,by=20) , q=1)

tune.win(x.valid, y.valid, seq(50,500,by=20) , q=1)

tune.win(x.test,y.test, seq(50,500,by=20) , q=1)

# Checking for multicollinearity in the test set with no lags 
# Multicollinearity seems to be an issue here

summary(lm(y.test[102:251] ~ x.test[101:250,]))


### Q5 ### 

# In this question, we do not use any S&P500 lags

# Function to find the sharpe ratio as before, but now using the principal components matrix instead of the returns matrix

pc.sharpe <- function(x, y, win = 150 , t0=250 , n.pc=2){
  
  # 'win' is the rolling window size to be used (Note that win must be less than t0)
  # 't0' is the starting point in the training set (same as how it is defined in the notes)
  # 'n.pc' is the number of principal components to be included as regressors
  
  if (win>t0){
    print('Window size must be smaller than t0!!')
  }
  
  n <- length(y)
  
  y.pred <- rep(0, n-t0) # Vector of predicted S&P returns
  y.true <- rep(0, n-t0) # Vector of actual S&P returns
  ret <- rep(0, n-t0) # Vector of returns
  position <- rep(1 , n-t0) # vector of signs to calculate our actual return 
  
  if (n.pc==1){
    x <- x[ , 1] # Select only the first principal component
    
  } else {
    x <- x[ , 1:n.pc] # We select only the first 'n.pc' principal components from the matrix x 
  }
  
  for (i in t0:n-1) { # (n-t0) iterations
    
    target <- y[(i-win+1): i] # Gives a vector of length win 
    
    if (n.pc==1){
      predictors <- x[(i-win):(i-1)] # This is no longer a matrix, it is a vector since only 1 PC is used, so we must slice in 1-dimension
      
    } else {
      predictors <- x[(i-win):(i-1) , ] 
    }
    
    if (n.pc==1){
      beta <- beta.ols(predictors,target)
      
    } else {
      beta <- beta.ols(predictors,target) # Get beta OLS solution of regressing Y(t) on X(t-1)
    }
    
    
    if (n.pc==1){
      pred <- sum( c(1, x[i]) * beta ) # This is the prediction
      
    } else {
      pred <- sum( c(1, x[i,]) * beta )
    }
    
    y.pred[i-t0+1] <- pred
    
    y.true[i-t0+1] <- y[i+1]
    
    if ( pred < 0){ # If predicted Y(t+1) is negative, add a position of -1 to the position vector
      position[i-t0+1] <- -1
    }
    
  }
  
  ret <- position * y.true # This is element-wise multiplication of the two vectors
  
  sharpe.ratio <- sqrt(250) * ( mean(ret) / sqrt(var(ret)) )
  
  sharpe.ratio
}



###################################################################

pc.sharpe2 <- function(x, y, win = 150 , t0=250 , n.pc=2){
  
  # 'win' is the rolling window size to be used (Note that win must be less than t0)
  # 't0' is the starting point in the training set (same as how it is defined in the notes)
  # 'n.pc' is the number of principal components to be included as regressors
  
  if (win>t0){
    print('Window size must be smaller than t0!!')
  }
  
  n <- length(y)
  
  y.pred <- rep(0, n-t0) # Vector of predicted S&P returns
  y.true <- rep(0, n-t0) # Vector of actual S&P returns
  ret <- rep(0, n-t0) # Vector of returns
  position <- rep(1 , n-t0) # vector of signs to calculate our actual return 
  
  x.pca <- prcomp(x, center = TRUE , scale=TRUE) # Performing PCA to the X data
  
  if (n.pc==1){
    x <- (x.pca$x)[ , 1] # Select only the first principal component
    weights <- (x.pca$rotation)[,1] # Selecting the relevant weights obtained through PCA
    
  } else {
    x <- (x.pca$x)[ , 1:n.pc] # We select only the first 'n.pc' principal components from the matrix x 
    weights <- (x.pca$rotation)[,1:n.pc]
  }
  
  for (i in t0:n-1) { # (n-t0) iterations
    
    target <- y[(i-win+1): i] # Gives a vector of length win 
    
    if (n.pc==1){
      predictors <- x[(i-win):(i-1)] # This is no longer a matrix, it is a vector since only 1 PC is used, so we must slice in 1-dimension
      
    } else {
      predictors <- x[(i-win):(i-1) , ] 
    }
    
    if (n.pc==1){
      beta <- beta.ols(predictors,target)
      beta <- as.matrix(beta[-1] , length(beta) , 1) # removing the intercept term
      beta <- weights %*% beta 
      
    } else {
      beta <- beta.ols(predictors,target)
      beta <- as.matrix(beta[-1] , length(beta) , 1) # removing the intercept term
      beta <- weights %*% beta
    }
    
    
    if (n.pc==1){
      pred <- sum( c(x[i]) * beta ) # This is the prediction
      
    } else {
      pred <- sum( c(x[i,]) * beta )
    }
    
    y.pred[i-t0+1] <- pred
    
    y.true[i-t0+1] <- y[i+1]
    
    if ( pred < 0){ # If predicted Y(t+1) is negative, add a position of -1 to the position vector
      position[i-t0+1] <- -1
    }
    
  }
  
  ret <- position * y.true # This is element-wise multiplication of the two vectors
  
  sharpe.ratio <- sqrt(250) * ( mean(ret) / sqrt(var(ret)) )
  
  sharpe.ratio
}



win.pc.sharpe2 <- function(x , y, window = seq(50,350,by=20) , n.pc=2 ){
  
  x.sharpe <- c(rep(0,length(window)))
  i <- 1 # Set as a counter
  
  
  for (d in window){
    x.sharpe[i] <- pc.sharpe2(x , y, win=d, t0=d+1 , n.pc=n.pc)
    i <- i + 1
  }
  s.mean <- mean(x.sharpe)
  
  plot(window,x.sharpe, pch=4, main=paste('Sharpe ratio using ',n.pc,' principal components'), xlab='Window length' , ylab='Sharpe ratio' )
  points(mean(window), s.mean, col = 4, pch = 19) # Plotting the mean sharpe ratio as a blue dot
  abline(lm(x.sharpe ~ window),col='red',lwd=1)
  
  list(sharpe = x.sharpe , mean = s.mean)
}

##########################################################################################################

# Function to compute the sharpe ratio for each window size, and with 'n.pc' principal components

win.pc.sharpe <- function(x , y, window = seq(50,350,by=20) , n.pc=2 ){
  
  x.sharpe <- c(rep(0,length(window)))
  i <- 1 # Set as a counter
  
  x <- prcomp(x, center = TRUE , scale=TRUE)$x  # Perform PCA and set X to the principal components matrix 
  
  for (d in window){
    x.sharpe[i] <- pc.sharpe(x , y, win=d, t0=d+1 , n.pc=n.pc)
    i <- i + 1
  }
  s.mean <- mean(x.sharpe)
  
  plot(window,x.sharpe, pch=4, main=paste('Sharpe ratio using ',n.pc,' principal components'), xlab='Window length' , ylab='Sharpe ratio' )
  points(mean(window), s.mean, col = 4, pch = 19) # Plotting the mean sharpe ratio as a blue dot
  abline(lm(x.sharpe ~ window),col='red',lwd=1)
  
  list(sharpe = x.sharpe , mean = s.mean)
}


### Results obtained using this method seem to be better in general compared to Q4

# Using 1 PC

win.pc.sharpe(x.train , y.train, seq(50,500,by=20), n.pc=1)

win.pc.sharpe(x.valid , y.valid, seq(50,500,by=20), n.pc=1)

win.pc.sharpe(x.test , y.test, seq(50,500,by=20), n.pc=1)

# Using 2 PCs

win.pc.sharpe(x.train , y.train, seq(50,500,by=20), n.pc=2)

win.pc.sharpe(x.valid , y.valid, seq(50,500,by=20), n.pc=2)

win.pc.sharpe(x.test , y.test, seq(50,500,by=20), n.pc=2)

