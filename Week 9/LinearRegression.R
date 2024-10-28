####################
#  Oct 16 Meeting  #
# Estimating slope #
####################

#LinearRegression
#Make your points
x <-c(1, 2, 3, 4)
y <-c(2, 3, 5, 7)

#Make your line 
model <- lm(y ~ x)
plot(x, y, pch = 19, col = "orchid")
abline(model, col = "orchid1", lwd = 2)

#Now lets assume we don't know the line of best fit
# Establish unknown variables
best_slope <- NA
best_intercept <- NA
#This value represents the min residual sum of squares(RSS), set to infinity to not limit the value/
min_rss <- Inf

slopes<- seq(-10, 10, by = 0.1)
intercepts<- seq(0, 10, by = 0.5)
like_matrix<- matrix(0, nrow = length(slopes), ncol = length(intercepts))

like<- function(y, x, slope, intercept, error=0.1){
 means<- intercept + slope * x
 probs<- dnorm(y, means, error, FALSE)
 ret <- prod(probs)
 return(ret)
}

max_like <- -Inf
# For loop through a different slopes and intercepts
# Varying slope from -5 to 5
### I want to keep it reasonable and measure every 0.1
for (slope in slopes) {        
  # Varying intercept from 0 to 10 
  #Y values range from 2-7 so I want to cover those and then some
  #I measure every 0.5 now to speed up the processes a litle bit
  for (intercept in intercepts) {
    # Calculate predicted values
    #Y = mx+b based on previous for loops
    predicted_y <- slope * x + intercept
    # Calculate residuals and RSS
    # sum of (y - yhat) squared
   ## rss <- sum((y - predicted_y)^2)
    mylike<- like(y, x, slope, intercept, 0.1)
    
    if (mylike > max_like){
      max_like <- mylike
      best_slope <- slope
      best_intercept <- intercept
    }
        
    # Check if this model is better
  #  if (rss < min_rss) {
   #   min_rss <- rss
    #  best_slope <- slope
    #  best_intercept <- intercept
    #}
  }
}

# Output the best slope and intercept from the if statement above
cat("Best Slope:", best_slope, "\n")
cat("Best Intercept:", best_intercept, "\n")
cat("Maximum Likelihood:", max_like, "\n")

# Plot the original data and the new best fit line
plot(x, y, pch = 19, col = "orchid")
abline(a = best_intercept, b = best_slope, col = "orchid4", lwd = 2) 


# Next try to use optim to do a proper max_like on your function
####Today 10-16 you did a grid search



###############################
#        Oct 23 Meeting       #
# Estimating slope with optim #
###############################

#Used the tutorial https://www.magesblog.com/post/2013-03-12-how-to-use-optim-in-r/
#Turn your points into data frame
dat=data.frame(x=c(1, 2, 3, 4), 
               y=c(2, 3, 5, 7))
#Goint to use RSS
min.RSS <- function(data, par) {
  with(data, sum((par[1] + par[2] * x - y)^2))
}
#Use optim
(result <- optim(par = c(0, 1), fn = min.RSS, data = dat))
#plot
plot(y ~ x, data = dat, main="Least square regression", pch = 19, col = "orchid")
abline(a = result$par[1], b = result$par[2], col = "orchid4")

#############################################
#Grace tries to do it on her own with chatgtp
#############################################
# Example data
x <- c(1, 2, 3, 4)
y <- c(2, 3, 5, 7)

# Define the negative log-posterior function
neg_log_posterior <- function(params, x, y) { #params is expected to be a vector containing the model parameters (intercept, slope, and log of the standard deviation)
  alpha <- params[1]  # Intercept
  beta <- params[2]   # Slope
  sigma <- exp(params[3])  # Sigma: The standard deviation is transformed using the exponential function to ensure it remains positive.
  
  # Likelihood
  likelihood <- sum(dnorm(y, mean = alpha + beta * x, sd = sigma, log = TRUE))  #What is the likelihood of the observed data under a normal distribution. 
  
  # Priors
  #Alpha and Beta are modeled with normal distributions centered at 0 with a standard deviation of 10.
  #Sigma uses a half-Cauchy distribution to model the prior for the standard deviation, ensuring it is always positive.
    #half-Cauchy distribution is used to keep values positive but extend essentially to inf.
  #prior_alpha <- dnorm(alpha, mean = 0, sd = 10, log = TRUE)
  #prior_beta <- dnorm(beta, mean = 0, sd = 10, log = TRUE)
  #prior_sigma <- dcauchy(sigma, location = 0, scale = 2, log = TRUE)
  
  # Total negative log-posterior
  return(- (likelihood)) #+ prior_alpha + prior_beta + prior_sigma))  #the function calculates the total negative log-posterior by summing the likelihood and the prior terms, then negating this sum
}

# Initial values for alpha, beta, and log(sigma)
init_params <- c(alpha = 0, beta = 0, log_sigma = log(1))

# Optimize the negative log-posterior
optim_result <- optim(init_params, neg_log_posterior, x = x, y = y)
  #The optim function starts with the initial parameters provided in init_params and uses 
  #optimization algorithms (default is the Nelder-Mead method) to adjust the parameters
  #It tries different parameters until the neg_log_posterior is minimized

# Extract results
alpha_est <- optim_result$par[1]
beta_est <- optim_result$par[2]
sigma_est <- exp(optim_result$par[3])  # Back-transform sigma
cat("Estimated Intercept:", alpha_est, "\n")
cat("Estimated Slope:", beta_est, "\n")
cat("Estimated Sigma:", sigma_est, "\n")

# Plot the original data
plot(x, y, pch = 19, col = "orchid")
abline(a = alpha_est, b = beta_est, col = "orchid4", lwd = 2)

#Questions: What is Nelder-Mead method? What are gradient-based optimization methods?

