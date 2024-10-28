#SIR bare bones
install.packages("deSolve")
library(deSolve)
sir_model <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  
  beta <- parameters[1]
  gamma <- parameters[2]
  
  dS <- -beta * S * I
  dI <- beta * S * I - gamma * I
  dR <- gamma * I
  
  return(list(c(dS, dI, dR)))
}

# Initial number of individuals
initial_state <- c(S = 990, I = 10, R = 0)  # 990 susceptible, 10 infected, 0 recovered

# Parameters
parameters <- c(beta = 0.3, gamma = 0.1)

# Time frame for the simulation
time_seq <- seq(0, 160, by = 1)  # from 0 to 160 days

# Run the SIR model
output <- ode(y = initial_state, times = time_seq, func = sir_model, parms = parameters)


# Convert output to a data frame
output <- as.data.frame(output)

# Basic plot
plot(output$time, output$S, type = "l", col = "blue", ylim = c(0, 1000), ylab = "Population", xlab = "Time (days)")
lines(output$time, output$I, col = "red")
lines(output$time, output$R, col = "green")
legend("right", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1)


#SIR with added regions
sir_model_multi_region <- function(time, state, parameters) {
  # Unpack states for each region
  S1 <- state[1]
  I1 <- state[2]
  R1 <- state[3]
  S2 <- state[4]
  I2 <- state[5]
  R2 <- state[6]
  
  # Unpack parameters for each region
  beta1 <- parameters[1]
  gamma1 <- parameters[2]
  
  beta2 <- parameters[3]
  gamma2 <- parameters[4]
  
  # SIR equations for region 1
  dS1 <- -beta1 * S1 * I1
  dI1 <- beta1 * S1 * I1 - gamma1 * I1
  dR1 <- gamma1 * I1
  
  # SIR equations for region 2
  dS2 <- -beta2 * S2 * I2
  dI2 <- beta2 * S2 * I2 - gamma2 * I2
  dR2 <- gamma2 * I2
  
  return(list(c(dS1, dI1, dR1, dS2, dI2, dR2)))
}

# Initial populations for each region
initial_state <- c(S1 = 990, I1 = 20, R1 = 0, 
                   S2 = 700, I2 = 5, R2 = 0)

# Parameters: beta and gamma for each region
parameters <- c(beta1 = 0.5, gamma1 = 0.1,  # Region 1-New york
                beta2 = 0.1, gamma2 = 0.05) # Region 2-Midwest

# Time frame for the simulation
time_seq <- seq(0, 160, by = 1)  # from 0 to 160 days

# Run the multi-region SIR model
output_multi <- ode(y = initial_state, times = time_seq, func = sir_model_multi_region, parms = parameters)

# Convert output to a data frame
output_multi <- as.data.frame(output_multi)

# Basic plot
plot(output_multi$time, output_multi$S1, type = "l", col = "blue", ylim = c(0, 1000), 
     ylab = "Population", xlab = "Time (days)", main = "SIR Model in Three Regions")
lines(output_multi$time, output_multi$I1, col = "red")
lines(output_multi$time, output_multi$R1, col = "green")

lines(output_multi$time, output_multi$S2, col = "darkblue")
lines(output_multi$time, output_multi$I2, col = "darkred")
lines(output_multi$time, output_multi$R2, col = "darkgreen")


legend("topright", legend = c("Susceptible 1", "Infected 1", "Recovered 1", 
                              "Susceptible 2", "Infected 2", "Recovered 2"),
       col = c("blue", "red", "green", "darkblue", "darkred", "darkgreen"),
       lty = c(1, 1, 1, 2, 2, 2), bty = "n")

#SIR model where species are moving from region 1 to region 2 as time increases
#As species move into region 2, the rate of infection will increase

sir_model_migration <- function(time, state, parameters) {
  # Unpack states for each region
  S1 <- state[1]
  I1 <- state[2]
  R1 <- state[3]
  S2 <- state[4]
  I2 <- state[5]
  R2 <- state[6]
  
  # Unpack parameters
  beta1 <- parameters[1]  # Infection rate for Region 1
  beta2 <- parameters[2]  # Infection rate for Region 2
  gamma <- parameters[3]   # Common recovery rate
  migration_rate <- parameters[4]  # Rate of migration from Region 1 to Region 2
  
  # SIR equations for Region 1
  dS1 <- -beta1 * S1 * I1 - migration_rate * S1
  dI1 <- beta1 * S1 * I1 - gamma * I1 - migration_rate * I1
  dR1 <- gamma * I1
  
  # SIR equations for Region 2
  dS2 <- migration_rate * S1 + migration_rate * I1 - beta2 * S2 * I2
  dI2 <- beta2 * S2 * I2 - gamma * I2
  dR2 <- gamma * I2
  
  return(list(c(dS1, dI1, dR1, dS2, dI2, dR2)))
}
# Initial populations for each region
initial_state <- c(S1 = 990, I1 = 10, R1 = 0, 
                   S2 = 980, I2 = 20, R2 = 0)

# Parameters: beta for each region, common gamma, and migration rate
parameters <- c(beta1 = 0.3, beta2 = 0.4, gamma = 0.1, migration_rate = 2)

# Time frame for the simulation
time_seq <- seq(0, 160, by = 1)  # from 0 to 160 days

# Run the SIR model with migration
output_migration <- ode(y = initial_state, times = time_seq, func = sir_model_migration, parms = parameters)

# Convert output to a data frame
output_migration <- as.data.frame(output_migration)

# Basic plot
plot(output_migration$time, output_migration$S1, type = "l", col = "blue", ylim = c(0, 1000), 
     ylab = "Population", xlab = "Time (days)", main = "SIR Model with Migration")
lines(output_migration$time, output_migration$I1, col = "red")
lines(output_migration$time, output_migration$R1, col = "green")

lines(output_migration$time, output_migration$S2, col = "blue", lty = 2)
lines(output_migration$time, output_migration$I2, col = "red", lty = 2)
lines(output_migration$time, output_migration$R2, col = "green", lty = 2)

legend("topright", legend = c("Susceptible 1", "Infected 1", "Recovered 1", 
                              "Susceptible 2", "Infected 2", "Recovered 2"),
       col = c("blue", "red", "green", "blue", "red", "green"),
       lty = c(1, 1, 1, 2, 2, 2), bty = "n")

