#################################
#  Birth Death Model Simulator  #
#################################
#          Grace Made           #
#################################

#Set seed allows for you to control randomness, stats tool
#Population is a vector
set.seed(0)

p_birth <- 0.1
p_death <- 0.05
p_no_event <- 1 - p_birth - p_death 
#Probability 
steps <- 100
population <- 1

for (i in 1:steps) {
  #If you still are alive draw again  
  if (tail(population, 1) > 0) {
    #Decide what happens, something or nothing    
    temp <- sample(x = c(-1, 0, 1),
                   size =1,
                   prob = c(p_death, p_no_event, p_birth))
    #Apply
    population <- c(population,
                    tail(population, 1) + temp)
    next
  }
  #Changed the code to be 0 if there is no members of the population
  else {population[i + 1] <- 0}
  #If you have nobody. can have no event of a birth event  
  # if(tail(population, 1) == 0){
  #   temp <- sample(x = c(0,1),
  #                  size=1,
  #                  prob = c(p_no_event / (p_no_event + p_birth),
  #                           p_birth / p_no_event + p_birth))
  #   population <- c(population,
  #                   tail(population, 1) + temp)
}
plot(population, type = "l")

##############################
#    BiSSE Model Simulator   #
##############################
#        Raymond Made        #
##############################

#Set seed allows for you to control randomness, stats tool
#Population is a vector
set.seed(0)

birth1 <- 0.1
birth2 <- 0.1
death1 <- 0.05
death2 <- 0.05
trans_1_2 <-0.5
trans_2_1 <-0.5
#Probability 
#Change steps to draw a waiting time
steps <- 100
curr_time <- 0
max_taxa <- 100
#Change to have binary states
population <- c(1,1)

#data frame containing info about taxa and past events and time
taxa_table_df <- data.frame(
  Taxa_1 = population[1],
  Taxa_2 = population[2],
  Time = curr_time,
  Event = 'NA'
)

while(curr_time <= steps && sum(population) <= max_taxa && sum(population)!=0){
  #find the birth of state 1
  birth_prob1 <- population[1] * birth1
  #find the birth of state 2
  birth_prob2 <- population[2] * birth2
  #find the death of state 1
  death_prob1 <- population[1] * death1
  #find the death of state 2
  death_prob2 <- population[2] * death2
  #Q-matrix
  #
  #        1   2 
  # q= 1 [-, rate]
  #   2 [rate, 1]
  trans_matrix <- matrix(c(0, population[1]*trans_1_2, population[2]*trans_2_1, 0), nrow = 2, byrow = TRUE)

  all_rates <- c(birth_prob1, birth_prob2, death_prob1, death_prob2, trans_matrix[1,2], trans_matrix[2,1])
  next_time<-  rexp(1, rate=1/sum(all_rates))
  print(next_time)
  curr_time=101
}


# What is next?
# Draw the what is going to happen - uniform




#####################################################
# Make a BiSSE model from scratch based on inference#
#####################################################
library(ape)
# Load your tree
tree <- read.tree("/Users/gracecoppinger/Desktop/RB_bisse_tutorial/anolis.phy")

# Load your data
data <- read.csv("/Users/gracecoppinger/Desktop/RB_bisse_tutorial/anolisDataAppended.csv", row.names = 1)

#Make ecomorphs into a state
ecomorph <- data[, "ecomorph"]
names(ecomorph) <- rownames(data)

# Create binary states: 1 for "TG", 0 for all others
states <- ifelse(ecomorph == "TG", 1, 0)

# Name the states based on species
names(states) <- rownames(data)

# Filter states to only include species present in the tree
states <- states[names(states) %in% tree$tip.label]

# Check the result
print(states)

# Define your parameters
lambda1 <- 0.10  # Speciation rate for state 1
lambda0 <- 0.05 # Speciation rate for state 0
mu1 <- 0.01     # Extinction rate for state 1
mu0 <- 0.02     # Extinction rate for state 0

#Create a transition matrix for the BiSSE model. Basically the Q matrix
trans_matrix <- matrix(c( lambda1, lambda0, mu1, mu0), nrow = 2, byrow = TRUE)
#assign names to the matrix
rownames(trans_matrix) <- c("State1", "State0")
colnames(trans_matrix) <- c("Speciation", "Extinction")

# Define a function to calculate the likelihood of the tree given the transition rates
calculate_likelihood <- function(tree, states, trans_matrix) {
  likelihood <- 1.0
#We start with a liklihood of 1 that gets multiplied by probabilities along the edges (connection between two nodes) of the tree.
 
   # For loop that moves through each edge in the tree
  for (i in 1:nrow(tree$edge)) {
    # ID the nodes
    child <- tree$edge[i, 2]  # Child node
    parent <- tree$edge[i, 1]  # Parent node
    
    # Check if parent and child states exist
    if (parent <= length(states) && child <= length(states)) {
      parent_state <- states[parent]
      child_state <- states[child]
      #Checks that both nodes have valid states. Retrieves the states for both the parent and child nodes.
      
      # Calculate the transition probabilities
      if (parent_state == child_state) {
        rate <- trans_matrix[parent_state + 1, 1]  # Speciation
      } else {
        rate <- trans_matrix[parent_state + 1, 2]  # Extinction
      }
      #Determines whether the transition is a speciation (same state) or extinction (different state) and selects the appropriate rate from the transition matrix.
      
      # Update the likelihood
      ####Updates the overall likelihood by multiplying it by the probability of the observed transition along the edge. The expression exp(-rate * tree$edge.length[i]) calculates the likelihood of surviving through the edge length given the rate.
      likelihood <- likelihood * exp(-rate * tree$edge.length[i])
    }
  }
  
  return(likelihood)
}

#Create a vector to store the four parameters in a single object. Allows you to pass all four parameters together as a single object
initial_params <- c(lambda1, lambda0, mu1, mu0)

# Define a function to minimize (negative log-likelihood)
neg_log_likelihood <- function(params) {
  lambda1 <- params[1]
  lambda0 <- params[2]
  mu1 <- params[3]
  mu0 <- params[4]
# Parameter Preparation: This part of the function is setting up the specific parameters that will be used later in the likelihood calculation. The function is likely intended to calculate the negative log-likelihood of a model given these parameters.
  #After extracting these parameters, the function would typically proceed to calculate the negative log-likelihood based on the model, the tree structure, and the observed data (e.g., states).
  # A lower value of NLL indicates a better fit of the model to the data.
  
  # Update transition matrix
  trans_matrix <- matrix(c(lambda1, lambda0, mu1, mu0), nrow = 2)
  
  return(-calculate_likelihood(tree, states, trans_matrix))
}
#calculates the likelihood of observing the states given the phylogenetic tree and the transition matrix.The negative sign in front of the likelihood indicates that you are returning the negative log-likelihood, which is a common approach when performing optimization (as mentioned previously).
#this code creates a transition matrix based on speciation and extinction rates, calculates the likelihood of the data given the model, and returns the negative of that likelihood for optimization purposes.

# Optimize parameters
result <- optim(initial_params, neg_log_likelihood)
#The primary goal of this line of code is to fit the BiSSE model to your data by finding the best estimates for the speciation and extinction rates (represented by the parameters in initial_params).

#Results
estimated_params <- result$par
print(estimated_params)

#Plot
library(ggplot2)
# Parameter names
param_names <- c("lambda1", "lambda0", "mu1", "mu0")
# Create a bar plot for estimated parameters
barplot(estimated_params, names.arg = param_names, col = "skyblue",
        main = "Estimated Parameters from BiSSE Model",
        ylab = "Rate") 

###################################
# Make a HiSSE model from scratch #
###################################
# Hypothesis: 
#There are hidden ecological states (e.g., habitat types, adaptive strategies) that influence 
#the evolution of ecomorphs. For instance, certain ecomorphs may experience different speciation 
#and extinction rates due to the hidden states.
#################################################################################################


#################################
#     BiSSE Model Simulator     #
#################################
#          Grace Made           #
#################################
#   Help from Raymond and Sean  #
#################################
set.seed(0)
#rates
birth1<- 0.1
birth2<- 0.2
death1<- 0.03
death2<- 0.03
q12<- 0.01
q21<- 0.01

#Preventions
stop_time <- 1000
curr_time <- 0
max_taxa <- 100

sim_BiSSE = function(birth1, birth2, death1, death2, q12, q21, stop_time, curr_time, max_taxa){
  
  #Change to have binary states
  population <- c(1,1)
  
  #Make a table
  taxa_table <- data.frame(
    Taxa_1 = population[1],
    Taxa_2 = population[2],
    Time = curr_time,
    Event = 'Start'
  )
  #The condition for this while loop prevents the current time from exceding steps, population exceding max taxa, and the population being below zero
  while (curr_time <= stop_time && sum(population) <= max_taxa && sum(population)!=0) {
    #find the prob of birth in state 1
    birth_prob1 <- population[1] * birth1
    #find the prob of birth in state 2
    birth_prob2 <- population[2] * birth2
    #find the prob of death in state 1
    death_prob1 <- population[1] * death1
    #find the prob of death in state 2
    death_prob2 <- population[2] * death2
    #Q matrix
    trans_matrix<- matrix(c(0, population[1]*q12, population[2]*q21, 0), nrow = 2, byrow = TRUE)
    #Make all the rate relative
    all_rates <- c(birth_prob1, birth_prob2, death_prob1, death_prob2, trans_matrix[1,2], trans_matrix[2,1])
    next_time <- rexp(1, rate=1/sum(all_rates))
    #Make it end at 100
    curr_time = curr_time + next_time
  
  
  relative <- cumsum(all_rates / sum(all_rates))
  
  #Draw an event
  draw <- runif(1)
  if (draw < relative[1]) {
    print("speciation state 1")
    population[1] = population[1] + 1
    event = "speciation state 1"
  } else if (draw < relative[2]) {
    print("speciation state 2")
    population[2] = population[2] + 1
    event = "speciation state 2"
  } else if (draw < relative[3]) {
    print("death state 1")
    population[1] = population[1] - 1
    event = "death state 1"
  } else if (draw < relative[4]) {
    print("death state 2")
    population[2] = population[2] - 1
    event = "death state 2"
  } else if (draw < relative[5]) {
    print("transition from state 1 to 2")
    population[1] = population[1] - 1
    population[2] = population[2] + 1
    event = "transition from state 1 to 2"
  } else if (draw < relative[6]) {
    print("transition from state 2 to 1")
    population[2] = population[2] - 1
    population[1] = population[1] + 1
    event = "transition from state 2 to 1"
  }
  
  
  taxa_table_row= c( population[1],
                     population[2],
                         curr_time,
                             event 
      )
  
  taxa_table = rbind(taxa_table, taxa_table_row)
  
  
  
  }
  return(taxa_table)

}

sim_BiSSE(birth1, birth2, death1, death2, q12, q21, stop_time, curr_time, max_taxa)
  
###################################
# Make a BiSSE model from scratch #
###################################
set.seed(0)
#rates
birth1<- 0.1
birth2<- 0.9
death1<- 0.5
death2<- 0.3
q12<- 0.4
q21<- 0.8
#Preventions
stop_time <- 1000
curr_time <- 0
max_taxa <- 100

sim_BiSSE = function(birth1, birth2, death1, death2, q12, q21, stop_time, curr_time, max_taxa){
  
  #Change to have binary states
  population <- c(1,1)
  
  # Running list of IDs
  taxa_id <- rep(list(1), sum(population))
  
  # Taxa Relationships contains the parent-daughter relationships between taxa that will be added as the process is simulated
  taxa_relationships <- data.frame(
    Taxa_ID = unlist(taxa_id) ,
    d_1 = rep(0, length(taxa_id)),
    d_2 = rep(0, length(taxa_id)),
    p = rep(0, length(taxa_id))
  )
  
  
  #Make a table
  taxa_table <- data.frame(
    Taxa_1 = population[1],
    Taxa_2 = population[2],
    Time = curr_time,
    Event = 'Start'
  )
  #The condition for this while loop prevents the current time from exceding steps, population exceding max taxa, and the population being below zero
  while (curr_time <= stop_time && sum(population) <= max_taxa && sum(population)!=0) {
    #find the prob of birth in state 1
    birth_prob1 <- population[1] * birth1
    #find the prob of birth in state 2
    birth_prob2 <- population[2] * birth2
    #find the prob of death in state 1
    death_prob1 <- population[1] * death1
    #find the prob of death in state 2
    death_prob2 <- population[2] * death2
    #Q matrix
    trans_matrix<- matrix(c(0, population[1]*q12, population[2]*q21, 0), nrow = 2, byrow = TRUE)
    #Make all the rate relative
    all_rates <- c(birth_prob1, birth_prob2, death_prob1, death_prob2, trans_matrix[1,2], trans_matrix[2,1])
    next_time <- rexp(1, rate=1/sum(all_rates))
    #Make it end at 100
    curr_time = curr_time + next_time
    
    
    relative <- cumsum(all_rates / sum(all_rates))
    
    #Draw an event
    draw <- runif(1)
    if (draw < relative[1]) {
      print("speciation state 1")
      population[1] = population[1] + 1
      event = "speciation state 1"
    } else if (draw < relative[2]) {
      print("speciation state 2")
      population[2] = population[2] + 1
      event = "speciation state 2"
    } else if (draw < relative[3]) {
      print("death state 1")
      population[1] = population[1] - 1
      event = "death state 1"
    } else if (draw < relative[4]) {
      print("death state 2")
      population[2] = population[2] - 1
      event = "death state 2"
    } else if (draw < relative[5]) {
      print("transition from state 1 to 2")
      population[1] = population[1] - 1
      population[2] = population[2] + 1
      event = "transition from state 1 to 2"
    } else if (draw < relative[6]) {
      print("transition from state 2 to 1")
      population[2] = population[2] - 1
      population[1] = population[1] + 1
      event = "transition from state 2 to 1"
    }
    
    
    taxa_table_row= c( population[1],
                       population[2],
                       curr_time,
                       event 
    )
    
    taxa_table = rbind(taxa_table, taxa_table_row)
    
    
    
  }
  return(taxa_table)
  
}

result <- sim_BiSSE(birth1, birth2, death1, death2, q12, q21, stop_time, curr_time, max_taxa)

#Trying to simulate a tree
# https://rdrr.io/cran/diversitree/man/simulate.html
library(diversitree)
library(geiger)

pars<- c(birth1, birth2, death1, death2, q12, q21)
set.seed(3)
phy <- tree.bisse(pars, max.taxa=4, x0=0, include.extinct = TRUE)
phy$tip.state
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)





#####################################
# Making this model more extensible #
#####################################
# SEAN                              #
#####################################

nstates = 2

birth_prior = rexp
death_prior = rexp



#If you have more births you can add them later
birthvec = rep(birth_prior(1), nstates)
deathvec = rep(death_prior(1), nstates)
populationvec = rep(x=0, nstates)
root_state = sample(x = 1:nstates, size = 1)

populationvec[sample(x = 1:nstates, size = 1)]

trans_matrix <- matrix(c(0, population[1]* q1_2, population[2]* q2_1, 0), nrow = 2, byrow = TRUE)






















