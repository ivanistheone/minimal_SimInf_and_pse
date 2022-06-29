#SIR model - simple for Ivan + wrapper functions for use with pse

library(SimInf)

## Define model - SIR model
transitions <- c("S -> S*beta*I/(S+I+R) -> I", 
                 "I -> gamma*I -> R")

compartments <- c("S", "I", "R") #number of model compartments, equivalent to Susceptible, Infected, Removed
n <- 1000 #number of simulations of the model
Psize <- 10000 #Population size
Isize <- 5 #Number of people initially infected


# Function that runs the SimInf model with parameters:
# - beta = transmission coefficient
# - gamma = infectious period (4 days^-1)
oneRun <- function (beta, gamma) {
  # initial starting conditions
  u0 <- data.frame(
    S=rep((Psize-Isize),n),
    I=rep(Isize,n),
    R=rep(0,n)
  )
  model <- mparse(transitions = transitions,
                  compartments = compartments,
                  gdata = c(beta = beta, gamma=gamma),
                  u0 = u0,
                  tspan = 1:400) #4oo days with each time step one day
  result <- run(model)
  return (result)
}


# Test run
result = oneRun(beta = 0.6, gamma = 0.25)


# Plot model result
plot.new()
par(mfrow=(c(1,1)))
plot(result)






#Now want to get the model results to talk with the pse package.. something like below 

library(pse)

LHS1 <- LHS(result, N=60, factors=2, repetitions=5, nboot=50)
