#SIR model - simple for Ivan + wrapper functions for use with pse

library(SimInf)

## Define model - SIR model
transitions <- c("S -> S*beta*I/(S+I+R) -> I", 
                 "I -> gamma*I -> R")

compartments <- c("S", "I", "R")
n <- 1000      # number of simulations of the SimInf model
Psize <- 10000 # Population size
Isize <- 5     # Number of people initially infected
tmax <- 400    # 4oo days with each time step one day


# Testing code to make sure basic wrapping function works...
# # Function that runs the SimInf model with parameters:
# # - beta = transmission coefficient
# # - gamma = infectious period (4 days^-1)
# # and returns the SimInf result (an object)
# oneRun <- function (beta, gamma) {
#   # initial starting conditions
#   u0 <- data.frame(
#     S=rep((Psize-Isize),n),
#     I=rep(Isize,n),
#     R=rep(0,n)
#   )
#   model <- mparse(transitions = transitions,
#                   compartments = compartments,
#                   gdata = c(beta = beta, gamma=gamma),
#                   u0 = u0,
#                   tspan = 1:tmax) 
#   result <- run(model)
#   return (result)
# }
# 
# 
# # Test run with parameters
# result = oneRun(beta = 0.6, gamma = 0.25)
# 
# 
# # Plot model result
# plot.new()
# par(mfrow=(c(1,1)))
# plot(result)


# Function that runs the SimInf model with parameters:
# - beta = transmission coefficient
# - gamma = infectious period (4 days^-1)
# and returns an array of the final values S and R at the end of tspan
oneRunFinalValues <- function (beta, gamma) {
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
                  tspan = 1:tmax) #4oo days with each time step one day
  result <- run(model)

  # NEW NEW NEW
  # Custom code to extract the final values of S, I, R at tmax
  
  # result.U = the result matrix with the number of individuals in each compartment in every node.
  # Integer matrix (NnodesNcomp × length(tspan)).
  # U[, j] contains the number of individuals in each compartment at tspan[j].
  # U[1:Ncomp, j] contains the number of individuals in each compartment in node 1 at tspan[j].
  # U[(Ncomp + 1):(2 * Ncomp), j] contains the number of individuals in each compartment in node 2 at tspan[j] etc.

  # Sample code that returns a single value = S at t=tmax (used for debugging)
  # SfinalVals <- result@U[seq(1, 3*n, by=3),tmax]
  # SfinalVal = mean(SfinalVals)
  # output = SfinalVal
  # return (output)

  # Extract the values from the result matrix U
  # result@U[,tmax] is a 3*n long array that contains 
  # [S,I,R, S,I,R, S,I,R, ...] values from the n runs
  # we use seq(1, 3*n, by=3) to get only the S values:
  SfinalVals <- result@U[seq(1, 3*n, by=3),tmax]
  SfinalVal  <- mean(SfinalVals)
  IfinalVals <- result@U[seq(2, 3*n, by=3),tmax]
  IfinalVal  <- mean(IfinalVals)
  RfinalVals <- result@U[seq(3, 3*n, by=3),tmax]
  RfinalVal  <- mean(RfinalVals)
  
  # Prepare array of returned values
  outputs <- array();
  outputs[1] <- SfinalVal
  outputs[2] <- RfinalVal
  # note: trying to include IfinalVal doesn't work out because it has zero variance
  #       and we get an error 
  #               Error in if (const(t, min(1e-08, mean(t, na.rm = TRUE)/1e+06))) { : 
  #               missing value where TRUE/FALSE needed
  #       from somewhere deep in the package.
  return (outputs)
}




# If your model accepts a single combination of values,
# it is easy to write a “wrapper” using the function mapply to your model.
# The function oneRun receives two parameters, corresponding to beta and gamma
# and returns a result single value corresponding to the final population.
# The function modelRun encapsulates this function, in a manner to receive a data.frame
# containing all parameter combinations and returning the results in one array.

modelRun <- function (my.data) {
  return (mapply(oneRunFinalValues, my.data[,1], my.data[,2]))
}





# Testing modelRun works 
# L <- as.data.frame(matrix(nrow=3, ncol=length(factors)))
# L[1,] <- list(0.6, 0.25)
# L[2,] <- list(0.62, 0.25)
# L[3,] <- list(0.63, 0.25)
# L
# 
# results <- modelRun(L)
# results
# 
# # results is of type list
# class(results)
# 
# # individual result
# res = results[[1]]







#Now want to get the model results to talk with the pse package.. something like below 

library(pse)


# SENSITIVITY ANALYSIS CHOICE
# How should we vary the parameters beta and gamma?
# Let's choose a uniform distribution of width 0.2 centered
# at the values beta <- 0.6 and gamma <- 0.25 which we used above
factors <-  c("beta",                       "gamma")
q <-        c("qunif",                      "qunif")
#             transmission coefficient      infectious period
#             beta <- 0.6                   gamma <- 0.25
q.arg <- list(list(min=0.5, max=0.7),       list(min=0.15, max=0.35) )



# SMALL TEST RUN WITH N=10 (to see what's going on)
LHS10 <- LHS(modelRun,
             factors=factors,
             N=10,
             q=q,
             q.arg=q.arg,
             res.names=c("S", "R"),
             nboot=50,
             repetitions=1)

# Looking around at parameters we used for each simulation
class(get.data(LHS10))
dim(get.data(LHS10))
get.data(LHS10)

# ... and results
class(get.results(LHS10))
dim(get.results(LHS10))
get.results(LHS10)

# try plots from the vignette pse_tutorial.pdf
# https://mran.microsoft.com/snapshot/2017-05-24/web/packages/pse/vignettes/pse_tutorial.pdf
plotecdf(LHS10, stack=TRUE)
plotscatter(LHS10, add.lm=FALSE)
plotprcc(LHS10)






# LARGE RUN WITH N=1000
# I keep repetitions=1 still, because 
LHS1K <- LHS(modelRun,
             factors=factors,
             N=1000,
             q=q,
             q.arg=q.arg,
             res.names=c("S", "R"),
             nboot=50,
             repetitions=1)

# Looking around at parameters we used for each simulation
class(get.data(LHS1K))
dim(get.data(LHS1K))
# get.data(LHS1K)

# ... and results
class(get.results(LHS1K))
dim(get.results(LHS1K))
# get.results(LHS1K)

# try plots from the vignette pse_tutorial.pdf
# https://mran.microsoft.com/snapshot/2017-05-24/web/packages/pse/vignettes/pse_tutorial.pdf
plotecdf(LHS1K, stack=TRUE)
plotscatter(LHS1K, add.lm=FALSE)
plotprcc(LHS1K)


