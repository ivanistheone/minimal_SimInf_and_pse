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



# The model that you wish to analyse must be formulated as an R function that:
# receives a data.frame, in which every column represent a different parameter,
# and every line represents a different combination of values for those parameters.
# The function must return an array with the same number of elements as there were lines
# in the original data frame, and each entry in the array should correspond to the result
# of running the model with the corresponding parameter combination.

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
  
  # relevant info:
  #    result.U = the result matrix with the number of individuals in each compartment in every node.
  #    Integer matrix (NnodesNcomp × length(tspan)).
  #    U[, j] contains the number of individuals in each compartment at tspan[j].
  #    U[1:Ncomp, j] contains the number of individuals in each compartment in node 1 at tspan[j].
  #    U[(Ncomp + 1):(2 * Ncomp), j] contains the number of individuals in each compartment in node 2 at tspan[j] etc.

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
  #       so to avoid this error we don't return the IfinalVal
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



# Now want to get the model results to talk with the pse package.. something like below 

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






# LARGE RUN WITH N=100
# I keep repetitions=1 still, because model 
LHS100 <- LHS(modelRun,
             factors=factors,
             N=100,
             q=q,
             q.arg=q.arg,
             res.names=c("S", "R"),
             nboot=50,
             repetitions=1)

# Looking around at parameters we used for each simulation
class(get.data(LHS100))
dim(get.data(LHS100))
# get.data(LHS100)

# ... and results
class(get.results(LHS100))
dim(get.results(LHS100))
# get.results(LHS100)

# try plots from the vignette pse_tutorial.pdf
# https://mran.microsoft.com/snapshot/2017-05-24/web/packages/pse/vignettes/pse_tutorial.pdf
plotecdf(LHS100, stack=TRUE)
plotscatter(LHS100, add.lm=FALSE)
plotprcc(LHS100)

