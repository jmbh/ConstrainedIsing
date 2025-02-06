# jonashaslbeck@protonmail.com; Feb 6th, 20256


source("IsingNegLogLik.R")
source("OptimConstrIsing.R")


## Generate Data
library(IsingSampler)
p <- 3
G <- matrix(0, p, p)
G[1,2] <- G[2,1] <- 1
th <- c(0.5, 0, 0)

set.seed(1)
X <- IsingSampler(n=1000, 
                  graph=G, 
                  thresholds = th, 
                  responses=c(-1, 1))

## Recover without constraints
estimate_ising_mle(X) # Looks quite good

## With constraints
J_fixed <- matrix(NA, p, p) # also possible, here no constraints
h_fixed <- c(1,1,1) # constrain all thresholds to 1
estimate_ising_mle(X, h_fixed = h_fixed, J_fixed = J_fixed) 

# fixing thresholds affects interaction-estimates







