Frontier_Combinations <- function(totals){
  # determine frontier
  CEmat <- totals %>% dplyr::select(!c(Combination, model))
  CEmat <- rbind(c(0, 0, 0, 0), CEmat)
  
  # initialize some variables
  costsCol <- which(colnames(CEmat)=="Costs")
  qalyCol  <- which(colnames(CEmat)=="QALYs")
  icerCol  <- which(colnames(CEmat)=="ICERs")
  
  # 1st screening: remove strategies with QALYs lower than the least costly strategy
  lowestcost  <- as.numeric(which.min(CEmat[ , costsCol]))
  CEmat       <- CEmat[(CEmat[ ,qalyCol] >= CEmat[lowestcost,qalyCol]), ]
  
  # 2nd screening: remove strategies with costs higher than the strategy with highest QALYs
  highestqaly <- as.numeric(which.max(CEmat[ , qalyCol]))
  CEmat       <- CEmat[(CEmat[ , costsCol] <= CEmat[highestqaly, costsCol]), ]
  
  # # 3rd screening: remove strategies with -ve ICERS
  # dominated <- CEmat[CEmat[icerCol] <0,]$number
  # CEmat     <- CEmat[!(CEmat$number %in% dominated),]
  
  # 4th screening: remove infinite ICERs
  CEmat <- CEmat[!is.infinite(CEmat$ICERs),]
  
  # Number of combinations remaining after screening
  numComb <- dim(CEmat)[1]
  
  # find WTP levels to test so that all strategies on frontier will be captured
  # this means testing on either side of all NMB intersections, which are just all the pairwise ICERs
  ICERmat <- matrix(1, numComb, numComb)
  for (i in 1:numComb) {
    indexStrat             <- matrix(1, numComb, 4)
    indexStrat[, costsCol] <- indexStrat[, costsCol] * CEmat[i, costsCol]
    indexStrat[, qalyCol]  <- indexStrat[, qalyCol] * CEmat[i, qalyCol]
    indexStrat[, icerCol]  <- indexStrat[, icerCol] * CEmat[i, icerCol]
    delCostQalys           <- CEmat - indexStrat
    ICERmat[, i]           <- delCostQalys[, icerCol]
  } 
  
  ICERvec <- numeric(numComb*numComb/2 - numComb)
  for (ij in 1:(numComb-1)){
    ICERvec[(((numComb-1) + (numComb-ij+1))*(ij-1)/2 + 1):(((numComb-1) + (numComb-ij))*ij/2)] <- ICERmat[ij, (ij+1):numComb]
  }
  
  intersections <- sort(unique(ICERvec))
  intersections <- intersections[is.finite(intersections)]
  maxWTP <- 100000000        # any positive value or Inf
  WTPtestPoints <- c(0, intersections[intersections >= 0], maxWTP)
  rm(delCostQalys); rm(ICERmat); rm(ICERvec); rm(intersections)
  
  # Find the strategy with the max NMB at each of the WTP test points
  indiciesOfMax <- numeric(length(WTPtestPoints))
  for (i in 1:length(WTPtestPoints)){
    romd <- which.max((WTPtestPoints[i] * CEmat[, qalyCol]) - CEmat[, costsCol])
    indiciesOfMax[i] <- CEmat[romd, 1]
  }
  frontier.v <- unique(indiciesOfMax)
  
  # Determine frontier
  CEthreshold <- 100000  #cost-effectiveness threshold
  
  if(frontier.v[1] == 0){frontier.v <- frontier.v[-1]}
  
  frontier.matrix  <- totals[frontier.v , c("QALYs", "Costs","ICERs")]
  frontier.matrix  <- cbind(frontier.v, frontier.matrix)
  colnames(frontier.matrix)[1] <- "combination"
  costsCol <-which(colnames(frontier.matrix)=="Costs")
  qalyCol <- which(colnames(frontier.matrix)=="QALYs")
  
  if (any(frontier.matrix$Costs < 0)){
    frontier.matrix$ICER_front[frontier.matrix$Costs < 0] <- "CS"
    NCS  <- sum(frontier.matrix$Costs < 0)  #number of cost savings
    if (NCS == nrow(frontier.matrix)){
      ocis       <- frontier.matrix$Strategy.ind[nrow(frontier.matrix)]
    } else{
      frontier.matrix$ICER_front[(NCS+1):nrow(frontier.matrix)] <-
        diff(frontier.matrix$Costs[NCS:nrow(frontier.matrix)]) / diff(frontier.matrix$QALYs[NCS:nrow(frontier.matrix)])
      icer <- as.numeric(frontier.matrix$ICER_front[(NCS+1):nrow(frontier.matrix)])  #ICER for non cost-saving strategies
      if (all(icer >= CEthreshold)){
        ocis       <- frontier.matrix$Strategy.ind[NCS]
      } else{
        ocis.ind   <- which.max(icer[icer < CEthreshold]) + NCS
        ocis       <- frontier.matrix$Strategy.ind[ocis.ind]
      }
    }
  } else{
    frontier.matrix$ICER_front[1] <- frontier.matrix$Costs[1] / frontier.matrix$QALYs[1]
    frontier.matrix$ICER_front[2:nrow(frontier.matrix)] <-
      diff(frontier.matrix$Costs) / diff(frontier.matrix$QALYs)
    ocis.ind   <- which.max(frontier.matrix$ICER_front[frontier.matrix$ICERs < CEthreshold])
    ocis       <- frontier.matrix$Strategy.ind[ocis.ind]
  }
  
  frontier.df <- as.data.frame(frontier.matrix)
  rm(indexStrat, highestqaly, i, ij, indiciesOfMax, lowestcost, maxWTP, 
     numComb, ocis, costsCol, qalyCol, romd, WTPtestPoints,
     frontier.matrix, CEmat)
  
  return(frontier.df)
}
