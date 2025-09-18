library(elo)
library(gtools)
library(heatsourcetools)
library(dplyr)
library(fitdistrplus)
library(beepr)
library(rstatix)
library(ggplot2)

#################################################################################
## Initalizing Main Variables
EloRange <- c(seq(0,2400,1))
x <- c(1/(1+(10^((EloRange-1200)/400))))# Xvalues for Prob
xs <- c(1/(1+(10^(((EloRange+1)-1200)/400))))
Teams <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T") # Teams
# <- c("A","B","C","D","E","F","G","H","I","J")
combos <- combos <- combn(c(1:length(Teams)), 2)
Order <- c(sample(c(1:length(combos[1,])), length(c(1:length(combos[1,]))), replace = FALSE))
N <- 20
TrueRatings <- rnorm(length(Teams), mean = 1200, sd = 200) # creates Referant elos
RatingResult <- list((pbeta(x,8,6)-pbeta(xs,8,6))) # Prior Elo Vals
Count <- 10# Prior Elo Weigth
TeamData <- as.data.frame(cbind(Teams,TrueRatings,Count,RatingResult))

diffs_mat <- outer(EloRange, EloRange, "-")
win_probs <- 1 / (1 + 10^(-diffs_mat/400))   

##################################################################################


#Generates distributions for combination
BetaUpdate <-function(N,S,C,VELO,HELO){
  alpha <- S+1 # Evidence Alpha
  beta <- 1+(N-S) #Evidence Beta
  HeroElo <- HELO # prior info about the team we wish to update (The Hero)
  prior <- VELO # prior info about the team the Hero is playing against (The Villan)
  posterior <- prior # just a filler to initialize the Elo range which Prior elo +Rating Differance will be placed into.
  data <- pbeta(x,alpha,beta)-pbeta(xs,alpha,beta) # Beta Dist obtained using alpha = 1+S, Beta = 1+(N-S),  probability of rating differance = (-400*(log((1/x)-1,10)))
  #Elo Shift calculation. Product of Sum (P(Villain Elo = X) * P(Villain Elo - Difference = Y)) for all (X-Y = Difference) = P(HeroElo = Y)
  
  W <- N/C

  for(i in 1:length(x)){ 

    db <- dbeta(win_probs[,i], alpha, beta)  # beta density evaluated at win_probs
    weighted_sum <- sum(data * prior[i] * db)
    posterior[i] <- weighted_sum
  }
  posterior<-posterior/sum(posterior) #standerdizze
  
  
  
  
  #Update Hero
  postAdj <- posterior
  HeroElo <-   (HeroElo+(W*postAdj))/sum((HeroElo+(W*postAdj))) # combine the Dists, if more rounds are contained in the Hero Elo then in the postAdj data they should be weighted acordingly.
  # ^ I.E If we have 20 rounds of info in the HeroElo table and only 10 in the postAdj table, then the formula would be  (((2/3)*(HeroElo[,2]))* ((1/3*postAdj[,2])))/sum(((2/3)*(HeroElo[,2]))* ((1/3*postAdj[,2])) 
  return(list(HeroElo))}

#Error Calculation
ErrorCalc <- function(K){
  print(RoundRecord$EloEstByRound[[K]] - TrueRatings)
  RMSE <- c()
  TrueRatingADJ <- TrueRatings + (1200 - mean(TrueRatings))
  EloEstsADJ <- RoundRecord$EloEstByRound[[K]] + (1200-mean(RoundRecord$EloEstByRound[[K]] ))
  TrueRatingsDif <- matrix(0,length(Teams),length(Teams))
  EloEstsDif <- matrix(0,length(Teams),length(Teams))
  for(i in 1:length(TrueRatingADJ)){
    for(j in 1:length(TrueRatingADJ)){
      TrueRatingsDif[i,j] <- TrueRatingADJ[i]-TrueRatingADJ[j]
      EloEstsDif[i,j] <- EloEstsADJ[i] - EloEstsADJ[j]
    }
  }
  
  EloEstsDif <- as.matrix(pull_lower_triangle(EloEstsDif, diagonal = TRUE))
  EloEstsDif <- as.double(as.vector(EloEstsDif)[which(as.vector(EloEstsDif) != "")])
  EloEstsDif <- EloEstsDif[which(EloEstsDif != 0)]
  
  TrueRatingsDif<- as.matrix(pull_lower_triangle(TrueRatingsDif, diagonal = TRUE))
  TrueRatingsDif <- as.double(as.vector(TrueRatingsDif)[which(as.vector(TrueRatingsDif) != "")])
  TrueRatingsDif <- TrueRatingsDif[which(TrueRatingsDif != 0)]
  
  ErrorMatrix <- EloEstsDif-TrueRatingsDif
  RMSE <- sqrt(mean(ErrorMatrix^2))
  MeanError <- mean(ErrorMatrix)
  print(paste("RMSE:", sqrt(mean(ErrorMatrix^2)), "ME:", MeanError))
  return(RMSE)
}
##################################################################################
TrueProb <- matrix(0,length(Teams),length(Teams))
IdealOutcomes <- matrix(0,length(Teams),length(Teams))
for(i in 1:length(Teams)){
  for(j in 1:length(Teams)){
    TrueProb[i,j] <- 1/(1+(10^((TrueRatings[j]-TrueRatings[i])/400))) # Calculates Win Rates
  }
}
for(i in 1:length(Teams)){
  for(j in 1:length(Teams)){
    #IdealOutcomes[i,j] <- round(N*(1/(1+(10^((TrueRatings[j]-TrueRatings[i])/400)))),0) # Calculates Ideal Outcomes
    IdealOutcomes[i,j] <- rbinom(1,N,(1/(1+(10^((TrueRatings[j]-TrueRatings[i])/400)))))
    IdealOutcomes[j,i] <- N-IdealOutcomes[i,j]
  }
}

print(TrueRatings)

###############################################################################
## Recording tool
RoundRecord <- data.frame(Round = c(1:5), EloEstByRound = c(1:5), ErrByRound = c(1:5), RMSE = c(1:5))
for( i in 1:length(RoundRecord$EloEstByRound)){
  RoundRecord$EloEstByRound[i] <- list(c(1:length(TrueRatings)))
  RoundRecord$ErrByRound[i] <- list(c(1:length(TrueRatings)))
}
ErrorRecord <- data.frame(Teams,Trial = 1 ,ErrorbyRound = c(1:length(Teams)))
ErrorRecord$ErrorbyRound[1:length(Teams)] <- list(c(1:length(Teams)))

for(k in 1){
  RunTime <- ((length(Teams)-1)*(length(Teams)))/2
  Trial <- 0
  ## Main Calc Loop
  
  for( i in Order){
    Start <- Sys.time()
    Home <- combos[1,i]
    Away <- combos[2,i]
    # Pass Info into BetaFunction
    HOLD <- BetaUpdate(20,IdealOutcomes[Away,Home],TeamData$Count[[Home]],TeamData$RatingResult[[Away]],TeamData$RatingResult[[Home]])
    HOLDAwy <- BetaUpdate(20,IdealOutcomes[Home,Away],TeamData$Count[[Away]],TeamData$RatingResult[[Home]],TeamData$RatingResult[[Away]])
    
    TeamData$RatingResult[[Home]] <- HOLD[[1]]
    TeamData$RatingResult[[Away]] <- HOLDAwy[[1]]
    TeamData$Count[[Home]] <- TeamData$Count[[Home]] +20
    TeamData$Count[[Away]] <-  TeamData$Count[[Away]] +20
    
    Trial <- Trial +1
    print(paste("Hero:",Teams[Home],"Villan:", Teams[Away], "Progress:", round((Trial/RunTime)*100,2), "%"))
    print(paste("True:",TrueRatings[Home] ,"Est:",sum(TeamData$RatingResult[[Home]]*EloRange)  ,"error:",sum(TeamData$RatingResult[[Home]]*EloRange) - TeamData$TrueRatings[[Home]]))
    
    plot(EloRange,TeamData$RatingResult[[Home]], type = 'l', main = paste("Hero:",Teams[Home],"Villan:",Teams[Away]))
    abline(v= sum(EloRange*TeamData$RatingResult[[Home]]), col = 'red' )
    abline(v = TrueRatings[Home])
    
    
    
    
    ErrorRecord$ErrorbyRound[[Home]][ErrorRecord$Trial[Home]] <- sum(TeamData$RatingResult[[Home]]*EloRange) -TeamData$TrueRatings[[Home]]
    ErrorRecord$Trial[Home] <- ErrorRecord$Trial[[Home]] + 1
    ErrorRecord$ErrorbyRound[[Away]][ErrorRecord$Trial[Away]] <- sum(TeamData$RatingResult[[Away]]*EloRange) -TeamData$TrueRatings[[Away]]
    ErrorRecord$Trial[Away] <- ErrorRecord$Trial[[Away]] + 1
    
    
    print(Sys.time()-Start)
  }
  
### Record Stages by Itteration
  RoundRecStore <- c()
  ErrRecStore <- c()
  for(i in 1:length(RoundRecord$EloEstByRound[[1]])){
    print(i)
    RoundRecStore[i] <- sum(EloRange*TeamData$RatingResult[[i]])
    ErrRecStore[i] <- sum(EloRange*TeamData$RatingResult[[i]]) - TrueRatings[i]
  }
  RoundRecord$EloEstByRound[[k]] <-RoundRecStore
  RoundRecord$ErrByRound[[k]] <- ErrRecStore - mean(ErrRecStore)
  RoundRecord$RMSE[k] <- sqrt(mean((ErrRecStore - mean(ErrRecStore))^2))
  Order <- c(sample(c(1:length(combos[1,])), length(c(1:length(combos[1,]))), replace = FALSE))
  
  ErrorCalc(k)
  
  beep(sound = 8)
}






##### Error Track
plot(c(1:(length(Teams)-1)),ErrorRecord$ErrorbyRound[[1]][1:19],type = 'l', ylim = c(min(unlist(as.vector(ErrorRecord$ErrorbyRound))),max(unlist(as.vector(ErrorRecord$ErrorbyRound)))))

for(i in 1:length(ErrorRecord$Teams)){
  lines(c(1:(length(Teams)-1)),ErrorRecord$ErrorbyRound[[i]][1:19],type = 'l')
}
abline(h = 0, col = "green")

RMSERound <- c()
RoundErr <- c()
for( j in 1:(length(ErrorRecord$ErrorbyRound[[1]])-1)){
  ErrJ <- c()
  for(i in 1:(length(ErrorRecord$Teams))){
    ErrJ[i] <- ErrorRecord$ErrorbyRound[[i]][j] 
  }
  print(sd(ErrJ))
  RMSERound[j] <-sqrt(mean((ErrJ -mean(ErrJ))^2))
  RoundErr[j] <- mean(ErrJ)
}

RMSEDat <- data.frame("Round" = c(1:(length(ErrorRecord$ErrorbyRound[[1]])-1)),"MeanErr" = RoundErr, "RoundVar" = RMSERound)


p1 <- ggplot(RMSEDat, aes(Round)) + geom_line(aes(y = c(sample(c(0),(length(ErrorRecord$ErrorbyRound[[1]])-1),replace = TRUE))), color = 'black', linewidth = 2)
p2<- p1 + geom_ribbon(aes(ymin = 1/(1+(10^((RoundVar)/400)))-0.5, ymax = 0.5-(1/(1+(10^((RoundVar)/400))))), color = "red", alpha = 0.2, fill = 'red') 
p3 <- p2 + geom_ribbon(aes(ymin =(1/(1+(10^((2*RoundVar)/400))))-0.5 , ymax = 0.5-(1/(1+(10^((2*RoundVar)/400))))), color = "red", alpha = 0.4, fill = 'red')
print(p3)
print((1/(1+10^(RMSEDat$RoundVar[(length(ErrorRecord$ErrorbyRound[[1]])-1)]/400)))-0.5)
