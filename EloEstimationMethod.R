library(elo)
library(gtools)
library(heatsourcetools)
library(dplyr)
library(fitdistrplus)
library(beepr)
library(rstatix)
###ErrorCountRelation
Teams <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T")
TrueRatings <- rnorm(length(Teams), mean = 1200, sd = 200) # creates Referant elos
EloRange <- c(seq(500,1900,1))
x <- c(1/(1+(10^((EloRange-1200)/400))))# Xvalues for Prob
xm1 <- c(1/(1+(10^(((EloRange+1)-1200)/400))))
 # Teams
#Teams <- c("A","B","C","D","E","F","G","H","I","J")
combos <- combos <- combn(c(1:length(Teams)), 2)
Order <- c(sample(c(1:length(combos[1,])), length(c(1:length(combos[1,]))), replace = FALSE))

#################################################################################
## Initalizing Main Variables
N <- 20
RatingResult <- list(cbind((-400*(log((1/x)-1,10)))+1200,(pbeta(x,8,6)-pbeta(xm1,8,6)))) # Prior Elo Vals
Count <- 7# Prior Elo Weigth
TeamData <- as.data.frame(cbind(Teams,TrueRatings,Count,RatingResult))
##################################################################################


#Generates distributions for combination
BetaUpdate <-function(N,S,C,CO,VELO,HELO){
  alpha <- S+1 # Evidence Alpha
  beta <- 1+(N-S) #Evidence Beta
  HeroElo <- HELO # prior info about the team we wish to update (The Hero)
  prior <- VELO # prior info about the team the Hero is playing against (The Villan)
  posterior <- prior # just a filler to initialize the Elo range which Prior elo +Rating Differance will be placed into.
  data <- data.frame("1"=(-400*(log((1/x)-1,10))),"X2" =(pbeta(x,alpha,beta)-pbeta(xm1,alpha,beta))) # Beta Dist obtained using alpha = 1+S, Beta = 1+(N-S),  probability of rating differance = (-400*(log((1/x)-1,10)))
  names(prior) = c("X1","X2")
  names(HeroElo) = c("X1","X2")
  names(posterior) = c("X1","X2")
  names(data) = c("X1","X2")
  

  
  
  #Elo Shift calculation. Product of Sum (P(Villain Elo = X) * P(Villain Elo - Difference = Y)) for all (X-Y = Difference) = P(HeroElo = Y)
  for(i in 1:length(x)){ 
    diffs <- prior[i,1] - prior[,1]  # vector of differences
    win_probs <- 1 / (1 + 10^(-diffs / 400))  # logistic win probabilities
    db <- dbeta(win_probs, alpha, beta)  # beta density evaluated at win_probs
    weighted_sum <- sum(data[,2] *prior[i,2] * db)
    posterior[i,2] <- weighted_sum
  }
  posterior[,2]<-posterior[,2]/sum(posterior[,2]) #standerdizze
  
  weight <- C/N
  weight2 <-1
  
  #Update Hero
  StanderdX <- data.frame(X1 = HeroElo[,1], X2 = sample(c(0), length(HeroElo[,1]), replace = TRUE)) #Reintroduce Elos of Intrest
  postAdj <- subset(near_join(StanderdX, posterior,by = "X1", tolerance = 20), X2!=0) # move posterior Data onto Elo Range
  HeroElo[,2] <- (((weight2*HeroElo[,2])+(weight*postAdj[,2])))/sum((weight2*HeroElo[,2])+(weight*postAdj[,2])) # combine the Dists, if more rounds are contained in the Hero Elo then in the postAdj data they should be weighted acordingly.
  # ^ I.E If we have 20 rounds of info in the HeroElo table and only 10 in the postAdj table, then the formula would be  (((2/3)*(HeroElo[,2]))* ((1/3*postAdj[,2])))/sum(((2/3)*(HeroElo[,2]))* ((1/3*postAdj[,2])) 
  return(HeroElo)}

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

#print(TrueRatings)

###############################################################################
## Recording tool
RoundRecord <- data.frame(Round = c(1:5), EloEstByRound = c(1:5))
for( i in 1:length(RoundRecord$EloEstByRound)){
  RoundRecord$EloEstByRound[i] <- list(c(1:length(TrueRatings)))
}
ErrorRecord <- data.frame(Teams,Trial = 1 ,ErrorbyRound = c(1:length(Teams)))
ErrorRecord$ErrorbyRound[1:length(Teams)] <- list(c(1:length(Teams)))


RunTime <- ((length(Teams)-1)*(length(Teams)))/2
Trial <- 0
## Main Calc Loop

  for( i in 1:length(Order)){
    Start <- Sys.time()
 
    
    # Pass Info into BetaFunction
    Hold <- BetaUpdate(20,IdealOutcomes[combos[2,Order[i]],combos[1,Order[i]]],TeamData$Count[[combos[1,Order[i]]]][1],TeamData$Count[[combos[2,Order[i]]]][1],as.data.frame(TeamData$RatingResult[[combos[1,Order[i]]]]),as.data.frame(TeamData$RatingResult[[combos[2,Order[i]]]])) #Update J
    TeamData$RatingResult[[combos[1,Order[i]]]] <- BetaUpdate(20,IdealOutcomes[combos[1,Order[i]],combos[2,Order[i]]],TeamData$Count[[combos[2,Order[i]]]][1],TeamData$Count[[combos[1,Order[i]]]][1],as.data.frame(TeamData$RatingResult[[combos[2,Order[i]]]]),as.data.frame(TeamData$RatingResult[[combos[1,Order[i]]]])) #Update i
    TeamData$Count[[combos[1,Order[i]]]][1] <- TeamData$Count[[combos[1,Order[i]]]][1] + 20 #Update Count i
    TeamData$Count[[combos[2,Order[i]]]][1] <- TeamData$Count[[combos[2,Order[i]]]][1] + 20 #update Count j
    TeamData$RatingResult[[combos[2,Order[i]]]] <- Hold # pass hold so it dosent influance i calc
    Trial <- Trial +1
    print(paste("Hero:",combos[1,Order[i]],"Villan:",combos[2,Order[i]], "Progress:", round((Trial/RunTime)*100,2), "%"))
    #print(paste("True:",TrueRatings[combos[1,Order[i]]] ,"Est:",sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2]) ,"error:",sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2])-TrueRatings[combos[1,Order[i]]]))
    
 
    ErrorRecord$ErrorbyRound[[combos[1,Order[i]]]][ErrorRecord$Trial[combos[1,Order[i]]]] <- sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2])
    ErrorRecord$Trial[combos[1,Order[i]]] <-  ErrorRecord$Trial[combos[1,Order[i]]] +1
    
    ErrorRecord$ErrorbyRound[[combos[2,Order[i]]]][ErrorRecord$Trial[combos[2,Order[i]]]] <- sum(TeamData$RatingResult[[combos[2,Order[i]]]][,1]*TeamData$RatingResult[[combos[2,Order[i]]]][,2])
    ErrorRecord$Trial[combos[2,Order[i]]] <-  ErrorRecord$Trial[combos[2,Order[i]]] +1
    
    #print(Sys.time()-Start)
    }

  ### Record Stages by Itteration

  
# Order <- c(sample(c(1:length(combos[1,])), length(c(1:length(combos[1,]))), replace = FALSE))


beep(sound = 8)


ErrorCalc <- function(RR){
  RMSE <- c()
  RMSERec <- c()
  TrueRatingADJ <- TrueRatings + (1200 - mean(TrueRatings)) #Normalized
  EloEstsADJ <- RR + (1200-mean(RR)) #Normalized
  #Generates Rating Differance Matrix
  TrueRatingsDif <- matrix(0,length(Teams),length(Teams))
  EloEstsDif <- matrix(0,length(Teams),length(Teams))
  #Inputes into rating Differance Matrix
  for(i in 1:length(TrueRatingADJ)){
    for(j in 1:length(TrueRatingADJ)){
      TrueRatingsDif[i,j] <- TrueRatingADJ[i]-TrueRatingADJ[j]
      EloEstsDif[i,j] <- EloEstsADJ[i] - EloEstsADJ[j]
    }
  }
  
  EloEstsDif <- as.matrix(pull_lower_triangle(EloEstsDif, diagonal = FALSE))
  
  #EloEstsDif <- as.double(as.vector(EloEstsDif)[which(as.vector(EloEstsDif) != "")])
  #EloEstsDif <- EloEstsDif[which(EloEstsDif != 0)]
  
  TrueRatingsDif<- as.matrix(pull_lower_triangle(TrueRatingsDif, diagonal = FALSE))
  
  #TrueRatingsDif <- as.double(as.vector(TrueRatingsDif)[which(as.vector(TrueRatingsDif) != "")])
  #TrueRatingsDif <- TrueRatingsDif[which(TrueRatingsDif != 0)]
  
  ErrorMatrix <- as.numeric(EloEstsDif[EloEstsDif != ""]) -as.numeric(TrueRatingsDif[TrueRatingsDif != ""])
  RMSE <- sqrt(mean(ErrorMatrix^2))
  MeanError <- mean(ErrorMatrix)
  print(RMSE)
  return(list(RMSE,ErrorMatrix))
}

##### Error Track
TrackRMSE <- c()
Variance <- c()
TrueElos <- c(c())
EstElos <- c(c())
for(i in 1:(length(ErrorRecord$Teams)-1)){
  EC <- ErrorCalc(sapply(ErrorRecord$ErrorbyRound, function(x){x[i]}))
  TrackRMSE[i] <- EC[[1]][1]
  Variance[i] <- sd(EC[[2]])

}


plot(1:19,TrackRMSE,type = 'l')
ErrorRedux <- lm(TrackRMSE~c(1:19))
lines(1:19,ErrorRedux$fitted.values,type = 'l', col='red')

Accuracy <- lm(EloEstsADJ~TrueRatingADJ)
plot(TrueRatingADJ,EloEstsADJ, xlab = "True Rating", ylab = "Estimated Rating", main = "Real vs Predicted Elos")
lines(TrueRatingADJ,Accuracy$fitted.values, col = "red")