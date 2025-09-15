library(elo)
library(gtools)
library(heatsourcetools)
library(dplyr)
library(fitdistrplus)
library(beepr)
library(rstatix)


set.seed(190022838)
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
RatingResult <- list(cbind((-400*(log((1/x)-1,10)))+1200,(pbeta(x,1,1)-pbeta(xs,1,1)))) # Prior Elo Vals
Count <- 8# Prior Elo Weigth
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
  data <- data.frame("1"=(-400*(log((1/x)-1,10))),"X2" =(pbeta(x,alpha,beta)-pbeta(xs,alpha,beta))) # Beta Dist obtained using alpha = 1+S, Beta = 1+(N-S),  probability of rating differance = (-400*(log((1/x)-1,10)))
  #Elo Shift calculation. Product of Sum (P(Villain Elo = X) * P(Villain Elo - Difference = Y)) for all (X-Y = Difference) = P(HeroElo = Y)
        

  B <- dbeta(win_probs, alpha, beta)    
  
  data_mat <- matrix(data[,2],length(VELO[,1]),length(VELO[,1]), byrow = TRUE)
  
  weighted <- (B) * t(data_mat)
  

  posterior[,2] <- rowSums(weighted)*(VELO[,2])
  posterior[,2]<-posterior[,2]/sum(posterior[,2]) #standerdizze
  
  
  #Update Hero
  postAdj <- posterior
  HeroElo[,2] <-   (((C/(C))*HeroElo[,2])+(((N/(C))*postAdj[,2])))/sum((((C/(C))*HeroElo[,2])+(((N/(C))*postAdj[,2])))) # combine the Dists, if more rounds are contained in the Hero Elo then in the postAdj data they should be weighted acordingly.
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

print(TrueRatings)

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
  Hold <- BetaUpdate(20,IdealOutcomes[combos[2,Order[i]],combos[1,Order[i]]],TeamData$Count[[combos[2,Order[i]]]][1],as.data.frame(TeamData$RatingResult[[combos[1,Order[i]]]]),as.data.frame(TeamData$RatingResult[[combos[2,Order[i]]]])) #Update J
  TeamData$RatingResult[[combos[1,Order[i]]]] <- BetaUpdate(20,IdealOutcomes[combos[1,Order[i]],combos[2,Order[i]]],TeamData$Count[[combos[1,Order[i]]]][1],as.data.frame(TeamData$RatingResult[[combos[2,Order[i]]]]),as.data.frame(TeamData$RatingResult[[combos[1,Order[i]]]])) #Update i
  TeamData$Count[[combos[1,Order[i]]]][1] <- TeamData$Count[[combos[1,Order[i]]]][1] + 20 #Update Count i
  TeamData$Count[[combos[2,Order[i]]]][1] <- TeamData$Count[[combos[2,Order[i]]]][1] + 20 #update Count j
  TeamData$RatingResult[[combos[2,Order[i]]]] <- Hold # pass hold so it dosent influance i calc
  ##Ploting To Track
  print(sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2]))
  plot(TeamData$RatingResult[[combos[1,Order[i]]]][,1],TeamData$RatingResult[[combos[1,Order[i]]]][,2], type = 'l', main = paste("Hero:",Teams[combos[1,Order[i]]],"Villan:",Teams[combos[2,Order[i]]]))
  abline(v= sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2]), col = 'red' )
  abline(v = TrueRatings[combos[1,Order[i]]])
  Trial <- Trial +1
  print(paste("Hero:",combos[1,Order[i]],"Villan:",combos[2,Order[i]], "Progress:", round((Trial/RunTime)*100,2), "%"))
  print(paste("True:",TrueRatings[combos[1,Order[i]]] ,"Est:",sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2]) ,"error:",sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2])-TrueRatings[combos[1,Order[i]]]))
  
  
  ErrorRecord$ErrorbyRound[[combos[1,Order[i]]]][ErrorRecord$Trial[combos[1,Order[i]]]] <- sum(TeamData$RatingResult[[combos[1,Order[i]]]][,1]*TeamData$RatingResult[[combos[1,Order[i]]]][,2])-TrueRatings[combos[1,Order[i]]]
  ErrorRecord$Trial[combos[1,Order[i]]] <-  ErrorRecord$Trial[combos[1,Order[i]]] +1
  
  ErrorRecord$ErrorbyRound[[combos[2,Order[i]]]][ErrorRecord$Trial[combos[2,Order[i]]]] <- sum(TeamData$RatingResult[[combos[2,Order[i]]]][,1]*TeamData$RatingResult[[combos[2,Order[i]]]][,2])-TrueRatings[combos[2,Order[i]]]
  ErrorRecord$Trial[combos[2,Order[i]]] <-  ErrorRecord$Trial[combos[2,Order[i]]] +1
  
  print(Sys.time()-Start)
}

### Record Stages by Itteration
RoundRecStore <- c()
for(i in 1:length(RoundRecord$EloEstByRound[[1]])){
  print(i)
  RoundRecStore[i] <- sum(TeamData$RatingResult[[i]][,1]*TeamData$RatingResult[[i]][,2])
}
RoundRecord$EloEstByRound[[1]] <-RoundRecStore

Order <- c(sample(c(1:length(combos[1,])), length(c(1:length(combos[1,]))), replace = FALSE))


beep(sound = 8)


print(RoundRecord$EloEstByRound[[1]] - TrueRatings)
RMSE <- c()
TrueRatingADJ <- TrueRatings + (1200 - mean(TrueRatings))
EloEstsADJ <- RoundRecord$EloEstByRound[[1]] + (1200-mean(RoundRecord$EloEstByRound[[1]] ))
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
plot(c(1000:1400),dnorm(c(1000:1400),mean=MeanError+1200, sd = RMSE))


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

RMSEDat <- data.frame("Round" = c(1:19),"MeanErr" = RoundErr, "RoundVar" = RMSERound)


p1 <- ggplot(RMSEDat, aes(Round)) + geom_line(aes(y = c(sample(c(0),19,replace = TRUE))), color = 'black', linewidth = 2)
p2<- p1 + geom_ribbon(aes(ymin = 1/(1+(10^((RoundVar)/400)))-0.5, ymax = 0.5-(1/(1+(10^((RoundVar)/400))))), color = "red", alpha = 0.2, fill = 'red') 
p3 <- p2 + geom_ribbon(aes(ymin =(1/(1+(10^((2*RoundVar)/400))))-0.5 , ymax = 0.5-(1/(1+(10^((2*RoundVar)/400))))), color = "red", alpha = 0.4, fill = 'red')
print(p3)
print((1/(1+10^(RMSEDat$RoundVar[19]/400)))-0.5)