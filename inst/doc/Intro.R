## -----------------------------------------------------------------------------
# Use HPL() to fit the 20-21 PL soccer data from url:https://www.football-data.co.uk/
library(SA23204181)
data(data)
attach(data)
# The data should be processed reserving 
# "HomeTeam,AwayTeam,FTHG,FTAG,FTR,MaxH,MaxD,MaxA,AvgH,AvgD,AvgA" columns of raw from url.
# size_train A number of training set.
# size_valid A number of prediction or validation set.
PL2021<- HPL(data,size_train,size_valid)
PL2021
# return shows prediction criteria: accurary,rps,likelihood-based criteria,profit.

