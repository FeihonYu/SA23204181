#' @title A hierarchical Poisson Linear model for Premier League.
#' @param data used data.
#' @param size_train A number of training set.
#' @param size_valid A number of prediction or validation set.
#' @return The summary of prediction and profit by gambling.
#' @import INLA
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' b<-HPL(data,size_train = 360,size_valid = 20)
#' }
#' @export
HPL<- function(data,size_train,size_valid){
 
  newdata= tail(data,n=size_valid)
  dat <- head(data,size_train)
  FTHG = c(data$FTHG[1:size_train])
  FTAG = c(data$FTAG[1:size_train])
  Att_HomeTeam =dat$HomeTeam
  Def_AwayTeam =dat$AwayTeam
  HomeTeam = dat$HomeTeam

  goal<- c(FTHG,FTAG)
  data1 <- data.frame(goal=goal,
                    Att_HomeTeam = factor(rep(Att_HomeTeam,2),levels =c("0", unique(data$HomeTeam))),
                    Def_AwayTeam = factor(rep(Def_AwayTeam,2),levels =c("0", unique(data$HomeTeam))),
                    athome = factor(c(rep(1,size_train),rep(0,size_train)) )
                    )

  formula1 <- goal ~ Att_HomeTeam + Def_AwayTeam + athome   
  
  model <- inla(formula=formula1, data = data1, family = "poisson"
           ,control.predictor = list(compute = TRUE) 
           ,control.compute=list(cpo=TRUE,#resulting object containsvalues for CPO and PIT, 
                                 return.marginals.predictor=TRUE#Posterior marginals
                                 ,config= TRUE) #计算抽样
           ,control.inla = list(strategy = "auto")
              )
  print(summary(model))
  
  size_sample<- 5000
  model.samples = inla.posterior.sample(size_sample, model)

  ht<-factor(newdata$HomeTeam,levels =c("0",unique(data$HomeTeam)))
  at<-factor(newdata$AwayTeam,levels =c("0",unique(data$HomeTeam)))
  htn<-as.numeric(ht)
  atn<-as.numeric(at)
  pred_probs<- matrix(0,nrow= size_valid,ncol = 3)
  for (j in 1:size_valid) {
    z1<-0
    z2<-0
    z3<-0
    i<-1  
    for (i in 1:size_sample) {
      latent<-model.samples[[i]]$latent
      length<-length(latent)
      coef<-latent[((length-41):length)]
      hg<-rpois(1,exp(coef[1]+coef[42]+coef[htn[j]+1]+coef[21+atn[j]]))
      ag<-rpois(1,exp(coef[1]+coef[atn[j]+1]+coef[21+htn[j]]))
      
      if(hg>ag) z3 <-z3 +1 
      else if(hg==ag) z2 <- z2+1
      else z1<- z1+1
    }
    pred_probs[j,]<-c(z1,z2,z3)/size_sample
    
  }
  print(pred_probs)
  
  obs_outcomes<-ComputeOut(newdata)
  
  rps<- RPS(pred_probs ,obs_outcomes)
  
  obs_results <- factor(sign(newdata$FTHG-newdata$FTAG) ,levels = c(-1,0,1))
  
  acc<- accuracy(pred_probs , as.numeric(obs_results))
  
  log_sum<-log_sum_max(pred_probs )
  
  sum<- sum_max(pred_probs )
  
  max_odds <-as.matrix (newdata[,c("MaxA","MaxD","MaxH")])
  
  avg_odds <-as.matrix (newdata[,c("AvgA","AvgD","AvgH")])
  
  ide_profit <- ideal_profit_sum(pred_probs  ,max_odds)
  
  act_profit <- actual_profit_sum (pred_probs  ,max_odds ,obs_outcomes)
  
  cat("\naccurary=",acc,"rps=",rps,"\n",sep=" ")
  cat("sum",sum,"logsum=",log_sum,"ideal odds profit",ide_profit," ",sep=" ")
  cat("actual profit", act_profit,sep = " ")
  return (data.frame(acc=acc,rps=rps,sum=sum,log_sum=log_sum
                     ,act_profit=act_profit,ide_profit=ide_profit
                                          )
          )
}
#' @title A dataset used for HPL.
#' @name size_train
#' @description This dataset is used to fit the HPL model
NULL
#' @title A dataset used for HPL.
#' @name size_valid
#' @description This dataset is used to fit the HPL model
NULL
#' @title A dataset used for HPL.
#' @name data
#' @description This dataset is used to fit the HPL model
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' b<-HPL(data,size_train = 360,size_valid = 20)
#' }
NULL
#' @import Rcpp 
#' @import Matrix
#' @import foreach
#' @import parallel
#' @import sp
#' @import INLA
#' @import microbenchmark
#' @import DAAG
#' @import boot
#' @import bootstrap
#' @import MASS
#' @import CDFt
#' @importFrom VGAM qlaplace
#' @importFrom stats rpois
#' @useDynLib SA23204181
NULL

