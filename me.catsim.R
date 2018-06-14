# R function implementing simulation study reported in Chapter 15 of Kreuter (2013) Paradata edited volume

# Author: Brady T. West
# Updated: February 2012

me.catsim <- function (beta.x.obs,beta.x.y,beta.obs.y,beta.x.p,beta.obs.p,fpr.r,fnr.r,fpr.nr,fnr.nr,popsize=10000,sampsize=500,nsamps=1000,myseed) 
{
   # browser()

   require(survey)
   require(rms)

   set.seed(myseed)

   # initialize matrix of simulation results
   simresults <- matrix(0,ncol=19,nrow=nsamps)

   # define population

   # initialize data vectors
   x <- c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000),rep(7,1000),rep(8,1000),rep(9,1000),rep(10,1000))
   y.t <- numeric(popsize)
   obs.t <- numeric(popsize)

   # generate true values for auxiliary variable being observed (function of x)
   pobs <- plogis(0.25 + beta.x.obs * x)
   uni <- runif(popsize)
   obs.t[uni <= pobs] <- 1
   obs.t[uni > pobs] <- 0

   # define binary outcome values
   py.t <- plogis(-0.5 + beta.x.y * x + beta.obs.y * obs.t)
   uni.y <- runif(popsize)
   y.t[uni.y <= py.t] <- 1
   y.t[uni.y > py.t] <- 0

   # save pseudo r-squared for reality check
   fit <- lrm(y.t ~ x + obs.t)
   r2 <- fit$stats[10]

   for (a in 1:nsamps)
   { # start of for loop

      simresults[a,1] <- mean(y.t)

      # draw sample from population 
      s <- sample(1:10000,sampsize)
      y.t.samp <- y.t[s]
      x.samp <- x[s]
      obs.t.samp <- obs.t[s]       
      obs.error <- numeric(sampsize)
      ri <- numeric(sampsize)

      # assign response indicators
      pi <- plogis(-1 + beta.x.p * x.samp + beta.obs.p * obs.t.samp)
      puni <- runif(sampsize)
      ri[puni <= pi] <- 1
      ri[puni > pi] <- 0

      # generate binary indicator measured with error according to specified error rates
      uni.error <- runif(sampsize)
      obs.error[ri == 1 & uni.error <= fnr.r & obs.t.samp == 1] <- 0
      obs.error[ri == 1 & uni.error <= fpr.r & obs.t.samp == 0] <- 1
      obs.error[ri == 1 & uni.error > fnr.r & obs.t.samp == 1] <- 1
      obs.error[ri == 1 & uni.error > fpr.r & obs.t.samp == 0] <- 0
      obs.error[ri == 0 & uni.error <= fnr.nr & obs.t.samp == 1] <- 0
      obs.error[ri == 0 & uni.error <= fpr.nr & obs.t.samp == 0] <- 1
      obs.error[ri == 0 & uni.error > fnr.nr & obs.t.samp == 1] <- 1
      obs.error[ri == 0 & uni.error > fpr.nr & obs.t.samp == 0] <- 0
 
      # save error rate to see what average error rates are
      simresults[a,19] <- length(obs.error[obs.error != obs.t.samp]) / sampsize

      d <- data.frame(y.t.samp,x.samp,obs.t.samp,obs.error,ri)

      # Complete Case Analysis

      simresults[a,2] <- mean(y.t.samp[ri == 1])
      mean0.ll <- mean(y.t.samp[ri == 1]) - 1.96*sd(y.t.samp[ri == 1])/sqrt(length(y.t.samp[ri == 1]))
      mean0.ul <- mean(y.t.samp[ri == 1]) + 1.96*sd(y.t.samp[ri == 1])/sqrt(length(y.t.samp[ri == 1]))
      if (mean0.ll <= simresults[a,1] && mean0.ul > simresults[a,1]) simresults[a,3] <- 1
      else simresults[a,3] <- 0
      simresults[a,4] <- (simresults[a,2] - simresults[a,1])^2 # squared error
      simresults[a,5] <- (simresults[a,2] - simresults[a,1]) # bias

      padj <- numeric(sampsize)
      padj2 <- numeric(sampsize)
      padj3 <- numeric(sampsize)

      # Weighting Class Analysis
      for (i in 1:sampsize)
      {
         for (j in 1:10)
	 {
            if (x.samp[i] == j & obs.error[i] == 1) padj[i] <- mean(ri[x.samp == j & obs.error == 1])
            else if (x.samp[i] == j & obs.error[i] == 0) padj[i] <- mean(ri[x.samp == j & obs.error == 0])
            if (x.samp[i] == j & obs.t.samp[i] == 1) padj2[i] <- mean(ri[x.samp == j & obs.t.samp == 1])
            else if (x.samp[i] == j & obs.t.samp[i] == 0) padj2[i] <- mean(ri[x.samp == j & obs.t.samp == 0]) 
            if (x.samp[i] == j) padj3[i] <- mean(ri[x.samp == j]) # xonly
         }
      }

      wc.adj <- 1 / padj
      wc.adj2 <- 1 / padj2
      wc.adj3 <- 1 / padj3

      temp.data <- data.frame(y.t.samp,x.samp,obs.t.samp,obs.error,ri,wc.adj,padj,wc.adj2,padj2)

      # Error-prone
      wc.data <- data.frame(y.t = y.t.samp[ri == 1], wc.adj = wc.adj[ri == 1])     
 
      wc.des <- svydesign(id=~1,weights=~wc.adj,data=wc.data)
      simresults[a,6] <- coef(svymean(~y.t,wc.des)) # Weighted R mean

      wcmean.ll <- confint(svymean(~y.t,wc.des))[1]
      wcmean.ul <- confint(svymean(~y.t,wc.des))[2]
      if (wcmean.ll <= simresults[a,1] && wcmean.ul > simresults[a,1]) simresults[a,7] <- 1
      else simresults[a,7] <- 0 
      simresults[a,8] <- (simresults[a,6] - simresults[a,1])^2 # squared error
      simresults[a,9] <- (simresults[a,6] - simresults[a,1]) # bias  

      # Error-free
      wc.data2 <- data.frame(y.t = y.t.samp[ri == 1], wc.adj2 = wc.adj2[ri == 1])     
 
      wc.des2 <- svydesign(id=~1,weights=~wc.adj2,data=wc.data2)
      simresults[a,11] <- coef(svymean(~y.t,wc.des2)) # Weighted R mean

      wcmean.ll <- confint(svymean(~y.t,wc.des2))[1]
      wcmean.ul <- confint(svymean(~y.t,wc.des2))[2]
      if (wcmean.ll <= simresults[a,1] && wcmean.ul > simresults[a,1]) simresults[a,12] <- 1
      else simresults[a,12] <- 0 
      simresults[a,13] <- (simresults[a,11] - simresults[a,1])^2 # squared error
      simresults[a,14] <- (simresults[a,11] - simresults[a,1]) # bias

      # XONLY
      wc.data3 <- data.frame(y.t = y.t.samp[ri == 1], wc.adj3 = wc.adj3[ri == 1])
      wc.des3 <- svydesign(id=~1,weights=~wc.adj3,data=wc.data3)
      simresults[a,15] <- coef(svymean(~y.t,wc.des3)) # Weighted R mean

      wcmean.ll <- confint(svymean(~y.t,wc.des3))[1]
      wcmean.ul <- confint(svymean(~y.t,wc.des3))[2]
      if (wcmean.ll <= simresults[a,1] && wcmean.ul > simresults[a,1]) simresults[a,16] <- 1
      else simresults[a,16] <- 0 
      simresults[a,17] <- (simresults[a,15] - simresults[a,1])^2 # squared error
      simresults[a,18] <- (simresults[a,15] - simresults[a,1]) # bias
 
      simresults[a,10] <- mean(ri)   

      cat("Sample", a, "\n")

   } # end of for loop for samples

   simresults2 <- matrix(0,ncol=13,nrow=1)

   # ratio of FPR/FNR for R to FPR/FNR for NR
   error.ratio <- (fpr.r/fnr.r) / (fpr.nr/fnr.nr)
  
   cat("True Mean of Y:"," ",mean(y.t),"\n")
   cat("Pseudo R-squared:"," ",r2,"\n")
   cat("Proportion of Pop with True Obs = 1:"," ",mean(obs.t),"\n")
   cat("Mean Response Rate:"," ",mean(simresults[,10]),"\n")
   cat("Mean Overall Error Rate:"," ",mean(simresults[,19]),"\n")
   cat("\n")
   simresults2[1,1] <- mean(simresults[,19])
   cat("Bias of CC mean estimate is:"," ",mean(simresults[,5]),"\n")
   simresults2[1,2] <- mean(simresults[,5]) # CC Bias
   cat("RMSE of CC mean is:"," ",sqrt(mean(simresults[,4])),"\n")
   simresults2[1,3] <- sqrt(mean(simresults[,4])) # CC RMSE
   cat("Coverage of mean is:"," ",mean(simresults[,3]),"\n")
   simresults2[1,4] <- mean(simresults[,3]) # CC Coverage
   cat("\n")
   cat("Bias of WCE mean estimate is:"," ",mean(simresults[,9]),"\n")
   simresults2[1,5] <- mean(simresults[,9]) # WCE Bias
   cat("RMSE of WCE mean is:"," ",sqrt(mean(simresults[,8])),"\n")
   simresults2[1,6] <- sqrt(mean(simresults[,8])) # WCE RMSE
   cat("Coverage of WCE mean is:"," ",mean(simresults[,7]),"\n")
   simresults2[1,7] <- mean(simresults[,7]) # WCE Coverage
   cat("\n")
   cat("Bias of WCT mean estimate is:"," ",mean(simresults[,14]),"\n")
   simresults2[1,8] <- mean(simresults[,14]) # WCT Bias
   cat("RMSE of WCT mean is:"," ",sqrt(mean(simresults[,13])),"\n")
   simresults2[1,9] <- sqrt(mean(simresults[,13])) # WCT RMSE
   cat("Coverage of WCT mean is:"," ",mean(simresults[,12]),"\n")
   simresults2[1,10] <- mean(simresults[,12]) # WCT coverage
   cat("\n")
   cat("Bias of XONLY mean estimate is:"," ",mean(simresults[,18]),"\n")
   simresults2[1,11] <- mean(simresults[,18]) # WCT Bias
   cat("RMSE of XONLY mean is:"," ",sqrt(mean(simresults[,17])),"\n")
   simresults2[1,12] <- sqrt(mean(simresults[,17])) # WCT RMSE
   cat("Coverage of XONLY mean is:"," ",mean(simresults[,16]),"\n")
   simresults2[1,13] <- mean(simresults[,16]) # WCT coverage
   cat("\n")

   return(simresults2) #simresults2
}

# 1A: ERR = 1, FPR > FNR for both R and NR, Bdp < 0, Bdy > 0 (constant neg. bias, bias increases as error increases, constant highest RMSE, good coverage)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.2,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.3,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.1,fpr.nr=0.5,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.6,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 1B: ERR = 1, FPR = FNR for both R and NR, Bdp < 0, Bdy > 0 (neg. bias increases as error increases, steadily increasing RMSE, good coverage that is slightly decreasing)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.4,fpr.nr=0.4,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.5,fpr.nr=0.5,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.6,fnr.r=0.6,fpr.nr=0.6,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 1C: ERR = 1, FPR < FNR for both R and NR, Bdp < 0, Bdy > 0 (neg. bias increases as error increases, increasing RMSE, decreasing coverage as error increases)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.1,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.1,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.1,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.5,fpr.nr=0.1,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.6,fpr.nr=0.1,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 2A: ERR = 1, FPR > FNR for both R and NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.2,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.3,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.1,fpr.nr=0.5,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.6,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.1,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 2B: ERR = 1, FPR = FNR for both R and NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.4,fpr.nr=0.4,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.5,fpr.nr=0.5,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.6,fnr.r=0.6,fpr.nr=0.6,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 2C: ERR = 1, FPR < FNR for both R and NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.1,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.1,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.1,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.5,fpr.nr=0.1,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.6,fpr.nr=0.1,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)


# 3A: FPR > FNR for both R and NR, Error rates greater for NR, Bdp < 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.1,fpr.nr=0.8,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5) 
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.9,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 3B: FPR = FNR for both R and NR, Error rates greater for NR, Bdp < 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.2,fpr.nr=0.5,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.3,fpr.nr=0.6,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.4,fpr.nr=0.7,fnr.nr=0.7,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.5,fpr.nr=0.8,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 3C: FNR > FPR for both R and NR, Error rates greater for NR, Bdp < 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.3,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.2,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.2,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5) 
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.9,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 4A: FPR > FNR for both R and NR, Error rates greater for NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.6,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.1,fpr.nr=0.8,fnr.nr=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5) 
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.1,fpr.nr=0.9,fnr.nr=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 4B: FPR = FNR for both R and NR, Error rates greater for NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.1,fpr.nr=0.4,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.2,fnr.r=0.2,fpr.nr=0.5,fnr.nr=0.5,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.3,fnr.r=0.3,fpr.nr=0.6,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.4,fnr.r=0.4,fpr.nr=0.7,fnr.nr=0.7,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.5,fnr.r=0.5,fpr.nr=0.8,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 4C: FNR > FPR for both R and NR, Error rates greater for NR, Bdp > 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.3,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.2,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.9,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.2,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=1,nsamps=1000,myseed=6) 
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 5A: FPR > FNR for both R and NR, Error rates greater for R, Bdp < 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.0,fnr.nr=0.0,fpr.r=0.0,fnr.r=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.2,fnr.nr=0.1,fpr.r=0.4,fnr.r=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.3,fnr.nr=0.1,fpr.r=0.6,fnr.r=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.2,fnr.nr=0.1,fpr.r=0.6,fnr.r=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.4,fnr.nr=0.1,fpr.r=0.8,fnr.r=0.2,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5) 
res5 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.3,fnr.nr=0.1,fpr.r=0.9,fnr.r=0.3,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# 5C: FNR > FPR for both R and NR, Error rates greater for R, Bdp < 0, Bdy > 0 (...)
res0 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.0,fnr.nr=0.0,fpr.r=0.0,fnr.r=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.1,fnr.nr=0.2,fpr.r=0.2,fnr.r=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.1,fnr.nr=0.2,fpr.r=0.3,fnr.r=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.1,fnr.nr=0.3,fpr.r=0.2,fnr.r=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.1,fnr.nr=0.3,fpr.r=0.3,fnr.r=0.9,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=0.1,fpr.nr=0.1,fnr.nr=0.4,fpr.r=0.2,fnr.r=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6) 
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

# Alternative correlations of X and D (Scenario 3)
res0 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.2,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.2,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4) 
res4 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.3,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5)
res5 <- me.catsim(beta.x.obs=-0.5,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.9,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)

res0 <- me.catsim(beta.x.obs=0.5,fpr.r=0.0,fnr.r=0.0,fpr.nr=0.0,fnr.nr=0.0,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=1) 
res1 <- me.catsim(beta.x.obs=0.5,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.2,fnr.nr=0.4,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=2)
res2 <- me.catsim(beta.x.obs=0.5,fpr.r=0.1,fnr.r=0.2,fpr.nr=0.3,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=3)
res3 <- me.catsim(beta.x.obs=0.5,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.2,fnr.nr=0.6,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=4)
res4 <- me.catsim(beta.x.obs=0.5,fpr.r=0.1,fnr.r=0.4,fpr.nr=0.2,fnr.nr=0.8,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=5) 
res5 <- me.catsim(beta.x.obs=0.5,fpr.r=0.1,fnr.r=0.3,fpr.nr=0.3,fnr.nr=0.9,beta.x.y=0.25,beta.obs.y=1,beta.x.p=0.25,beta.obs.p=-1,nsamps=1000,myseed=6)
res <- rbind(res0,res1,res2,res3,res4,res5)
resd <- data.frame(res)



# Plotting functions

par(mfrow=c(3,1))

plot(resd$X1, resd$X2, ylim=c(-0.05,0.05), type='l', xlab="Overall Average Error Rate", ylab = "Empirical Bias")
lines(resd$X1, resd$X5, lty=2)
lines(resd$X1, resd$X8, lty=3)
#lines(resd$X1, resd$X11, lty=4)
legend(locator(1),c("CC","WCE","WCT"),lty=c(1,2,3),cex=0.7)
#legend(locator(1),c("CC","WCE","WCT","XONLY"),lty=c(1,2,3,4),cex=0.7)

plot(resd$X1, resd$X3, ylim=c(0.02,0.05), type='l', xlab="Overall Average Error Rate", ylab = "Empirical RMSE")
lines(resd$X1, resd$X6, lty=2)
lines(resd$X1, resd$X9, lty=3)
#lines(resd$X1, resd$X12, lty=4)
legend(locator(1),c("CC","WCE","WCT"),lty=c(1,2,3),cex=0.7)
#legend(locator(1),c("CC","WCE","WCT","XONLY"),lty=c(1,2,3,4),cex=0.7)

plot(resd$X1, resd$X4, ylim=c(0.6,1.0), type='l', xlab="Overall Average Error Rate", ylab = "95% CI Coverage")
lines(resd$X1, resd$X7, lty=2)
lines(resd$X1, resd$X10, lty=3)
#lines(resd$X1, resd$X13, lty=4)
legend(locator(1),c("CC","WCE","WCT"),lty=c(1,2,3),cex=0.7)
#legend(locator(1),c("CC","WCE","WCT","XONLY"),lty=c(1,2,3,4),cex=0.7)





