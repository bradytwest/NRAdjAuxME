# Simulation from West and Little (2013), using PMM-MI approach based on a PMM

trivar.sim.pmm.mi <- function(nsamps,sampsize,pi1,rho,myseed,m=100)
{
   require(mnormt)

   vcov.0 <- cbind(c(1,rho,0.25),c(rho,1,0.5),c(0.25,0.5,1))
   #meanvec.0 <- c(1.1,1,9.5) # for rho = 0.9
   meanvec.0 <- c(1.4,1,10.5) # for rho = 0.6

   vcov.1 <- cbind(c(1,rho,0.25),c(rho,1,0.5),c(0.25,0.5,1))
   #meanvec.1 <- c(2,2,10) # for rho = 0.9
   meanvec.1 <- c(2,2,11) # for rho = 0.6
   
   results <- matrix(0,nrow=nsamps,ncol=14)
   set.seed(41279)

   for (s in 1:nsamps)
   {

      # simulate samples from models for each pattern, using value of pi1 from original model
    
      n0 <- rbinom(1,sampsize,1-pi1)
      n1 <- sampsize-n0
      trivar.0 <- rmnorm(n0, meanvec.0, vcov.0)
      trivar.1 <- rmnorm(n1, meanvec.1, vcov.1)

      # check marginal means defined by PMM

      results[s,1] <- (1-pi1)*meanvec.0[1] + pi1*meanvec.1[1]
      results[s,10] <- n0/sampsize # response rate
      results[s,11] <- (1-pi1)*meanvec.0[2] + pi1*meanvec.1[2]
      results[s,12] <- (1-pi1)*meanvec.0[3] + pi1*meanvec.1[3]

        # Draws for Bayesian inference

       s11.c.0.draws <- numeric(m)
       b1c.c.0.draws <- numeric(m)
       s11.c.1.draws <- numeric(m)
       b1c.c.1.draws <- numeric(m)
       s11.2c.draws <- numeric(m)
       s33.2c.draws <- numeric(m)
       b12.2c.draws <- numeric(m)
       b10.2c.draws <- numeric(m)
       b32.2c.draws <- numeric(m)
       b30.2c.draws <- numeric(m)

       # For Draws of ML estimates
       b2c.c.1.draws <- numeric(m)
       b3c.c.1.draws <- numeric(m)
       s12.c.1.draws <- numeric(m)
       s22.c.1.draws <- numeric(m)
       s13.c.1.draws <- numeric(m)
       s23.c.1.draws <- numeric(m)
       s33.c.1.draws <- numeric(m)

       # For Draws of Regression Parameters in Conditional Distributions, for missing pattern imputations
       b2c.1c.1.draws <- numeric(m)
       b21.1c.1.draws <- numeric(m)
       s22.1c.1.draws <- numeric(m)
       b3c.21c.1.draws <- numeric(m)
       b31.21c.1.draws <- numeric(m)
       b32.21c.1.draws <- numeric(m)
       s33.21c.1.draws <- numeric(m)

       # Regress X1 on 1 for complete cases

       fit1c.0 <- lm(trivar.0[,1] ~ 1)
       b1c.c.0 <- summary(fit1c.0)$coef[,1]
       sigma11.c.0 <- summary(fit1c.0)$sigma^2
       ssc1.c.0 <- t(model.matrix(fit1c.0)) %*% model.matrix(fit1c.0)

       # Regress X1 on 1 for incomplete cases

       fit1c.1 <- lm(trivar.1[,1] ~ 1)
       b1c.c.1 <- summary(fit1c.1)$coef[,1]
       sigma11.c.1 <- summary(fit1c.1)$sigma^2
       ssc1.c.1 <- t(model.matrix(fit1c.1)) %*% model.matrix(fit1c.1)

       # Regress X1 on X2 and 1 for complete cases

       fit12c.0 <- lm(trivar.0[,1] ~ trivar.0[,2])
       b12.2c <- summary(fit12c.0)$coef[2,1]
       sigma11.2c.0 <- summary(fit12c.0)$sigma^2

       # Regress X3 on X2 and 1 for complete cases

       fit32c.0 <- lm(trivar.0[,3] ~ trivar.0[,2])
       b32.2c <- summary(fit32c.0)$coef[2,1]
       sigma33.2c.0 <- summary(fit32c.0)$sigma^2

       # Regress X2 on 1 for complete cases
 
       fit2c.0 <- lm(trivar.0[,2] ~ 1)
       b2c.c.0 <- summary(fit2c.0)$coef[,1]
       sigma22.c.0 <- summary(fit2c.0)$sigma^2

       # Regress X3 on 1 for complete cases
 
       fit3c.0 <- lm(trivar.0[,3] ~ 1)
       b3c.c.0 <- summary(fit3c.0)$coef[,1]
       sigma33.c.0 <- summary(fit3c.0)$sigma^2

       # Estimated covariances conditional on C for complete cases
      
       s12.c.0 <- cov(fit1c.0$resid,fit2c.0$resid)
       s13.c.0 <- cov(fit1c.0$resid,fit3c.0$resid)
       s23.c.0 <- cov(fit2c.0$resid,fit3c.0$resid)

       # Means for complete cases

       mu1.0 <- mean(trivar.0[,1])
       mu2.0 <- mean(trivar.0[,2])
       mu3.0 <- mean(trivar.0[,3])

       # create vectors to hold estimates from each MI analysis

       x2meanimpvec <- numeric(m)
       x2varimpvec <- numeric(m)
       x3meanimpvec <- numeric(m)
       x3varimpvec <- numeric(m)

       for (d in 1:m) 
       {

          while (s22.1c.1.draws[d] <= 0 | s33.21c.1.draws[d] <= 0) # make sure draws of residual variances are positive
          { 

          # sequence of draws described in paper

          s11.c.0.draws[d] <- n0*sigma11.c.0 / rchisq(1,n0-1)
          b1c.c.0.draws[d] <- rmnorm(1, b1c.c.0, solve(ssc1.c.0)*s11.c.0.draws[d])
          s11.c.1.draws[d] <- n1*sigma11.c.1 / rchisq(1,n1-1)
          b1c.c.1.draws[d] <- rmnorm(1, b1c.c.1, solve(ssc1.c.1)*s11.c.1.draws[d])
          s11.2c.draws[d] <- n0*sigma11.2c.0 / rchisq(1,n0-1)
          s33.2c.draws[d] <- n0*sigma33.2c.0 / rchisq(1,n0-1)
          b12.2c.draws[d] <- rnorm(1, b12.2c, sqrt(s11.2c.draws[d] / (n0*sigma22.c.0)))
          b10.2c.draws[d] <- rnorm(1, mu1.0 - b12.2c.draws[d]*mu2.0, sqrt(s11.2c.draws[d] / n0))
          b32.2c.draws[d] <- rnorm(1, b32.2c, sqrt(s33.2c.draws[d] / (n0*sigma22.c.0)))
          b30.2c.draws[d] <- rnorm(1, mu3.0 - b32.2c.draws[d]*mu2.0, sqrt(s33.2c.draws[d] / n0))

          # compute draws of ML estimates based on draws above
          b2c.c.1.draws[d] <- b2c.c.0 + (b1c.c.1.draws[d] - b1c.c.0.draws[d]) / b12.2c.draws[d]
          b3c.c.1.draws[d] <- b3c.c.0 + (b1c.c.1.draws[d] - b1c.c.0.draws[d]) * (b32.2c.draws[d] / b12.2c.draws[d])
          s12.c.1.draws[d] <- s12.c.0 + (s11.c.1.draws[d] - s11.c.0.draws[d]) / b12.2c.draws[d]
          s22.c.1.draws[d] <- sigma22.c.0 + (s11.c.1.draws[d] - s11.c.0.draws[d]) / b12.2c.draws[d]^2
          s13.c.1.draws[d] <- s13.c.0 + b32.2c.draws[d] * (s11.c.1.draws[d] - s11.c.0.draws[d]) / b12.2c.draws[d]       
          s23.c.1.draws[d] <- s23.c.0 + b32.2c.draws[d] * (s11.c.1.draws[d] - s11.c.0.draws[d]) / b12.2c.draws[d]^2 
          s33.c.1.draws[d] <- sigma33.c.0 + b32.2c.draws[d]^2 * (s11.c.1.draws[d] - s11.c.0.draws[d]) / b12.2c.draws[d]^2 

          # compute draws of parameters in conditional distributions for imputations
	  b2c.1c.1.draws[d] <- b2c.c.1.draws[d] - b1c.c.1.draws[d] * (s12.c.1.draws[d] / s11.c.1.draws[d])
          b21.1c.1.draws[d] <- s12.c.1.draws[d] / s11.c.1.draws[d]
          s22.1c.1.draws[d] <- s22.c.1.draws[d] - s12.c.1.draws[d]^2 / s11.c.1.draws[d]

          b3c.21c.1.draws[d] <- (b3c.c.1.draws[d] - b1c.c.1.draws[d] * (s13.c.1.draws[d] / s11.c.1.draws[d])) -
	  	(b2c.c.1.draws[d] - b1c.c.1.draws[d] * (s12.c.1.draws[d] / s11.c.1.draws[d])) *
	  	(s23.c.1.draws[d] - s12.c.1.draws[d]*s13.c.1.draws[d] / s11.c.1.draws[d]) / (s22.c.1.draws[d] - s12.c.1.draws[d]^2 / s11.c.1.draws[d])

          b31.21c.1.draws[d] <- (s13.c.1.draws[d] / s11.c.1.draws[d]) - (s12.c.1.draws[d] / s11.c.1.draws[d]) *
	  	(s23.c.1.draws[d] - s12.c.1.draws[d]*s13.c.1.draws[d] / s11.c.1.draws[d]) / (s22.c.1.draws[d] - s12.c.1.draws[d]^2 / s11.c.1.draws[d])

          b32.21c.1.draws[d] <- (s23.c.1.draws[d] - s12.c.1.draws[d]*s13.c.1.draws[d] / s11.c.1.draws[d]) / 
	  	(s22.c.1.draws[d] - s12.c.1.draws[d]^2 / s11.c.1.draws[d]) 	
          
	  s33.21c.1.draws[d] <- (s33.c.1.draws[d] - s13.c.1.draws[d]^2 / s11.c.1.draws[d]) -
	  	(s23.c.1.draws[d] - s12.c.1.draws[d]*s13.c.1.draws[d] / s11.c.1.draws[d])^2 / (s22.c.1.draws[d] - s12.c.1.draws[d]^2 / s11.c.1.draws[d])

 	  } # end of while loop

          # impute missing values based on conditional distributions
                       
          trivar.1[,2] <- rnorm(length(trivar.1[,2]), mean = model.matrix(fit1c.1) %*% b2c.1c.1.draws[d] + b21.1c.1.draws[d] * trivar.1[,1], 
	  	sd = sqrt(s22.1c.1.draws[d]))

	  trivar.1[,3] <- rnorm(length(trivar.1[,3]), mean = model.matrix(fit1c.1) %*% b3c.21c.1.draws[d] + b31.21c.1.draws[d] * trivar.1[,1] +
	  	b32.21c.1.draws[d] * trivar.1[,2], sd = sqrt(s33.21c.1.draws[d])) 

         # stack data sets
         trivar.imp <- rbind(trivar.0,trivar.1)
         trivar.imp.data <- data.frame(trivar.imp)
         
         x2meanimpvec[d] <- mean(trivar.imp.data$X2)
         x2varimpvec[d] <- var(trivar.imp.data$X2) / length(trivar.imp.data$X2)
         x3meanimpvec[d] <- mean(trivar.imp.data$X3)
         x3varimpvec[d] <- var(trivar.imp.data$X3) / length(trivar.imp.data$X3)

       } # end of multiple loop
      
      # MI Inference

      results[s,2] <- mean(x2meanimpvec) 
      results[s,3] <- mean(x3meanimpvec) 

      # MI RMSE Calculations

      results[s,4] <- (results[s,2] - results[s,11])^2 
      results[s,5] <- (results[s,3] - results[s,12])^2 

      # For MI confidence interval coverage and width

      mi.var.x2 <- mean(x2varimpvec) + (1 + 1/m) * var(x2meanimpvec)
      mi.var.x3 <- mean(x3varimpvec) + (1 + 1/m) * var(x3meanimpvec)
      df.x2 <- (m-1) / ((1 + 1/m) * var(x2meanimpvec) / mi.var.x2)^2
      df.x3 <- (m-1) / ((1 + 1/m) * var(x3meanimpvec) / mi.var.x3)^2

      x2mean.ll <- mean(x2meanimpvec) - qt(0.975,df.x2)*sqrt(mi.var.x2)
      x2mean.ul <- mean(x2meanimpvec) + qt(0.975,df.x2)*sqrt(mi.var.x2)
      x3mean.ll <- mean(x3meanimpvec) - qt(0.975,df.x3)*sqrt(mi.var.x3)
      x3mean.ul <- mean(x3meanimpvec) + qt(0.975,df.x3)*sqrt(mi.var.x3)

      if (x2mean.ll <= results[s,11] && x2mean.ul > results[s,11]) results[s,6] <- 1
      else results[s,6] <- 0
      results[s,7] <- x2mean.ul - x2mean.ll

      if (x3mean.ll <= results[s,12] && x3mean.ul > results[s,12]) results[s,8] <- 1
      else results[s,8] <- 0
      results[s,9] <- x3mean.ul - x3mean.ll

      results[s,13] <- (1 + 1/m) * var(x2meanimpvec) / mi.var.x2 # FMI for X2 mean
      results[s,14] <- (1 + 1/m) * var(x3meanimpvec) / mi.var.x3 # FMI for X3 mean

      cat("Sample:"," ",s,"\n")

   } # end of for loop for simulations

   cat("Mean response rate is:"," ",mean(results[,10]),"\n")

   cat("True mean of X2 is:"," ",mean(results[,11]),"\n")
   cat("\n")

   cat("Mean of PMM-MI estimates of X2 mean for respondents is:"," ",mean(results[,2]),"\n")
   cat("RelBias of PMM-MI estimates of X2 mean is:"," ",((mean(results[,2]) - mean(results[,11]))/mean(results[,11]))*100,"\n")
   cat("RMSE of PMM-MI estimates of X2 mean is:"," ",sqrt(mean(results[,4])),"\n")
   cat("95% CI coverage for PMM-MI estimate of X2 mean is:"," ",mean(results[,6]),"\n")
   cat("95% CI width for PMM-MI estimate of X2 mean is:"," ",mean(results[,7]),"\n")
   cat("FMI for mean of X2 is:"," ",mean(results[,13]),"\n")
   cat("\n")

   cat("True mean of X3 is:"," ",mean(results[,12]),"\n")
   cat("\n")

   cat("Mean of PMM-MI estimates of X3 mean for respondents is:"," ",mean(results[,3]),"\n")
   cat("RelBias of PMM-MI estimates of X3 mean is:"," ",((mean(results[,3]) - mean(results[,12]))/mean(results[,12]))*100,"\n")
   cat("RMSE of PMM-MI estimates of X3 mean is:"," ",sqrt(mean(results[,5])),"\n")
   cat("95% CI coverage for PMM-MI estimate of X3 mean is:"," ",mean(results[,8]),"\n")
   cat("95% CI width for PMM-MI estimate of X3 mean is:"," ",mean(results[,9]),"\n")
   cat("FMI for mean of X3 is:"," ",mean(results[,14]),"\n")
   cat("\n")

   cat("Mean of Marginal means of X1 in full data set simulated from PMM is:"," ",mean(results[,1]),"\n")
   cat("Mean of Marginal means of X2 in full data set simulated from PMM is:"," ",mean(results[,11]),"\n")
   cat("Mean of Marginal means of X3 in full data set simulated from PMM is:"," ",mean(results[,12]),"\n")   
   cat("\n")

   return(results)
}

tres <- trivar.sim.pmm.mi(nsamps=1000,sampsize=1000,pi1=0.50,rho=0.9,myseed=41279,m=100)
tres <- trivar.sim.pmm.mi(nsamps=1000,sampsize=1000,pi1=0.25,rho=0.9,myseed=41279,m=100)

# CHANGE PARAMETERS
tres <- trivar.sim.pmm.mi(nsamps=1000,sampsize=1000,pi1=0.50,rho=0.6,myseed=41279,m=100)
tres <- trivar.sim.pmm.mi(nsamps=1000,sampsize=1000,pi1=0.25,rho=0.6,myseed=41279,m=100)



