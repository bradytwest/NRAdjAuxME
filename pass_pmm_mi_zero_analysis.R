###########################
# zero.mi.analysis function
###########################

# Author: Brady T. West
# Updated: January 2013
# Details: see West and Little (2013)

# aux.data has n rows and 6 columns (last three are strata, PSUs, weights),
# response has n indicators of response, and m is the number of multiples (draws).
# six columns are X1, X2, X3, strata, PSU, and weight (order matters); no (zero) auxiliary variables

zero.mi.analysis <- function(aux.data, response, m=100)
{
       # load survey package

       require(survey)

       # load package for multivariate normal draws

       require(mnormt)

       # read in data for two patterns

       trivar <- aux.data
       trivar.0 <- trivar[response == 1,]
       trivar.1 <- trivar[response == 0,]

       # sample sizes

       n0 <- length(trivar.0[,1]) # also defined as r in the paper (number of respondents)
       n1 <- length(trivar.1[,1]) # also defined as n - r in the paper (number of nonrespondents)

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

          # set survey features, apply survey mean to complete data (no missing data allowed on design features)

          trivar.imp.data2 <- trivar.imp.data[!is.na(trivar.imp.data[,6]) &  !is.na(trivar.imp.data[,5]) & !is.na(trivar.imp.data[,4]),]
          svyd <- svydesign(strata=~trivar.imp.data2[,4], id=~trivar.imp.data2[,5], weights=~trivar.imp.data2[,6] , data=trivar.imp.data2, nest=T, na.rm=T)

          # save estimate and variance for MI analysis

          x2meanimpvec[d] <- coef(svymean(~trivar.imp.data2[,2],svyd,na.rm=T))
          x3meanimpvec[d] <- coef(svymean(~trivar.imp.data2[,3],svyd,na.rm=T))
          x2varimpvec[d] <- SE(svymean(~trivar.imp.data2[,2],svyd,na.rm=T))^2
          x3varimpvec[d] <- SE(svymean(~trivar.imp.data2[,3],svyd,na.rm=T))^2
       }

       # MI Inference

       x2mean <- mean(x2meanimpvec)
       pred.mean <- exp(x2mean)
       cat("PMM-MI estimate of mean of X2 is:"," ",pred.mean,"\n")

       # For MI confidence interval coverage and width

       mi.var.x2 <- mean(x2varimpvec) + (1 + 1/m) * var(x2meanimpvec)
       mi.var.x3 <- mean(x3varimpvec) + (1 + 1/m) * var(x3meanimpvec)
       df.x2 <- (m - 1) / ((1 + 1/m) * var(x2meanimpvec) / mi.var.x2)^2
       df.x3 <- (m - 1) / ((1 + 1/m) * var(x3meanimpvec) / mi.var.x3)^2

       x2mean.ll <- mean(x2meanimpvec) - qt(0.975,df.x2)*sqrt(mi.var.x2)
       pred.mean.ll <- exp(x2mean.ll)
       x2mean.ul <- mean(x2meanimpvec) + qt(0.975,df.x2)*sqrt(mi.var.x2)
       pred.mean.ul <- exp(x2mean.ul)
       cat("95% PMM-MI CI for mean of X2 is:"," ",pred.mean.ll,",",pred.mean.ul,"\n")
       cat("95% PMM-MI CI Width is:"," ",pred.mean.ul - pred.mean.ll,"\n")
       cat("FMI for mean of X2 is:"," ",(1+1/m)*var(x2meanimpvec) / mi.var.x2,"\n")
       cat("\n")

       x3mean <- mean(x3meanimpvec)
       pred.mean3 <- exp(x3mean)
       cat("PMM-MI estimate of mean of X3 is:"," ",pred.mean3,"\n")

       x3mean.ll <- mean(x3meanimpvec) - qt(0.975,df.x3)*sqrt(mi.var.x3)
       pred.mean3.ll <- exp(x3mean.ll)
       x3mean.ul <- mean(x3meanimpvec) + qt(0.975,df.x3)*sqrt(mi.var.x3)
       pred.mean3.ul <- exp(x3mean.ul)
       cat("95% PMM-MI CI for mean of X3 is:"," ",pred.mean3.ll,",",pred.mean3.ul,"\n")
       cat("95% PMM-MI CI Width is:"," ",pred.mean3.ul - pred.mean3.ll,"\n")
       cat("FMI for mean of X3 is:"," ",(1+1/m)*var(x3meanimpvec) / mi.var.x3,"\n")
       cat("\n")

       imp.results <- data.frame(cbind(x2meanimpvec,x2varimpvec,x3meanimpvec,x3varimpvec))
       return(imp.results)

}

# Example (no auxiliary variables)
trivar <- read.csv("J:\\Dissertation\\Paper 3\\passdata.csv",h=T)
trivar$proxy <- trivar[,1] / 10000
trivar2 <- data.frame(cbind(trivar[,1],trivar[,3],trivar[,4],trivar[,7],trivar[,6],trivar[,5]))
resp <- trivar[,2]

zero.mi.analysis(trivar2,resp,m=100)

