# Other general simulation studies from West and Little (2013), based on a PMM

trivar.sim.pmm <- function(nsamps,sampsize,pi1,rho)
{
   require(mnormt)
   require(survey)
   require(mi)

   vcov.0 <- cbind(c(1,rho,0.25),c(rho,1,0.5),c(0.25,0.5,1))
   #meanvec.0 <- c(1.1,1,9.5) # for rho = 0.9
   meanvec.0 <- c(1.4,1,10.5) # for rho = 0.6

   vcov.1 <- cbind(c(1,rho,0.25),c(rho,1,0.5),c(0.25,0.5,1))
   #meanvec.1 <- c(2,2,10) # for rho = 0.9
   meanvec.1 <- c(2,2,11) # for rho = 0.6
   
   results <- matrix(0,nrow=nsamps,ncol=42)
   set.seed(41279)

   for (s in 1:nsamps)
   {

      # simulate samples from models for each pattern, using value of pi1 from original model
    
      n0 <- rbinom(1,sampsize,1-pi1)
      n1 <- sampsize-n0
      trivar.0 <- rmnorm(n0, meanvec.0, vcov.0)
      trivar.1 <- rmnorm(n1, meanvec.1, vcov.1)

      # prep for logistic regression (predict probability of RESPONSE with X1 for GW approach)

      trivar.stack <- rbind(trivar.0,trivar.1)  
      resp <- c(rep(1,n0),rep(0,n1)) # 1 for respondents (n0), 0 for nonrespondents (n1)
      resplogit <- glm(resp ~ trivar.stack[,1], family=binomial(link="logit"))
      nrweights <- 1 / resplogit$fit[resp == 1]

      # prep for MI procedure

      trivar.stack2 <- data.frame(trivar.stack)
      trivar.stack2$X2[resp == 0] <- NA
      trivar.stack2$X3[resp == 0] <- NA 

      # check marginal means defined by PMM

      results[s,9] <- (1-pi1)*meanvec.0[1] + pi1*meanvec.1[1]
      results[s,10] <- (1-pi1)*meanvec.0[2] + pi1*meanvec.1[2]
      results[s,11] <- (1-pi1)*meanvec.0[3] + pi1*meanvec.1[3]

      # check covariances of variables in respondent pattern

      results[s,14] <- vcov.0[1,3]
      results[s,15] <- vcov.0[2,3]
      results[s,16] <- vcov.0[1,2]

      # Compute CC estimates based on simulated data from two patterns

      results[s,1] <- mean(trivar.0[,1]) # R mean of X1
      results[s,2] <- mean(trivar.0[,2]) # R mean of X2
      results[s,3] <- mean(trivar.0[,3]) # R mean of X3
      results[s,4] <- n0/sampsize # response rate

      # CC RMSE Calculations

      results[s,5] <- (results[s,2] - results[s,10])^2 # For RMSE R mean X2
      results[s,6] <- (results[s,3] - results[s,11])^2 # For RMSE R mean X3

      # For CC interval coverage and width

      x2mean.ll <- results[s,2] - qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)
      x2mean.ul <- results[s,2] + qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)
      if (x2mean.ll <= results[s,10] && x2mean.ul > results[s,10]) results[s,7] <- 1
      else results[s,7] <- 0
      results[s,12] <- 2*qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)

      x3mean.ll <- results[s,3] - qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)
      x3mean.ul <- results[s,3] + qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)
      if (x3mean.ll <= results[s,11] && x3mean.ul > results[s,11]) results[s,8] <- 1
      else results[s,8] <- 0
      results[s,13] <- 2*qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)

      # Compute GW estimates based on simulated data from response pattern
     
      gwdes <- svydesign(id=~1,weights=~nrweights,data=data.frame(trivar.0))
      results[s,17] <- coef(svymean(~X2,gwdes)) # Weighted R mean of X2
      results[s,18] <- coef(svymean(~X3,gwdes)) # Weighted R mean of X3

      # GW RMSE Calculations

      results[s,19] <- (results[s,17] - results[s,10])^2 # For RMSE Weighted R mean X2
      results[s,20] <- (results[s,18] - results[s,11])^2 # For RMSE Weighted R mean X3

      # For GW interval coverage and width

      x2meanGW.ll <- confint(svymean(~X2,gwdes))[1]
      x2meanGW.ul <- confint(svymean(~X2,gwdes))[2]
      if (x2meanGW.ll <= results[s,10] && x2meanGW.ul > results[s,10]) results[s,21] <- 1
      else results[s,21] <- 0
      results[s,22] <- confint(svymean(~X2,gwdes))[2] - confint(svymean(~X2,gwdes))[1]

      x3meanGW.ll <- confint(svymean(~X3,gwdes))[1]
      x3meanGW.ul <- confint(svymean(~X3,gwdes))[2]
      if (x3meanGW.ll <= results[s,11] && x3meanGW.ul > results[s,11]) results[s,23] <- 1
      else results[s,23] <- 0
      results[s,24] <- confint(svymean(~X3,gwdes))[2] - confint(svymean(~X3,gwdes))[1]

      # Draws for Bayesian inference

      pi0draws <- numeric(1000)
      pi1draws <- numeric(1000)
      s22draws <- numeric(1000)
      s11.1.draws <- numeric(1000)
      s11.2.0.draws <- numeric(1000)
      b12.2.draws <- numeric(1000)
      b10.2.draws <- numeric(1000)
      mu.2.0.draws <- numeric(1000)
      mu.1.1.draws <- numeric(1000)
      s33.0.draws <- numeric(1000)
      s33.2.0.draws <- numeric(1000)
      b32.2.draws <- numeric(1000)
      b30.2.draws <- numeric(1000)
      mu2est.draws <- numeric(1000)
      mu3est.draws <- numeric(1000)

      # Regress X1 on X2 for complete cases

      fit12 <- lm(trivar.0[,1] ~ trivar.0[,2])
      b10.2 <- summary(fit12)$coef[1,1]
      b12.2 <- summary(fit12)$coef[2,1]
      b12.2.se <- summary(fit12)$coef[2,2]
      sigma11.2 <- summary(fit12)$sigma^2

      # Regress X3 on X2 for complete cases

      fit32 <- lm(trivar.0[,3] ~ trivar.0[,2])
      b30.2 <- summary(fit32)$coef[1,1]
      b30.2.se <- summary(fit32)$coef[1,2]
      b32.2 <- summary(fit32)$coef[2,1]
      b32.2.se <- summary(fit32)$coef[2,2]
      sigma33.2 <- summary(fit32)$sigma^2

      for (d in 1:1000) # sequence of 12 draws described in paper
      {
         pi0draws[d] <- rbeta(1, n0 + 0.5, n1 + 0.5)
         pi1draws[d] <- 1 - pi0draws[d]
         s22draws[d] <- n0*var(trivar.0[,2]) / rchisq(1, n0-1)
         while (s11.1.draws[d] <= s11.2.0.draws[d])
         {
            s11.1.draws[d] <- n1*var(trivar.1[,1]) / rchisq(1, n1-1)
            s11.2.0.draws[d] <- n0*sigma11.2 / rchisq(1, n0-2)
         }
         b12.2.draws[d] <- rnorm(1, b12.2, sqrt(s11.2.0.draws[d] / (n0*var(trivar.0[,2]))))
         b10.2.draws[d] <- rnorm(1, mean(trivar.0[,1]) - b12.2.draws[d]*mean(trivar.0[,2]), sqrt(s11.2.0.draws[d] / n0))
         mu.2.0.draws[d] <- rnorm(1, mean(trivar.0[,2]), sqrt(s22draws[d] / n0))
         mu.1.1.draws[d] <- rnorm(1, mean(trivar.1[,1]), sqrt(s11.1.draws[d] / n1))
         while (s33.0.draws[d] <= s33.2.0.draws[d])
         {
            s33.0.draws[d] <- n0*var(trivar.0[,3]) / rchisq(1, n0-1)
            s33.2.0.draws[d] <- n0*sigma33.2 / rchisq(1, n0-2)
         }
         b32.2.draws[d] <- rnorm(1, b32.2, sqrt(s33.2.0.draws[d] / (n0*var(trivar.0[,3]))))
         b30.2.draws[d] <- rnorm(1, mean(trivar.0[,3]) - b32.2.draws[d]*mean(trivar.0[,2]), sqrt(s33.2.0.draws[d] / n0))

         mu2est.draws[d] <- pi1draws[d]*((mu.1.1.draws[d] - b10.2.draws[d])/b12.2.draws[d]) + pi0draws[d]*mu.2.0.draws[d]
         mu3est.draws[d] <- b30.2.draws[d] + b32.2.draws[d]*mu2est.draws[d]         
      }
 
      # Compute PMM estimates based on simulated data from posterior distributions
     
      results[s,25] <- median(mu2est.draws) # posterior median of X2 means
      results[s,26] <- median(mu3est.draws) # posterior median of X3 means

      # PMM RMSE Calculations

      results[s,27] <- (results[s,25] - results[s,10])^2 
      results[s,28] <- (results[s,26] - results[s,11])^2 

      # For PMM interval coverage and width

      x2meanPMM.ll <- quantile(mu2est.draws,0.025)
      x2meanPMM.ul <- quantile(mu2est.draws,0.975)
      if (x2meanPMM.ll <= results[s,10] && x2meanPMM.ul > results[s,10]) results[s,29] <- 1
      else results[s,29] <- 0
      results[s,30] <- x2meanPMM.ul - x2meanPMM.ll

      x3meanPMM.ll <- quantile(mu3est.draws,0.025)
      x3meanPMM.ul <- quantile(mu3est.draws,0.975)
      if (x3meanPMM.ll <= results[s,11] && x3meanPMM.ul > results[s,11]) results[s,31] <- 1
      else results[s,31] <- 0
      results[s,32] <- x3meanPMM.ul - x3meanPMM.ll

################ NEW MULTIPLE IMPUTATION ADDITIONS

      # Regress X2 on X1 and intercept for CC

      fit2 <- lm(trivar.0[,2] ~ trivar.0[,1])
      b2 <- summary(fit2)$coef[,1]
      sigma2 <- summary(fit2)$sigma^2
      mat2 <- t(model.matrix(fit2)) %*% model.matrix(fit2)

      # Regress X3 on X1 and intercept for CC

      fit3 <- lm(trivar.0[,3] ~ trivar.0[,1])
      b3 <- summary(fit3)$coef[,1]
      sigma3 <- summary(fit3)$sigma^2
      mat3 <- t(model.matrix(fit3)) %*% model.matrix(fit3)

      # Initialize vectors 
      s2draws <- numeric(100)
      b2draws <- matrix(0,nrow=2,ncol=100)
      s3draws <- numeric(100)
      b3draws <- matrix(0,nrow=2,ncol=100)
      x2meanimpvec <- numeric(100)
      x2varimpvec <- numeric(100)
      x3meanimpvec <- numeric(100)
      x3varimpvec <- numeric(100)

      for (d in 1:100)
      {
         s2draws[d] <- n0*sigma2 / rchisq(1,n0-2)
         b2draws[,d] <- rmnorm(1, b2, solve(mat2)*s2draws[d])
         s3draws[d] <- n0*sigma3 / rchisq(1,n0-2)
         b3draws[,d] <- rmnorm(1, b3, solve(mat3)*s3draws[d]) 
         
	 # impute missing values based on conditional distributions
                       
         trivar.1[,2] <- rnorm(length(trivar.1[,2]), mean = b2draws[1,d] + b2draws[2,d] * trivar.1[,1], 
	  	sd = sqrt(s2draws[d]))

	 trivar.1[,3] <- rnorm(length(trivar.1[,3]), mean = b3draws[1,d] + b3draws[2,d] * trivar.1[,1], 
		sd = sqrt(s3draws[d])) 

         # stack data sets
         trivar.imp <- rbind(trivar.0,trivar.1)
         trivar.imp.data <- data.frame(trivar.imp)    

         x2meanimpvec[d] <- mean(trivar.imp.data[,2])
         x2varimpvec[d] <- var(trivar.imp.data[,2]) / length(trivar.imp.data[,2])
         x3meanimpvec[d] <- mean(trivar.imp.data[,3])
         x3varimpvec[d] <- var(trivar.imp.data[,3]) / length(trivar.imp.data[,3])
      }

      # Compute MI estimates of two means based on m = 100 imputations
     
      results[s,33] <- mean(x2meanimpvec) 
      results[s,34] <- mean(x3meanimpvec) 

      # MI RMSE Calculations

      results[s,35] <- (results[s,33] - results[s,10])^2 
      results[s,36] <- (results[s,34] - results[s,11])^2 

      # For MI confidence interval coverage and width

      mi.var.x2 <- mean(x2varimpvec) + (101/100) * var(x2meanimpvec)
      mi.var.x3 <- mean(x3varimpvec) + (101/100) * var(x3meanimpvec)
      df.x2 <- 99 / ((101/100) * var(x2meanimpvec) / mi.var.x2)^2
      df.x3 <- 99 / ((101/100) * var(x3meanimpvec) / mi.var.x3)^2
      results[s,41] <- ((101/100) * var(x2meanimpvec)) / ((101/100) * var(x2meanimpvec) + mean(x2varimpvec))
      results[s,42] <- ((101/100) * var(x3meanimpvec)) / ((101/100) * var(x3meanimpvec) + mean(x3varimpvec))

      x2mean.ll <- mean(x2meanimpvec) - qt(0.975,df.x2)*sqrt(mi.var.x2)
      x2mean.ul <- mean(x2meanimpvec) + qt(0.975,df.x2)*sqrt(mi.var.x2)
      x3mean.ll <- mean(x3meanimpvec) - qt(0.975,df.x3)*sqrt(mi.var.x3)
      x3mean.ul <- mean(x3meanimpvec) + qt(0.975,df.x3)*sqrt(mi.var.x3)

      if (x2mean.ll <= results[s,10] && x2mean.ul > results[s,10]) results[s,37] <- 1
      else results[s,37] <- 0
      results[s,38] <- x2mean.ul - x2mean.ll

      if (x3mean.ll <= results[s,11] && x3mean.ul > results[s,11]) results[s,39] <- 1
      else results[s,39] <- 0
      results[s,40] <- x3mean.ul - x3mean.ll
      cat("Sample:"," ",s,"\n")

   } # end of for loop for simulations

   cat("Mean response rate is:"," ",mean(results[,4]),"\n")

   cat("True mean of X2 is:"," ",mean(results[,10]),"\n")
   cat("\n")

   cat("Mean of CC estimates of X2 mean for respondents is:"," ",mean(results[,2]),"\n")
   cat("RelBias of CC estimates of X2 mean is:"," ",((mean(results[,2]) - mean(results[,10]))/mean(results[,10]))*100,"\n")
   cat("RMSE of CC estimates of X2 mean is:"," ",sqrt(mean(results[,5])),"\n")
   cat("95% CI coverage for CC estimate of X2 mean is:"," ",mean(results[,7]),"\n")
   cat("95% CI width for CC estimate of X2 mean is:"," ",mean(results[,12]),"\n")
   cat("\n")

   cat("Mean of GW estimates of X2 mean for respondents is:"," ",mean(results[,17]),"\n")
   cat("RelBias of GW estimates of X2 mean is:"," ",((mean(results[,17]) - mean(results[,10]))/mean(results[,10]))*100,"\n")
   cat("RMSE of GW estimates of X2 mean is:"," ",sqrt(mean(results[,19])),"\n")
   cat("95% CI coverage for GW estimate of X2 mean is:"," ",mean(results[,21]),"\n")
   cat("95% CI width for GW estimate of X2 mean is:"," ",mean(results[,22]),"\n")
   cat("\n")

   cat("Mean of PMM estimates of X2 mean for respondents is:"," ",mean(results[,25]),"\n")
   cat("RelBias of PMM estimates of X2 mean is:"," ",((mean(results[,25]) - mean(results[,10]))/mean(results[,10]))*100,"\n")
   cat("RMSE of PMM estimates of X2 mean is:"," ",sqrt(mean(results[,27])),"\n")
   cat("95% CI coverage for PMM estimate of X2 mean is:"," ",mean(results[,29]),"\n")
   cat("95% CI width for PMM estimate of X2 mean is:"," ",mean(results[,30]),"\n")
   cat("\n")

   cat("Mean of MI estimates of X2 mean for respondents is:"," ",mean(results[,33]),"\n")
   cat("RelBias of MI estimates of X2 mean is:"," ",((mean(results[,33]) - mean(results[,10]))/mean(results[,10]))*100,"\n")
   cat("RMSE of MI estimates of X2 mean is:"," ",sqrt(mean(results[,35])),"\n")
   cat("Average FMI for X2 mean is:"," ",mean(results[,41]),"\n")
   cat("95% CI coverage for MI estimate of X2 mean is:"," ",mean(results[,37]),"\n")
   cat("95% CI width for MI estimate of X2 mean is:"," ",mean(results[,38]),"\n")
   cat("\n")

   cat("True mean of X3 is:"," ",mean(results[,11]),"\n")
   cat("\n")

   cat("Mean of CC estimates of X3 mean for respondents is:"," ",mean(results[,3]),"\n")
   cat("RelBias of CC estimates of X3 mean is:"," ",((mean(results[,3]) - mean(results[,11]))/mean(results[,11]))*100,"\n")
   cat("RMSE of CC estimates of X3 mean is:"," ",sqrt(mean(results[,6])),"\n")
   cat("95% CI coverage for CC estimate of X3 mean is:"," ",mean(results[,8]),"\n")
   cat("95% CI width for CC estimate of X3 mean is:"," ",mean(results[,13]),"\n")
   cat("\n")

   cat("Mean of GW estimates of X3 mean for respondents is:"," ",mean(results[,18]),"\n")
   cat("RelBias of GW estimates of X3 mean is:"," ",((mean(results[,18]) - mean(results[,11]))/mean(results[,11]))*100,"\n")
   cat("RMSE of GW estimates of X3 mean is:"," ",sqrt(mean(results[,20])),"\n")
   cat("95% CI coverage for GW estimate of X3 mean is:"," ",mean(results[,23]),"\n")
   cat("95% CI width for GW estimate of X3 mean is:"," ",mean(results[,24]),"\n")
   cat("\n")

   cat("Mean of PMM estimates of X3 mean for respondents is:"," ",mean(results[,26]),"\n")
   cat("RelBias of PMM estimates of X3 mean is:"," ",((mean(results[,26]) - mean(results[,11]))/mean(results[,11]))*100,"\n")
   cat("RMSE of PMM estimates of X3 mean is:"," ",sqrt(mean(results[,28])),"\n")
   cat("95% CI coverage for PMM estimate of X3 mean is:"," ",mean(results[,31]),"\n")
   cat("95% CI width for PMM estimate of X3 mean is:"," ",mean(results[,32]),"\n")
   cat("\n")

   cat("Mean of MI estimates of X3 mean for respondents is:"," ",mean(results[,34]),"\n")
   cat("RelBias of MI estimates of X3 mean is:"," ",((mean(results[,34]) - mean(results[,11]))/mean(results[,11]))*100,"\n")
   cat("RMSE of MI estimates of X3 mean is:"," ",sqrt(mean(results[,36])),"\n")
   cat("Average FMI for X3 mean is:"," ",mean(results[,42]),"\n")
   cat("95% CI coverage for MI estimate of X3 mean is:"," ",mean(results[,39]),"\n")
   cat("95% CI width for MI estimate of X3 mean is:"," ",mean(results[,40]),"\n")
   cat("\n")

   cat("Mean of Marginal means of X1 in full data set simulated from PMM is:"," ",mean(results[,9]),"\n")
   cat("Mean of Marginal means of X2 in full data set simulated from PMM is:"," ",mean(results[,10]),"\n")
   cat("Mean of Marginal means of X3 in full data set simulated from PMM is:"," ",mean(results[,11]),"\n")   
   cat("\n")

   return(results)
}

tres <- trivar.sim.pmm(nsamps=1000,sampsize=1000,pi1=0.50,rho=0.9)
tres <- trivar.sim.pmm(nsamps=1000,sampsize=1000,pi1=0.25,rho=0.9)

# CHANGE PARAMETERS
tres <- trivar.sim.pmm(nsamps=1000,sampsize=1000,pi1=0.50,rho=0.6)
tres <- trivar.sim.pmm(nsamps=1000,sampsize=1000,pi1=0.25,rho=0.6)




# X2 Bias, High RR

plot(c(2,2,2,1,1,1,0.5,0.5,0.5,0,0,0,2,2,2,1,1,1,0.5,0.5,0.5,0,0,0),
c(-0.01,8.86,35.14,0.13,5.10,22.44,0.12,3.17,14.59,0.07,0.03,-0.15,-0.05,18.79,26.09,-0.27,15.82,22.95,-0.21,7.85,11.54,-0.36,-0.35,-0.35),
xlab=expression(beta),ylab=expression(paste("Relative Empirical Bias of   ",hat(mu)[2],"  (%)")),pch=".",main=expression(paste("Relative Empirical Bias Comparison,  ",hat(mu)[2])))

lines(c(0,0.5,1,2),c(0.07,0.12,0.13,-0.01),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(-0.36,-0.21,-0.27,-0.05),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.03,3.17,5.10,8.86),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(-0.35,7.85,15.82,18.79),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(-0.15,14.59,22.44,35.14),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(-0.35,11.54,22.95,26.09),lty=2,col="blue") #CC 0.6
legend('topleft',c("PMM-High ME","PMM-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"),lty=c(2,1,2,1,2,1),col=c("black","black","red","red","blue","blue"))

# X3 Bias, High RR

plot(c(2,2,2,1,1,1,0.5,0.5,0.5,0,0,0,2,2,2,1,1,1,0.5,0.5,0.5,0,0,0),
c(0.01,1.27,2.11,0.02,0.74,1.18,0.01,0.43,0.78,0.01,-0.01,-0.01,-0.03,0.98,1.25,-0.01,0.88,1.16,0.01,0.46,0.62,-0.01,-0.01,-0.01),
xlab=expression(beta),ylab=expression(paste("Relative Empirical Bias of   ",hat(mu)[3],"  (%)")),pch=".",main=expression(paste("Relative Empirical Bias Comparison,  ",hat(mu)[3])))
lines(c(0,0.5,1,2),c(0.01,0.01,0.02,0.01),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(-0.01,0.01,-0.01,-0.03),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(-0.01,0.43,0.74,1.27),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(-0.01,0.46,0.88,0.98),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(-0.01,0.78,1.18,2.11),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(-0.01,0.62,1.16,1.25),lty=2,col="blue") #CC 0.6
legend('topleft',c("PMM-High ME","PMM-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"),lty=c(2,1,2,1,2,1),col=c("black","black","red","red","blue","blue"))


