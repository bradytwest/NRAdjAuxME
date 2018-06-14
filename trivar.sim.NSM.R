# Other general simulation studies from West and Little (2013), based on a normal selection model

trivar.sim.nsm <- function(nsamps,sampsize,rrate,beta,rho)
{
   require(mnormt)
   require(survey)
   require(mi)

   # Step 1 Model
   vcov <- cbind(c(1,rho,0.25),c(rho,1,0.5),c(0.25,0.5,1))
   meanvec <- c(1,1,10)
   
   results <- matrix(0,nrow=nsamps,ncol=42)
   set.seed(41279)

   for (s in 1:nsamps)
   {

      # simulate sample from normal selection model

      trivar <- rmnorm(sampsize, meanvec, vcov)

      resp <- numeric(sampsize)
      pi <- plogis(qlogis(rrate)+beta*trivar[,2]) 
      uni <- runif(sampsize)
      for (j in 1:sampsize)
      {
         if (uni[j] <= pi[j]) resp[j] <- 1
         else resp[j] <- 0
      }
    
      n0 <- sum(resp)
      n1 <- sampsize - sum(resp)
      trivar.0 <- trivar[(resp == 1),]
      trivar.1 <- trivar[(resp == 0),]

      # prep for logistic regression (predict probability of RESPONSE with X1 for GW approach)

      resplogit <- glm(resp ~ trivar[,1], family=binomial(link="logit"))
      nrweights <- 1 / resplogit$fit[resp == 1]

      # prep for MI procedure

      trivar.stack2 <- data.frame(trivar)
      trivar.stack2$X2[resp == 0] <- NA
      trivar.stack2$X3[resp == 0] <- NA 

      # Compute CC estimates based on simulated data 

      results[s,1] <- mean(trivar.0[,1]) # R mean of X1
      results[s,2] <- mean(trivar.0[,2]) # R mean of X2
      results[s,3] <- mean(trivar.0[,3]) # R mean of X3
      results[s,4] <- n0/sampsize # response rate

      # CC RMSE Calculations

      results[s,5] <- (results[s,2] - 1)^2 # For RMSE R mean X2
      results[s,6] <- (results[s,3] - 10)^2 # For RMSE R mean X3

      # For CC interval coverage and width

      x2mean.ll <- results[s,2] - qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)
      x2mean.ul <- results[s,2] + qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)
      if (x2mean.ll <= 1 && x2mean.ul > 1) results[s,7] <- 1
      else results[s,7] <- 0
      results[s,12] <- 2*qt(0.975,n0-1)*sd(trivar.0[,2])/sqrt(n0)

      x3mean.ll <- results[s,3] - qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)
      x3mean.ul <- results[s,3] + qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)
      if (x3mean.ll <= 10 && x3mean.ul > 10) results[s,8] <- 1
      else results[s,8] <- 0
      results[s,13] <- 2*qt(0.975,n0-1)*sd(trivar.0[,3])/sqrt(n0)

      # Compute GW estimates based on simulated data from response pattern
     
      gwdes <- svydesign(id=~1,weights=~nrweights,data=data.frame(trivar.0))
      results[s,17] <- coef(svymean(~X2,gwdes)) # Weighted R mean of X2
      results[s,18] <- coef(svymean(~X3,gwdes)) # Weighted R mean of X3

      # GW RMSE Calculations

      results[s,19] <- (results[s,17] - 1)^2 # For RMSE Weighted R mean X2
      results[s,20] <- (results[s,18] - 10)^2 # For RMSE Weighted R mean X3

      # For GW interval coverage and width

      x2meanGW.ll <- confint(svymean(~X2,gwdes))[1]
      x2meanGW.ul <- confint(svymean(~X2,gwdes))[2]
      if (x2meanGW.ll <= 1 && x2meanGW.ul > 1) results[s,21] <- 1
      else results[s,21] <- 0
      results[s,22] <- confint(svymean(~X2,gwdes))[2] - confint(svymean(~X2,gwdes))[1]

      x3meanGW.ll <- confint(svymean(~X3,gwdes))[1]
      x3meanGW.ul <- confint(svymean(~X3,gwdes))[2]
      if (x3meanGW.ll <= 10 && x3meanGW.ul > 10) results[s,23] <- 1
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

      results[s,27] <- (results[s,25] - 1)^2 
      results[s,28] <- (results[s,26] - 10)^2 

      # For PMM interval coverage and width

      x2meanPMM.ll <- quantile(mu2est.draws,0.025)
      x2meanPMM.ul <- quantile(mu2est.draws,0.975)
      if (x2meanPMM.ll <= 1 && x2meanPMM.ul > 1) results[s,29] <- 1
      else results[s,29] <- 0
      results[s,30] <- x2meanPMM.ul - x2meanPMM.ll

      x3meanPMM.ll <- quantile(mu3est.draws,0.025)
      x3meanPMM.ul <- quantile(mu3est.draws,0.975)
      if (x3meanPMM.ll <= 10 && x3meanPMM.ul > 10) results[s,31] <- 1
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


      # Use mi() procedure to perform MI and compute MI estimates

      #info <- mi.info(trivar.stack2)
      #imp <- mi(trivar.stack2, info, run.past.converge=TRUE, check.coef.convergence=TRUE, n.imp=20)
      #imp <- mi(trivar.stack2, info, n.imp=20, n.iter=6, add.noise=noise.control(post.run.iter=0))

      # allow MI procedure to converge
      #while(converged(imp,check="data") == FALSE)
      #{
         #imp <- mi(imp, run.past.convergence = TRUE, n.iter = 20)
      #}
     
      # Compute MI estimates based on simulated data from two patterns

      #x2meanimpvec <- numeric(20)
      #x2varimpvec <- numeric(20)
      #x3meanimpvec <- numeric(20)
      #x3varimpvec <- numeric(20)
      #for (m in 1:20)
      #{
         #impvals <- mi.data.frame(imp,m)
         #x2meanimpvec[m] <- mean(impvals$X2)
         #x2varimpvec[m] <- var(impvals$X2) / length(impvals$X2)
         #x3meanimpvec[m] <- mean(impvals$X3)
         #x3varimpvec[m] <- var(impvals$X3) / length(impvals$X3)
      #}

      # Compute MI estimates of two means based on m = 100 imputations
     
      results[s,33] <- mean(x2meanimpvec) 
      results[s,34] <- mean(x3meanimpvec) 

      # MI RMSE Calculations

      results[s,35] <- (results[s,33] - 1)^2 
      results[s,36] <- (results[s,34] - 10)^2 

      # For MI confidence interval coverage and width

      mi.var.x2 <- mean(x2varimpvec) + (101/100) * var(x2meanimpvec)
      mi.var.x3 <- mean(x3varimpvec) + (101/100) * var(x3meanimpvec)
      results[s,41] <- ((101/100) * var(x2meanimpvec)) / ((101/100) * var(x2meanimpvec) + mean(x2varimpvec))
      results[s,42] <- ((101/100) * var(x3meanimpvec)) / ((101/100) * var(x3meanimpvec) + mean(x3varimpvec))
      df.x2 <- 99 / ((101/100) * var(x2meanimpvec) / mi.var.x2)^2
      df.x3 <- 99 / ((101/100) * var(x3meanimpvec) / mi.var.x3)^2

      x2mean.ll <- mean(x2meanimpvec) - qt(0.975,df.x2)*sqrt(mi.var.x2)
      x2mean.ul <- mean(x2meanimpvec) + qt(0.975,df.x2)*sqrt(mi.var.x2)
      x3mean.ll <- mean(x3meanimpvec) - qt(0.975,df.x3)*sqrt(mi.var.x3)
      x3mean.ul <- mean(x3meanimpvec) + qt(0.975,df.x3)*sqrt(mi.var.x3)

      if (x2mean.ll <= 1 && x2mean.ul > 1) results[s,37] <- 1
      else results[s,37] <- 0
      results[s,38] <- x2mean.ul - x2mean.ll

      if (x3mean.ll <= 10 && x3mean.ul > 10) results[s,39] <- 1
      else results[s,39] <- 0
      results[s,40] <- x3mean.ul - x3mean.ll
      cat("Sample:"," ",s,"\n")

   } # end of for loop for simulations

   cat("Mean response rate is:"," ",mean(results[,4]),"\n")

   cat("True mean of X2 is:"," ",1,"\n")
   cat("\n")

   cat("Mean of CC estimates of X2 mean for respondents is:"," ",mean(results[,2]),"\n")
   cat("RelBias of CC estimates of X2 mean is:"," ",((mean(results[,2]) - 1)/1)*100,"\n")
   cat("RMSE of CC estimates of X2 mean is:"," ",sqrt(mean(results[,5])),"\n")
   cat("95% CI coverage for CC estimate of X2 mean is:"," ",mean(results[,7]),"\n")
   cat("95% CI width for CC estimate of X2 mean is:"," ",mean(results[,12]),"\n")
   cat("\n")

   cat("Mean of GW estimates of X2 mean for respondents is:"," ",mean(results[,17]),"\n")
   cat("RelBias of GW estimates of X2 mean is:"," ",((mean(results[,17]) - 1)/1)*100,"\n")
   cat("RMSE of GW estimates of X2 mean is:"," ",sqrt(mean(results[,19])),"\n")
   cat("95% CI coverage for GW estimate of X2 mean is:"," ",mean(results[,21]),"\n")
   cat("95% CI width for GW estimate of X2 mean is:"," ",mean(results[,22]),"\n")
   cat("\n")

   cat("Mean of PMM estimates of X2 mean for respondents is:"," ",mean(results[,25]),"\n")
   cat("RelBias of PMM estimates of X2 mean is:"," ",((mean(results[,25]) - 1)/1)*100,"\n")
   cat("RMSE of PMM estimates of X2 mean is:"," ",sqrt(mean(results[,27])),"\n")
   cat("95% CI coverage for PMM estimate of X2 mean is:"," ",mean(results[,29]),"\n")
   cat("95% CI width for PMM estimate of X2 mean is:"," ",mean(results[,30]),"\n")
   cat("\n")

   cat("Mean of MI estimates of X2 mean for respondents is:"," ",mean(results[,33]),"\n")
   cat("RelBias of MI estimates of X2 mean is:"," ",((mean(results[,33]) - 1)/1)*100,"\n")
   cat("RMSE of MI estimates of X2 mean is:"," ",sqrt(mean(results[,35])),"\n")
   cat("Average FMI for X2 mean is:"," ",mean(results[,41]),"\n")
   cat("95% CI coverage for MI estimate of X2 mean is:"," ",mean(results[,37]),"\n")
   cat("95% CI width for MI estimate of X2 mean is:"," ",mean(results[,38]),"\n")
   cat("\n")

   cat("True mean of X3 is:"," ",10,"\n")
   cat("\n")

   cat("Mean of CC estimates of X3 mean for respondents is:"," ",mean(results[,3]),"\n")
   cat("RelBias of CC estimates of X3 mean is:"," ",((mean(results[,3]) - 10)/10)*100,"\n")
   cat("RMSE of CC estimates of X3 mean is:"," ",sqrt(mean(results[,6])),"\n")
   cat("95% CI coverage for CC estimate of X3 mean is:"," ",mean(results[,8]),"\n")
   cat("95% CI width for CC estimate of X3 mean is:"," ",mean(results[,13]),"\n")
   cat("\n")

   cat("Mean of GW estimates of X3 mean for respondents is:"," ",mean(results[,18]),"\n")
   cat("RelBias of GW estimates of X3 mean is:"," ",((mean(results[,18]) - 10)/10)*100,"\n")
   cat("RMSE of GW estimates of X3 mean is:"," ",sqrt(mean(results[,20])),"\n")
   cat("95% CI coverage for GW estimate of X3 mean is:"," ",mean(results[,23]),"\n")
   cat("95% CI width for GW estimate of X3 mean is:"," ",mean(results[,24]),"\n")
   cat("\n")

   cat("Mean of PMM estimates of X3 mean for respondents is:"," ",mean(results[,26]),"\n")
   cat("RelBias of PMM estimates of X3 mean is:"," ",((mean(results[,26]) - 10)/10)*100,"\n")
   cat("RMSE of PMM estimates of X3 mean is:"," ",sqrt(mean(results[,28])),"\n")
   cat("95% CI coverage for PMM estimate of X3 mean is:"," ",mean(results[,31]),"\n")
   cat("95% CI width for PMM estimate of X3 mean is:"," ",mean(results[,32]),"\n")
   cat("\n")

   cat("Mean of MI estimates of X3 mean for respondents is:"," ",mean(results[,34]),"\n")
   cat("RelBias of MI estimates of X3 mean is:"," ",((mean(results[,34]) - 10)/10)*100,"\n")
   cat("RMSE of MI estimates of X3 mean is:"," ",sqrt(mean(results[,36])),"\n")
   cat("Average FMI for X3 mean is:"," ",mean(results[,42]),"\n")
   cat("95% CI coverage for MI estimate of X3 mean is:"," ",mean(results[,39]),"\n")
   cat("95% CI width for MI estimate of X3 mean is:"," ",mean(results[,40]),"\n")
   cat("\n")

   return(results)
}

tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=2,rho=0.9)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=1,rho=0.9)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=0,rho=0.9)

tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=2,rho=0.6)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=1,rho=0.6)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.5,beta=0,rho=0.6)

tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=2,rho=0.9)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=1,rho=0.9)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=0,rho=0.9)

tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=2,rho=0.6)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=1,rho=0.6)
tres <- trivar.sim.nsm(nsamps=1000,sampsize=1000,rrate=0.268945,beta=0,rho=0.6)


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

# X2 RMSE, High RR

plot(c(2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0,2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0),
c(0.035,0.090,0.096,0.340,0.032,0.060,0.059,0.223,0.033,0.044,0.045,0.149,0.035,0.035,0.034,0.045,0.040,0.186,0.187,0.258,0.044,0.159,0.158,0.226,0.048,0.089,0.087,0.122,0.053,0.042,0.040,0.044),
xlab=expression(beta),ylab=expression(paste("RMSE of   ",hat(mu)[2])),pch=".",main=expression(paste("RMSE Comparison,  ",hat(mu)[2])))

lines(c(0,0.5,1,2),c(0.035,0.033,0.032,0.035),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.053,0.048,0.044,0.040),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.035,0.044,0.060,0.090),lty=1,col="green") #MI 0.9
lines(c(0,0.5,1,2),c(0.042,0.089,0.159,0.186),lty=2,col="green") #MI 0.6
lines(c(0,0.5,1,2),c(0.034,0.045,0.059,0.096),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.040,0.087,0.158,0.187),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.045,0.149,0.223,0.340),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.044,0.122,0.226,0.258),lty=2,col="blue") #CC 0.6
legend('topleft',c("PMM-High ME","PMM-Low ME","MI-High ME","MI-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"), lty=c(2,1,2,1,2,1,2,1), col=c("black","black","green","green","red","red","blue","blue"), cex=0.8)

# X3 RMSE, High RR

plot(c(2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0,2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0),
c(0.041,0.137,0.134,0.213,0.037,0.084,0.082,0.122,0.040,0.060,0.058,0.087,0.045,0.045,0.044,0.045,0.039,0.105,0.103,0.129,0.040,0.095,0.095,0.122,0.044,0.060,0.061,0.074,0.048,0.047,0.043,0.044),
xlab=expression(beta),ylab=expression(paste("RMSE of   ",hat(mu)[3])),pch=".",main=expression(paste("RMSE Comparison,  ",hat(mu)[3])))

lines(c(0,0.5,1,2),c(0.045,0.040,0.037,0.041),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.048,0.044,0.040,0.039),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.045,0.060,0.084,0.137),lty=1,col="green") #MI 0.9
lines(c(0,0.5,1,2),c(0.047,0.060,0.095,0.103),lty=2,col="green") #MI 0.6
lines(c(0,0.5,1,2),c(0.044,0.058,0.082,0.134),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.043,0.061,0.095,0.103),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.045,0.087,0.122,0.213),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.044,0.074,0.122,0.129),lty=2,col="blue") #CC 0.6
legend('topleft',c("PMM-High ME","PMM-Low ME","MI-High ME","MI-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"), lty=c(2,1,2,1,2,1,2,1),col=c("black","black","green","green","red","red","blue","blue"),cex=0.8)

# X2 Coverage, High RR

plot(c(2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0,2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0),
c(0.947,0.592,0.591,0.000,0.951,0.642,0.794,0.000,0.954,0.863,0.920,0.041,0.954,0.950,0.990,0.945,0.958,0.001,0.000,0.000,0.947,0.010,0.006,0.000,0.950,0.451,0.480,0.163,0.946,0.949,0.973,0.957),
xlab=expression(beta),ylab=expression(paste("Coverage of 95% CI for   ",mu[2])),pch=".",main=expression(paste("95% CI Coverage Comparison, ", mu[2])))

lines(c(0,0.5,1,2),c(0.954,0.954,0.951,0.947),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.946,0.950,0.947,0.958),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.950,0.863,0.642,0.292),lty=1,col="green") #MI 0.9
lines(c(0,0.5,1,2),c(0.949,0.451,0.010,0.001),lty=2,col="green") #MI 0.6
lines(c(0,0.5,1,2),c(0.990,0.920,0.794,0.591),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.973,0.480,0.006,0.000),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.945,0.041,0.000,0.000),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.957,0.163,0.000,0.000),lty=2,col="blue") #CC 0.6
legend(locator(1),c("PMM-High ME","PMM-Low ME","MI-High ME","MI-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"), lty=c(2,1,2,1,2,1,2,1),col=c("black","black","green","green","red","red","blue","blue"),cex=0.8)

# X3 Coverage, High RR

plot(c(2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0,2,2,2,2,1,1,1,1,0.5,0.5,0.5,0.5,0,0,0,0),
c(0.939,0.143,0.183,0.000,0.937,0.535,0.510,0.094,0.940,0.819,0.798,0.512,0.929,0.956,0.954,0.950,0.952,0.225,0.211,0.041,0.947,0.390,0.344,0.118,0.937,0.802,0.794,0.647,0.949,0.944,0.952,0.945),
xlab=expression(beta),ylab=expression(paste("Coverage of 95% CI for   ",mu[3])),pch=".",main=expression(paste("95% CI Coverage Comparison, ",mu[3])))

lines(c(0,0.5,1,2),c(0.929,0.940,0.937,0.939),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.949,0.937,0.947,0.952),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.956,0.819,0.535,0.143),lty=1,col="green") #MI 0.9
lines(c(0,0.5,1,2),c(0.944,0.802,0.390,0.225),lty=2,col="green") #MI 0.6
lines(c(0,0.5,1,2),c(0.954,0.798,0.510,0.183),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.952,0.794,0.344,0.211),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.950,0.512,0.094,0.000),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.945,0.647,0.118,0.041),lty=2,col="blue") #CC 0.6
legend(locator(1),c("PMM-High ME","PMM-Low ME","MI-High ME","MI-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"), lty=c(2,1,2,1,2,1,2,1),col=c("black","black","green","green","red","red","blue","blue"),cex=0.8)

# X2 Width, High RR

plot(c(2,2,2,1,1,1,0.5,0.5,0.5,0,0,0,2,2,2,1,1,1,0.5,0.5,0.5,0,0,0),
c(0.140,0.207,0.127,0.127,0.154,0.134,0.131,0.158,0.152,0.138,0.177,0.177,0.162,0.125,0.119,0.175,0.140,0.136,0.190,0.154,0.153,0.211,0.176,0.177),
xlab=expression(beta),ylab=expression(paste("Mean Width of 95% CIs for   ",mu[2])),pch=".",main=expression(paste("95% CI Mean Width Comparison,  ",mu[2])))

lines(c(0,0.5,1,2),c(0.138,0.131,0.127,0.140),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.211,0.190,0.175,0.162),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.177,0.158,0.154,0.207),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.176,0.154,0.140,0.125),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.177,0.152,0.134,0.127),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.177,0.153,0.136,0.119),lty=2,col="blue") #CC 0.6
legend(locator(1),c("PMM-High ME","PMM-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"),lty=c(2,1,2,1,2,1),col=c("black","black","red","red","blue","blue"),cex=0.8)

# X3 Width, High RR

plot(c(2,2,2,1,1,1,0.5,0.5,0.5,0,0,0,2,2,2,1,1,1,0.5,0.5,0.5,0,0,0),
c(0.160,0.175,0.150,0.144,0.148,0.143,0.151,0.157,0.155,0.168,0.176,0.177,0.147,0.136,0.133,0.158,0.146,0.144,0.166,0.154,0.154,0.186,0.176,0.177),
xlab=expression(beta),ylab=expression(paste("Mean Width of 95% CIs for   ",mu[3])),pch=".",main=expression(paste("95% CI Mean Width Comparison,  ",mu[3])))

lines(c(0,0.5,1,2),c(0.168,0.151,0.144,0.160),lty=1,col="black") #PMM 0.9
lines(c(0,0.5,1,2),c(0.186,0.166,0.158,0.147),lty=2,col="black") #PMM 0.6
lines(c(0,0.5,1,2),c(0.176,0.157,0.148,0.175),lty=1,col="red") #GW 0.9
lines(c(0,0.5,1,2),c(0.176,0.154,0.146,0.136),lty=2,col="red") #GW 0.6
lines(c(0,0.5,1,2),c(0.177,0.155,0.143,0.150),lty=1,col="blue") #CC 0.9
lines(c(0,0.5,1,2),c(0.177,0.154,0.144,0.133),lty=2,col="blue") #CC 0.6
legend(locator(1),c("PMM-High ME","PMM-Low ME","GW-High ME","GW-Low ME","CC-High ME","CC-Low ME"),lty=c(2,1,2,1,2,1),col=c("black","black","red","red","blue","blue"))


