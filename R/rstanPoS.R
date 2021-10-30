#'@export find_weibull
#'@export rdsamples_generator
#'@export data_looktime
#'@export prediction_samples
#'@export data_studyend
#'@export interim
#'@export secondstage_samples
#'@export data_studyend2
#'@export final
#'@import doParallel
#'@import parallel
#'@import methods

find_weibull <- function(PFS_1_time, PFS_1_probability, PFS_2_time,
                         PFS_2_probability) {
  x1=PFS_1_time; x2=PFS_2_time; p1=1-PFS_1_probability; p2=1-PFS_2_probability
  shape <- (log(-log(1-p2)) - log(-log(1-p1))) / (log(x2)-log(x1))
  scale <- x1 / (-log(1-p1))^(1/shape)
  median_survival <- scale*log(2)^(1/shape)
  output <- list()
  output$Shape <- shape
  output$Scale <- scale
  output$Median <- median_survival
  output$PFS1 <- 1-pweibull(q=PFS_1_time,shape=shape,scale=scale)
  output$PFS2 <- 1-pweibull(q=PFS_2_time,shape=shape,scale=scale)
  output$PFS1time <- PFS_1_time
  output$PFS2time <- PFS_2_time

  class(output) <- "find_weibull"
  print(output)
}

#'@export
plot.find_weibull<- function(x, ...) {

  par(mfrow=c(1,1))
  x_seq <- seq(0,24, by = 0.1)
  dw <- 1- pweibull(x_seq,x$Shape, x$Scale)
  plot(x_seq,dw,xlab = "Time in Months",
       ylab = "Survival Probability", main = paste("Weibull Survival Curve with\nshape =", round(x$Shape,2),"and scale =", round(x$Scale,2)), type = "l",...)
}

a <- find_weibull(PFS_1_time = 9, PFS_1_probability = 0.127, PFS_2_time = 4,
                  PFS_2_probability = 0.4)
plot(a)
#Generate the data
rdsamples_generator <- function(n1,tau,k,omega1,kappa,theta){

  #generate the accrual time
  e <- runif(n1,0,omega1)
  #generate time from e till event
  y <- rweibull(n1,kappa,theta)
  # e + y time of event
  y_e <- e + y
  #generate time from e till administrative censoring
  z <- rep(max(e) + k,n1)

  #generate indicator variable for event = 1
  delta <- (y_e<=z)*1

  #combine all the generated data into one data frame
  dataSurv <- data.frame(e,y,y_e,z,delta)

  return(dataSurv)
}

#A FUNCTION IS NEEDED WHICH TRANSFORMS THE NO EVENTS INTO EVENT TIMES OF 1000!!!

#From the generated data only take everything up to k
data_looktime <- function(dataSurv,k,till,to,to2) {
  #RoundUp Function - only required if data is simulated
  roundUp <- function(x,till,to,to2)
  {
    if(x<till) {
      to*(x%/%to + as.logical(x%%to))}
    else(to2*(x%/%to2 + as.logical(x%%to2)))
  }
  interim <- max(dataSurv$e)+k
  output <- list()
  #active people (no event, no censoring, no loss to follow up at looktime)
  output$m1 <- nrow(dataSurv[dataSurv$y_e>interim,])
  #events out of n1-m1 patients
  output$n1event <- as.numeric(
    nrow(dataSurv[dataSurv$y_e<=interim,]))
  #Survival times up to looktime
  output$y_looktime <- dataSurv[dataSurv$y_e<=interim,]$y
  for(i in 1:length(output$y_looktime)) {
    output$y_looktime_round[i] <- roundUp(output$y_looktime[i],till,to,to2) }  #Censoring time b/c no event until looktime
  output$z_looktime <- interim - dataSurv[dataSurv$y_e>interim,]$e
  if(output$m1 > 0) {
    for(i in 1:length(output$z_looktime)) {
      output$z_looktime_round[i] <- roundUp(output$z_looktime[i],till,to,to2)}}
  else{output$z_looktime_round <- output$z_looktime}
  return(output)
}


#for each gamma j we will predict survival data at study end
prediction_samples <- function(dataSurv,n,n1,m1,tau,omega2,kappa,theta,k) {
  interim <- max(dataSurv$e)+k
  #generate the accrual time
  e1 <- dataSurv[dataSurv$y_e>interim,]$e
  e2 <- runif(n-n1,max(e1) + k,max(e1) + k + omega2)
  e <- c(e1,e2)

  #generate time from e1 till event consid. the time spent in trial already
  y1 <- rweibull(m1,kappa,theta) + (interim - dataSurv[dataSurv$y_e>interim,]$e)
  #generate time from e2 till event
  y2 <- rweibull(n-n1,kappa,theta)
  y <- c(y1,y2)

  #e1 + y1 time of event
  y_e1 <- e1 + y1
  #e2 + y2 time of event
  y_e2 <- e2 + y2
  y_e <- c(y_e1,y_e2)


  #generate time from e till administrative censoring
  z <- rep(max(e)+tau,n-n1+m1)

  #generate indicator variable for event = 1
  delta <- (y_e<=z)*1

  #combine all the generated data into one data frame
  dataSurv_predicted <- data.frame(e,y,y_e,z,delta)

  return(dataSurv_predicted)
}

#Generate the predicted survival data at studyend
data_studyend <- function(dataSurv,prediction,tau,k,till,to,to2) {
  #RoundUp Function - only required if data is simulated
  roundUp <- function(x,till,to,to2)
  {
    if(x<till) {
      to*(x%/%to + as.logical(x%%to))}
    else(to2*(x%/%to2 + as.logical(x%%to2)))
  }
  final <- max(prediction$e)+tau
  interim <- max(dataSurv$e)+k
  output <- list()
  #events out of m = m1 + n-n1 patients (m1 patients from first stage and second stage patients)
  output$mevent <- as.numeric(
    nrow(prediction[prediction$y_e<=final,]))
  #aC out of m patients
  output$ncens <- as.numeric(
    nrow(prediction[prediction$y_e>final,]))
  #events out of n patients
  output$nevent <- output$mevent + as.numeric(
    nrow(dataSurv[dataSurv$y_e<=interim,]))
  #Survival time
  output$y_studyend <- c(dataSurv[dataSurv$y_e<=interim,]$y,prediction[prediction$y_e<=final,]$y)
  for(i in 1:length(output$y_studyend)) {
    output$y_studyend_round[i] <- roundUp(output$y_studyend[i],till,to,to2) }  #Censoring time
  output$z_studyend <- final - prediction[prediction$y_e>final,]$e
  if(output$ncens > 0) {
    for(i in 1:length(output$z_studyend)) {
      output$z_studyend_round[i] <- roundUp(output$z_studyend[i],till,to,to2) }}
  else{output$z_studyend_round <- output$z_studyend}
  return(output)
}

interim <- function(n.cores,rddata,discrete,data,n,n1,p0,p2_0,k,tau,omega1,omega2,kappa,theta,d,eta_star,lambda,mu,sigma,a,b,till,to,to2) {
  require(doParallel)
  require(parallel)
  PoS <- c()
  decision <- c()

  if(rddata == 1) {
    #Data generated
    dataSurv <- rdsamples_generator(n1,tau,k,omega1,kappa,theta)
  }
  else{dataSurv <- data
  kappa <- find_weibull(tau, p0, k,
                        p2_0)$Shape
  theta <- find_weibull(tau, p0, k,
                        p2_0)$Scale
  }

  #Generate data only up to looktime
  output_looktime <- data_looktime(dataSurv,k,till,to,to2)

  if (output_looktime$m1 > 0) {
    #Make a list for Sampling
    if(discrete == 0) {
      survivaldata_looktime <- list("n1" = n1, "m1"=output_looktime$m1,
                                    "n1event"=output_looktime$n1event,
                                    "y_looktime"=
                                      as.array(output_looktime$y_looktime),
                                    "z_looktime" =
                                      as.array(output_looktime$z_looktime), "mu" = mu,
                                    "sigma" = sigma, "a" = a, "b" = b)
    }

    else{survivaldata_looktime <- list("n1" = n1, "m1"=output_looktime$m1,
                                       "n1event"=output_looktime$n1event,
                                       "y_looktime"=
                                         as.array(output_looktime$y_looktime_round),
                                       "z_looktime" =
                                         as.array(output_looktime$z_looktime_round), "mu" = mu,
                                       "sigma" = sigma, "a" = a, "b" = b)
    }

    #Sampling
    samples_looktime <- rstan::sampling(stanmodels$looktime,
                                        data = survivaldata_looktime,
                                        iter = 10000, chains = 4, cores = 6,
                                        refresh = 0,
                                        control = list(adapt_delta = 0.999))

    #######Second Step#######
    #Extract Parameter Values
    sequence_theta <- rstan::extract(samples_looktime,'theta')$theta
    sequence_kappa <- rstan::extract(samples_looktime,'kappa')$kappa

    obs_theta <- data.frame(median(sequence_theta), quantile(sequence_theta,0.05),
                            quantile(sequence_theta,0.95))
    obs_kappa <- data.frame(median(sequence_kappa), quantile(sequence_kappa,0.05),
                            quantile(sequence_kappa,0.95))

    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    result <- foreach (r=1:d, .combine = rbind) %dopar% {
      #ensures random gamma j is taken
      random_index <- sample(x=length(sequence_theta), size=1)

      #gives us the a,e,x,closs,z,y,delta for gammaj
      prediction <- prediction_samples(dataSurv,n,n1,output_looktime$m1,tau,omega2,
                                       sequence_kappa[random_index],
                                       sequence_theta[random_index],k)
      #combines observed and predicted data for each gamma j
      output_studyend  <- data_studyend(dataSurv,prediction, tau, k,till,to,to2)
      #the list for rstan model
      if(discrete == 0) {
        survivaldata_studyend  <- list("ncens"=output_studyend$ncens,
                                       "nevent"=output_studyend$nevent,
                                       "y_studyend"=
                                         as.array(output_studyend$y_studyend),
                                       "z_studyend" =
                                         as.array(output_studyend$z_studyend),"mu" = mu,
                                       "sigma" = sigma, "a" = a, "b" = b)
      }
      else{survivaldata_studyend  <- list("ncens"=output_studyend$ncens,
                                          "nevent"=output_studyend$nevent,
                                          "y_studyend"=
                                            as.array(output_studyend$y_studyend_round),
                                          "z_studyend" =
                                            as.array(output_studyend$z_studyend_round),"mu" = mu,
                                          "sigma" = sigma, "a" = a, "b" = b)
      }

      #Estimation predictive posterior distribution
      samples_studyend  <- rstan::sampling(stanmodels$studyend,
                                           data = survivaldata_studyend,
                                           iter = 10000, chains = 4, cores = 6,
                                           control = list(adapt_delta = 0.999))

      #Extract Parameter Values
      sequence_theta2 <- rstan::extract(samples_studyend,'theta')$theta
      sequence_kappa2 <- rstan::extract(samples_studyend,'kappa')$kappa
      #Get the quantile with which the 12months Posterior is reached -->Prob.
      p0_calc <- 1 - pweibull(tau, shape = sequence_kappa2, scale = sequence_theta2)


      pred_theta <- c(median(sequence_theta2))
      pred_kappa <- c(median(sequence_kappa2))


      #for each of the d predictions for a specific gammaj we look if the gamma
      #j is greather than the gamma under the Alternative
      prob <- length(which(p0_calc > p0)) /
        length(sequence_theta2)

      if (prob > eta_star) {success <- 1} else {success <- 0}
      return(data.frame("success" = success, "pred_theta" = pred_theta, "pred_kappa"=pred_kappa))

    }
    stopCluster(cl)
    all_pred_theta <- data.frame(median(result$pred_theta),
                                 quantile(result$pred_theta,0.05),
                                 quantile(result$pred_theta,0.95))
    all_pred_kappa <- data.frame(median(result$pred_kappa),
                                 quantile(result$pred_kappa,0.05),
                                 quantile(result$pred_kappa,0.95))


    PoS <- sum(result$success)/d
    if(PoS>lambda) {decision <- "Trial continues"} else{
      decision <- "Stop for futility"}}
  else {PoS <- "NONE, as there are no patients left w/o event"}




  colnames(all_pred_theta) <- c("Median" , "5%_Quantile" , "95%_Quantile")
  colnames(obs_theta) <- c("Median" , "5%_Quantile" , "95%_Quantile")
  colnames(all_pred_kappa) <- c("Median" , "5%_Quantile" , "95%_Quantile")
  colnames(obs_kappa) <- c("Median" , "5%_Quantile" , "95%_Quantile")


  output <- list()
  output$PoS <- PoS
  output$decision <- decision
  output$kappa <- kappa
  output$theta <- theta
  output$p0 <- p0
  output$p2_0 <- p2_0
  output$obs_kappa <- obs_kappa
  output$obs_theta <- obs_theta
  output$all_pred_kappa <- all_pred_kappa
  output$all_pred_theta <- all_pred_theta
  output$dataSurv <- dataSurv


  #if(rddata == 1) {output$dataSurv <- dataSurv}

  #define the S3 class

  class(output) <- "interim"
  print(output)
}

#'@export
print.interim <- function(x, ...) {
  cat("\nResults:\n")
  results <- data.frame(x$p0,
                        x$p2_0,
                        x$kappa,
                        x$theta,
                        x$PoS,
                        x$decision)
  colnames(results) <- c("Ass. Survival Probability EP1", "Ass. Survival Probability EP2", "Ass. Shape", "Ass. Scale","PoS","Decision")
  print(results)
  cat("\nSummary of Posterior Values at looktime:\n")
  rownames(x$obs_kappa) <- "Shape"
  rownames(x$obs_theta) <- "Scale"
  print(x$obs_kappa)
  cat("\n")
  print(x$obs_theta)
  cat("\nSummary of Posterior Values at studyend:\n")
  rownames(x$all_pred_kappa) <- "Shape"
  rownames(x$all_pred_theta) <- "Scale"
  print(x$all_pred_kappa)
  cat("\n")
  print(x$all_pred_theta)
  invisible(x)
}

#second stage samples generated
secondstage_samples <- function(data,n,n1,tau,omega2,kappa,theta,k) {

  #generate the accrual time
  e <- runif(n-n1, max(data$e) + k,max(data$e) + k + omega2)

  #generate time from e till event
  y <- rweibull(n-n1,kappa,theta)

  #e1 + y1 time of event
  y_e <- y + e

  #generate time from e till administrative censoring
  z <- rep(max(e)+tau,n-n1)

  #generate indicator variable for event = 1
  delta <- (y_e<=z)*1

  #combine all the generated data into one data frame
  dataSurv <- data.frame(e,y,y_e,z,delta)

  return(dataSurv)

}


data_studyend2 <- function(dataSurv,tau,till,to,to2) {
  #RoundUp Function - only required if data is simulated
  roundUp <- function(x,till,to,to2)
  {
    if(x<till) {
      to*(x%/%to + as.logical(x%%to))}
    else(to2*(x%/%to2 + as.logical(x%%to2)))
  }
  final <- max(dataSurv$e) + tau
  output <- list()
  #events out of n patients
  output$nevent <- as.numeric(
    nrow(dataSurv[dataSurv$y_e<=final,]))
  #aC out of n patients
  output$ncens <- as.numeric(
    nrow(dataSurv[dataSurv$y_e>final,]))
  #Survival time
  output$y_studyend <- dataSurv[dataSurv$y_e<=final,]$y
  for(i in 1:length(output$y_studyend)) {
    output$y_studyend_round[i] <- roundUp(output$y_studyend[i],till,to,to2) }
  #Censoring time
  output$z_studyend <- final - dataSurv[dataSurv$y_e>final,]$e
  if(output$ncens > 0) {
    for(i in 1:length(output$z_studyend)) {
      output$z_studyend_round[i] <- roundUp(output$z_studyend[i],till,to,to2) }}
  else{output$z_studyend_round <- output$z_studyend}
  return(output)
}

final <- function(rddata,discrete,interim,data,n,n1,p0,p2_0,k,tau,omega2,kappa,theta,eta_star,mu,sigma,a,b,till,to,to2) {
  PProb <- c()
  decision <- c()
  success <- c()
  if(rddata == 1) {
    dataSurv_first <- interim$dataSurv
    #Data generated in the second stage
    dataSurv_second <- secondstage_samples(dataSurv_first,n,n1,tau,omega2,kappa,theta,k)
    dataSurv_total <- rbind(dataSurv_first,dataSurv_second)
  }
  else{dataSurv_total <- data}

  output_studyend <- data_studyend2(dataSurv_total,tau,till,to,to2)
  #Make a list for Sampling
  if(discrete == 0) {
    survivaldata_studyend  <- list("ncens"=output_studyend$ncens,
                                   "nevent"=output_studyend$nevent,
                                   "y_studyend"=
                                     as.array(output_studyend$y_studyend),
                                   "z_studyend" =
                                     as.array(output_studyend$z_studyend),"mu" = mu,
                                   "sigma" = sigma, "a" = a, "b" = b)
  }

  else{survivaldata_studyend  <- list("ncens"=output_studyend$ncens,
                                      "nevent"=output_studyend$nevent,
                                      "y_studyend"=
                                        as.array(output_studyend$y_studyend_round),
                                      "z_studyend" =
                                        as.array(output_studyend$z_studyend_round),"mu" = mu,
                                      "sigma" = sigma, "a" = a, "b" = b)
  }

  #Sampling
  samples_studyend <- rstan::sampling(stanmodels$studyend,
                                      data = survivaldata_studyend,
                                      iter = 10000, chains = 4, cores = 6,
                                      refresh = 0,
                                      control = list(adapt_delta = 0.9999))

  #Extract Parameter Values
  sequence_theta <- rstan::extract(samples_studyend,'theta')$theta
  sequence_kappa <- rstan::extract(samples_studyend,'kappa')$kappa

  obs_theta <- data.frame(median(sequence_theta), quantile(sequence_theta,0.05),
                          quantile(sequence_theta,0.95))
  obs_kappa <- data.frame(median(sequence_kappa), quantile(sequence_kappa,0.05),
                          quantile(sequence_kappa,0.95))

  p0_calc <- 1 - pweibull(tau, shape = sequence_kappa, scale = sequence_theta)

  PProb <- length(which(p0_calc > p0)) / length(sequence_theta)

  if (PProb > eta_star) {success <- 1} else {success <- 0}
  if(success == 1) {decision <- "Reject H0"} else{
    decision <- "Do not reject H0"}

  output <- list()
  output$PProb <- PProb
  output$decision <- decision
  output$p0 <- p0
  output$p2_0 <- p2_0
  output$kappa <- kappa
  output$theta <- theta
  output$sampling <- samples_studyend
  output$obs_kappa <- obs_kappa
  output$obs_theta <- obs_theta

  #define the S3 class

  class(output) <- "final"
  print(output)
}
#'@export
print.final <- function(x, ...) {
  cat("\nResults:\n")
  results <- data.frame(x$p0,
                        x$p2_0,
                        x$kappa,
                        x$theta,
                        x$PProb,
                        x$decision)
  colnames(results) <- c("Ass. Survival Probability EP1","Ass. Survival Probability EP2", "Ass. Shape", "Ass. Scale","Posterior Probability","Decision")
  print(results)
  cat("\nRStan Sampling Output:\n")
  print(x$sampling)
  cat("\nSummary of Posterior Values at studyend:\n")
  colnames(x$obs_kappa) <- c("Median" , "5%_Quantile" , "95%_Quantile")
  colnames(x$obs_theta) <- c("Median" , "5%_Quantile" , "95%_Quantile")
  rownames(x$obs_kappa) <- "Shape"
  rownames(x$obs_theta) <- "Scale"
  print(x$obs_kappa)
  cat("\n")
  print(x$obs_theta)

  invisible(x)
}

