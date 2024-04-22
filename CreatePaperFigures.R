############################### Running this file will recreate all figures shown in the paper ###########################

# load libraries
library(ggplot2)
library(PIST)
library(fitdistrplus)
library(devtools)
######### GENERAL ###########

means_AUC <- determine.means()

means <- means_AUC[[2]]
AUC <- means_AUC[[1]]

# Create Concept figure

means.figure <- means[c(4,5,6)]
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))

#############################################################
##################### Comparison of all methods #############
#############################################################

# Part A
plot.1.prev.1 <- create.figure(N_tests = 10000, N = 50, means = means, AUC = AUC, true_prev = 0.01)
plot.1.prev.1

# Part B
plot.4 <- create.figure(N_tests = 10000, N = 50, means = means, AUC = AUC, true_prev = 0.04)
plot.4

# Part C
plot.8 <- create.figure(N_tests = 10000, N = 50, means = means, AUC = AUC, true_prev = 0.08)
plot.8

#############################################################
#################### power analysis ###############################
############################################################

# load datafile from cluster
setwd("/Users/judit/pist/")
power.data <- read.table('power_analysis_fullData_4.txt')
colnames(power.data) <- c('NULL','AUC','Sample.size', 'true.prev','estimate','CI')
power.data[,2] <- round(power.data[,2],2)

#power.data <- power.data[power.data$true.prev==0.1,]

sample.sizes <- c(1,100, 500, seq(1000, 10000, 1000), 12500, 15000)
AUC <- sort(unique(power.data[,2]))[2:length(unique(power.data[,2]))]

power <- matrix(0, nrow = length(sample.sizes), ncol = length(AUC) )

for (i in seq(1, length(sample.sizes),1)){
  for (j in seq(1, length(AUC),1)){
    data.now <- power.data[ ( power.data[,2]==AUC[j] & power.data[,3]==sample.sizes[i] ),  ]
    #print(data.now)
    lower <-  data.now[1,4] - data.now[1,4]*0.25
    upper <- data.now[1,4] + data.now[1,4]*0.25
    tot <- dim(data.now)[1]
    print(tot)
    power[i,j] <- sum( (data.now[,5]>lower & data.now[,5]<upper) & ( data.now[,4]>(data.now[,5]-2*data.now[,6]) & data.now[,4]<(data.now[,5]+2*data.now[,6]))  )/tot
  }
}

# plot power graphs

cols = c("violetred","springgreen3","peru", "dodgerblue" , "orange", "pink","cyan")
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))
plot(sample.sizes,  power[,1], xlab="Number of individuals tested", ylab = "Power", type = 'l',
     lwd = 3, cex.axis = 2, cex.lab = 2, col = cols[1], ylim = c(0,1))

for (i in seq(2, length(AUC), 1)){
  lines(sample.sizes, power[,i], col = cols[i], lwd = 3)
}
lines(c(0, 20000), c(0.9, 0.9), lwd = 2.5, lty = 2)

legend("bottomright", c(as.character(AUC), 'Power level of 0.9'), col = c(cols,'black'), lwd = c(rep(3,7),2),
       lty = c(rep(1,7),2),title = "AUC-ROC")

# Figure with minimal number of samples for a power of 0.9
min_sample_size <- rep(0, length(AUC))

for (i in seq(1, length(AUC),1)){
  min_sample_size[i] <- sample.sizes[ power[,i]>0.9 ][1]
}

plot(AUC, min_sample_size, ylab='Minimal required number of ind. tested',
     xlab = 'Area under the ROC-curve (AUC-ROC)', lwd = 3, col = 'red', type = 'l', cex.lab=2, cex.axis = 2)


#############################################################
#################### Temporal trend #########################
#############################################################

mean <- means[6] # take mean such that AUC = 0.95

prevs <- c(0.015, 0.15)
N <- 50
N_tests <- 5000

results_value <- matrix(0, nrow = length(prevs), ncol = N)
results_naive_classical <- matrix(0, nrow = length(prevs), ncol = N)
results_naive_spec <- matrix(0, nrow = length(prevs), ncol = N)
results_naive_classical_ncor <- matrix(0, nrow = length(prevs), ncol = N)
results_naive_spec_ncor <- matrix(0, nrow = length(prevs), ncol = N)

threshold.clas = classical.threshold(mean = mean)
threshold.spec = classical.threshold(mean = mean, method = 'specificity')

for (i in seq(1, length(prevs), 1)){

  for (j in seq(1,N,1)){
    # simulate data given those distributions
    data.serosurvey <- sim.serosurvey( N_tests = N_tests, mean = mean, prev = prevs[i])

    # analyze data and estimate prev
    output <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean = mean, type = 'likelihood')
    results_value[i,j] <- output[[1]]

    # analyze data and estimate prev
    output.clas <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean=mean, type = 'naive',
                                 threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = TRUE)
    output.spec <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean=mean, type = 'naive',
                                 threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = TRUE)

    # what is the probability that our estimated prevalence is equal to the true prevalence?

    results_naive_classical[i,j] <- output.clas[[1]]
    results_naive_spec[i,j] <- output.spec[[1]]

    output.clas_ncor <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean=mean, type = 'naive',
                                 threshold = threshold.clas[1], TNR = threshold.clas[3], TPR = threshold.clas[2], correct = FALSE)
    output.spec_ncor <- analyze.serosurvey(data.serosurvey = data.serosurvey, mean=mean, type = 'naive',
                                 threshold = threshold.spec[1], TNR = threshold.spec[3], TPR = threshold.spec[2], correct = FALSE)

    # what is the probability that our estimated prevalence is equal to the true prevalence?

    results_naive_classical_ncor[i,j] <- output.clas_ncor[[1]]
    results_naive_spec_ncor[i,j] <- output.spec_ncor[[1]]
  }

}

increase_likelihood <- results_value[2,]/results_value[1,]
increase_classical <- results_naive_classical[2,]/results_naive_classical[1,]
increase_spec <- results_naive_spec[2,]/results_naive_spec[1,]
increase_classical_ncor <- results_naive_classical_ncor[2,]/results_naive_classical_ncor[1,]
increase_spec_ncor <- results_naive_spec_ncor[2,]/results_naive_spec_ncor[1,]

# create dataframe
increase_all <- matrix(0, nrow = 5*N, ncol = 2 )
increase_all[1:N,1] <- increase_likelihood
increase_all[(N+1):(2*N),1] <- increase_classical
increase_all[(2*N+1):(3*N),1] <- increase_spec
increase_all[(3*N+1):(4*N),1] <- increase_classical_ncor
increase_all[(4*N+1):(5*N),1] <- increase_spec_ncor

increase_all[1:N,2] <- 'likelihood'
increase_all[(N+1):(2*N),2] <- 'max Youden corrected'
increase_all[(2*N+1):(3*N),2] <- 'high specificity corrected'
increase_all[(3*N+1):(4*N),2] <- 'max Youden'
increase_all[(4*N+1):(5*N),2] <- 'high specificity'

colnames(increase_all) <- c('value','type')

increase_all <- as.data.frame(increase_all)
increase_all[,1] <- as.numeric(paste(increase_all[,1]))
increase_all[,2] <-factor(increase_all[,2], c('max Youden', 'max Youden corrected', 'high specificity', 'high specificity corrected', 'likelihood'))

ggplot()+ylim(0, 30)+
  geom_boxplot(data = increase_all, aes(x = type, y = value), colour = "dodgerblue4", fill = "dodgerblue", alpha = 0.2,outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=TRUE ) +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),axis.text = element_text(size=20),axis.title = element_text(size=22), text = element_text(size=22),legend.position="bottom",
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),   axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5))+
  geom_hline(aes(yintercept = 10, color = "True increase"), size = 1)+
  scale_colour_manual("", values=c("True increase" = "red"),
                      aesthetics = c("colour"))  + ylab('Estimated fold increase')



#############################################################
#################### Sampling cases and controls #################
#############################################################

# TAKE AUC = 0.9, 0.95 and 0.98
examples.AUC <- c(2.4924925, 4.4644645 , 5.3953954) # 0.3703704  1.0710711  1.7617618  2.4924925  3.3333333  4.4644645  5.3953954 10.0000000

# Sample from the true distributions
N_tests <- 10000
real.prev <- 0.08
X <- c(75, 100, 150, 200, 300, 500, 1000)# number of sampled people from both distributions
N <- 100

results.ll <- estimate.CC.distr(true.mean = examples.AUC[1], sample.sizes = X, name="0.9",N_tests = N_tests,  N = N, real.prev = real.prev)
results.ll.2 <- estimate.CC.distr(true.mean = examples.AUC[2], sample.sizes = X, name="0.9",N_tests = N_tests,  N = N, real.prev = real.prev)
results.ll.3 <- estimate.CC.distr(true.mean = examples.AUC[3], sample.sizes = X, name="0.9",N_tests = N_tests,  N = N, real.prev = real.prev)

all.results <- rbind(results.ll, results.ll.2, results.ll.3)

all.results <- as.data.frame(all.results)
colnames( all.results ) <- c('Value','resample', 'type')
all.results[,2] <- factor(all.results[,2], c( "75", "100", "150", "200","300","500","1000"))
all.results[,1] <- as.numeric(paste((all.results[,1])))

ggplot() +
  geom_violin(data=all.results, aes(x = resample, y = Value,colour=type, fill=type), alpha = .2, size=0.3)+
  theme(axis.text = element_text(size=24),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom",
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_colour_manual("AUC-ROC", values = c("0.9"="violetred3","0.95" = "dodgerblue" , "0.98"="springgreen3",
                                            "True prevalence" = "red" ),
                      aesthetics = c("colour","fill"))+
  geom_hline(aes(yintercept = real.prev, color = "True prevalence"), size = 1 )+
  ylab("Estimated seroprevalence") + xlab("Number of samples from controls and cases")


# Sampling conceptual figure

means.figure <- means[c(3,5,7)]
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))

X = 100

titer.TN.sampled <- rgamma(X, shape = 1, scale = 1)# sample(titers, size = X, prob = titer.TN[,2], replace = TRUE)
titer.TP.sampled <- rgamma(X, shape = means.figure[1], scale = 1)  #sample(titers, size = X, prob = titer.TP[,2], replace = TRUE)

fit.TP <- fitdist(titer.TP.sampled, distr = "gamma")
fit.TN <- fitdist(titer.TN.sampled, distr = "gamma", fix.arg = list(shape = 1))

titer.sampled <- c(titer.TN.sampled, titer.TP.sampled)
titer.sampled <- cbind(as.numeric(paste(titer.sampled)), (c(rep("Controls", X), rep("Cases",X))))

titer.sampled <- as.data.frame(titer.sampled)
titer.sampled[,1]<- as.numeric(paste(titer.sampled[,1]))
colnames(titer.sampled) <- c('value', 'name')

titers <- seq(0,10,0.01)
titer.TP <- cbind(seq(0,10,0.01), dgamma(seq(0,10,0.01), shape = means.figure[1], scale =1)/sum(dgamma(seq(0,10,0.01), shape = means.figure[1], scale =1)))
titer.TN <- cbind(seq(0,10,0.01), dgamma(seq(0,10,0.01), shape = 1, scale =1)/sum(dgamma(seq(0,10,0.01), shape = 1, scale =1)))

TP.df <- as.data.frame(titer.TP)
TP.df[,2] <- TP.df[,2]
TN.df <- as.data.frame(titer.TN)
TN.df[,2] <- TN.df[,2]

#estimated curve
titer.TP.est <- cbind(seq(0,10,0.01), dgamma(seq(0,10,0.01), shape = fit.TP[[1]][1], scale =fit.TP[[1]][1])/sum(dgamma(seq(0,10,0.01), shape = fit.TP[[1]][1], scale =fit.TP[[1]][1])))
titer.TN.est <- cbind(seq(0,10,0.01), dgamma(seq(0,10,0.01), shape = 1, scale =fit.TN[[1]][1])/sum(dgamma(seq(0,10,0.01), shape = 1, scale =fit.TN[[1]][1])))

TP.df.est <- as.data.frame(titer.TP.est)
TP.df.est[,2] <- TP.df.est[,2]
TN.df.est <- as.data.frame(titer.TN.est)
TN.df.est[,2] <- TN.df.est[,2]

ggplot() +
  geom_histogram(data = subset(titer.sampled, name =="Controls"), aes(x=value, y = 0.01*..density..) ,
                 fill = "grey42", alpha = 0.4, binwidth = 0.1)+
  geom_histogram(data = subset(titer.sampled, name =="Cases"), aes(x=value, y = 0.01*..density..) ,
                 fill = "orange", alpha = 0.4, binwidth = 0.1)+
  geom_line(data = TP.df, aes(x =titers, y = V2, colour = "Cases"), size = 2) +
  geom_line(data = TN.df, aes(x =titers, y = V2, colour = "Controls"), size = 2) +
  geom_line(data = TP.df.est, aes(x =titers, y = V2, colour = "Cases"), size = 1.5, alpha = 1, linetype = "dashed") +
  geom_line(data = TN.df.est, aes(x =titers, y = V2, colour = "Controls"), size = 1.5, alpha = 1, linetype = "dashed") +
  theme(axis.text = element_text(size=24),axis.title = element_text(size=24), text = element_text(size=24),legend.position="bottom",
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+
  ylab("Frequency") +
  xlab("Raw serological test output") +
  scale_colour_manual("", values=c("Controls" = "grey42", "Cases"="orange"),
                      aesthetics = c("colour"))

#############################################################
#################### Sens and spec for cutoff based methods #################
#############################################################

# sensitivity and specificity of various methods and AUC values

sens_spec <- rep(0, nrow = length(means))
spec_spec <- rep(0, nrow = length(means))
sens_max <- rep(0, nrow = length(means))
spec_max <- rep(0, nrow = length(means))

titers <- seq(0,20,0.1)
titer.TN <- cbind(titers, 1000*exp(-1*titers)/sum(1000*exp(-1*titers)))

for (i in seq(1, length(means), 1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means[i], rate = 1))

  # use as a threshold the crossing between both curves --> different choises would have a different effect!
  threshold.clas = classical.threshold(mean = means[i])
  threshold.spec = classical.threshold(mean = means[i], method = 'specificity')

  sens_spec[i] <- threshold.spec[2]
  spec_spec[i] <- threshold.spec[3]
  sens_max[i] <- threshold.clas[2]
  spec_max[i] <- threshold.clas[3]

}

plot(AUC, sens_spec,   ylim = c(0,1), col = "violetred3", lwd = 4, pch= 19, xlab = "Area under ROC-curve (AUC-ROC)", ylab = "Sensitivity or specificity", cex.axis = 2, cex.lab = 2)
points(AUC, sens_max, col =  "dodgerblue", lwd = 4, pch = 19)
points(AUC, spec_spec, col = "violetred3", lwd = 4, pch = 17)
points(AUC, spec_max, col =  "dodgerblue", lwd = 4, pch = 17)
legend("bottomright", c("Sensitivity","Specificity","High specificity","Max Youden"), lty = c(0,0,1,1), col = c("black","black","violetred3","dodgerblue"),
       pch = c(19,17,26,26,26), lwd = 3)


#############################################################
########### Create bootstrap figure #########################
#############################################################

create.figure.bootstrap(means = means, AUC = AUC )


