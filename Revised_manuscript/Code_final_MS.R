# load libraries
library(ggplot2)
library(SeroPrev)
library(fitdistrplus)
library(devtools)
library(gridBase)
library(msm)
library(truncdist)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)
library(rstan)

library(tidyverse)
library(cowplot)
library(gridGraphics)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

######### GENERAL ###########

# I fix the scale parameter of the case data such that I obtain AUC values for the ROC curve between the control
# and case data which are realistic for the serological assays

means <- c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287) # can be determined with: determine.means()
AUC <- c(0.801,0.849,0.9,0.95,0.975,0.990) # can be determined with: determine.means()

means.figure <- means[c(4,5,6)]
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))

# the variable code determines whether all code is ran (TRUE) or whether previously ran data is used for the
# generation of the figures
code = FALSE

# Set working directory
setwd('/Users/judith/pist/Revised_manuscript')

#####################################################
# Figure 1 -- Overview of comparison of all methods
#####################################################

if (code){

  # run the bayesian models
  RunBayesian_models(N=50, true_prevs=c(0.01,0.04,0.08), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
             AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990), file.save = 'data.bayesian.Fig1.txt', rstan.model = "/Users/judith/testing/Revised_manuscript/cutoff_based.stan")

  # run the non-bayesian models
  RunNonBayesian_models(N=1, true_prevs=c(0.01,0.04,0.08), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
                                    AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990), file.save = 'data.nonbayesian.Fig1.txt')

  all <- rbind(read_delim(file = 'Data/data.bayesian.Fig1.txt', col_names = FALSE, delim = " "), read_delim(file = 'Data/Data_fig1_complete.txt', col_names = FALSE, delim = " "))
  colnames(all) <- c('row', 'value',  'AUC', 'True prevalence',     'type', 'mean',  'lower CI', 'upper CI')
  write_delim(all, file = 'Data/Data_fig1_complete_coded.txt')

}

data.F1 <- as.data.frame(read_delim(file = 'Data/Data_fig1_complete.txt', col_names = TRUE, delim = " "))[,2:8] # make sure that this file corresponds with the file where you have saved the data before.
data.F1$AUC[data.F1$AUC==0.849] <- 0.85

new_range_AUC <- c(0.80, 0.85, 0.90, 0.95, 0.975, 0.99)

theme_set(theme_bw())

# Create correct dataframe
#est = as_tibble(data.F1)
est = subset(data.F1, data.F1$AUC %in% new_range_AUC )
est$AUC = factor(est$AUC, levels= new_range_AUC)
est$size.interval = est$`upper CI`-est$`lower CI`
est$prev2 <- paste0("True value = ",est$`True prevalence`*100,"%")

g1 = ggplot(est) +
  ggtitle('A')+
  geom_violin(aes(x=AUC,y=value,colour=type, fill = type), width = 1, alpha = .2, size=0.7, position = position_dodge(0.8) ) +
  stat_summary(fun=median, aes(y = value, x = AUC, fill = type, colour = type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept=`True prevalence`),linetype=2,size=.5,alpha=.6) +
  facet_grid(prev2~., scales = "free_y") +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  scale_colour_manual('',values=c("High specificity - RG" =  "seagreen",
                                  "Max Youden - RG" =  "steelblue",
                                  "High specificity - Bayesian" = "cyan",
                                  'Max Youden - Bayesian' = "violetred3",
                                  'Mixture model' = "orange"), aesthetics = c('colour',"fill")) +
  labs(x="AUC-ROC",y="Estimated cumulative incidence") +
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), text = element_text(size=20),
        legend.position="bottom", title = element_text(size = 26, face='bold'),
        panel.grid.major =element_blank(), panel.background = element_blank())+
  guides(colour=guide_legend(nrow=2,ncol = 3, byrow=TRUE))
g2 = ggplot(est) +
  ggtitle('B')+
  geom_violin(aes(x=AUC,y=size.interval,colour=type, fill = type), width = 1, alpha = .2, size=0.7, position = position_dodge(0.8) ) +
  stat_summary(fun=median, aes(y = size.interval, x = AUC, fill = type, colour = type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  #geom_line(aes(x=AUC,y=half_95CI,group=type,colour=type), size = 1.5) +
  facet_grid(prev2~.) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  scale_colour_manual('', values=c("High specificity - RG" =  "seagreen",
                                   "Max Youden - RG" =  "steelblue",
                                   "High specificity - Bayesian" = "cyan",
                                   'Max Youden - Bayesian' = "violetred3",
                                   'Mixture model' = "orange")) +
  #coord_cartesian(ylim=c(0,0.055)) +
  labs(x="AUC-ROC",y="Size of 95% uncertainty interval")+
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), text = element_text(size=20),
        legend.position="bottom", title = element_text(size = 26, face='bold'),
        panel.grid.major =element_blank(), panel.background = element_blank())+
  guides(colour=guide_legend(nrow=2,ncol = 3, byrow=TRUE))
mylegend<-g_legend(g1)

pdf(file = '/Users/judith/testing/Figures/Revised/Overview_all_methods_new.pdf', width = 16, height = 18)
all <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                g2 + theme(legend.position="none"),
                                nrow=1), mylegend, nrow = 2,heights=c(10,1))
dev.off()



#####################################################
# Figure 2 -- Temporal trend in estimating cumulative incidence
#####################################################

if (code){

  # run the bayesian models
  RunBayesian_models(N=50, true_prevs=c(0.015,0.15), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
                     AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990), file.save = 'data.bayesian.Fig2.txt', rstan.model = "/Users/judith/testing/Revised_manuscript/cutoff_based.stan")

  # run the non-bayesian models
  RunNonBayesian_models(N=1, true_prevs=c(0.015,0.15), means = c(2.322645, 2.719439, 3.314629, 4.306613, 5.320641, 6.643287),
                        AUCs = c(0.801,0.849,0.9,0.95,0.975,0.990), file.save = 'data.nonbayesian.Fig2.txt' )

  all.2 <- rbind(read_delim(file = 'Data/data.bayesian.Fig2.txt', col_names = FALSE, delim = " "), read_delim(file = 'Data/Data_fig2_complete.txt', col_names = FALSE, delim = " "))
  colnames(all.2) <- c('row', 'value',  'AUC', 'True prevalence',     'type', 'mean',  'lower CI', 'upper CI')
  write_delim(all.2, file = 'Data/Data_fig2_complete_coded.txt')
}

#Load data that is produced
data.F2 <- as.data.frame(read_delim(file = 'Revised_manuscript/Data/Data_fig2_complete.txt', col_names = TRUE, delim = " "))[,2:8] # make sure that this file corresponds with the file where you have saved the data before.
data.F2$AUC[data.F2$AUC==0.849] <- 0.85

fig2_results <- data.frame(matrix(0, nrow = 350, ncol = 2))
colnames(fig2_results) <- c('value', 'type')
fig2_results$value <- data.F2$value[data.F2$`true prevalence`==0.15]/data.F2$value[data.F2$`true prevalence`==0.015]
fig2_results$type <- data.F2$type[data.F2$`true prevalence`==0.15]
est_f2 = as_tibble(fig2_results)

est_f2$type <- factor(est_f2$type, levels =c("High specificity - uncorrected",
                                             "Max Youden - uncorrected",
                                             "High specificity - RG",
                                             "Max Youden - RG",
                                             "High specificity - Bayesian",
                                             'Max Youden - Bayesian',
                                             'Mixture model' ) )


pdf(file = '/Users/judith/testing/Figures/Revised/Fold_increase.pdf', width = 10, height = 12)
ggplot(est_f2) +
  geom_violin(aes(x = type, y=value,colour=type, fill = type), width = 1, alpha = .2, size=0.7, position = position_dodge(0.8) ) +
  stat_summary(fun=median, aes(y = value, x = type, fill = type, colour = type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept=10),linetype=2,size=.5,alpha=.6) +
  scale_colour_manual('',values=c("High specificity - uncorrected" =  "dodgerblue",
                                  "Max Youden - uncorrected" =  "maroon",
                                  "High specificity - RG" =  "seagreen",
                                  "Max Youden - RG" =  "steelblue",
                                  "High specificity - Bayesian" = "cyan",
                                  'Max Youden - Bayesian' = "violetred3",
                                  'Mixture model' = "orange"), aesthetics = c('colour',"fill")) +
  labs(x="",y="Estimated fold increase") +
  ylim(0, 26)+
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), text = element_text(size=20),
        legend.position="none", title = element_text(size = 26, face='bold'),
        axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5))
dev.off()

#####################################################
# Figure 3 -- Power figure, impact of the size of the serosurvey
#####################################################
# the data to create the power figure has been run on a cluster,
# as it requires many simulations, here we just load the data and create the figure

power.data.1 <- read.delim('power_analysis_revised_full.txt', fill =T,sep = " ", quote=" ", header = FALSE)
power.data.2 <- read.delim('power_analysis_revised_full_2.txt', fill =T , sep = ' ', quote = " ", header= FALSE)
power.data <- rbind(power.data.1, power.data.2)[,2:7]

# filter out lines with problems
power.data[,1] <- as.numeric(power.data[,1])
power.data[,2] <- as.numeric(power.data[,2])
power.data[,3] <- as.numeric(power.data[,3])
power.data[,4] <- as.numeric(power.data[,4])
power.data[,5] <- as.numeric(power.data[,5])
power.data[,6] <- as.numeric(power.data[,6])
power.data <- power.data[!is.na(power.data[,1]),]
power.data <- power.data[!is.na(power.data[,2]),]
power.data <- power.data[!is.na(power.data[,3]),]
power.data <- power.data[!is.na(power.data[,4]),]
power.data <- power.data[!is.na(power.data[,5]),]
power.data <- power.data[!is.na(power.data[,6]),]

colnames(power.data) <- c('AUC','Sample.size', 'true.prev','Validation.size','estimate', 'CI')
power.data[,1] <- round(power.data[,1],3)

# select only the 5000 individuals in testing size
power.data_fig3 <- power.data[power.data$Validation.size==5000&(power.data$AUC %in% c(0.801,0.849,0.9,0.95,0.975,0.98,0.99)),]

sample.sizes <- c(1,100, 500, seq(1000, 10000, 1000))
AUC <- c(0.801,0.849,0.9,0.95,0.975,0.990)

power <- matrix(0, nrow = length(sample.sizes), ncol = length(AUC) )

for (i in seq(2, length(sample.sizes),1)){
  for (j in seq(1, length(AUC),1)){
    data.now <- power.data_fig3[ ( power.data_fig3[,1]==AUC[j] & power.data_fig3[,2]==sample.sizes[i] ),  ]
    lower <-  data.now[1,3] - data.now[1,3]*0.25
    upper <- data.now[1,3] + data.now[1,3]*0.25
    tot <- dim(data.now)[1]
    print(tot)
    power[i,j] <- sum( (data.now[,5]>lower & data.now[,5]<upper) & ( data.now[,3]>(data.now[,6]-2*data.now[,6]) & data.now[,3]<(data.now[,5]+2*data.now[,6]))  )/tot
  }
}

# Check to which AUC value the minimal requirement is valid
power.data.additional <- read.delim('/Users/judith/testing/Revised_manuscript/Data/power_analysis_revised_additional.txt',
                           fill =T , sep = ' ', quote = " ", header= FALSE)[,2:7]
colnames(power.data.additional) <- c('AUC','serosize','true.prev','val.size','upper','lower')
AUC.tested <- sort(unique(power.data.additional$AUC))
sample.sizes.extra <- sort(unique(power.data.additional$serosize))
power.add <- matrix(0, nrow = length(AUC.tested), ncol = length(sample.sizes))

for ( j in seq(1, length(sample.sizes.extra),1)){
  for (i in seq(1, length(AUC.tested),1)){
    data.now <- power.data.additional[ ( power.data.additional[,1]==AUC.tested[i] ) & (power.data.additional$serosize==sample.sizes.extra[j]),  ]
    lower <-  data.now[3,3] - data.now[3,3]*0.25
    upper <- data.now[3,3] + data.now[3,3]*0.25
    tot <- dim(data.now)[1]
    power.add[i, j] <- sum( (data.now[,5]>lower & data.now[,5]<upper) & ( data.now[,3]>(data.now[,6]-2*data.now[,6]) & data.now[,3]<(data.now[,5]+2*data.now[,6])), na.rm = T )/tot
  }
}

colnames(power.add) <- sample.sizes.extra
rownames(power.add) <- AUC.tested

# plot power graphs
pdf(file = '/Users/judith/testing/Figures/Revised/Power_plot_serosurvey_size_5000.pdf', width = 16, height = 10)

par(mfrow=c(1,2),oma=c(1.5,2.5,1,4),mar=c(4,5,2,2))
cols = c("violetred","springgreen3","peru", "dodgerblue" , "orange", "pink")
plot(sample.sizes,  power[,1], xlab="Number of individuals in serosurvey", ylab = "Power", type = 'l',
     lwd = 3, cex.axis = 2, cex.lab = 2, col = cols[1], ylim = c(0,1))
mtext('A', side = 3, line = 1, at = -3, cex = 2, outer =F, font = 2)
for (i in seq(2, length(AUC), 1)){
  lines(sample.sizes, power[,i], col = cols[i], lwd = 3)
}
lines(c(0, 20000), c(0.9, 0.9), lwd = 2.5, lty = 2)

legend("bottomright", c('0.8','0.85','0.9','0.95','0.975','0.99', 'Power level of 0.9'), col = c(cols,'black'), lwd = c(rep(3,6),2),
       lty = c(rep(1,6),2),title = "AUC-ROC")

# Figure with minimal number of samples for a power of 0.9
min_sample_size <- rep(0, length(AUC))
min_sample_size.extra <- rep(0, length(AUC.tested))

for (i in seq(1, length(AUC),1)){
  min_sample_size[i] <- sample.sizes[ power[,i]>0.9 ][1]
}

for (i in seq(1, length(AUC.tested)-1,1)){
  min_sample_size.extra[i] <- sample.sizes.extra[ power.add[i,]>0.9 ][1]
}
min_sample_size.extra[length(AUC.tested)] <- 3000

plot(AUC, (min_sample_size), ylab='Minimal required number of ind. in serosurvey',
     xlab = 'Area under the ROC-curve (AUC-ROC)', lwd = 3, col = 'red', type = 'l', cex.lab=2, cex.axis = 2,
     ylim = c(1,100000), log = 'y')
lines(AUC.tested, min_sample_size.extra, lwd = 3, col = 'red')
mtext('B', side = 3, cex = 2, line = 1, at = 0.75, font = 2)
rect(0.8, 1, 0.85, 100000, border = "grey", col = "grey")

dev.off()


#####################################################
# Figure 4 -- Power figure, impact of the size of the validation set
#####################################################
# Panel A, the conceptual figure
means.figure <- means[c(3,5,6)]
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))

X = 150

titer.TN.sampled <- rgamma(X, shape = 1, scale = 1)# sample(titers, size = X, prob = titer.TN[,2], replace = TRUE)
titer.TP.sampled <- rgamma(X, shape = means.figure[3], scale = 1)  #sample(titers, size = X, prob = titer.TP[,2], replace = TRUE)

pdfs <- sim.validation.data(control.data = titer.TN.sampled , case.data = titer.TP.sampled)

titer.sampled <- c(titer.TN.sampled, titer.TP.sampled)
titer.sampled <- cbind(as.numeric(paste(titer.sampled)), (c(rep("Controls", X), rep("Cases",X))))

titer.sampled <- as.data.frame(titer.sampled)
titer.sampled[,1]<- as.numeric(paste(titer.sampled[,1]))
colnames(titer.sampled) <- c('value', 'name')

titers <- seq(0,20,0.01)
titer.TP <- cbind(seq(0,20,0.01), dgamma(seq(0,20,0.01), shape = means.figure[3], scale =1)/sum(dgamma(seq(0,20,0.01), shape = means.figure[1], scale =1))*100)
titer.TN <- cbind(seq(0,20,0.01), dgamma(seq(0,20,0.01), shape = 1, scale =1)/sum(dgamma(seq(0,20,0.01), shape = 1, scale =1))*100)

TP.df <- as.data.frame(titer.TP)
TP.df[,2] <- TP.df[,2]
TN.df <- as.data.frame(titer.TN)
TN.df[,2] <- TN.df[,2]

#estimated curve
titer.TP.est <- cbind(seq(0,20,0.01), pdfs[[1]](titers))
titer.TN.est <- cbind(seq(0,20,0.01), pdfs[[2]](titers))

TP.df.est <- as.data.frame(titer.TP.est)
TP.df.est[,2] <- TP.df.est[,2]
TN.df.est <- as.data.frame(titer.TN.est)
TN.df.est[,2] <- TN.df.est[,2]

P3_A <-  ggplot() +
  geom_histogram(data = subset(titer.sampled, name =="Controls"), aes(x=value, y = ..density..,
                                                                      fill = "Controls") , alpha = 0.4)+
  geom_histogram(data = subset(titer.sampled, name =="Cases"), aes(x=value, y = ..density..,
                                                                   fill = "Cases") , alpha = 0.4)+
  geom_line(data = TP.df, aes(x =titers, y = V2, linetype = "True distribution"), colour = "orange", size = 2) +
  geom_line(data = TN.df, aes(x =titers, y = V2, linetype = "True distribution"), colour = "grey42", size = 2) +
  geom_line(data = TP.df.est, aes(x =titers, y = V2, linetype = "Estimated distribution"), colour = "orange", size = 1.5, alpha = 1) +
  geom_line(data = TN.df.est, aes(x =titers, y = V2, linetype = "Estimated distribution"), colour = "grey42", size = 1.5, alpha = 1) +
  theme(axis.text = element_text(size=20),axis.title = element_text(size=20), text = element_text(size=20),legend.position="bottom",
        panel.grid.major =element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0))+
  ylab("Frequency") +
  xlab("Raw serological test output") +
  scale_fill_manual("", values=c("Controls" = "grey42", "Cases"="orange"),
                    aesthetics = c("fill")) +
  scale_linetype_manual("",values=c("True distribution"=1,"Estimated distribution"=4))+
  ggtitle("A")

# Panel B
#### Effect of sample size on Estimated cumulative incidence

examples.AUC <- means[4:6]
AUCs <- AUC[4:6]

# Sample from the true distributions
N_tests <- 10000
real.prev <- 0.08
sample.sizes <- c(100, 200, 300, 500, 1000, 5000)# number of sampled people from both distributions
N <- 50

results_value <- matrix(0, nrow = length(sample.sizes)*length(examples.AUC), ncol = N)
results_sd <- matrix(0, nrow = length(sample.sizes)*length(examples.AUC), ncol = N)

#simulate results for likelihood method
for (k in seq(1, length(examples.AUC),1)){
  for (i in seq(1, length(sample.sizes), 1)){

    for (j in seq(1,N,1)){

      # Simulate validation data
      control.data <- rgamma(sample.sizes[i], shape = 1, scale = 1)
      case.data <- rgamma(sample.sizes[i], shape = examples.AUC[k], scale = 1)

      # simulate data given those distributions
      data.serosurvey <- sim.serosurvey( N_tests = N_tests,  mean = examples.AUC[k], prev=real.prev)

      ############### likelihood-based method #######################
      # simulate validation data and create pdf for both distributions from there
      pdfs <- sim.validation.data(control.data = control.data, case.data = case.data, nn = sample.sizes[i], np = sample.sizes[i])

      # analyze data and estimate prev
      output <- analyze.serosurvey.emperical(data.serosurvey = data.serosurvey, pdf.case.smooth = pdfs[[1]],
                                             pdf.control.smooth = pdfs[[2]], type = 'likelihood')

      # what is the probability that our estimated prevalence is equal to the true prevalence?
      results_value[(k-1)*length(sample.sizes)+i,j] <- output[[1]]
      results_sd[(k-1)*length(sample.sizes)+i,j] <- output[[2]]

    }
  }
}

results.ll <- matrix(0, nrow = length(sample.sizes)*length(examples.AUC)*N, ncol = 3)
for(k in seq(1, length(examples.AUC),1)){
  for (i in seq(1, length(sample.sizes), 1)){
    for (j in seq(1, N,1 )){
      results.ll[ ((k-1)*(length(sample.sizes)*N) + (i-1)*N + j) ,1] <- results_value[(k-1)*length(sample.sizes)+i,j]
      results.ll[ ((k-1)*(length(sample.sizes)*N) + (i-1)*N + j) ,2] <- round(AUCs[k], 3)
      results.ll[ ((k-1)*(length(sample.sizes)*N) + (i-1)*N + j) ,3] <- sample.sizes[i]
    }
  }
}

results.ll <- as.data.frame(results.ll)
colnames( results.ll ) <- c('Value','AUC', 'Sample.size')
results.ll$Sample.size <- factor(results.ll$Sample.size, c('100','200', '300', '500', '1000', '5000'))
results.ll$AUC <- factor(results.ll$AUC, c('0.95','0.975','0.99'))

cols = c("violetred","springgreen3","peru", "dodgerblue" , "orange", "pink","cyan")

P3_B <- ggplot() +
  geom_violin(data=results.ll, aes(x = Sample.size, y = Value, colour=AUC, fill=AUC), alpha = .2, size=0.3, position = position_dodge(0.7))+
  ggtitle('B')+
  theme(axis.text = element_text(size=20),axis.title = element_text(size=20), text = element_text(size=20),legend.position="bottom",
        panel.grid.major =element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),
        plot.title = element_text(hjust = 0))+
  scale_colour_manual("AUC-ROC", values = c("0.95"=cols[4],"0.975" = cols[5] , "0.99"=cols[6]),
                      aesthetics = c("colour", 'fill'))+
  geom_hline(aes(yintercept = real.prev), color = "red", size = 1 )+
  ylab("Estimated cumulative incidence") + xlab("Number of controls and cases in validation data")

# panel C
# add power graph for sample.size of validation data

sample.sizes.validation <- sort(unique(power.data$Validation.size))
sample.sizes.serosurvey <- sort(unique(power.data$Sample.size))

power.10000 <- power(power.data, sample.sizes.validation, AUC, serosurvey.size=10000)

par(mfrow=c(2,2),oma=c(1.5,2.5,1,4),mar=c(4,5,2,2))
plot.new()
plot.new()
mtext('A', side = 3, line = 1, at = -3, cex = 2, outer =F)
plot(sample.sizes.validation,  power.10000[,1], xlab="Number of controls and cases in validation data", ylab = "Power", type = 'l',
     lwd = 3, cex.axis = 2, cex.lab = 2, col = cols[1], ylim = c(0,1), xlim = c(0,5000))
mtext('C', side = 3, line = 1, at = -3, cex = 2, outer =F)
for (i in seq(2, length(AUC), 1)){
  lines(sample.sizes.validation, power.10000[,i], col = cols[i], lwd = 3)
}
lines(c(0, 20000), c(0.9, 0.9), lwd = 2.5, lty = 2)
legend("bottomright", c('0.8','0.85','0.9','0.95','0.975','0.99', 'Power level of 0.9'), col = c(cols,'black'), lwd = c(rep(3,7),2),
       lty = c(rep(1,7),2),title = "AUC-ROC")

# panel D

min_serosurvey <- matrix(0, nrow = length(sample.sizes.validation), ncol = length(AUC))

for ( i in seq(1, length(sample.sizes.validation),1)){
  power.now <- power(power.data, sample.sizes.serosurvey, AUC, validation.size=sample.sizes.validation[i])
  for ( j in seq(1, length(AUC),1)){
    min_serosurvey[i,j] <- sample.sizes.serosurvey[power.now[,j]>0.9 ][1]
  }
}

# make plot
plot(sample.sizes.validation,  min_serosurvey[,3], xlab="Number of controls and cases in validation data", ylab = "Minimal required number of ind. in serosurvey", type = 'l',
     lwd = 3, cex.axis = 2, cex.lab = 2, col = cols[3], ylim =c(min(min_serosurvey, na.rm = T), max(min_serosurvey, na.rm = T)), xlim = c(0,5000))
mtext('D', side = 3, line = 1, at = -3, cex = 2, outer =F)
for (i in seq(4, length(AUC), 1)){
  lines(sample.sizes.validation, min_serosurvey[,i], col = cols[i], lwd = 3)
}
legend("topright", c(as.character(AUC[3:7])), col = c(cols[3:7]), lwd = c(rep(3,5)),
       lty = c(rep(1,5)),title = "AUC-ROC")

ptwo <- recordPlot()

# now add the first ggplot element which is going to take
# up the left two quadrants
vp <- viewport(height = unit(1,"npc"), width=unit(0.5, "npc"),
               just = c("left","top" ),  y = 1, x = 0.5)

# add the second ggplot element in the bottom right quadrant
vp_B <- viewport(height = unit(0.48,"npc"), width=unit(0.45, "npc"), just = c("left","top"), y = 1, x = 0.48)
vp_A <- viewport(height = unit(0.48,"npc"), width=unit(0.45, "npc"), just = c("left","top"), y = 1, x = 0.02)

pdf(file = '/Users/judith/testing/Figures/Revised/Power_plot_validation_size.pdf', width = 16, height = 16)
print(ptwo, vp = vp)
print(P3_B, vp = vp_B)
print(P3_A, vp = vp_A)
dev.off()


#####################################################
# Figure 5 -- Concept of asymptomatic/severe cases
#####################################################

pdf("/Users/judith/testing/Figures/Revised/detect-discrepancy-concept.pdf", height=7, width=11)

n <- 1e4
cases <- rgamma(n,40,5)
milds <- rgamma(n,20,5)
controls <- rgamma(n,1,1)

par(mfcol=c(2,2)) -> op
hist(controls, col=rgb(0,0,0,alpha=1), breaks=30, xlim=c(0,15), ylim=c(0,n/2), freq=T, xlab="Quantitative serological measure",
     main = '')
title("Test validation data with severe cases and controls", line = -0.5)
hist(cases, col=rgb(1,0,0,alpha=0.3), breaks=30, add=T)
mtext('A', side = 3, line = 1, at = c(0), cex = 2, outer =F, font = 2)
legend(12,n/4,legend=c("controls","severe cases"), col=c(rgb(0,0,0,alpha=1), rgb(1,0,0,alpha=0.3)), pch=16, bty="n")
##
hist(c(controls,cases[1:round(0.15*n)]), col=rgb(0,0,0,alpha=0.3), breaks=30, xlim=c(0,15), ylim=c(0,n/2),
     freq=T, main="", xlab="Quantitative serological measure")
title("Serosurvey data with 15% severe cases" , line = -0.5)
mtext('B', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
#hist(controls, col=rgb(0,0,0,alpha=0.3), breaks=30, add=T)
##
##
hist(controls, col=rgb(0,0,0,alpha=1), breaks=30, xlim=c(0,15), ylim=c(0,n/2), freq=T,
     main="", xlab="Quantitative serological measure")
hist(cases, col=rgb(1,0,0,alpha=0.3), breaks=30, add=T)
hist(milds, col=rgb(0,0,1,alpha=0.3), breaks=30, add=T)
hist(controls, col=rgb(0,0,0,alpha=0.3), breaks=30, add=T)
title("Test validation data with severe cases and controls", line = -0.5)
mtext('C', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
legend(10,n/4,legend=c("controls","severe cases","asymptomatic cases"), col=c(rgb(0,0,0,alpha=1), rgb(1,0,0,alpha=0.3), rgb(0,0,1,alpha=0.3)), pch=16, bty="n")
##
hist(c(controls, milds[1:round(0.10*n)], cases[1:round(0.05*n)]), col=rgb(0,0,0,alpha=0.3), breaks=30, xlim=c(0,15), ylim=c(1,n/2), freq=T,
     main="", xlab="Quantitative serological measure")
mtext('D', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
title("Serosurvey data with 10% asym. and 5% severe cases", line = -0.5)
#hist(cases[1:round(0.05*n)], col=rgb(1,0,0,alpha=0.3), breaks=30, add=T)
#hist(milds[1:round(0.10*n)], col=rgb(0,0,1,alpha=0.3), breaks=30, add=T)
par(op)
rm(op,cases,controls,milds,n)

dev.off()

#####################################################
# Figure 6 -- Extended method for estimating asymptomatic cases
#####################################################

if (code){
  N = 50
  # Take the distribution of the severe case sera such that the AUC-ROC for the severe cases is 0.975
  mean.S <- 12 # to have an AUC-ROC between control and severe case sera to be 1
  #all <- determine.means(desired.AUC = c(0.7, 0.75,0.8,0.85,0.9,0.95), shape.controls = 12, rate.controls = 1, min.mean = 1, max.mean = 12 )

  means.asym <- c(9.575150, 8.935872, 8.252505, 7.480962, 6.599198, 5.364729) # all[[2]]
  AUC.asym <- c(0.6991440, 0.7494334, 0.7997224, 0.8506663, 0.8996441, 0.9499151) # all[[1]]

  # Panel A
  results.MS.l.b <- matrix(0, nrow = length(means.asym), ncol = N)
  results.MS.l.s <- matrix(0, nrow = length(means.asym), ncol = N)

  true.prev <- 0.08
  p.M <- 0.2

  for (j in seq(1,length(means.asym),1)){

    for (i in seq(1,N,1)){
      # Simulate data from both
      data.serosurvey <- simulate.serosurvey.MS(N_tests = 1e4, mean.M = means.asym[j], mean.S = mean.S, prev = true.prev, p.M = p.M)

      # Simulate validation data
      control.data <- rgamma(10000, shape = 1, scale = 1)
      case.data <- rgamma(10000, shape = mean.S, scale = 1)
      get.asym <- sum( rbinom(10000,size = 1, prob = p.M))
      case.data.b <- c(rgamma(get.asym, shape = means.asym[j], scale = 1), rgamma(10000-get.asym, shape = mean.S, scale = 1))

      ############### likelihood-based method #######################
      # simulate validation data and create pdf for both distributions from there
      pdfs <- sim.validation.data(control.data = control.data, case.data = case.data )

      pdfs.b <- sim.validation.data(control.data = control.data, case.data = case.data.b )

      # analyze with both distributions
      #likelihood
      results.MS.l.b[j,i] <- analyze.serosurvey.AS.emp(data.serosurvey = data.serosurvey, prev.mild = p.M*0.1, pdf.control.smooth=pdfs.b[[2]], pdf.case.smooth=pdfs.b[[1]],
                                                       type = 'likelihood', mean.M = means.asym[j], threshold = NULL, TNR = NULL, TPR = NULL, correct = TRUE, complex = FALSE)[1]

      # test for bias when we only analyze with severe cases
      # likelihood
      results.MS.l.s[j,i] <- analyze.serosurvey.AS.emp(data.serosurvey = data.serosurvey, prev.mild = p.M*0.1, pdf.control.smooth=pdfs[[2]], pdf.case.smooth=pdfs[[1]],
                                                       type = 'likelihood', mean.M = means.asym[j], threshold = NULL, TNR = NULL, TPR = NULL, correct = TRUE, complex = FALSE)[1]
    }
  }

  # PART 2
  # Can we, given the distribution of the survey results say if there were to few/no mild cases tested.
  # more general, is the validation set of individuals representative of the total variation in cases?

  methods <- c("Mixture model", "likelihood Complex")

  # rethink on how to code this, because likelihood has changed... Use combined severe and mild distribution?
  results.MS.complex <- matrix(0, nrow = length(means.asym), ncol = N)
  results.MS.complex.pM <- matrix(0, nrow = length(means.asym), ncol = N)
  results.MS.complex.mM <- matrix(0, nrow = length(means.asym), ncol = N)

  ll.values.complex <- matrix(0, nrow = length(means.asym), ncol = N)

  # do this for varying levels of mean.M --> determine various levels of overlap I would like to show.
  for (j in seq(1, length(means.asym),1)){
    for (i in seq(1,N,1)){
      # Simulate data from both
      data.serosurvey <- simulate.serosurvey.MS(N_tests = 1e4, mean.M = means.asym[j], mean.S = mean.S, prev = true.prev, p.M = p.M)

      # Simulate validation data
      control.data <- rgamma(5000, shape = 1, scale = 1)
      case.data <- rgamma(5000, shape = mean.S, scale = 1)

      # simulate validation data and create pdf for both distributions from there
      pdfs <- sim.validation.data(control.data = control.data, case.data = case.data )

      # analyze with both distributions
      # likelihood
      result2 <- tryCatch( {result2 <- analyze.serosurvey.AS.emp(data.serosurvey = data.serosurvey, prev.mild = 0.4*0.1, pdf.control.smooth=pdfs[[2]],
                                                                 pdf.case.smooth=pdfs[[1]],
                                                                 type = 'likelihood', mean.M = means.asym[j], threshold = NULL, TNR = NULL, TPR = NULL, correct = TRUE,
                                                                 complex = TRUE, case.data=case.data)},
                           error = function(e){i <<- i-1
                           return(rep(NA,7))})

      results.MS.complex[j,i] <- result2[1]
      results.MS.complex.pM[j,i] <- result2[2]
      results.MS.complex.mM[j,i] <- result2[3]
      ll.values.complex[j,i] <- result2[7]
    }
  }

  df.results.MS.complex <-  matrix(0, nrow = length(means.asym)*N, ncol = 3)
  df.results.MS.complex[,3] <- "Total cum. inc."

  df.results.MS.complex.pM <-  matrix(0, nrow = length(means.asym)*N, ncol = 3)
  df.results.MS.complex.pM[,3] <- "Cum. inc. of asym. cases"

  df.results.MS.complex.pS <-  matrix(0, nrow = length(means.asym)*N, ncol = 3)
  df.results.MS.complex.pS[,3] <- "Cum. inc. of severe cases"

  df.results.MS.l.s <-  matrix(0, nrow = length(means.asym)*N, ncol = 3)
  df.results.MS.l.s[,3] <- "Mixture model -- only sev. dist."
  df.results.MS.l.b <-  matrix(0, nrow = length(means.asym)*N, ncol = 3)
  df.results.MS.l.b[,3] <- "Mixture model -- both dist."

  for (i in seq(1, length(means.asym), 1)){
    for (j in seq(1, N,1 )){
      df.results.MS.complex[ ((i-1)*N + j) ,1] <- results.MS.complex[i,j]
      df.results.MS.complex[ ((i-1)*N + j) ,2] <- round( AUC.asym[i] , 2)

      df.results.MS.complex.pM[ ((i-1)*N + j) ,1] <- results.MS.complex.pM[i,j]
      df.results.MS.complex.pM[ ((i-1)*N + j) ,2] <- round( AUC.asym[i] , 2)

      df.results.MS.complex.pS[ ((i-1)*N + j) ,1] <- results.MS.complex[i,j] - results.MS.complex.pM[i,j]
      df.results.MS.complex.pS[ ((i-1)*N + j) ,2] <- round( AUC.asym[i] , 2)

      df.results.MS.l.s[ ((i-1)*N + j) ,1] <- results.MS.l.s[i,j]
      df.results.MS.l.s[ ((i-1)*N + j) ,2] <- round(AUC.asym[i],2)

      df.results.MS.l.b[ ((i-1)*N + j) ,1] <- results.MS.l.b[i,j]
      df.results.MS.l.b[ ((i-1)*N + j) ,2] <- round(AUC.asym[i],2)
    }
  }

  df.plot.1 <- as.data.frame(rbind(df.results.MS.l.b, df.results.MS.l.s))
  colnames( df.plot.1  ) <- c('Value','means.asym', 'Type')
  df.plot.1[,2] <- factor(df.plot.1[,2])
  df.plot.1[,3] <- as.factor(df.plot.1[,3])
  df.plot.1[,1] <- as.numeric(paste((df.plot.1[,1])))

  df.plot.2 <- as.data.frame(rbind(df.results.MS.complex,df.results.MS.complex.pM, df.results.MS.complex.pS))
  colnames( df.plot.2  ) <- c('Value','means.asym', 'Type')
  df.plot.2[,2] <- factor(df.plot.2[,2])
  df.plot.2[,3] <- as.factor(df.plot.2[,3])
  df.plot.2[,1] <- as.numeric(paste((df.plot.2[,1])))

  write.table(df.plot.1, 'data_fig_6A.txt' )
  write.table(df.plot.2, 'data_fig_6B.txt' )
}

df.plot.1<- as.data.frame(read.table('data_fig_6A.txt' ))
df.plot.2 <- as.data.frame(read.table('data_fig_6B.txt' ))

df.plot.1$means.asym <- as.factor(df.plot.1$means.asym)
df.plot.2$means.asym <- as.factor(df.plot.2$means.asym)

# plot 1

P5_A <- ggplot(data=df.plot.1) +
  geom_violin( aes(x = means.asym, y = Value, fill = Type, colour = Type), width = 1, alpha = .2, size=0.3)+
  stat_summary(fun=median, aes(y = Value, x = means.asym, colour = Type ), geom="point", shape=19, size=2, position = position_dodge(1)) +
  theme(axis.text = element_text(size=20),axis.title = element_text(size=20), text = element_text(size=20),legend.position="bottom",
        panel.grid.major =element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0))+
  geom_hline(aes(yintercept = true.prev, col = "True cum. inc."), size = 1 , lty = 2)+
  ggtitle('A')+
  ylab("Estimated cumulative incidence")+
  xlab("AUC-ROC between asym. & sev. case distributions")+
  scale_fill_manual("", values=c("True cum. inc." = "red", "Mixture model -- both dist." = "mediumpurple" ,
                                 "Mixture model -- only sev. dist."="cyan"),
                    aesthetics = c("fill"), breaks =c("Mixture model -- both dist." ,  "Mixture model -- only sev. dist."))+
  scale_colour_manual("", values=c("True cum. inc." = "red", "Mixture model -- both dist." = "mediumpurple" ,
                                   "Mixture model -- only sev. dist."="cyan"),
                      aesthetics = c("colour"), breaks = c("True cum. inc."))

# plot 2
P5_B <- ggplot(data=df.plot.2) +
  ylim(0,0.13)+
  geom_violin( aes(x = means.asym, y = Value, fill = Type, col=Type), width = 1, alpha = .2, size=0.3)+
  stat_summary(fun=median, aes(y = Value, x = means.asym, colour = Type), geom="point", shape=19, size=2, position = position_dodge(1)) +
  theme(axis.text = element_text(size=20),axis.title = element_text(size=20), text = element_text(size=20),legend.position="bottom",
        panel.grid.major =element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0))+
  geom_hline(aes(yintercept = true.prev*p.M, color ="Cum. inc. of asym. cases"), size = 1 )+
  geom_hline(aes(yintercept = true.prev*(1-p.M), color = "Cum. inc. of severe cases"), size = 1 )+
  geom_hline(aes(yintercept = true.prev), color = "orange", size = 1 )+
  ggtitle('B')+
  xlab("AUC-ROC between asym. & sev. case distributions")+
  ylab("Estimated cumulative incidence")+
  scale_colour_manual("", values = c("Cum. inc. of asym. cases"="orchid1",
                                     "Cum. inc. of severe cases"="blueviolet",
                                     "Total cum. inc." = "orange"),
                      aesthetics = c("colour","fill"))


pdf(file = '/Users/judith/testing/Figures/Revised/Asym_severe.pdf', width = 10.5, height = 16)
grid.arrange( arrangeGrob(P5_A, P5_B, nrow = 2, ncol = 1) )
dev.off()



#############################################################
# Figure 7 -- Conceptual Figure
#############################################################

pdf(file = '/Users/judith/testing/Figures/Revised/Conceptual_figure.pdf', width = 16, height = 16)

means.figure <- means[c(4,5,6)]
par(mfrow=c(2,2),oma=c(1.5,2.5,1,4),mar=c(4,5,2,2))

titers <- seq(0,20,0.01)
titer.TN <- cbind(titers, dgamma(titers, shape = 1, rate = 1)/sum(dgamma(titers, shape = 1, rate = 1)))

# part A
plot(titers[1:1000], titer.TN[1:1000,2], col = "grey42", type = 'l', lwd = 5, cex.lab = 2, cex.axis = 2,
     xlab = "Quantitative test measures", ylab = 'Frequency')

cols <- c("springgreen3", "dodgerblue", "orange")

for (i in seq(1, length(means.figure),1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[i], rate = 1)/sum(dgamma(titers, shape = means.figure[i], rate = 1)))
  lines(titers[1:1000], titer.TP[1:1000,2],col = cols[i], lwd = 5)
}
mtext('A', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
legend('topright', c("Controls", "Cases -  AUC-ROC = 0.85", "Cases -  AUC-ROC = 0.9", "Cases -  AUC-ROC = 0.95")
       , col = c("grey42", cols), lwd = 5, cex = 1.5)


# part B

plot(c(0,1), c(0,1), type = 'l', lty = 2, lwd = 2.5,
     xlab = "False positive rate", ylab = "True positive rate", cex.lab=2, cex.axis =2)
mtext('B', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)

for (j in seq(1, length(means.figure),1)){
  titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[j], rate = 1)/sum(dgamma(titers, shape = means.figure[j], rate = 1)))

  n <- length(titer.TN[,1])

  FPR = rep(0, n)
  TPR = rep(0, n)

  for (i in seq(1,n,1)){
    TN <- sum(titer.TN[1:i,2])
    TP <- sum(titer.TP[(i:length(titer.TN[,1])),2])
    FN <- sum(titer.TP[1:i,2])
    FP <- sum(titer.TN[(i:length(titer.TN[,1])),2])
    TPR[i] <- TP/(TP+FN)
    FPR[i] <- FP/(TN+FP)
  }

  TPR <- rev(TPR)
  FPR <- rev(FPR)

  lines(FPR, TPR, col = cols[j], lwd = 5)
}

# Part C
plot(titers[1:1000], titer.TN[1:1000,2], col = "grey42", type = 'l', lwd = 5, cex.lab = 2, cex.axis = 2,
     xlab = "Quantitative test measures", ylab = 'Frequency')
mtext('C', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[3], rate = 1)/sum( dgamma(titers, shape = means.figure[3], rate = 1)))
lines(titers[1:1000], titer.TP[1:1000,2], col = cols[3], lwd = 5)

# add threshold values:

youden <- classical.threshold(mean = means.figure[3], method = 'youden')[1]
specificity <- classical.threshold(mean = means.figure[3], method = 'specificity')[1]

linecolors = c("turquoise3","maroon")

lines(rep(youden,2) , c(0,0.01), lty = 3, lwd = 3, col = linecolors[1])
lines(rep(specificity,2) , c(0,0.01), lty = 4, lwd = 3, col = linecolors[2])

legend('topright', c("Controls", "Cases -  AUC-ROC = 0.95", "Max Youden", "High specificity")
       , col = c("grey42", cols[3], linecolors), lwd = c(5,5,3,3,3,3) , cex = 1.5, lty=c(1,1,3,4,5,6))

# Part D

plot(c(0,1), c(0,1), type = 'l', lty = 2, lwd = 2.5,
     xlab = "False positive rate", ylab = "True positive rate", cex.lab=2, cex.axis =2)
mtext('D', side = 3, line = 1, at = 0, cex = 2, outer =F, font = 2)
titer.TP <- cbind(titers, dgamma(titers, shape = means.figure[3], rate = 1))

n <- length(titer.TN[,1])

FPR = rep(0, n)
TPR = rep(0, n)

for (i in seq(1,n,1)){
  TN <- sum(titer.TN[1:i,2])
  TP <- sum(titer.TP[(i:length(titer.TN[,1])),2])
  FN <- sum(titer.TP[1:i,2])
  FP <- sum(titer.TN[(i:length(titer.TN[,1])),2])
  TPR[i] <- TP/(TP+FN)
  FPR[i] <- FP/(TN+FP)
}

TPR <- rev(TPR)
FPR <- rev(FPR)

lines(FPR, TPR, col = cols[3], lwd = 5)

# add points for the various thresholds
k <- seq(1,n,1)[abs(titers-youden) == min(abs(titers-youden))]
TN <- sum(titer.TN[1:k,2])
TP <- sum(titer.TP[(k:length(titer.TN[,1])),2])
FN <- sum(titer.TP[1:k,2])
FP <- sum(titer.TN[(k:length(titer.TN[,1])),2])
TPR.youden <- TP/(TP+FN)
FPR.youden <- FP/(TN+FP)

points(FPR.youden, TPR.youden, lwd = 5, col = linecolors[1])

k <- seq(1,n,1)[abs(titers-specificity) == min(abs(titers-specificity))]
TN <- sum(titer.TN[1:k,2])
TP <- sum(titer.TP[(k:length(titer.TN[,1])),2])
FN <- sum(titer.TP[1:k,2])
FP <- sum(titer.TN[(k:length(titer.TN[,1])),2])
TPR.specificity <- TP/(TP+FN)
FPR.specificity <- FP/(TN+FP)

points(FPR.specificity, TPR.specificity, lwd = 5, col = linecolors[2])

legend('bottomright', c("Max Youden", "High specificity")
       , col = c(linecolors), lwd = c(5,5) , cex = 1.5, lty=c(0,0), pch=c(1,1))

dev.off()

#####################################################
# Figure S1 -- Sensitivity and Specificity corresponding to the AUC-ROC range
#####################################################

# Plot for sensitivity/specificity
means = c(  2.322645,  2.719439,  3.314629,  4.306613,  5.320641, 6.643287)
AUCs = c( 0.8, 0.85, 0.9, 0.95, 0.975,  0.99)

sensitivities <- matrix(0, nrow = 2, ncol = length(means))
specificities <- matrix(0, nrow = 2, ncol = length(means))

for (i in seq(1, length(means), 1)){
  control.data <- rgamma(10000, shape = 1, scale = 1)
  case.data <- rgamma(10000, shape = means[i], scale = 1)
  max.youden.emp <- observed.threshold(control.data, case.data, method = 'youden')
  high.sens.emp <- observed.threshold(control.data, case.data, method = 'specificity')

  sensitivities[1,i] <- max.youden.emp[[2]]
  specificities[1,i] <- max.youden.emp[[3]]

  sensitivities[2,i] <- high.sens.emp[[2]]
  specificities[2,i] <- high.sens.emp[[3]]
}

pdf(file = '/Users/judith/testing/Figures/Revised/sens_spec.pdf', width = 10, height = 10)
par(mfrow=c(1,1), mar=c(5.1,5.5,4.1,2.1))
plot(AUCs, sensitivities[1,], ylab = 'Sensitivity/Specificity', xlab='Area under ROC-curve (AUC-ROC)', xaxt = "n",
     cex = 2, cex.axis = 2, cex.lab = 2, col = 'blue', pch = 16, ylim = c(0,1))
points(AUCs, sensitivities[2,], col = 'red', pch = 16 , cex = 2)
points(AUCs, specificities[1,], ylab = 'Sensitivity/Specificity', xlab='Area under ROC-curve (AUC-ROC)', xaxt = "n",
       cex = 2, cex.axis = 2, cex.lab = 2, col = 'blue', pch = 5, ylim = c(0,1))
points(AUCs, specificities[2,], col = 'red', pch =5 , cex = 2)
axis(1, at=AUCs, labels=AUCs, cex = 2,cex = 2, cex.axis = 2, cex.lab = 2)
legend('bottomright', col = c('black', 'black', 'blue', 'red'), pch = c(16, 5, NA,NA ), lty = c(NA, NA, 1,1), lwd =c(NA, NA, 3,3),
       legend = c('Sensitivity', 'Specificity', 'Max Youden', 'High Specificity'))
dev.off()


#####################################################
# Figure S2 -- Overview of uncorrected methods
#####################################################

if (code){

}

est_NC <- as.data.frame(read.table('data_S1.txt'))
est_NC$prev2 <- paste0("True value = ",est_NC$`True.prevalence`*100,"%")
est_NC$AUC <- as.factor(est_NC$AUC)

g1_NC = ggplot(est_NC) +
  ggtitle('A')+
  geom_violin(aes(x=AUC,y=value,colour=type, fill = type), width = 1, alpha = .2, size=0.7, position = position_dodge(0.8) ) +
  stat_summary(fun=median, aes(y = value, x = AUC, fill = type, colour = type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept=True.prevalence),linetype=2,size=.5,alpha=.6) +
  facet_grid(prev2~., scales = "free_y") +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  scale_colour_manual('',values=c("High specificity - uncorrected" =  "dodgerblue",
                                  "Max Youden - uncorrected" =  "maroon"), aesthetics = c('colour',"fill")) +
  labs(x="AUC-ROC",y="Estimated cumulative incidence") +
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), text = element_text(size=20),
        legend.position="bottom", title = element_text(size = 26, face='bold'),
        panel.grid.major =element_blank(), panel.background = element_blank())+
  guides(colour=guide_legend(nrow=2,ncol = 3, byrow=TRUE))
g2_NC = ggplot(est_NC) +
  ggtitle('B')+
  geom_violin(aes(x=AUC,y=size.interval,colour=type, fill = type), width = 1, alpha = .2, size=0.7, position = position_dodge(0.8) ) +
  stat_summary(fun=median, aes(y = size.interval, x = AUC, fill = type, colour = type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  #geom_line(aes(x=AUC,y=half_95CI,group=type,colour=type), size = 1.5) +
  facet_grid(prev2~.) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0)) +
  scale_colour_manual('',values=c("High specificity - uncorrected" =  "dodgerblue",
                                  "Max Youden - uncorrected" =  "maroon"), aesthetics = c('colour',"fill")) +
  #coord_cartesian(ylim=c(0,0.055)) +
  labs(x="AUC-ROC",y="Size of 95% confidence interval")+
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), text = element_text(size=20),
        legend.position="bottom", title = element_text(size = 26, face='bold'),
        panel.grid.major =element_blank(), panel.background = element_blank())+
  guides(colour=guide_legend(nrow=2,ncol = 3, byrow=TRUE))
mylegend<-g_legend(g1_NC)
#plot_grid(g1+theme(legend.position = "none"),g2+theme(legend.position = "none"),ncol=2,nrow=1,rel_widths=c(1,1),mylegend)

pdf(file = '/Users/judith/testing/Figures/Revised/Overview_all_methods_new_NC.pdf', width = 16, height = 18)
all_NC <- grid.arrange(arrangeGrob(g1_NC + theme(legend.position="none"),
                                   g2_NC + theme(legend.position="none"),
                                   nrow=1), mylegend, nrow = 2,heights=c(10,1))
dev.off()

#####################################################
# Figure S3 -- Comparing of bootstrap intevals and variancen in point estimates
#####################################################

ll.all.bootstrap <- create.figure.bootstrap(means = means, AUC = AUCs )

pdf(file = '/Users/judith/testing/Figures/Revised/comparison_bootstrap_MM.pdf', width = 16, height = 16)
ggplot(data=ll.all.bootstrap) +
  geom_violin( aes(x = AUC, y = Value, fill=Type,colour=Type), width = 1, alpha = .2, size=0.3, position = position_dodge(0.8) )+
  theme(axis.text = element_text(size=26),axis.title = element_text(size=26), text = element_text(size=26),legend.position="bottom",
        panel.grid.major =element_line(colour = "grey"), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(), legend.background = element_blank(),plot.title = element_text(hjust = 0.5))+
  stat_summary(fun=median, aes(y = Value, x = AUC, fill = Type, colour = Type), geom="point", shape=19, size=2,
               position = position_dodge(width = 0.8)) +
  ylab("Deviation from mean cum. inc. est.") + xlab("Area under ROC-curve (AUC - ROC)") +
  scale_colour_manual("", values = c("High specificity bootstrap"="violetred3","High specificity"="orange","Max Youden bootstrap" = "dodgerblue" ,
                                     "Max Youden" = "darkviolet","Mixture model bootstrap"="springgreen3" ,"Mixture model"="cyan"),
                      aesthetics = c("colour","fill"))
dev.off()

############# TO DO ################

# add final functions to basic functions file
# make saving and retreating data work
# Add R-package description
# test code
# check figrure S3




