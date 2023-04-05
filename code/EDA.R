#### setup ####
rm(list=ls())
#library(tidyverse)
library(dplyr)
library(ggplot2)

## Set working directory to the project folder 
setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Q2 2021-22 Spring/ADAP/code/LCA-ADAP/code/poLCA-master/poLCA-master/")
## Go to R terminal (next to console); type R CMD SHLIB ./src/poLCA.c

## Back in console type: You should not need this step as my folder already contains it
## Please let me know if this causing issues

dyn.load("src/poLCA.dll")

## Source all functions from R folder
setwd("./R/")
files.sources = list.files()
sapply(files.sources, source)

# set you working directory to location where main_script and data is saved
setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Q2 2021-22 Spring/ADAP/code/LCA-ADAP/")

# Load all the data from the folder
survey_dat <- read.csv('data/2011-2018_nhats_data_lca_sub.csv') 

#### Cleaning self-report-sample ####
self_report_sample <- survey_dat %>% filter(proxyres5_sp==0,dementia5_cd==1)
self_report_sample <- self_report_sample[complete.cases(self_report_sample), ]

#reordering values to start from 1 to max value instead of 0
self_report_sample$phq5_oad <- self_report_sample$phq5_oad+1
self_report_sample$gad5_oad <- self_report_sample$gad5_oad+1 
self_report_sample$exercise5_sp <- self_report_sample$exercise5_sp+1
self_report_sample$infse_sp5_con <- self_report_sample$infse_sp5_con+1
self_report_sample$fse_sp5_con <- self_report_sample$fse_sp5_con+1
names(self_report_sample)[names(self_report_sample) == 'w5varstrat'] <- 'strata'
names(self_report_sample)[names(self_report_sample) == 'w5varunit'] <- 'cluster'

# Select relevant columns
self_report_outcomes <- self_report_sample %>% dplyr::select(ordfallyr15,
                                                             fof5_cat, 
                                                             slfhealth5_sp, 
                                                             pain5_sp, 
                                                             phq5_oad, 
                                                             gad5_oad,
                                                             exercise5_sp,
                                                             infse_sp5_con,
                                                             fse_sp5_con,
                                                             fatigue_sp5,
                                                             strata,cluster)
# Summary for table/ Figure 1
self_outcomes_summary <- apply(self_report_outcomes, 2, table)
self_outcomes_perc <- apply(self_report_outcomes,2, function(x) round(table(x)/length(x)*100,1))

# Normalize weight to sum to number of observations
w_self <- self_report_sample$w5anfinwgt0/sum(self_report_sample$w5anfinwgt0)*nrow(self_report_sample)

y_self <- self_report_outcomes

w_self_df <- as.data.frame(w_self)
ggplot(w_self_df, aes(x=w_self)) + geom_histogram(bins = 50)+theme_bw() + 
  xlab("sampling weights")


#### Self report analysis ####
f <- cbind(ordfallyr15, fof5_cat, slfhealth5_sp,pain5_sp,phq5_oad,
           gad5_oad,exercise5_sp,infse_sp5_con,fse_sp5_con,fatigue_sp5) ~ 1

BIC.best.self <- Inf
AIC.k.self <- matrix(0, nrow = 1, ncol = 6)
BIC.k.self <- matrix(0, nrow = 1, ncol = 6)
# Loop to run over different number of clusters
for(k in 1:6){
  mod.k <- poLCA.ordinal(f, data = y_self,weights = w_self, nclass = k, 
                         nrep = 20, calc.se = FALSE,
                         maxiter = 3000, tol = 1e-5,
                         lt = 0.001,seed=1234)
  AIC.k.self[1,k] <- mod.k$aic
  BIC.k.self[1,k] <- mod.k$bic
  # If increasing #clusters decreases BIC by more than 0.2% then 
  # accept new #clusters and corresponding model as the best model
  if(((mod.k$bic/BIC.best.self-1)*100 < -0.2) | k==1){
    BIC.best.self <- mod.k$bic
    mod.best.self <- mod.k
    k.best.self <- k
  }
}

# Plot AIC and BIC 
plot(AIC.k.self[1, ], main = "AIC vs #Classes", xlab = "Classes", ylab = "AIC", type = "b")
plot(BIC.k.self[1, ], main = "BIC vs #Classes", xlab = "Classes", ylab = "BIC", type = "b")


#### Computationally Intensive se calculation for self-report ####
# Calculate se for the best model fit 
mod.best.self <- poLCA.ordinal(f, data = y_self,weights = w_self, 
                           nclass = k.best.self, 
                           nrep = 20, calc.se = TRUE,
                           maxiter = 3000, tol = 1e-5,
                           lt = 0.001)

# Doing jittered plots faceted by predicted latent class
library(tidyr)
y.wide <- cbind(self_report_sample$spid,y_self)
var.names <- c("id","Fallhist","Fearfall","health","pains","Depression",
               "Anxiety","Exercise","InformalSoc","FormalSoc","Fatigue","strata","cluster")
colnames(y.wide) <- var.names
y.long <- gather(y.wide, "Variable","Response", -id )

y.wide.ord <- cbind(y.wide,mod.best.self$predclass)
colnames(y.wide.ord)[ncol(y.wide.ord)] <- "Predclass"
y.long.ord <- gather(y.wide.ord, "Variable","Response",Fallhist:Fatigue)
y.long.ord$Variable <- factor(y.long.ord$Variable, 
                              levels = c("Fallhist","Fearfall","health",
                                         "pains","Fatigue","Depression",
                                         "Anxiety","Exercise","InformalSoc",
                                         "FormalSoc"))
p <- ggplot(data = y.long.ord, aes(x = Variable, y = Response, group = id))
p + 
  geom_line(alpha=1/10, position = position_jitter(width = 0, height = .1)) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  facet_grid(. ~ Predclass)


# Plot for class conditioned expected score
class.cond.exp.score <- matrix(NA,ncol = length(mod.best.self$probs),
                               nrow = k.best.self)
colnames(class.cond.exp.score) <- var.names[2:11]

rownames(class.cond.exp.score) <- paste("Class ",1:k.best.self," (",round(mod.best.self$P*100,0),"%)",sep="")

for(i in 1:ncol(class.cond.exp.score)){
  max.outcome <- ncol(mod.best.self$probs[[i]])
  class.cond.exp.score[,i] <- (mod.best.self$probs[[i]] %*% 
                                  (1:max.outcome))/max.outcome*100
}

df_cces <- as.data.frame(class.cond.exp.score)
df_cces <- df_cces%>% mutate(Class=rownames(df_cces))

df_cces_melt <- reshape2::melt(df_cces,id.vars="Class")
df_cces_melt$variable <- factor(df_cces_melt$variable, 
                                levels=c("Fallhist","Fearfall",
                                         "health","pains","Fatigue",
                                         "Depression","Anxiety",
                                         "Exercise","InformalSoc",
                                         "FormalSoc"))

ggplot(df_cces_melt, aes(x = variable, y = value)) + 
  geom_line(aes(color = Class, group = Class),size=1.2)+
  scale_color_discrete(name = "Class", labels = c("Class 1: High FRO and poor QoL",
                                                 "Class 2: Low FRO and good QoL", 
                                                 "Class 3: Intermediate FRO and QoL"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, 
                                   hjust=1), 
        legend.text=element_text(size=7))+
  xlab("Variable")+
  ylab("Class conditioned expected score (%)")


#### Cleaning proxy-report-sample ####
proxy_report_sample <- survey_dat %>% filter(proxyres5_sp==1, dementia5_cd==1)
proxy_report_sample <- proxy_report_sample[ ,-8]
proxy_report_sample <- proxy_report_sample[complete.cases(proxy_report_sample), ]

#reordering values to start from 1 to max value instead of 0
proxy_report_sample$phq5_oad <- proxy_report_sample$phq5_oad+1
proxy_report_sample$gad5_oad <- proxy_report_sample$gad5_oad+1 
proxy_report_sample$exercise5_sp <- proxy_report_sample$exercise5_sp+1
proxy_report_sample$infse_sp5_con <- proxy_report_sample$infse_sp5_con+1
proxy_report_sample$fse_sp5_con <- proxy_report_sample$fse_sp5_con+1
names(proxy_report_sample)[names(proxy_report_sample) == 'w5varstrat'] <- 'strata'
names(proxy_report_sample)[names(proxy_report_sample) == 'w5varunit'] <- 'cluster'

# Select relevant columns
proxy_report_outcomes <- proxy_report_sample %>% dplyr::select(ordfallyr15,
                                                      fof5_cat, 
                                                      slfhealth5_sp, 
                                                      pain5_sp, 
                                                      phq5_oad, 
                                                      gad5_oad,
                                                      exercise5_sp,
                                                      infse_sp5_con,
                                                      fse_sp5_con,
                                                      fatigue_sp5,
                                                      strata,cluster)
# Summary for table/ Figure 1
proxy_outcomes_summary <- apply(proxy_report_outcomes, 2, table)
proxy_outcomes_perc <- apply(proxy_report_outcomes,2, function(x) round(table(x)/length(x)*100,1))

# Normalize weight to sum to number of observations
w_proxy <- proxy_report_sample$w5anfinwgt0/sum(proxy_report_sample$w5anfinwgt0)*nrow(proxy_report_sample)
y_proxy <- proxy_report_outcomes


#### Creating joint-report-sample ####
joint_report_sample <- rbind(self_report_sample[ ,-8],proxy_report_sample)
joint_report_outcomes <- joint_report_sample %>% dplyr::select(ordfallyr15,
                                                               fof5_cat, 
                                                               slfhealth5_sp, 
                                                               pain5_sp, 
                                                               phq5_oad, 
                                                               gad5_oad,
                                                               exercise5_sp,
                                                               infse_sp5_con,
                                                               fse_sp5_con,
                                                               fatigue_sp5,
                                                               strata,cluster)
w_joint <- joint_report_sample$w5anfinwgt0/sum(joint_report_sample$w5anfinwgt0)*nrow(joint_report_sample)
y_joint <- joint_report_outcomes

#### Self+ proxy joint report analysis ####
f <- cbind(ordfallyr15, fof5_cat, slfhealth5_sp,pain5_sp,phq5_oad,
           gad5_oad,exercise5_sp,infse_sp5_con,fse_sp5_con,fatigue_sp5) ~ 1

BIC.best.joint <- Inf
AIC.k.joint <- matrix(NA, nrow = 1, ncol = 6)
BIC.k.joint <- matrix(NA, nrow = 1, ncol = 6)
for(k in 1:6){
  mod.k.joint <- poLCA.ordinal(f, data = y_joint,weights = w_joint, nclass = k, 
                               nrep = 25, calc.se = FALSE,
                               maxiter = 3000, tol = 1e-5,lt=0.001,
                               seed=1234)
  AIC.k.joint[1,k] <- mod.k.joint$aic
  BIC.k.joint[1,k] <- mod.k.joint$bic
  if(((mod.k.joint$bic/BIC.best.joint-1)*100 < -0.2) | k==1){
    BIC.best.joint <- mod.k.joint$bic
    mod.best.joint <- mod.k.joint
    k.best.joint <- k
  }
}

#### Computationally Intensive se calculation for joint-report ####
# Calculate se for the best model fit 
mod.best.joint <- poLCA.ordinal(f, data = y_joint,weights = w_joint, 
                                nclass = k.best.joint, 
                             nrep = 25, calc.se = TRUE,
                             maxiter = 3000, tol = 1e-5,lt=0.001,
                             seed=1234)


plot(AIC.k.joint[1, ], main = "AIC vs #Classes", xlab = "Classes", ylab = "AIC", type = "b")
plot(BIC.k.joint[1, ], main = "BIC vs #Classes", xlab = "Classes", ylab = "BIC", type = "b")


# Doing jittered plots faceted by predicted latent class
library(tidyr)
y.wide <- cbind(joint_report_sample$spid,y_joint)
var.names <- c("id","Fallhist","Fearfall","health","pains","Depression",
               "Anxiety","Exercise","InformalSoc","FormalSoc","Fatigue","strata","cluster")
colnames(y.wide) <- var.names
y.long <- gather(y.wide, "Variable","Response", -id )

y.wide.ord <- cbind(y.wide,mod.best.joint$predclass)
colnames(y.wide.ord)[ncol(y.wide.ord)] <- "Predclass"
y.long.ord <- gather(y.wide.ord, "Variable","Response",Fallhist:Fatigue)
y.long.ord$Variable <- factor(y.long.ord$Variable, 
                              levels = c("Fallhist","Fearfall","health",
                                         "pains","Fatigue","Depression",
                                         "Anxiety","Exercise","InformalSoc",
                                         "FormalSoc"))
p <- ggplot(data = y.long.ord, aes(x = Variable, y = Response, group = id))
p + 
  geom_line(alpha=1/10, position = position_jitter(width = 0, height = .1)) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  facet_grid(. ~ Predclass)


# Plot for class conditioned expected score
class.cond.exp.score <- matrix(NA,ncol = length(mod.best.joint$probs),
                               nrow = k.best.joint)
colnames(class.cond.exp.score) <- var.names[2:11]
rownames(class.cond.exp.score) <- paste("Class ",1:k.best.joint,sep="")

for(i in 1:ncol(class.cond.exp.score)){
  max.outcome <- ncol(mod.best.joint$probs[[i]])
  class.cond.exp.score[,i] <- (mod.best.joint$probs[[i]] %*% 
                                 (1:max.outcome))/max.outcome*100
}

df_cces <- as.data.frame(class.cond.exp.score)
df_cces <- df_cces%>% mutate(Class=rownames(df_cces))

df_cces_melt <- reshape2::melt(df_cces,id.vars="Class")
df_cces_melt$variable <- factor(df_cces_melt$variable, 
                                levels=c("Fallhist","Fearfall",
                                         "health","pains","Fatigue",
                                         "Depression","Anxiety",
                                         "Exercise","InformalSoc",
                                         "FormalSoc"))

ggplot(df_cces_melt, aes(x = variable, y = value)) + 
  geom_line(aes(color = Class, group = Class),size=1.2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("Variable")+
  ylab("Class conditioned expected score (%)")
