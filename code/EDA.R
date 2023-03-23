rm(list=ls())
#library(tidyverse)
library(dplyr)
library(ggplot2)
setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Q2 2021-22 Spring/ADAP/code/LCA-ADAP/")
survey_dat <- read.csv('data/2011-2018_nhats_data_lca_sub.csv') 
#survey_dat <- survey_dat[complete.cases(survey_dat), ]


self_report_sample <- survey_dat %>% filter(proxyres5_sp==0,dementia5_cd==1)
proxy_report_sample <- survey_dat %>% filter(proxyres5_sp==1, dementia5_cd==1)


proxy_report_sample <- proxy_report_sample[ ,-8]
proxy_report_sample <- proxy_report_sample[complete.cases(proxy_report_sample), ]
self_report_sample <- self_report_sample[complete.cases(self_report_sample), ]

#reordering values to start from 1 to max value instead of 0
self_report_sample$phq5_oad <- self_report_sample$phq5_oad+1
self_report_sample$gad5_oad <- self_report_sample$gad5_oad+1 
self_report_sample$exercise5_sp <- self_report_sample$exercise5_sp+1
self_report_sample$infse_sp5_con <- self_report_sample$infse_sp5_con+1
self_report_sample$fse_sp5_con <- self_report_sample$fse_sp5_con+1

#reordering values to start from 1 to max value instead of 0
proxy_report_sample$phq5_oad <- proxy_report_sample$phq5_oad+1
proxy_report_sample$gad5_oad <- proxy_report_sample$gad5_oad+1 
proxy_report_sample$exercise5_sp <- proxy_report_sample$exercise5_sp+1
proxy_report_sample$infse_sp5_con <- proxy_report_sample$infse_sp5_con+1
proxy_report_sample$fse_sp5_con <- proxy_report_sample$fse_sp5_con+1

self_report_outcomes <- self_report_sample %>% select(ordfallyr15,
                                                      fof5_cat, 
                                                      slfhealth5_sp, 
                                                      pain5_sp, 
                                                      phq5_oad, 
                                                      gad5_oad,
                                                      exercise5_sp,
                                                      infse_sp5_con,
                                                      fse_sp5_con,
                                                      fatigue_sp5)
self_outcomes_summary <- apply(self_report_outcomes, 2, table)
self_outcomes_perc <- apply(self_report_outcomes,2, function(x) round(table(x)/length(x)*100,1))
w_self <- self_report_sample$w5anfinwgt0/sum(self_report_sample$w5anfinwgt0)*nrow(self_report_sample)


proxy_report_outcomes <- proxy_report_sample %>% select(ordfallyr15,
                                                      fof5_cat, 
                                                      slfhealth5_sp, 
                                                      pain5_sp, 
                                                      phq5_oad, 
                                                      gad5_oad,
                                                      exercise5_sp,
                                                      infse_sp5_con,
                                                      fse_sp5_con,
                                                      fatigue_sp5)
proxy_outcomes_summary <- apply(proxy_report_outcomes, 2, table)
proxy_outcomes_perc <- apply(proxy_report_outcomes,2, function(x) round(table(x)/length(x)*100,1))
w_proxy <- proxy_report_sample$w5anfinwgt0/sum(proxy_report_sample$w5anfinwgt0)*nrow(proxy_report_sample)


setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Q2 2021-22 Spring/ADAP/code/LCA-ADAP/code/poLCA-master/poLCA-master/")
dyn.load("src/poLCA.dll")


## Source all functions from R folder

setwd("./R/")
files.sources = list.files()
sapply(files.sources, source)

# set you working directory to location where main_script and data is saved
setwd("C:/Users/Anupreet Porwal/Dropbox/Academics UW Seattle/Q2 2021-22 Spring/ADAP/code/LCA-ADAP/")

y_self <- self_report_outcomes
y_proxy <- proxy_report_outcomes

#### Self report analysis ####
set.seed(12)
f <- cbind(ordfallyr15, fof5_cat, slfhealth5_sp,pain5_sp,phq5_oad,
           gad5_oad,exercise5_sp,infse_sp5_con,fse_sp5_con,fatigue_sp5) ~ 1

BIC.best.ordinal <- BIC.best.polca <- Inf
AIC.k <- matrix(0, nrow = 2, ncol = 6)
BIC.k <- matrix(0, nrow = 2, ncol = 6)
chisq.k <- matrix(0, nrow = 2, ncol = 6)
for(k in 1:6){
  mod.k <- poLCA.ordinal(f, data = y_self,weights = w_self, nclass = k, 
                         nrep = 20, calc.se = FALSE,
                         maxiter = 3000, tol = 1e-4)
  AIC.k[1,k] <- mod.k$aic
  chisq.k[1,k] <- mod.k$Chisq
  BIC.k[1,k] <- mod.k$bic
  if(mod.k$bic<BIC.best.ordinal){
    BIC.best.ordinal <- mod.k$bic
    mod.best.ord <- mod.k
  }
  mod.k.polca <- poLCA(f, data = y_self, nclass = k, nrep = 20)
  AIC.k[2,k] <- mod.k.polca$aic
  chisq.k[2,k] <- mod.k.polca$Chisq
  BIC.k[2,k] <- mod.k.polca$bic
  if(mod.k.polca$bic<BIC.best.polca){
    BIC.best.polca <- mod.k.polca$bic
    mod.best.polca <- mod.k.polca
  }
  
}

plot(AIC.k[1, ], main = "AIC vs K", xlab = "Classes", ylab = "AIC", type = "b")
plot(chisq.k[1, ], main = "Chi-Squared vs K", xlab = "Classes", ylab = "Chi-Squared", type = "b")
plot(BIC.k[1, ], main = "BIC vs R", xlab = "Classes", ylab = "BIC", type = "b")

library(tidyr)
y.wide <- cbind(self_report_sample$spid,y_self)
var.names <- c("id","Fallhist","Fearfall","health","pains","Depression",
               "Anxiety","Exercise","InformalSoc","FormalSoc","Fatigue")
colnames(y.wide) <- var.names
y.long <- gather(y.wide, "Variable","Response", -id )


# Doing jittered plots faceted by predicted latent class for both our model and 
# LCA by poLCA package 
y.wide.ord <- cbind(y.wide,mod.best.ord$predclass)
colnames(y.wide.ord)[ncol(y.wide.ord)] <- "Predclass"
y.long.ord <- gather(y.wide.ord, "Variable","Response",Fallhist:Fatigue)
p <- ggplot(data = y.long.ord, aes(x = Variable, y = Response, group = id))
p + 
  geom_line(alpha=1/10, position = position_jitter(width = 0, height = .1)) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  facet_grid(. ~ Predclass)





#### proxy report analysis ####
set.seed(12)
f <- cbind(ordfallyr15, fof5_cat, slfhealth5_sp,pain5_sp,phq5_oad,
           gad5_oad,exercise5_sp,infse_sp5_con,fse_sp5_con,fatigue_sp5) ~ 1

BIC.best.ordinal.proxy <- BIC.best.polca.proxy <- Inf
AIC.k.proxy <- matrix(0, nrow = 2, ncol = 6)
BIC.k.proxy <- matrix(0, nrow = 2, ncol = 6)
chisq.k.proxy <- matrix(0, nrow = 2, ncol = 6)
for(k in 1:6){
  mod.k.proxy <- poLCA.ordinal(f, data = y_proxy,weights = w_proxy, nclass = k, 
                         nrep = 20, calc.se = FALSE,
                         maxiter = 3000, tol = 1e-4)
  AIC.k.proxy[1,k] <- mod.k.proxy$aic
  chisq.k.proxy[1,k] <- mod.k.proxy$Chisq
  BIC.k.proxy[1,k] <- mod.k.proxy$bic
  if(mod.k.proxy$bic<BIC.best.ordinal.proxy){
    BIC.best.ordinal.proxy <- mod.k.proxy$bic
    mod.best.ord.proxy <- mod.k.proxy
  }
  mod.k.polca.proxy <- poLCA(f, data = y_proxy, nclass = k, nrep = 20)
  AIC.k.proxy[2,k] <- mod.k.polca.proxy$aic
  chisq.k.proxy[2,k] <- mod.k.polca.proxy$Chisq
  BIC.k.proxy[2,k] <- mod.k.polca.proxy$bic
  if(mod.k.polca.proxy$bic<BIC.best.polca.proxy){
    BIC.best.polca.proxy <- mod.k.polca.proxy$bic
    mod.best.polca.proxy <- mod.k.polca.proxy
  }
  
}


y.wide.proxy <- cbind(proxy_report_sample$spid,y_proxy)
var.names <- c("id","Fallhist","Fearfall","health","pains","Depression",
               "Anxiety","Exercise","InformalSoc","FormalSoc","Fatigue")
colnames(y.wide.proxy) <- var.names
y.long.proxy <- gather(y.wide.proxy, "Variable","Response", -id )


# Doing jittered plots faceted by predicted latent class for both our model and 
# LCA by poLCA package 
y.wide.ord.proxy <- cbind(y.wide.proxy,mod.best.ord.proxy$predclass)
colnames(y.wide.ord.proxy)[ncol(y.wide.ord.proxy)] <- "Predclass"
y.long.ord.proxy <- gather(y.wide.ord.proxy, "Variable","Response",Fallhist:Fatigue)
p.proxy <- ggplot(data = y.long.ord.proxy, aes(x = Variable, y = Response, group = id))
p.proxy + 
  geom_line(alpha=3/10, position = position_jitter(width = 0, height = .1)) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  facet_grid(. ~ Predclass)



#### Self+ proxy joint report analysis ####
joint_report_sample <- rbind(self_report_sample[ ,-8],proxy_report_sample)

joint_report_outcomes <- joint_report_sample %>% select(ordfallyr15,
                                                        fof5_cat, 
                                                        slfhealth5_sp, 
                                                        pain5_sp, 
                                                        phq5_oad, 
                                                        gad5_oad,
                                                        exercise5_sp,
                                                        infse_sp5_con,
                                                        fse_sp5_con,
                                                        fatigue_sp5)
w_joint <- joint_report_sample$w5anfinwgt0/sum(joint_report_sample$w5anfinwgt0)*nrow(joint_report_sample)
y_joint <- joint_report_outcomes
set.seed(12)
f <- cbind(ordfallyr15, fof5_cat, slfhealth5_sp,pain5_sp,phq5_oad,
           gad5_oad,exercise5_sp,infse_sp5_con,fse_sp5_con,fatigue_sp5) ~ 1

BIC.best.ordinal.joint <- BIC.best.polca.joint <- Inf
AIC.k.joint <- matrix(0, nrow = 2, ncol = 6)
BIC.k.joint <- matrix(0, nrow = 2, ncol = 6)
chisq.k.joint <- matrix(0, nrow = 2, ncol = 6)
for(k in 1:6){
  mod.k.joint <- poLCA.ordinal(f, data = y_joint,weights = w_joint, nclass = k, 
                               nrep = 20, calc.se = FALSE,
                               maxiter = 3000, tol = 1e-4)
  AIC.k.joint[1,k] <- mod.k.joint$aic
  chisq.k.joint[1,k] <- mod.k.joint$Chisq
  BIC.k.joint[1,k] <- mod.k.joint$bic
  if(mod.k.joint$bic<BIC.best.ordinal.joint){
    BIC.best.ordinal.joint <- mod.k.joint$bic
    mod.best.ord.joint <- mod.k.joint
  }
  # mod.k.polca.joint <- poLCA(f, data = y_joint, nclass = k, nrep = 20)
  # AIC.k.joint[2,k] <- mod.k.polca.joint$aic
  # chisq.k.joint[2,k] <- mod.k.polca.joint$Chisq
  # BIC.k.joint[2,k] <- mod.k.polca.joint$bic
  # if(mod.k.polca.joint$bic<BIC.best.polca.joint){
  #   BIC.best.polca.joint <- mod.k.polca.joint$bic
  #   mod.best.polca.joint <- mod.k.polca.joint
  # }
  
}

library(tidyr)
y.wide.joint <- cbind(joint_report_sample$spid,
                      joint_report_sample$proxyres5_sp,
                      y_joint)
var.names <- c("id","proxy","Fallhist","Fearfall","health","pains","Depression",
               "Anxiety","Exercise","InformalSoc","FormalSoc","Fatigue")
colnames(y.wide.joint) <- var.names
y.long.joint <- gather(y.wide.joint, "Variable","Response", -c(id,proxy) )


# Doing jittered plots faceted by predicted latent class for both our model and 
# LCA by poLCA package 
y.wide.ord.joint <- cbind(y.wide.joint,mod.best.ord.joint$predclass)
colnames(y.wide.ord.joint)[ncol(y.wide.ord.joint)] <- "Predclass"
y.long.ord.joint <- gather(y.wide.ord.joint, "Variable","Response",Fallhist:Fatigue)
y.long.ord.joint$proxy <- factor(y.long.ord.joint$proxy)
p.joint <- ggplot(data = y.long.ord.joint, aes(x = Variable, y = Response, group = id,color= proxy))
p.joint + 
  geom_line(alpha=3/10, position = position_jitter(width = 0, height = .1)) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  facet_grid(. ~ Predclass)





