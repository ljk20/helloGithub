library(foreign)
library(optmatch)
library(RItools)
library(xtable)

rm(list=ls())

setwd("~/Documents/Full_Fine/Data")
data <- read.dta("phe_chi_match_old.dta", warn.missing.labels=FALSE)

#Read in Additional Home Grown Functions
setwd("~/Documents/AIR Projects/Full Matching/DOS Match Functions")
source("smahal.R")
source("addcaliper.R")

###############################Pair Match###############################

#Estimate p-score
p <- glm(treat ~ mj_use + drink30d + drink5 + drugb_past30 + partner_num 
         +sex_yesno + sex_know + coc_know + drugb_past30_type 
         + sex_safe + dec_making + frpl + health_know + black 
         + mulrace + white + hispanic + stu_gender 
         + healthy_eating + active_num + dis_504 + dis_e 
         + dis_ld + mj_use_miss + drink30d_miss + drink5_miss 
         + partner_num_miss + sex_yesno_miss + sex_know_miss 
         + coc_know_miss + sex_safe_miss + dec_making_miss 
         + frpl_miss + health_know_miss + healthy_eating_miss 
         + active_num_miss, data=data, family=binomial)$fitted.values


#Caliper setting - Std Dev of P Score Divided By 2 Rule of Thumb From Rosenbaum
cal.set <- sd(p)/2

#Covaraiates for Mahal Distance
X <- cbind(data$mj_use, data$drink30d, data$drink5, data$drugb_past30, data$partner_num, data$sex_yesno, data$sex_know, data$coc_know, data$drugb_past30_type, data$sex_safe, data$dec_making, data$frpl, data$health_know, data$black, data$mulrace, data$white, data$hispanic, data$stu_gender, data$healthy_eating, data$active_num, data$dis_504, data$dis_e, data$dis_ld, data$mj_use_miss, data$drink30d_miss, data$drink5_miss, data$partner_num_miss, data$sex_yesno_miss, data$sex_know_miss, data$coc_know_miss, data$sex_safe_miss, data$dec_making_miss, data$health_know_miss, data$healthy_eating_miss, data$active_num_miss)

#Create Rank Based Mahal Distance Matrix
dmat <- smahal(data$treat, X)

#Add a caliper based on the p-score 
dmat <- addcaliper(dmat, data$treat, p, caliper=.01)

#Note: All the steps above are also needed for the full match as well. 
#Next bloc of code is for pair matching - full match is below.

#Pair Match
pm <- pairmatch(dmat, controls=1, data=data)
summary(pm)
stratumStructure(pm)

#Check Balance - Using Function from RITools
bal.pm <- xBalance(treat ~ black + mulrace + white + hispanic 
                   + stu_gender + dis_504 + dis_e + dis_ld + frpl 
                   + mj_use + drink30d + drink5 + drugb_past30 + drugb_past30_type + partner_num + sex_yesno + sex_know + coc_know + sex_safe + dec_making + health_know +  healthy_eating + active_num + mj_use_miss + drink30d_miss + drink5_miss + partner_num_miss + sex_yesno_miss + sex_know_miss + coc_know_miss + sex_safe_miss + dec_making_miss + health_know_miss + healthy_eating_miss + active_num_miss, data = data, report= c("std.diffs", "adj.means", "z.scores"), strata = data.frame(original = factor("none"), pm))

names <- c("African American 1/0", "Multi-Racial 1/0", "White 1/0", 
           "Hispanic 1/0","Female 1/0", "Disability type 1 1/0", "Disability type 2 1/0", "Disability type 3 1/0", "Free or reduced price lunch 1/0", "Marijuana use 1/0", "Drunk in past 30 days 1/0", "5 or more drinks in past 30 days",  "Drug use past 30 days",
           "Type of drugs used", "Number of sexual partners",  "Ever had sex 1/0", "Understand cause of pregnancy 1/0", "Can obtain contraception 1/0", "Perception of sex safety", "Decision-making skill", "Knowledge of healthy eating", "Number of times eating healthy", "Number of days physically active", "Marijuana missing 1/0", "Drinking 30 missing 1/0", "Drink 5 missing 1/0", "Sex partners missing 1/0", "Had sex missing 1/0", "Pregnancy Missing 1/0", "Contraception missing 1/0",
           "Sex safe missing 1/0",  "Decision-making missing", "Eating knowledge missing",  "Healthy eating missing 1/0", "Active missing 1/0")

#print(bal.pm, show.pvals=TRUE)

tab <- xtable(bal.pm)
unmatch <- cbind(round(tab[,1], 3),  round(tab[,2], 3), round(tab[,3], 3), (1 - round(pnorm(abs(tab[,4])), 3))*2)
colnames(unmatch) <- c("Control Mean", "Treat Mean", "Std Diff", "p-val")
rownames(unmatch) <- names
unmatch
xtable(unmatch, digits=3)

setwd("~/Documents/Full_Fine/Analysis/pvals")
save(unmatch, file="pmatch.Rdata")


match <- cbind(round(tab[,6], 3),  round(tab[,7], 3), round(tab[,8], 3), (1 - round(pnorm(abs(tab[,9])), 3))*2)
colnames(match) <- c("Control Mean", "Treat Mean", "Std Diff", "p-val")
rownames(match) <- names
xtable(match, digits=3)


##########################Propensity Score Plots###############################
#Now Format the Matched Data - Split the Data By Matched Groups
data$pscore <- p
#Create Matched Strata Indicators
matchresult.transformed <-  as.numeric(pm)
matchresult.transformed[is.na(matchresult.transformed)] <- 0 
matchresult.transformed

data$match_grp <- matchresult.transformed
match.data <- data[data$match_grp!=0,]


#Propensity Score Comparisons
#Boxplot of Propensity Scores
setwd("~/Documents/Full_Fine/Drafts")
pdf("box-pscore.pdf", width=7, height=4, onefile=FALSE, paper="special")
par(mfrow=c(1,2))
boxplot(p ~ data$treat, ,names=c("Control","Treatment"), ylab="Propensity Score")

boxplot(match.data$pscore ~ match.data$treat,names=c("Control","Treatment"), ylab="Propensity Score")
dev.off()

mean(as.numeric(match.data$pscore))


##-------------------------------------------------------------
## Variable Ratio Match
##-------------------------------------------------------------
fm.em.2 <- fullmatch(dmat, min.controls=1, max.controls=10, data=data)

#Summarize the Strata Structure
summary(fm.em.2)
stratumStructure(fm.em.2)

length(tapply(data$treat, fm.em.2, sum))

bal.em <- xBalance(treat ~ black + mulrace + white + hispanic + stu_gender + dis_504 + dis_e + dis_ld + frpl + mj_use + drink30d + drink5 + drugb_past30 + drugb_past30_type + partner_num + sex_yesno + sex_know + coc_know + sex_safe + dec_making + health_know +  healthy_eating + active_num + mj_use_miss + drink30d_miss + drink5_miss + partner_num_miss + sex_yesno_miss + sex_know_miss + coc_know_miss + sex_safe_miss + dec_making_miss + health_know_miss + healthy_eating_miss + active_num_miss, data = data, report= c("std.diffs", "adj.means", "z.scores"), strata = data.frame(original = factor("none"), fm.em.2))
#print(bal.full.em.2, show.pvals=TRUE)


tab.em <- xtable(bal.em)
match.em <- cbind(round(tab.em[,6], 3),  round(tab.em[,7], 3), round(tab.em[,8], 3), (1 - round(pnorm(abs(tab.em[,9])), 3))*2)
colnames(match.em) <- c("Control Mean", "Treat Mean", "Std Diff", "p-val")
rownames(match.em) <- names
match.em
xtable(match.em, digits=3)



