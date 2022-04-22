#QDM-013: Analyzing host-adapted microbiomes and their effects on quinoa and BVM growth

##Experimental design: 7x2 factorial CRD with 6 reps per treatment 
#(6 reps water, 3 biological reps for each soil treatment and 2 technical reps per biological rep)
#Treatments:
#Host inoculated: Quinoa, BVM
#Microbiome slurry/origin: Chenopod (BVM), Amaranth, Tomato, steriled Chenopod, sterilized Amaranth,
#sterilized tomato, sterile water

##Analyzed as a RCBD design with experimental replicates being the blocking factor

##data files needed: In Box Plant Pathology Lab/Quinoa DM project/PhD Experiments/QDM-013/pass1_dat

library(car)
library(agricolae)
library(emmeans)
library(ggplot2)
library(ggpubr)

#analyzing passages individually, all organized as two factor completely randomized designs

#Question 1: Does the presence of microbes affect Chenopodium species growth differently?

#####PASSAGE 1 

#make sure everything is numeric
pass1_dat$Num_leaves<-as.numeric(pass1_dat$Num_leaves)
pass1_dat$Tot_biomass<-as.numeric(pass1_dat$Tot_biomass)
pass1_dat$Ag_biomass<-as.numeric(pass1_dat$Ag_biomass)
pass1_dat$Bg_biomass<-as.numeric(pass1_dat$Bg_biomass)

####Below ground biomass models 

#model with exp_rep as a blocking factor
pass1_dat<-na.omit(pass1_dat) #look into this - make sure this command is doing what it is supposed to
pass1_mod<-lm(Bg_biomass ~ Host*Slurry + Exp_rep, pass1_dat)
#need Anova() function because we have an unbalanced data in experimental reps 2 and 3 
Anova(pass1_mod) #slurry approaching significance p = 0.06235

#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test (this is assessing experiments separately so there is no blocking factor)
#resids x preds plot, generate resids for SW test:
pass1_dat$resids<-residuals(pass1_mod)
pass1_dat$preds<-predict(pass1_mod)
pass1_dat$sq_preds<-pass1_dat$preds^2
plot(resids ~ preds, data = pass1_dat) # looks like a funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(pass1_dat$resids) # p-value = 2.328e-10


#levene's test - Tests for homogeneity of variance
#library(car)
leveneTest(Bg_biomass ~ Treatment, data = pass1_dat) # p-value = 0.05167

####Above ground biomass models
#model with exp_rep as a blocking factor
pass1_abmod<-lm(Ag_biomass ~ Host*Slurry + Exp_rep, pass1_dat)
Anova(pass1_abmod) #High effect of Exp_rep p = 2.594e-06


#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test (this is assessing experiments separately so there is no blocking factor)

#resids x preds plot, generate resids for SW test:
pass1_dat$abresids<-residuals(pass1_abmod)
pass1_dat$abpreds<-predict(pass1_abmod)
pass1_dat$absq_preds<-pass1_dat$abpreds^2
plot(abresids ~ abpreds, data = pass1_dat) # looks like funnel


#shapiro wilk - Tests for normality of residuals
shapiro.test(pass1_dat$abresids) # p-value = 0.0002523
#levene's test - Tests for homogeneity of variance
#library(car)
leveneTest(Ag_biomass ~ Treatment, data = pass1_dat) # p-value = 0.3396

#NEED TO TRANSFORM
pass1_dat$trans_ag_biomass<-sqrt(0.5+pass1_dat$Ag_biomass)

trans_pass1_abmod<-lm(trans_ag_biomass ~ Host*Slurry + Exp_rep, pass1_dat)
Anova(trans_pass1_abmod) #High effect of Exp_rep p = 2.923e-06

alt_trans_abmod_p1<-lm(trans_ag_biomass ~ Treatment.vague + Exp_rep, pass1_dat)
Anova(alt_trans_abmod_p1)


#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test (this is assessing experiments separately so there is no blocking factor)

#resids x preds plot, generate resids for SW test:
pass1_dat$trans_abresids<-residuals(trans_pass1_abmod)
pass1_dat$trans_abpreds<-predict(trans_pass1_abmod)
pass1_dat$trans_absq_preds<-pass1_dat$trans_abpreds^2
plot(trans_abresids ~ trans_abpreds, data = pass1_dat) # looks like funnel

#shapiro wilk - Tests for normality of residuals
shapiro.test(pass1_dat$trans_abresids) # p-value = 0.0005635
#levene's test - Tests for homogeneity of variance
#library(car)
leveneTest(trans_ag_biomass ~ Treatment, data = pass1_dat) # p-value = 0.3472


####Total biomass models
#model with exp_rep as a blocking factor
pass1_tbmod<-lm(Tot_biomass ~ Treatment + Exp_rep, pass1_dat)
pass1_tbmod<-lm(Tot_biomass ~ Host*Slurry + Exp_rep, pass1_dat)
Anova(pass1_tbmod) #Significant effect of slurry (p = 0.03496), 
tuk<-HSD.test(pass1_tbmod, "Slurry")
lsd<-LSD.test(pass1_tbmod, "Slurry") #separating means - but no pattern that makes sense

#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test 

#resids x preds plot, generate resids for SW test:
pass1_dat$tbresids<-residuals(pass1_tbmod)
pass1_dat$tbpreds<-predict(pass1_tbmod)
pass1_dat$tbsq_preds<-pass1_dat$tbpreds^2
plot(tbresids ~ tbpreds, data = pass1_dat) # looks a little funnely


#shapiro wilk - Tests for normality of residuals
shapiro.test(pass1_dat$tbresids) # p-value = 2.163e-08
leveneTest(Tot_biomass ~ Treatment, data = pass1_dat) # p-value = 0.05024




#Number of leaves models - experimental reps were significantly different

#model with exp_rep as a blocking factor
pass1_leaves_mod<-lm(Num_leaves ~ Host*Slurry + Exp_rep, pass1_dat)
Anova(pass1_leaves_mod) #significant effect of host (p = 0.0001743) and Exp_rep (p = 1.092e-07)
tuk<-HSD.test(pass1_leaves_mod, "Host")

#need to check if these models fit ANOVA assumptions: shapiro wilks, levene's test (this is assessing experiments separately so there is no blocking factor)

#resids x preds plot, generate resids for SW test:
pass1_dat$resids<-residuals(pass1_leaves_mod)
pass1_dat$preds<-predict(pass1_leaves_mod)
pass1_dat$sq_preds<-pass1_dat$preds^2
plot(resids ~ preds, data = pass1_dat) # 

#shapiro wilk - Tests for normality of residuals
shapiro.test(pass1_dat$resids) # p-value = 0.001341
leveneTest(Num_leaves ~ Treatment, data = pass1_dat) # p-value = 0.1147


#ANCOVA on height data using Height5 (final height) as Y and Height1 (initial height before addition of slurries) as the covariable X
#using pass1_dat csv file for this - but for linear regression I think I need the pass1_heightdat file - they are formatted differently
pass1_dat<-na.omit(pass1_dat)
regression<-lm(Height5 ~ Height1, pass1_dat)
plot(Height5 ~ Height1, pass1_dat)
anova(regression)
summary(regression) #R2 = 0.642 , adjusted R2 = 0.6405, so Height1 accounts for roughly 64% of variation in the experiment
#Height1 estimate from regression summary = 0.8510




#contrasts
#Q1 = comparing microbial slurries to sterilized slurries
#Q2 = comparing chenopod microbial slurries to nonchenopod microbial slurries
#Q3 = comparing sterile water to sterilized soil slurries
#Q4 = comparing A300 plants inoculated with Chenopod slurry vs Faro plants inoculated with Chenopod slurry
#Q5 = comparing A300 plants inoculated with any microbes vs Faro plants inoculated with any microbes
#Q6 = comparing A300 plants inoculated with non-Chenopod slurries vs Faro plants inoculated with non-Chenopod slurries

pass1_dat$Treatment<-as.factor(pass1_dat$Treatment)
#contrasts: add in an H2O vs sterilized slurries question - see if there's a difference, or maybe do H2O vs everything else
contrastmatrix<-cbind(c(-1,-1,+1,+1,+1,0,-1,-1,-1,+1,+1,+1,0,-1),
                      c(-1,+2,0,0,0,0,-1,-1,+2,0,0,0,0,-1),
                      c(0,0,-1,-1,-1,+3,0,0,0,-1,-1,-1,+3,0), 
                      c(0,-1,0,0,0,0,0,0,+1,0,0,0,0,0),
                      c(-1,-1,0,0,0,0,-1,+1,+1,0,0,0,0,+1),
                      c(-1,0,0,0,-1,0,0,+1,0,0,0,0,0,+1))
contrasts(pass1_dat$Treatment)<-contrastmatrix

contrast_bb_mod<-aov(Bg_biomass ~ Treatment+Exp_rep, pass1_dat) #sterile water vs slurries, p = 0.07, 
summary(contrast_bb_mod, split = list(Treatment = list("Microbial vs Sterile" =1, "Microbial chenopod vs microbial non-Chenopod"=2, "Sterile water vs sterile slurries"=3, "A300-Chenopod vs Faro-Chenopod"=4, "A300 microbial vs Faro microbial"=5, "A300 nonChen vs Faro nonChen"=6)))

contrast_tb_mod<-aov(Tot_biomass ~ Treatment+Exp_rep, pass1_dat) #sterile vs slurries p = 0.07, try sterile water vs sterilized slurries
summary(contrast_tb_mod, split = list(Treatment = list("Microbial vs Sterile" =1, "Microbial chenopod vs microbial non-Chenopod"=2, "Sterile water vs sterile slurries"=3, "A300-Chenopod vs Faro-Chenopod"=4, "A300 microbial vs Faro microbial"=5, "A300 nonChen vs Faro nonChen"=6)))

contrast_ab_mod<-aov(trans_ag_biomass ~ Treatment+Exp_rep, pass1_dat) #significant difference between Microbial chenopod vs Microbial non-chenopod slurries
summary(contrast_ab_mod, split = list(Treatment = list("Microbial vs Sterile" =1, "Microbial chenopod vs microbial non-Chenopod"=2, "Sterile water vs sterile slurries"=3, "A300-Chenopod vs Faro-Chenopod"=4, "A300 microbial vs Faro microbial"=5, "A300 nonChen vs Faro nonChen"=6)))

contrast_leaves_mod<-aov(Num_leaves ~ Treatment+Exp_rep, pass1_dat) #approaching significance between microbial chenopod vs non-chenopod (p=0.06)
summary(contrast_leaves_mod, split = list(Treatment = list("Microbial vs Sterile" =1, "Microbial chenopod vs microbial non-Chenopod"=2, "Sterile water vs sterile slurries"=3, "A300-Chenopod vs Faro-Chenopod"=4, "A300 microbial vs Faro microbial"=5, "A300 nonChen vs Faro nonChen"=6)))

pass1_abmod_vague<-lm(Ag_biomass ~ Treatment.vague + Exp_rep, pass1_dat)
Anova(pass1_abmod_vague)
tuk<-HSD.test(pass1_abmod_vague, "Treatment.vague")




#####SUBSETTING DATA SETS BY HOST AND ANALYZING THAT WAY

##passage 1
subp1_bvm_dat<-subset(pass1_dat, Host == "A300")
subp1_q_dat<-subset(pass1_dat, Host == "Faro")

p1bvm_bbmod<-lm(Bg_biomass ~ Slurry + Exp_rep, subp1_bvm_dat)
#need Anova() function because we have an unbalanced data in experimental reps 2 and 3 
Anova(p1bvm_bbmod) #slurry p = 0.04
tuk<-HSD.test(p1bvm_bbmod, "Slurry") 
lsd<-LSD.test(p1bvm_bbmod, "Slurry") #separation but no pattern

subp1_bvm_dat$Slurry<-as.factor(subp1_bvm_dat$Slurry)
subp2_q_dat$Slurry<-as.factor(subp1_q_dat$Slurry)
#contrasts: 1) microbial vs sterilized slurries, 2)microbial chen vs microbial nonchen, 3)h2o vs everyting
contrastmatrix<-cbind(c(-1,-1,+1,+1,+1,0,-1),
                      c(0,+1,0,-1,0,0,0),
                      c(+1,+1,+1,+1,+1,-6,+1))

contrasts(subp1_bvm_dat$Slurry)<-contrastmatrix
contrasts(subp2_q_dat$Slurry)<-contrastmatrix
contrast_p1bvm_bbmod<-aov(Bg_biomass ~ Treatment+Exp_rep, subp1_bvm_dat) #microbial vs sterile p = 0.08
summary(contrast_p1bvm_bbmod, split = list(Treatment = list("Microbial vs Sterilized" =1, "Microbial Chen vs Sterile Chen"=2, "H2O vs everything"=3)))


p1q_bbmod<-lm(Bg_biomass ~ Slurry + Exp_rep, subp1_q_dat)
Anova(p1q_bbmod) #slurry p = 0.22

p1bvm_tbmod<-lm(Tot_biomass ~ Slurry + Exp_rep, subp1_bvm_dat)
Anova(p1bvm_tbmod) #slurry p = 0.05, exp_rep p = 0.04
tuk<-HSD.test(p1bvm_tbmod, "Slurry")
lsd<-LSD.test(p1bvm_tbmod, "Slurry")

contrast_p1bvm_tbmod<-aov(Tot_biomass ~ Treatment+Exp_rep, subp1_bvm_dat) #microbial vs sterile p = 0.08
summary(contrast_p1bvm_tbmod, split = list(Treatment = list("Microbial vs Sterilized" =1, "Microbial Chen vs Sterile Chen"=2, "H2O vs everything"=3)))



p1q_tbmod<-lm(Tot_biomass ~ Slurry + Exp_rep, subp1_q_dat)
Anova(p1q_tbmod) #slurry p = 0.1

p1bvm_abmod<-lm(Ag_biomass ~ Slurry + Exp_rep, subp1_bvm_dat)
Anova(p1bvm_abmod) #slurry p = 0.47, exp_rep p = 0.0003

p1q_abmod<-lm(Ag_biomass ~ Slurry + Exp_rep, subp1_q_dat)
Anova(p1q_abmod) #slurry p = 0.07

p1bvm_nlmod<-lm(Num_leaves ~ Slurry + Exp_rep, subp1_bvm_dat)
Anova(p1bvm_nlmod) #slurry p = 0.15, exp_rep p = 0.02

p1q_nlmod<-lm(Num_leaves ~ Slurry + Exp_rep, subp1_q_dat)
Anova(p1q_nlmod) #slurry p = 0.8 exp_rep p < 0.0001




#####March 29, 2022 - talked to anissa and she thinks I should subset my data and non include
##sterilized slurries, run that anova and do a dunnets

#subsetting bvm data and quinoa data to get sterile water, microbial amaranth, tomato, and chenopod treatments
sub2p1_bvm<-subset(subp1_bvm_dat, Slurry == "sterile water")
sub2p1_bvm_ch<-subset(subp1_bvm_dat, Slurry == "Chenopod")
sub2p1_bvm_am<-subset(subp1_bvm_dat, Slurry == "Amaranth")
sub2p1_bvm_to<-subset(subp1_bvm_dat, Slurry == "Tomato")

sub3p1_bvm<-rbind(sub2p1_bvm, sub2p1_bvm_am, sub2p1_bvm_ch, sub2p1_bvm_to)

sub3p1_bvm_abmod<-lm(Ag_biomass ~ Slurry + Exp_rep, sub3p1_bvm)
Anova(sub3p1_bvm_abmod)

sub3p1_bvm_bgmod<-lm(Bg_biomass ~ Slurry + Exp_rep, sub3p1_bvm)
Anova(sub3p1_bvm_bgmod)

sub3p1_bvm_tbmod<-lm(Tot_biomass ~ Slurry + Exp_rep, sub3p1_bvm)
Anova(sub3p1_bvm_tbmod)

sub3p1_bvm_nlmod<-lm(Num_leaves ~ Slurry + Exp_rep, sub3p1_bvm)
Anova(sub3p1_bvm_nlmod)


sub2p1_q<-subset(subp1_q_dat, Slurry == "sterile water")
sub2p1_q_ch<-subset(subp1_q_dat, Slurry == "Chenopod")
sub2p1_q_am<-subset(subp1_q_dat, Slurry == "Amaranth")
sub2p1_q_to<-subset(subp1_q_dat, Slurry == "Tomato")

sub3p1_q<-rbind(sub2p1_q, sub2p1_q_am, sub2p1_q_ch, sub2p1_q_to)

sub3p1_q_abmod<-lm(Ag_biomass ~ Slurry + Exp_rep, sub3p1_q)
Anova(sub3p1_q_abmod)

sub3p1_q_bgmod<-lm(Bg_biomass ~ Slurry + Exp_rep, sub3p1_q)
Anova(sub3p1_q_bgmod)

ggplot(sub3p1_q, aes(x=Treatment, y=Bg_biomass)) +
  geom_boxplot() +
  xlab("Treatment") + 
  ylab("Belowground biomass (g)")

sub3p1_q_tbmod<-lm(Tot_biomass ~ Slurry + Exp_rep, sub3p1_q)
Anova(sub3p1_q_tbmod) #slurry p = 0.09

ggplot(sub3p1_q, aes(x=Treatment, y=Tot_biomass)) +
  geom_boxplot() +
  xlab("Treatment") + 
  ylab("Total biomass (g)")

sub3p1_q_nlmod<-lm(Num_leaves ~ Slurry + Exp_rep, sub3p1_q)
Anova(sub3p1_q_nlmod)

ggplot(sub3p1_q, aes(x=Treatment, y=Num_leaves)) +
  geom_boxplot() +
  xlab("Treatment") + 
  ylab("Number of leaves")


# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



