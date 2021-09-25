
##### Data processing and analysis code for Generalist Evolution project - species lifespan data (100 time step output)
##### Code by Jonathan R. Morris
##### Latest update 9-23-2021
##### Used for data processing, analysis, and figure preparation for manuscript accepted to Scientific Reports, Sept. 2021 (J. R. Morris, K. A. Allhoff, F. S. Valdovinos)

###Load required packages
library(Rmisc)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(plotrix)
library(pryr)
library(rsq)
library(mgcv)

###__________________________________________________________________________________________________________###
###__________________________________________________________________________________________________________###
###_____________________________________________DATA PROCESSING______________________________________________###
###__________________________________________________________________________________________________________###

#####**************************************SET DATA ANALYSIS CONDITIONS************************************#####

REPSIZE<-100                        #####Use this to SET NUMBER OF REPS (can set to lower values to run portion of reps)
IS<-0.1                             #####Use this to SET INCREMENT SIZE of parameter sweep
IL<-1:50                            #####Use this to SET INCREMENT LENGTH of parameter sweep
FOLDER<-"Z_SWEEP_DATA_FINAL"        #####Use this to SET FOLDER LOCATION for simulation data

#####******************************************************************************************************#####

#################################################################################################################
###-------------------------------DATA PROCESSING: IMPORT AND PROCESS DATA------------------------------------###

sumGSLSEXTALL<-list()
sumMIFRALL<-list()
sumMIFRALLmed<-list()
sumSLOPEALL<-list()
  for (k in 1:REPSIZE) {
    HTEMP <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",k,"/z_*/lifespans.csv")), read.csv)           ###NOTE: folder names may need to be updated below in order to run code
    sumGSEXT<-list()
    sumMIFR<-list()
    sumMIFRmed<-list()
    slo<-c()
     for(i in 1:length(HTEMP)) {
       HTEMP[[i]]$extinction[HTEMP[[i]]$extinction == Inf] <- 25000000                                       ###Extant species at end of simulation are set to extinct to avoid NAs
       HTEMP[[i]]$lifespan<-HTEMP[[i]]$extinction-HTEMP[[i]]$introduction                                    ###Calculates LIFESPANS for all species
       HTEMP[[i]]$typeEXT<-factor(ifelse(HTEMP[[i]]$feedingRange>=.39,                                       ###Defines GENERALISTS/SPECIALISTS (JUST EXTREMES)
       "Generalists", ifelse(HTEMP[[i]]$feedingRange>.32, 
                             "Neither", "Specialists")))   
       HTEMP[[i]]<-subset(HTEMP[[i]], introduction > 50000)                                                  ###Removes TRANSIENT DYNAMICS
       HTEMP[[i]]$lifespan[HTEMP[[i]]$lifespan < 101] <- NA                                                  ###Removes species that go extinct instantly (NON-Survivor Species)
       HTEMP[[i]]<-HTEMP[[i]][complete.cases(HTEMP[[i]]), ]                                                  ###Removes NAs
       sumGSEXT[[i]]<-summarySE(HTEMP[[i]], measurevar="lifespan", 
                                groupvars=c("typeEXT"))    %>%  mutate(z=(i*IS))                             ###Creates LIFESPAN summary for GENERALISTS vs. SPECIALISTS (EXTREME BINS)
       sumMIFR[[i]]<-summarySE(HTEMP[[i]], measurevar="feedingRange", 
                               groupvars=c("speciesType"))    %>%  mutate(z=(i*IS))                          ###Creates FEEDING RANGE summary for MUTANTS vs. INVADERS
       sumMIFRmed[[i]]<-HTEMP[[i]] %>% group_by(speciesType) %>% 
         summarise(median=median(feedingRange),n=n()) %>% mutate(z=(i*IS))                                   ###Creates MEDIAN FEEDING RANGE summary for MUTANTS vs. INVADERS
       slo[i]<-round((summary(glm(log10(HTEMP[[i]]$lifespan)~log10(HTEMP[[i]]$feedingRange), 
                                  family=Gamma(link="log")))$coefficients[2]), digits=4)                     ###Calculates slope of LIFESPAN TO FEEDINGRANGE regressions for all simulations
     }
  if(k==17) {                                                                   ###Select replicate to use for individual lifespan slope figures
      for (i in 1:length(HTEMP)) {
        HTEMP[[i]]$zvalue<-i*.1
      }
      SLOPEEXAMPLE<-bind_rows(HTEMP) 
    }
  rm(HTEMP)
  sumGSLSEXTALL[[k]]<-bind_rows(sumGSEXT)
  sumMIFRALL[[k]]<-bind_rows(sumMIFR)
  sumMIFRALLmed[[k]]<-bind_rows(sumMIFRmed)
  sumSLOPEALL[[k]]<-data.frame(slope=slo, z=IL*IS)
}
GSLSEXTALL<-subset(bind_rows(sumGSLSEXTALL), typeEXT!="Neither")                ###GENERALIST/SPECIALIST LIFESPAN data
MIFRALL<-subset(bind_rows(sumMIFRALL), speciesType!="Ancestor")                 ###MUTANT/INVADER FEEDING RANGE (mean) data
MIFRALLmed<-subset(bind_rows(sumMIFRALLmed), speciesType!="Ancestor")           ###MUTANT/INVADER FEEDING RANGE (median) data
SLOPEALL<-bind_rows(sumSLOPEALL)                                                ###LIFESPAN SLOPE data

rm(list=ls(pattern="^sum"))                                                     ###Deletes unneeded files


#################################################################################################################
###--------------------------------DATA PROCESSING: EXPORT/LOAD-----------------------------------------------###

###Export summary data frames as CSV (to avoid long run time of data processing)
# write.csv(GSLSEXTALL, file = "~/GSLSEXTALL.csv", row.names = FALSE)
# write.csv(MIFRALL, file = "~/MIFRALL.csv", row.names = FALSE)
# write.csv(MIFRALLmed, file = "~/MIFRALLmed.csv", row.names = FALSE)
# write.csv(SLOPEALL, file = "~/SLOPEALL.csv", row.names = FALSE)
# write.csv(SLOPEEXAMPLE, file = "~/SLOPEEXAMPLE.csv", row.names = FALSE)

###Import summary data frames as CSV (to avoid long run time of data processing)
GSLSEXTALL<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/GSLSEXTALL.csv")
MIFRALL<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/MIFRALL.csv")
MIFRALLmed<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/MIFRALLmed.csv")
SLOPEALL<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/SLOPEALL.csv")
SLOPEEXAMPLE<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/SLOPEEXAMPLE.csv")

#####______________________________________________________________________________________________________#####
###__________________________________________________________________________________________________________###
###__________________________________________STASTICIAL ANALYSIS_____________________________________________###
###__________________________________________________________________________________________________________###

###GAM (Gamma(loglink)) for mean feeding range (mutants only)
gamMEANFR<-gam(feedingRange~s(z), family=Gamma(link="log"), data=subset(MIFRALL, speciesType=="Mutant"))
summary(gamMEANFR)
rsq(gamMEANFR)

###GAM (Gamma(loglink)) for median feeding range (mutants only)
gamMEDFR<-gam(median~s(z), family=Gamma(link="log"), data=subset(MIFRALLmed, speciesType=="Mutant"))
summary(gamMEDFR)
rsq(gamMEDFR)

###GAM for lifespan slopes
gamSLOPE<-gam(slope~s(z), data=SLOPEALL)
summary(gamSLOPE)
rsq(gamSLOPE)

###GLM (Poisson) for relative lifespan data (with both reference points) (Extreme bins)
GSLSEXTALL$typeEXT<-factor(GSLSEXTALL$typeEXT,levels=c("Specialists", "Generalists"))
glmtest<-glm(lifespan~z*typeEXT, family=poisson, data=GSLSEXTALL)
summary(glmtest)
rsq(glmtest)

GSLSEXTALL$typeEXT<-factor(GSLSEXTALL$typeEXT,levels=c("Generalists", "Specialists"))
glmtest2<-glm(lifespan~z*typeEXT, family=poisson, data=GSLSEXTALL)
summary(glmtest2)
rsq(glmtest2)

###CREATES NEW DATA SETS WITH STATS FOR PLOTTING
SLOPEALL<-cbind(SLOPEALL,predict.gam(gamSLOPE,SLOPEALL,type="response",se.fit=TRUE))
MIFRALL<-cbind(subset(MIFRALL, speciesType=="Mutant"),predict.gam(gamMEANFR,subset(MIFRALL, speciesType=="Mutant"),type="response",se.fit=TRUE))
MIFRALLmed<-cbind(subset(MIFRALLmed, speciesType=="Mutant"),predict.gam(gamMEDFR,subset(MIFRALLmed, speciesType=="Mutant"),type="response",se.fit=TRUE))
GSLSEXTALL<-cbind(GSLSEXTALL,predict.glm(glmtest,GSLSEXTALL,type="response",se.fit=TRUE))

#####______________________________________________________________________________________________________#####
###__________________________________________________________________________________________________________###
###________________________________________________FIGURES___________________________________________________###
###__________________________________________________________________________________________________________###

###############################################################################################
###------------------------LIFESPAN OF GENERALISTS VS SPECIALISTS---------------------------###

##FIGURE: Mean Lifespan by Species Type (ONLY EXTREMES: Uses cut-offs of .39 and .32) (SE vanishingly small)
p1<-ggplot(GSLSEXTALL, aes(x=z, y=lifespan, group=typeEXT, colour=typeEXT, ribbon=typeEXT)) +
  geom_point(shape=19, size=2, stroke=0, alpha=.2) +
  geom_smooth(method="glm",formula=y~x, method.args=list(family="poisson"), size=.5, colour="red", se=TRUE, fill="red") +
  scale_colour_manual(values=c("black", "blue"), name="type") +
  scale_y_continuous("Mean species lifespan", labels = scales::scientific, limits=c(0,2000000), breaks=seq(0,2000000,by=500000)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position=c(.78,.88))+
  guides(colour = guide_legend(override.aes = list(alpha = .75)))

ggsave("GENERALISTS vs. SPECIALISTS Mean Lifespan by Z (EXTREMES).pdf", p1, path="~/", width=4, height=4)

###############################################################################################
###----------------------------------LIFESPAN TO FEEDING RANGE------------------------------###

###FIGURE: Lifespan slopes (coefficient of FEEDING RANGE to LIFESPAN regression) of simulations
px<-ggplot(SLOPEALL, aes(x=z, y=slope))+
  geom_point(shape=19, size=2, stroke=0, alpha=.2) +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_y_continuous("Lifespan slope",  breaks=seq(-3,0,.5)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold"))+
  ggtitle("(a)")

###FIGURE: Example plot of all lifespan slope regressions for each z value (from one parameter sweep replicate)
py<-ggplot(SLOPEEXAMPLE, aes(x=feedingRange, y=lifespan, group=zvalue))+
  geom_smooth(method="glm",formula=y~x, se=FALSE, fullrange=FALSE, method.args=list(family=Gamma(link="log")), size=1, alpha=.5, aes(colour=zvalue)) +
  scale_color_gradient(low = "red", high = "blue")+
  scale_y_log10("Species lifespan") +
  scale_x_log10(expression(paste("Feeding range (", italic("s"), ")"))) +
  theme_classic() +
  theme(legend.position=c(.85,.75),
        plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold"))+
  labs(color=expression(paste(italic("z"))))+
  ggtitle("(b)")

multislop<-grid.arrange(px, py, nrow=1)
ggsave("Lifespan Multiplot.pdf", multislop, path="~/", width=8, height=4)

###############################################################################################
###--------------------------------------FEEDING RANGE--------------------------------------###

###FIGURE: Mean FEEDING RANGE Plot (JUST MUTANTS)
p1<-ggplot(subset(MIFRALL, speciesType=="Mutant"), aes(x=z, y=feedingRange, group=speciesType, colour=speciesType)) + 
  geom_point(shape=19, size=2, stroke=0, alpha=.2) +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_colour_manual(values=c("black", "blue"), name="type") +
  scale_y_continuous(expression(paste("Mean feeding range (", italic("s"), ")")), limits=c(.305,.37), breaks=seq(.3,.37,.01)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none",
        plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold"))+
  ggtitle("(a)")

###FIGURE: MEDIAN FEEDING RANGE Plot (JUST MUTANTS)
p2<-ggplot(subset(MIFRALLmed, speciesType=="Mutant"), aes(x=z, y=median, group=speciesType, colour=speciesType)) + 
  geom_point(shape=19, size=2, stroke=0, alpha=.2) +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_colour_manual(values=c("black", "blue"), name="type") +
  scale_y_continuous(expression(paste("Median feeding range (", italic("s"), ")")), limits=c(.3,.33), breaks=seq(.3,.33,.01)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="none")

ggsave("Median FEEDING RANGE by Z (JUST MUTANTS).pdf", p2, path="~/", width=4, height=4)


