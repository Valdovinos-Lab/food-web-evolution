
##### Data processing and analysis code for Generalist Evolution project - community data (50000 timestep output) and turnover data (10000 timestep output)
##### Code by Jonathan R. Morris
##### Latest update 9-23-2021
##### Used for data processing, analysis, and figure preparation for manuscript submitted to Scientific Reports, Sept. 2021 (J. R. Morris, K. A. Allhoff, F. S. Valdovinos)

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
###________________________________________IMPORT AND TEST DATA______________________________________________###
###__________________________________________________________________________________________________________###

#####**************************************SET DATA ANALYSIS CONDITIONS************************************#####

REPSIZE<-100                     #####Use this to SET NUMBER OF REPS (can set to lower values to run portion of reps)
IS<-0.1                          #####Use this to SET INCREMENT SIZE of parameter sweep
FOLDER<-"Z_SWEEP_DATA_FINAL"     #####Use this to FOLDER LOCATION for simulation data

#####******************************************************************************************************#####

############################################################################################################
##--------------------------------IMPORT DATA FILES ACROSS PARAMETER SWEEP--------------------------------##

###NOTE: File location will need to be changed to run code below

###RES files### (resource biomass and other general data)
#Nested for loop to import replicated data (deals with parameter sweeps nested within replicate folders)
RES<- list() 
for (i in 1:REPSIZE) {
  RES[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_res.txt")), read.table)
}

###SPECIES TURNOVER files###
TO<- list() 
for (i in 1:REPSIZE) {
  TO[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/turnover.txt")), read.table)
}

###FEEDING RANGE files###  ***Manually set length to 100 (max species possible in Korinna's model)***
S<- list() 
for (i in 1:REPSIZE) {
  S[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_s.txt")), read.table,
                   col.names = paste0("V",seq_len(100)), fill=TRUE, header=FALSE)
}

###FEEDING CENTER files###  ***Manually set length to 100 (max species possible in Korinna's model)***
C<- list() 
for (i in 1:REPSIZE) {
  C[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_c.txt")), read.table,
                   col.names = paste0("V",seq_len(100)), fill=TRUE, header=FALSE)
}

###BODY SIZE files###  ***Manually set length to 100 (max species possible in Korinna's model)***
BS<- list() 
for (i in 1:REPSIZE) {
  BS[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_m.txt")), read.table,
                    col.names = paste0("V",seq_len(100)), fill=TRUE, header=FALSE)
}

###POP BIOMASS files###  ***Manually set length to 100 (max species possible in Korinna's model)***
PB<- list() 
for (i in 1:REPSIZE) {
  PB[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_pop.txt")), read.table,
                    col.names = paste0("V",seq_len(100)), fill=TRUE, header=FALSE)
}

###TROPHIC POSITION files###  ***Manually set length to 100 (max species possible in Korinna's model)***
TP<- list() 
for (i in 1:REPSIZE) {
  TP[[i]] <- lapply(Sys.glob(paste0("~/",FOLDER,"/R_",i,"/z_*/EvoInv_tp.txt")), read.table,
                    col.names = paste0("V",seq_len(100)), fill=TRUE, header=FALSE)
}


############################################################################################################
##--------------------------------------TEST FOR SIMULATION CRASHES---------------------------------------##

##########Loops through all files, creates dataset of all run lengths, tests for CRASHED runs
TOSrep<-list()
for (k in 1:length(TO)) {
  TOS<-list()
  for (i in 1:length(TO[[1]])) {
    TOS[[i]]<-NROW(TO[[k]][[i]])
  }
  TOS<-data.frame(TOS)
  colnames(TOS)<-c(1:length(TO[[1]])*.1)
  TOSrep[[k]]<-TOS
}

TOStest<-bind_rows(TOSrep)  
which(TOStest!=2500, arr.ind=TRUE)
crashed<-data.frame(which(TOStest!=2500, arr.ind=TRUE))       ###Use to record locations of crashed runs


############################################################################################################
##--------------------------------------TEST FOR PARAMETER VALUES-----------------------------------------##

###PAR files### (parameter record to check that import worked correctly)
PAR<-list()
ZPARTEST<-list()
PPARTEST<-list()
SSPARTEST<-list()
for (k in 1:REPSIZE) {
  ZPAR<-c()
  PPAR<-c()
  SSPAR<-c()
  PAR[[k]]<-lapply(Sys.glob(paste0("~/",FOLDER,"/R_",k,"/z_*/parameters.txt")),
                   read.table, fill=TRUE, header=FALSE, stringsAsFactors=F)
  for (i in 1:length(PAR[[1]])) {
    ZPAR[i]<-PAR[[k]][[i]][7,2]        ### z parameter
    PPAR[i]<-PAR[[k]][[i]][8,2]        ### p parameter
    SSPAR[i]<-PAR[[k]][[i]][25,1]      ### seed set parameter
  }
  ZPARTEST[[k]]<-as.numeric(ZPAR)                     ### Set as numeric to avoid being stored as factor
  PPARTEST[[k]]<-PPAR
  SSPARTEST[[k]]<-SSPAR
}
ZPARTEST<-data.frame(ZPARTEST); rownames(ZPARTEST)<-c((1:length(PAR[[1]]))); colnames(ZPARTEST)<-paste("R",c((1:REPSIZE)),sep="")
PPARTEST<-data.frame(PPARTEST); rownames(PPARTEST)<-c((1:length(PAR[[1]]))); colnames(PPARTEST)<-paste("R",c((1:REPSIZE)),sep="")
SSPARTEST<-data.frame(SSPARTEST); rownames(SSPARTEST)<-c((1:length(PAR[[1]]))); colnames(SSPARTEST)<-paste("R",c((1:REPSIZE)),sep="")

###SEED SET CHECK: records which have altered seed sets and stores in list
seedsetcheck<-list()
for(i in 1:REPSIZE) {
  seedsetcheck[[i]]<-separate(SSPARTEST,colnames(SSPARTEST[i]),c("stuff", "seed"))$seed
}
SSPARTEST<-data.frame(seedsetcheck); colnames(SSPARTEST)<-paste("R",c((1:REPSIZE)),sep="")
SEEDSUBS<-list()
for(i in 1:REPSIZE){
    SEEDSUBS[[i]]<-which(SSPARTEST[,i]!=i, arr.ind=TRUE)
}

###Z PARAMETER CHECK: Records where Z parameter values are incorrect
parmcheck<-c()
for(i in 1:nrow(ZPARTEST)) {
  parmcheck[i]<-isTRUE((max(ZPARTEST[i,])-(i*.1)<.00000001) & (min(ZPARTEST[i,])-(i*.1)<.00000001))   ###Note: cannot make these equal to 0 because values are slightly off
}

###Provides warning if crashed runs detected
if(length(matrix(which(TOStest!=2500)))==0){
  print("NO CRASHES")
} else {
  print("!!!!!ERROR: CRASHED RUNS!!!!!")
}

###Provides warning if z parameters are wrong
if(length(which(parmcheck!=TRUE))==0){
  print("Z PARAMETER SWEEP CORRECT") 
} else {
  print("!!!!!ERROR: Z PARAMTER SWEEP!!!!!")
}

###Provides warning if p parameters are wrong
if(length(matrix(which(PPARTEST==0.20)))==0){
  print("P PARAMETER CORRECT")
} else {
  print("!!!!ERROR: P PARAMTER!!!!!")
}


############################################################################################################
##-------------------------MANIPULATE DATA TO EXTRACT INSTANTANEOUS NEW SPECIES---------------------------##

###Function to remove last value (invader/mutant) from rows###
removelastvalue <- function(x) head(x[!is.na(x)], -1)

#Removes last value from FEEDING RANGE data
SCrep<-list()
for (k in 1:length(RES)) {
  SC<-list()
  for (i in 1:length(RES[[1]])) {                              #Nested for loop: cycles through all rows in all data sets, removes last value, then remakes data frame
    SC[[i]]<-apply(S[[k]][[i]], 1, removelastvalue)            #Goes through all rows in data frame and removes last value before NAs (transposes format)
    for (l in 1:501) {                                         #For loop runs through columns (original rows) and sets length to 100
      length(SC[[i]][[l]]) <-100
    }
    SC[[i]]<-as.data.frame(SC[[i]])                            #Remake as data frame
    colnames(SC[[i]])<-c(1:100)                                #Rename columns
    SC[[i]]<-as.data.frame(t(SC[[i]]))                         #Transpose back to original format
  }
  SCrep[[k]]<-SC
}
rm(S)                                                          ###Deletes file

#Removes last value from BODY SIZE data
BSCrep<-list()
for (k in 1:length(RES)) {
  BSC<-list()
  for (i in 1:length(RES[[1]])) {
    BSC[[i]]<-apply(BS[[k]][[i]], 1, removelastvalue)
    for (l in 1:501) {
      length(BSC[[i]][[l]]) <-100
    }
    BSC[[i]]<-as.data.frame(BSC[[i]])
    colnames(BSC[[i]])<-c(1:100)  
    BSC[[i]]<-as.data.frame(t(BSC[[i]]))
  }
 BSCrep[[k]]<-BSC
}
rm(BS)

#Removes last value from FEEDING CENTER data
CCrep<-list()
for (k in 1:length(RES)) {
  CC<-list()
  for (i in 1:length(RES[[1]])) {
    CC[[i]]<-apply(C[[k]][[i]], 1, removelastvalue)
    for (l in 1:501) {
      length(CC[[i]][[l]]) <-100
    }
    CC[[i]]<-as.data.frame(CC[[i]])
    colnames(CC[[i]])<-c(1:100)  
    CC[[i]]<-as.data.frame(t(CC[[i]]))
  }
  CCrep[[k]]<-CC
}
rm(C)

#Removes last value from POP BIOMASS data
PBCrep<-list()
for (k in 1:length(RES)) {
  PBC<-list()
  for (i in 1:length(RES[[1]])) {
    PBC[[i]]<-apply(PB[[k]][[i]], 1, removelastvalue)
    for (l in 1:501) {
      length(PBC[[i]][[l]]) <-100
    }
    PBC[[i]]<-as.data.frame(PBC[[i]])
    colnames(PBC[[i]])<-c(1:100)  
    PBC[[i]]<-as.data.frame(t(PBC[[i]]))
  }
  
 PBCrep[[k]]<-PBC
}
rm(PB)

#Removes last value from TROPHIC POSITION data
TPCrep<-list()
for (k in 1:length(RES)) {
  TPC<-list()
  for (i in 1:length(RES[[1]])) {
    TPC[[i]]<-apply(TP[[k]][[i]], 1, removelastvalue)
    for (l in 1:501) {
      length(TPC[[i]][[l]]) <-100
    }
    TPC[[i]]<-as.data.frame(TPC[[i]])
    colnames(TPC[[i]])<-c(1:100)  
    TPC[[i]]<-as.data.frame(t(TPC[[i]]))
  }
  
  TPCrep[[k]]<-TPC
}
rm(TP)

###__________________________________________________________________________________________________________###
###__________________________________________________________________________________________________________###
###__________________________________________DATA PROCESSING_________________________________________________###
###__________________________________________________________________________________________________________###


#################################################################################################################
###---------------------------DATA PROCESSING: TRAITS and POPULATION METRICS----------------------------------###

###Creates dataframe of summary data for all TRAITS and consumer POPULATION metrics by Z and REPLICATE

CBsum<-list()                            ###Community Biomass
ALL<-list()                              ###ALL traits/communities (for REALIZED FR)
ALLTROP<-list()                          ###ALL trophic position, misc.
for (k in 1:length(RES)) {
  alltrait<-list()
  allcombm<-list()
  for (i in 1:length(RES[[1]])) {
    na.omit(cbind(melt((BSCrep[[k]][[i]][-c(3)]), id.vars=c("V1", "V2")), melt((CCrep[[k]][[i]][-c(1:3)]), id.vars=c()),
                  melt((SCrep[[k]][[i]][-c(1:3)]), id.vars=c()), melt((PBCrep[[k]][[i]][-c(1:3)]), id.vars=c()),
                  melt((TPCrep[[k]][[i]][-c(1:3)]), id.vars=c()))[-c(3,5,7,9,11)] %>%
              rename(specrich=V1, time=V2, bodysize=value, feedingcenter=value.1, 
                     feedingrange=value.2, popbiomass=value.3, trophicpos=value.4) %>%
              mutate(specrich=specrich-1, popdensity=popbiomass/bodysize, z=i*IS, trophicposCLEAN=trophicpos,
                     tfeedingcenter=10^feedingcenter, tbodysize=10^bodysize))->alltrait[[i]]                          ###Calculates population density, transformed traits, cleaned trophic pos
    alltrait[[i]]$trophicposCLEAN[alltrait[[i]]$trophicposCLEAN<1]<-NA                                                ###Cleans trophic position (values less than 1)
    alltrait[[i]]$trophicposCLEAN[alltrait[[i]]$trophicposCLEAN>10]<-NA                                               ###Cleans trophic position (values greater than 10)
    data.frame(combm=tapply(alltrait[[i]]$popbiomass, alltrait[[i]]$time, FUN=sum, na.rm=TRUE),
               z=i*IS, time=unique(alltrait[[i]]$time),                                                               ###Calculates community biomass
               maxtrop=tapply(alltrait[[i]]$trophicposCLEAN, alltrait[[i]]$time, FUN=max, na.rm=TRUE),                ###Calculates max trophic position for each time output
               meantrop=tapply(alltrait[[i]]$trophicposCLEAN, alltrait[[i]]$time, FUN=mean, na.rm=TRUE),              ###Calculates mean trophic position for each time output
               alltrait[[i]] %>% 
                 group_by(time) %>% 
                 summarise(nugen=sum(feedingrange>.35), nuall=sum(feedingrange>=.30)) %>% 
                 mutate(propgen=nugen/nuall)                                                                          ###Generates proportion of generalists. CHANGE threshold if needed
               )->allcombm[[i]]
  }
  allcombmbind<-bind_rows(allcombm)
  allcombmbind$maxtrop[allcombmbind$maxtrop == -Inf] <- NA                                                            ###Cleans max trophic position
  alltraitbind<-bind_rows(alltrait)
  CBsum[[k]]<-summarySE(allcombmbind, measurevar="combm", groupvars=c("z"), na.rm=TRUE)
  ALL[[k]]<-alltraitbind
  ALLTROP[[k]]<-allcombmbind
}
CBsum<-bind_rows(CBsum) %>% mutate(replicate=rep(1:REPSIZE, each=length(RES[[1]])))

###Adds replicate labels to ALL TRAIT data and converts list to single file
for (k in 1:REPSIZE) {
ALL[[k]]<-ALL[[k]] %>% mutate(replicate=rep(k, each=nrow(ALL[[k]])))
}
ALLBIND<-bind_rows(ALL)

###Adds replicate labels to ALL TROPOS data and converts list to single file
ALLTROPt<-list()
for (k in 1:REPSIZE) {
  ALLTROPt[[k]]<-ALLTROP[[k]] %>% mutate(replicate=rep(k, each=nrow(ALLTROP[[k]])))
}
ALLTROP<-bind_rows(ALLTROPt)


#################################################################################################################
###----------------------DATA PROCESSING: SPECIES RICHNESS & RESOURCE BIOMASS---------------------------------###

###Creates dataframe of summary data for SPECIES RICHNESS and RESOURCE BIOMASS metrics by Z and REPLICATE

SRsum<-list()
RBMsum<-list()
for (k in 1:length(RES)) {
  allRES<-list()
  for (i in 1:length(RES[[1]])) {
    data.frame(cbind(RES[[k]][[i]]$V1, RES[[k]][[i]]$V3, RES[[k]][[i]]$V10))  %>%
      rename(time=X1, speciesrich=X2, resourcebiomass=X3) %>%
      mutate(z=i*IS, speciesrichext=speciesrich-1) -> allRES[[i]]
  }
  allRESbind<-bind_rows(allRES)
  allRESbind<-subset(allRESbind, time!=0)                                                                 ###Removes time 0 (ancestor community)
  SRsum[[k]]<-summarySE(allRESbind, measurevar="speciesrichext", groupvars=c("z"), na.rm=TRUE)            ###Uses only extant species in SR data (remove instantaneous invader)
  RBMsum[[k]]<-summarySE(allRESbind, measurevar="resourcebiomass", groupvars=c("z"), na.rm=TRUE)
}
SRsum<-bind_rows(SRsum) %>% mutate(replicate=rep(1:REPSIZE, each=length(RES[[1]])))
RBMsum<-bind_rows(RBMsum) %>% mutate(replicate=rep(1:REPSIZE, each=length(RES[[1]])))


#################################################################################################################
###--------------------------------DATA PROCESSING: SPECIES TURNOVER------------------------------------------###

###Creates dataframe of summary data for SPECIES TURNOVER metrics by Z and REPLICATE

TOSsum<-list()
for (k in 1:length(TO)) {
  allTOS<-list()
  for (i in 1:length(TO[[1]])) {
    allTOS[[i]]<-data.frame(turnover=1-((TO[[k]][[i]]$V2)/(TO[[k]][[i]]$V3)), overlap=(TO[[k]][[i]]$V5), z=i*IS)
  }
  allTOSbind<-bind_rows(allTOS)
  TOSsum[[k]]<-summarySE(allTOSbind, measurevar="turnover", groupvars=c("z"), na.rm=TRUE)
}
TOSsum<-bind_rows(TOSsum) %>% mutate(replicate=rep(1:REPSIZE, each=length(RES[[1]])))

#################################################################################################################
#################################################################################################################
###-------------------------------DATA PROCESSING: REALIZED FEEDING RANGE-------------------------------------###

#ALLBIND<-subset(ALLBIND, replicate==1)                                 ###Subset REPLICATES for testing

###Calculates maximum attack rate for each species
ALLBIND$maxattack<-((ALLBIND$tbodysize^.75))*(1/(ALLBIND$feedingrange*(sqrt(2*pi))))*
  (exp(-((((log10(ALLBIND$tfeedingcenter)-log10(ALLBIND$tfeedingcenter)))^2)/
           (2*ALLBIND$feedingrange*ALLBIND$feedingrange))))
FINALFR<-list()
length(FINALFR)<-length(unique(ALLBIND$replicate))
for (v in unique(ALLBIND$replicate)) {
  alltraitbind<-subset(ALLBIND, replicate==v)                     ###Subsets by replicate
  bindedtemp<-list()
  length(bindedtemp)<-length(unique(alltraitbind$z))
  for (g in unique(alltraitbind$z)) {
    tempz<-subset(alltraitbind, z==g)                             ###Subsets by z value
    tempstore<-list()
    length(tempstore)<-length(unique(tempz$time))
    for (t in unique(tempz$time)) {
      temp<-subset(tempz, time==t)                          ###Subsets each unique community by TIME OUTPUT
      tempvecpreysize<-temp$tbodysize                       ###Makes vector of all specices body size in community
      temp$numberofsp<-length(tempvecpreysize)              ###Number of species in community
      temp$allmean<-mean(tempvecpreysize)                   ###Mean body size of community
      temp$allsd<-sd(tempvecpreysize)                       ###SD body size of community
      temp$allmax<-max(tempvecpreysize)                     ###Max body size in community
      temp$allmin<-min(tempvecpreysize)                     ###Min body size in community
      fullcommsize<-c(tempvecpreysize, 1)                   ###Adding in basal resource with body size 1
      for (i in 1:nrow(temp)) {
        ###Calculates attack rate for each species on all species in community (including basal resource)
        fullcomm<-((temp[i,]$tbodysize^.75))*(1/(temp[i,]$feedingrange*(sqrt(2*pi))))*
          (exp(-((((log10(temp[i,]$tfeedingcenter)-log10((c(temp$tbodysize,1)))))^2)/
                   (2*temp[i,]$feedingrange*temp[i,]$feedingrange))))
        thresh<-fullcomm/temp[i,]$maxattack               ###Calculates proportion of realized attack rate to max attack rate for each resource (prey & basal resourse)
        fullcommsizeR<-fullcommsize[!(thresh)<.1]         ###Removes prey below THRESHOLD from community body size vector
        temp$numberofspR[i]<-length(fullcommsizeR)        ###Number of species consumed by focal species (above threshold)
        temp$allmeanR[i]<-mean(fullcommsizeR)             ###Mean body size of species consumed (above threshold)
        temp$allsdR[i]<-sd(fullcommsizeR)                 ###SD body size of species consumed (above threshold)
        temp$allmaxR[i]<-max(fullcommsizeR)               ###MAX species body size consumed (above threshold)
        temp$allminR[i]<-min(fullcommsizeR)               ###MIN species body size consumed (above threshold)
      }
      tempstore[[t/50000]]<-temp
    }
    bindedtemp[[g*10]]<-bind_rows(tempstore)
  }
  FINALFR[[v]]<-bind_rows(bindedtemp)
}
bindedtempt<-bind_rows(bindedtemp)                                 ###For comparison check for one replicate
FINALFRBIND<-bind_rows(FINALFR) %>%
  mutate(bredth=allmaxR-allminR, totalspec=numberofsp+1,           ###Calculates feeding breath and proportion of community consumed (1 added for basal resource)
         perccom=numberofspR/totalspec)

###Summarize data for realized FR and diet breadth across all replicates
REALFRSUM<-summarySE(FINALFRBIND, measurevar="perccom", groupvars=c("z","replicate"), na.rm=TRUE)


#################################################################################################################
###------------------------------DATA PROCESSING: CLEANING/EXPORTING/LOADING----------------------------------###

# rm(list=setdiff(ls(),c(ls(pattern="sum"),"SLOPEall", "alltrait", "alltraitbind", 
#                        "allRESbind", "allcombmbind", "REPSIZE", "IS", "RES")))           ###DELETES unnecessary objects

###Subset for time series figure
ALLBINDsub<-subset(ALLBIND[ALLBIND$replicate %in% seq(11,15,1),],
                   z==0.1 | z==0.2 | z==0.5 | z==1.0 | z==2.0 | z==5.0)[c(2,3,9,13)] 
ALLTROPsub<-subset(ALLTROP[ALLTROP$replicate %in% seq(11,15,1),],
                   z==0.1 | z==0.2 | z==0.5 | z==1.0 | z==2.0 | z==5.0)[c(2,3,4,10)] 

###Export summary data frames as CSV (to avoid long run time of data processing)
# write.csv(CBsum, file = "~/CBsum.csv", row.names = FALSE)
# write.csv(RBMsum, file = "~/RBMsum.csv", row.names = FALSE)
# write.csv(TOSsum, file = "~/TOSsum.csv", row.names = FALSE)
# write.csv(ALLBIND, file = "~/ALLBIND.csv", row.names = FALSE)
# write.csv(ALLTROP, file = "~/ALLTROP.csv", row.names = FALSE)
# write.csv(ALLBINDsub, file = "~/ALLBINDsub.csv", row.names = FALSE)
# write.csv(ALLTROPsub, file = "~/ALLTROPsub.csv", row.names = FALSE)
# write.csv(FINALFRBIND, file = "~/FINALFRBIND.csv", row.names = FALSE)
# write.csv(subset(FINALFRBIND, z==5.0 & replicate==2), file = "~/FINALFRBINDSAMPLEgl.csv", row.names = FALSE)
# write.csv(REALFRSUM, file = "~/REALFRSUMgl.csv", row.names = FALSE)

###Import summary data frames as CSV (to avoid long run time of data processing)
ALLBINDsub<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/ALLBINDsub.csv")
ALLTROPsub<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/ALLTROPsub.csv")
CBsum<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/CBsum.csv")
RBMsum<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/RBMsum.csv")
TOSsum<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/TOSsum.csv")
FINALFRBINDSAMPLE<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/FINALFRBINDSAMPLE.csv")
REALFRSUM<-read.csv("~/Z_SWEEP_DATA_SUMMARY_FINAL/REALFRSUM.csv")

#####______________________________________________________________________________________________________#####
###__________________________________________________________________________________________________________###
###__________________________________________STASTICIAL ANALYSIS_____________________________________________###
###__________________________________________________________________________________________________________###

###VARIATION METRIC 
gamTURNOVER<-gam(turnover~s(z), family=Gamma(link=log), data=TOSsum)
summary(gamTURNOVER)
rsq(gamTURNOVER)
gamRESOURCE<-gam(sd~s(z), family=Gamma(link=log), data=RBMsum)
summary(gamRESOURCE)
rsq(gamRESOURCE)
gamCOMBIOMASS<-gam(sd~s(z), family=Gamma(link=log), data=CBsum)
summary(gamCOMBIOMASS)
rsq(gamCOMBIOMASS)

###REALIZED REEDING RANGE STATS
glmREALIZEDFR<-glm(perccom~z, family=quasibinomial(link=logit), data=REALFRSUM)
summary(glmREALIZEDFR)
rsq(glmREALIZEDFR)

###MODIFIES DATA SETS WITH STATS FOR PLOTTING
TOSsum<-cbind(TOSsum,predict.gam(gamTURNOVER,TOSsum,type="response",se.fit=TRUE))
RBMsum<-cbind(RBMsum,predict.gam(gamRESOURCE,RBMsum,type="response",se.fit=TRUE))
CBsum<-cbind(CBsum,predict.gam(gamCOMBIOMASS,CBsum,type="response",se.fit=TRUE))
REALFRSUM<-cbind(REALFRSUM,predict.glm(glmREALIZEDFR,REALFRSUM,type="response",se.fit=TRUE))

#####______________________________________________________________________________________________________#####
###__________________________________________________________________________________________________________###
###________________________________________________FIGURES___________________________________________________###
###__________________________________________________________________________________________________________###


############################################################################################################
###-------------------------------------------SPECIES TURNOVER-------------------------------------------###

###FIGURE: Mean Species Turnover Across Z Sweep
MSTP<-ggplot(TOSsum, aes(x=z, y=turnover)) + 
  geom_point(shape=19, size=2, stroke=0, alpha=.2, fill="black") +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_y_continuous("Mean species turnover", limits=c(.0,.4), breaks=seq(.0,.4,by=.08)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold")) +
  ggtitle("(c)")

############################################################################################################
###-----------------------------------------RESOURCE BIOMASS---------------------------------------------###

###FIGURE: SD Resource Biomass Plot Across Z Sweep
RBSD<-ggplot(RBMsum, aes(x=z, y=sd)) + 
  geom_point(shape=19, size=2, stroke=0, alpha=.2, fill="black") +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_y_continuous("Basal resource biomass SD",  limits=c(0,.16), breaks=seq(0,.16,by=.04)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold")) +
  ggtitle("(a)")

############################################################################################################
###--------------------------------------------COMMUNITY BIOMASS-----------------------------------------###

###FIGURE: SD Community Biomass Plot Across Z Sweep
CBSD<-ggplot(CBsum, aes(x=z, y=sd)) + 
  geom_point(shape=19, size=2, stroke=0, alpha=.2, fill="black") +
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_y_continuous("Community biomass SD", limits=c(.02,.4), breaks=seq(.02,.4,by=.06)) +
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold")) +
  ggtitle("(b)")

######---------------------------Variation Figure----------------------------#####

summary<-grid.arrange(arrangeGrob(RBSD, CBSD, MSTP, heights=c(.25,.25,.3), layout_matrix=cbind(c(1,2,3))))
ggsave("Variation Figure.pdf", summary, path="~/", width=4, height=9)

############################################################################################################
###------------------------------------REALIZED FEEDING RANGE------------------------------------------#####

###FIGURE: Fundamental FR vs. Realized FR for example simulation
fundvsreal<-ggplot(FINALFRBINDSAMPLE, aes(x=feedingrange, y=perccom))+
  geom_point(shape=19, size=2, stroke=0, alpha=.2)+   
  geom_smooth(method="glm",formula=y~x, method.args=list(family=quasibinomial(link="logit")), size=.5, colour="red", se=TRUE, fill="red") +
  scale_y_continuous("Prop. of community consumed", limits=c(0,1))+
  scale_x_continuous(expression(paste("Fundamental feeding range (", italic("s"), ")")), limits=c(.3,.4))+
  theme_classic()

ggsave("Realized Feeding Range vs. Fundamental.pdf", fundvsreal, path="~/", width=4, height=4)

###FIGURE: Mean realized FR across all simulations
propcom<-ggplot(REALFRSUM, aes(x=z, y=perccom))+
  geom_point(shape=19, size=2, stroke=0, alpha=.2)+   
  geom_ribbon(aes(x=z, ymin=fit-se.fit*1.96, ymax=fit+se.fit*1.96), alpha=0.3, fill="red", size=0, colour=NA) +
  geom_line(aes(x=z, y=fit), colour="red", size=.5) +
  scale_y_continuous("Mean prop. of community consumed", limits=c(.22,.385), breaks=seq(.22,.38,.04))+
  scale_x_continuous(expression(paste("Invader strangeness (", italic("z"), ")"))) +
  theme_classic()+
  theme(plot.title = element_text(margin=margin(-3.5,0), size=11, hjust = -.16, face="bold"))+
  ggtitle("(b)")

###Combined figure fundamental and realized feeding range (NOTE: must use data from ANALYSIS_ZSWEEP_HIST to generate combined figure) 
combifr<-grid.arrange(p1, propcom, nrow=1)
ggsave("Feeding Range (Fundamental and Realized) by z.pdf", combifr, path="~/", width=8, height=4)

############################################################################################################
###------------------------------BODY SIZE TIME SERIES (W/ MAX TP)-------------------------------------#####

###FIGURE: Body size across time (multiple z, multiple replicates)
zsample<-c(.1,.2,.5,2,5)
a<-list()
for(k in zsample) {
  for(i in 11:15) {
    a[[1+length(a)]]<-ggplot(subset(ALLBINDsub, replicate==i & z==k), aes(x=time, y=bodysize)) +    
      geom_point(size=1.5, shape=19, stroke=0, alpha=.1) +
      geom_line(data=subset(ALLTROPsub, replicate==i & z==k), aes(x=time, y=maxtrop), size=.6, color="green", alpha=.7) +
      scale_y_continuous(limits=c(0,6)) +
      scale_x_continuous(breaks=seq(0,25000000,8000000))+
      theme_classic() +
      theme(axis.title=element_blank())
  }   
}
p<-list()
for(k in zsample) {
  for(i in 11:15) {
    p[[1+length(p)]]<-ggplot(subset(ALLBINDsub, replicate==i & z==k), aes(x=time, y=bodysize)) +    
      geom_point(size=1.5, shape=19, stroke=0, alpha=.1) +
      geom_line(data=subset(ALLTROPsub, replicate==i & z==k), aes(x=time, y=maxtrop), size=.6, color="green", alpha=.7) +
      scale_y_continuous(limits=c(0,6)) +
      scale_x_continuous(breaks=seq(0,25000000,8000000))+
      theme_classic() +
      theme(plot.title=element_text(margin=margin(-3.5,0), size=12, hjust=.5, face=2),
            axis.title=element_blank())+
      ggtitle(paste0("z ", k))
  }   
}
greencurvesimple1<-grid.arrange(p[[1]],p[[6]],p[[11]],p[[16]],p[[21]],
                                a[[2]],a[[7]],a[[12]],a[[17]],a[[22]],
                                a[[3]],a[[8]],a[[13]],a[[18]],a[[23]],
                                a[[4]],a[[9]],a[[14]],a[[19]],a[[24]],
                                a[[5]],a[[10]],a[[15]],a[[20]],a[[25]],
                                nrow=5, bottom="Time", left="log(Body size)")

ggsave("Body Size by Time with Max Trophic Position.pdf", greencurvesimple1, path="~/", width = 12, height = 12)



