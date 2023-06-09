# ####################DRAFT HAWKES NEEDS TO EDIT###############
# #PPPPPXXXXXX### library(tidyverse)
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #########################################################################################################################
# ##########################################>>>>>>>>>>>>>>>>>      2010     <<<<<<<<<<<<<<<<<<#############################
# #########################################################################################################################
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2010=convert.inp('ALLDATA.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
#  
# #Provides a summary of data within the dataframe
# summary(PN2010)
# PN2010[1:199,]
# 
# #Privides a count of the varibles within the group (Hatch vs Wild)
# summary(PN2010$origin)
# 
# ####INFORMATIVE PLOTS FOR ANALYSIS #########
# 
# #Cumulative Plot (ECDF plot) of Length by Origin 
# par(pty='s')  # force the plot to be square before we start
# 
# plot(ecdf(PN2010$fl[PN2010$origin=="hatchery"]),
#      xlim=c(135,225),
#      xlab="Fork Length (mm)",
#      ylab="Cumulative Proportion",
#      pch=4,
#      main="Fork Length of Tagged Smolts By Origin",
#      col="blue")
# lines(ecdf(PN2010$fl[PN2010$origin=="wild"]), pch=4,
#       col="red") 
#       
# 
# legend('bottomright', 
#        legend=c("Hatch","Wild"),  # text in the legend
#        col=c("blue","red"),  # point colors
#        pch=4)  # specify the point type to be a square
# 
# #Two-Sample Kolmogorov-Smirnov (KS) Test - in support of the ECDF plot 
# ks.test(PN2010$fl[PN2010$origin=="hatchery"],PN2010$fl[PN2010$origin=="wild"])
# 
# #Library needs to be called up to provide summary of data per origin
# library(plyr)
# library(dplyr)
# #Summary Data for FL
# PN2010 %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(fl,na.rm=TRUE),1),
#                    sd=round(sd(fl,na.rm=TRUE),1),
#                    min=min(fl,na.rm=TRUE),
#                    fQuant=quantile(fl,0.25,na.rm=TRUE),
#                    median=quantile(fl,0.5,na.rm=TRUE), 
#                    tQuant=quantile(fl,0.75,na.rm=TRUE),
#                    max=max(fl,na.rm=TRUE))
# 
# #Summary Data for tagb - could also use summary(PN2010), but data omits stdev
# PN2010 %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(tagb,na.rm=TRUE),3),
#                    sd=round(sd(tagb,na.rm=TRUE),3),
#                    min=min(tagb,na.rm=TRUE),
#                    fQuant=quantile(tagb,0.25,na.rm=TRUE),
#                    median=quantile(tagb,0.5,na.rm=TRUE), 
#                    tQuant=quantile(tagb,0.75,na.rm=TRUE),
#                    max=max(tagb,na.rm=TRUE))
# 
#  ###########CJS Model Work#############################WORK TAKEN FROM LEAK WORKSHOP
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2010=convert.inp('2010_H_W.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
# ##### code to get km intervals###
# ##Allsites#RKMs<-c(44.00, 43.23, 42.23, 40.75, 38.85, 34.60, 31.94, 30.12, 28.33, 26.68, 24.54, 23.61, 22.98, 22.35, 19.64, 17.61, 15.91, 14.72, 13.25, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# #10 Sites- Like G.Goulette Archive
# RKMs<-c(44.00, 38.85, 31.94, 23.61, 15.91, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# 
# RKM2<-rep(NA,length(RKMs)-1)
# for (i in 1:length(RKM2)){
#   RKM2[i]<-RKMs[i]-RKMs[i+1]
# }
# RKM2
# ##############
# #STEP 4 create process object: include: model, groups name, time intervals
# #PN2010.process=process.data(PN2010, time.intervals= RKM2, model= "CJS", groups="origin")
# PN2010.process=process.data(PN2010, time.intervals= RKM2, model= "CJS", groups="origin")
# #STEP 5: create design matrix, remove.unused = T (as all fish were released in the same site)
# #PN2010.ddl=make.design.data(PN2010.process,remove.unused = T)
# PN2010.ddl=make.design.data(PN2010.process,remove.unused = T)
# 
# #STEP 6: create models 
# #time means reach (i.e which interval in being evaluated), a model that includes time explores if there are differences in Phi or P, between different reaches (intervals)
# # a model that includes time, gives a value for every interval
# 
# multi.models= function(){
#   Phi.dot=list(formula=~1) #this is a single survival value (one estimate for the whole river)
#   Phi.time=list(formula=~time) #difference in surv between reaches (time)
#   Phi.time.origin.fl=list(formula=~time+origin+fl)  
#   Phi.time.origin=list(formula=~time+origin)
#   Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   p.time.origin=list(formula=~time+origin)
#   
#   cml=create.model.list("CJS")
#   
#   results=mark.wrapper(cml,data=PN2010.process,ddl=PN2010.ddl,adjust=F)
# }
# 
# PN2010.results=multi.models() 
# 
# #run this to delete junk files
# cleanup(ask=F)
# #This is to get the AIC table
# PN2010.table<-as.data.frame(model.table(PN2010.results))
# 
# #This is to get a specific model:
# #to get the model NAG2017.results$Phi$name.model$p$name.model$results$real
# ###in general model with lowest AIC is explored
# modeltoexplore<-PN2010.results$Phi.time.origin.p.time$results$real #everytime a $ is written, a menu comes out with all the models
# 
# #to save survival results into an excel file:
# write.csv(modeltoexplore,"resultsPN2010.csv")
# 
# ###to explore a covariate
# minfl=min(PN2010$fl)
# maxfl=max(PN2010$fl)
# fl.values=minfl+(0:30)*(maxfl-minfl)/30
# fl.values
# 
# Phibyfl=covariate.predictions(PN2010.results$Phi.time.origin.tagb.fl.p.time.origin,data = data.frame(fl=fl.values),indices=c(1:66))
# 
# 
# 
# #########################################################################################################################
# ##########################################>>>>>>>>>>>>>>>>>      2011     <<<<<<<<<<<<<<<<<<#############################
# #########################################################################################################################
# 
# 
# library(RMark)
# 
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2011=convert.inp('2011_H_W.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
# 
# #Provides a summary of data within the dataframe
# summary(PN2011)
# PN2011[1:203,]
# 
# #Privides a count of the varibles within the group (Hatch vs Wild)
# summary(PN2011$origin)
# 
# ####INFORMATIVE PLOTS FOR ANALYSIS #########
# 
# #Cumulative Plot (ECDF plot) of Length by Origin 
# par(pty='s')  # force the plot to be square before we start
# 
# plot(ecdf(PN2011$fl[PN2011$origin=="hatchery"]),
#      xlim=c(135,225),
#      xlab="Fork Length (mm)",
#      ylab="Cumulative Proportion",
#      pch=4,
#      main="Fork Length of Tagged Smolts By Origin",
#      col="blue")
# lines(ecdf(PN2011$fl[PN2011$origin=="wild"]),pch=4,
#       col="red")
# 
# legend('bottomright', 
#        legend=c("Hatch","Wild"),  # text in the legend
#        col=c("blue","red"),  # point colors
#        pch=4)  # specify the point type to be a square
# 
# #Two-Sample Kolmogorov-Smirnov (KS) Test - in support of the ECDF plot 
# ks.test(PN2011$fl[PN2011$origin=="hatchery"],PN2011$fl[PN2011$origin=="wild"])
# 
# 
# #Library needs to be called up to provide summary of data per origin
# library(plyr)
# library(dplyr)
# #Summary Data for FL
# PN2011 %>%
# dplyr::group_by(origin)%>%
# dplyr::summarize(n=n(),
#           mean=round(mean(fl,na.rm=TRUE),1),
#           sd=round(sd(fl,na.rm=TRUE),1),
#           min=min(fl,na.rm=TRUE),
#           fQuant=quantile(fl,0.25,na.rm=TRUE),
#           median=quantile(fl,0.5,na.rm=TRUE), 
#           tQuant=quantile(fl,0.75,na.rm=TRUE),
#           max=max(fl,na.rm=TRUE))
# 
# #Summary Data for tagb
# PN2011 %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(tagb,na.rm=TRUE),1),
#                    sd=round(sd(tagb,na.rm=TRUE),1),
#                    min=min(tagb,na.rm=TRUE),
#                    fQuant=quantile(tagb,0.25,na.rm=TRUE),
#                    median=quantile(tagb,0.5,na.rm=TRUE), 
#                    tQuant=quantile(tagb,0.75,na.rm=TRUE),
#                    max=max(tagb,na.rm=TRUE))
# 
# 
# ###########CJS Model Work############################
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2011=convert.inp('2011_H_W.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
# ##### code to get km intervals###
# ##Allsites#RKMs<-c(44.00, 43.23, 42.23, 40.75, 38.85, 34.60, 31.94, 30.12, 28.33, 26.68, 24.54, 23.61, 22.98, 22.35, 19.64, 17.61, 15.91, 14.72, 13.25, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# #10 Sites- Like G.Goulette Archive
# RKMs<-c(44.00, 38.85, 31.94, 23.61, 15.91, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# 
# RKM2<-rep(NA,length(RKMs)-1)
# for (i in 1:length(RKM2)){
#   RKM2[i]<-RKMs[i]-RKMs[i+1]
# }
# RKM2
# ##############
# #STEP 4 create process object: include: model, groups name, time intervals
# #PN2010.process=process.data(PN2010, time.intervals= RKM2, model= "CJS", groups="origin")
# PN2011.process=process.data(PN2011, time.intervals= RKM2, model= "CJS", groups="origin")
# #STEP 5: create design matrix, remove.unused = T (as all fish were released in the same site)
# #PN2010.ddl=make.design.data(PN2010.process,remove.unused = T)
# PN2011.ddl=make.design.data(PN2011.process,remove.unused = T)
# 
# #STEP 6: create models 
# #time means reach (i.e which interval in being evaluated), a model that includes time explores if there are differences in Phi or P, between different reaches (intervals)
# # a model that includes time, gives a value for every interval
# 
# multi.models= function(){
#   Phi.dot=list(formula=~1) #this is a single survival value (one estimate for the whole river)
#   Phi.time=list(formula=~time) #difference in surv between reaches (time)
#   Phi.time.origin.fl=list(formula=~time+origin+fl)  
#   Phi.time.origin=list(formula=~time+origin)
#   Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   p.time.origin=list(formula=~time+origin)
#   
#   cml=create.model.list("CJS")
#   
#   results=mark.wrapper(cml,data=PN2011.process,ddl=PN2011.ddl,adjust=F)
# }
# 
# PN2011.results=multi.models() 
# 
# #run this to delete junk files
# cleanup(ask=F)
# #This is to get the AIC table
# PN2011.table<-as.data.frame(model.table(PN2011.results))
# 
# #This is to get a specific model:
# #to get the model NAG2017.results$Phi$name.model$p$name.model$results$real
# ###in general model with lowest AIC is explored
# modeltoexplore<-PN2011.results$Phi.time.origin.p.time$results$real#everytime a $ is written, a menu comes out with all the models
# 
# #to save survival results into an excel file:
# write.csv(modeltoexplore,"resultsPN2011.csv")
# 
# ###to explore a covariate
# minfl=min(PN2011$fl)
# maxfl=max(PN2011$fl)
# fl.values=minfl+(0:30)*(maxfl-minfl)/30
# fl.values
# 
# Phibyfl=covariate.predictions(PN2011.results$Phi.time.origin.tagb.fl.p.time.origin,data = data.frame(fl=fl.values),indices=c(1:66))
# 
# 
# 
# 
# #########################################################################################################################
# ##########################################>>>>>>>>>>>>>>>>>      2012     <<<<<<<<<<<<<<<<<<#############################
# #########################################################################################################################
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2012=convert.inp('2012_H_W.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
# 
# #Provides a summary of data within the dataframe
# summary(PN2012)
# PN2012[1:202,]
# 
# #Privides a count of the varibles within the group (Hatch vs Wild)
# summary(PN2012$origin)
# 
# ####INFORMATIVE PLOTS FOR ANALYSIS #########
# 
# #Cumulative Plot (ECDF plot) of Length by Origin
# par(pty='s')  # force the plot to be square before we start
# 
# plot(ecdf(PN2012$fl[PN2012$origin=="hatchery"]),
#      xlim=c(140,250),
#      xlab="Fork Length (mm)",
#      ylab="Cumulative Proportion",
#      pch-4,
#      main="Fork Length of Tagged Smolts By Origin",
#      col="blue")
# lines(ecdf(PN2012$fl[PN2012$origin=="wild"]), pch=4,
#       col="red")
# 
# legend('bottomright', 
#        legend=c("Hatch","Wild"),  # text in the legend
#        col=c("blue","red"),  # point colors
#        pch=4)  # specify the point type to be a square
# 
# #Two-Sample Kolmogorov-Smirnov (KS) Test - in support of the ECDF plot 
# ks.test(PN2012$fl[PN2012$origin=="hatchery"],PN2012$fl[PN2012$origin=="wild"])
# 
# #Library needs to be called up to provide summary of data per origin
# library(plyr)
# library(dplyr)
# #Summary Data for FL
# PN2012 %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(fl,na.rm=TRUE),1),
#                    sd=round(sd(fl,na.rm=TRUE),1),
#                    min=min(fl,na.rm=TRUE),
#                    fQuant=quantile(fl,0.25,na.rm=TRUE),
#                    median=quantile(fl,0.5,na.rm=TRUE), 
#                    tQuant=quantile(fl,0.75,na.rm=TRUE),
#                    max=max(fl,na.rm=TRUE))
# 
# #Summary Data for tagb
# PN2012 %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(tagb,na.rm=TRUE),3),
#                    sd=round(sd(tagb,na.rm=TRUE),3),
#                    min=min(tagb,na.rm=TRUE),
#                    fQuant=quantile(tagb,0.25,na.rm=TRUE),
#                    median=quantile(tagb,0.5,na.rm=TRUE), 
#                    tQuant=quantile(tagb,0.75,na.rm=TRUE),
#                    max=max(tagb,na.rm=TRUE))
# 
# 
# 
# ###########CJS Model Work############################
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN2012=convert.inp('2012_H_W.inp', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb"))
# ##### code to get km intervals###
# ##Allsites#RKMs<-c(44.00, 43.23, 42.23, 40.75, 38.85, 34.60, 31.94, 30.12, 28.33, 26.68, 24.54, 23.61, 22.98, 22.35, 19.64, 17.61, 15.91, 14.72, 13.25, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# #10 Sites- Like G.Goulette Archive
# RKMs<-c(44.00, 38.85, 31.94, 23.61, 15.91, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# 
# RKM2<-rep(NA,length(RKMs)-1)
# for (i in 1:length(RKM2)){
#   RKM2[i]<-RKMs[i]-RKMs[i+1]
# }
# RKM2
# ##############
# #STEP 4 create process object: include: model, groups name, time intervals
# #PN2010.process=process.data(PN2010, time.intervals= RKM2, model= "CJS", groups="origin")
# PN2012.process=process.data(PN2012, time.intervals= RKM2, model= "CJS", groups="origin")
# #STEP 5: create design matrix, remove.unused = T (as all fish were released in the same site)
# #PN2010.ddl=make.design.data(PN2010.process,remove.unused = T)
# PN2012.ddl=make.design.data(PN2012.process,remove.unused = T)
# 
# #STEP 6: create models 
# #time means reach (i.e which interval in being evaluated), a model that includes time explores if there are differences in Phi or P, between different reaches (intervals)
# # a model that includes time, gives a value for every interval
# 
# 
# multi.models= function(){
#   Phi.dot=list(formula=~1) #this is a single survival value (one estimate for the whole river)
#   Phi.time=list(formula=~time) #difference in surv between reaches (time)
#   Phi.time.origin.fl=list(formula=~time+origin+fl)  
#   Phi.time.origin=list(formula=~time+origin)
#   Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   #Phi.time.origin.tagb.fl=list(formula=~time+origin+tagb+fl)
#   
#   p.time.origin=list(formula=~time+origin)
#   
#   cml=create.model.list("CJS")
#   
#   results=mark.wrapper(cml,data=PN2012.process,ddl=PN2012.ddl,adjust=F)
# }
# 
# PN2012.results=multi.models() 
# 
# #run this to delete junk files
# cleanup(ask=F)
# #This is to get the AIC table
# PN2012.table<-as.data.frame(model.table(PN2012.results))
# 
# #This is to get a specific model:
# #to get the model NAG2017.results$Phi$name.model$p$name.model$results$real
# ###in general model with lowest AIC is explored
# modeltoexplore<-PN2012.results$Phi.time.origin.p.time$results$real#everytime a $ is written, a menu comes out with all the models
# 
# #to save survival results into an excel file:
# write.csv(modeltoexplore,"resultsPN2012.csv")
# 
# ###to explore a covariate
# minfl=min(PN2012$fl)
# maxfl=max(PN2012$fl)
# fl.values=minfl+(0:30)*(maxfl-minfl)/30
# fl.values
# 
# Phibyfl=covariate.predictions(PN2012.results$Phi.time.origin.tagb.fl.p.time.origin,data = data.frame(fl=fl.values),indices=c(1:66))
# 
# #########################################################################################################################
# ##########################################>>>>>>>>>>>>>>>>>      All Years    <<<<<<<<<<<<<<<<<<#############################
# #########################################################################################################################
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
# PN_AllYears=convert.inp('AllYears_H_W', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb", "year"))
# 
# #Provides a summary of data within the dataframe
# summary(PN_AllYears)
# PN_AllYears[1:604,]
# 
# #Privides a count of the varibles within the group (Hatch vs Wild)
# summary(PN_AllYears$origin)
# 
# ####INFORMATIVE PLOTS FOR ANALYSIS #########
# 
# #Cumulative Plot (ECDF plot) of Length by Origin
# par(pty='s')  # force the plot to be square before we start
# 
# plot(ecdf(PN_AllYears$fl[PN_AllYears$origin=="hatchery"]),
#      xlim=c(140,250),
#      xlab="Fork Length (mm)",
#      ylab="Cumulative Proportion",
#      pch=4,
#      main="Fork Length of Tagged Smolts By Origin",
#      col="blue")
# lines(ecdf(PN_AllYears$fl[PN_AllYears$origin=="wild"]), pch=4,
#       col="red")
# 
# legend('bottomright', 
#        legend=c("Hatch","Wild"),  # text in the legend
#        col=c("blue","red"),  # point colors
#        pch=4)  # specify the point type to be a square
# 
# #Two-Sample Kolmogorov-Smirnov (KS) Test - in support of the ECDF plot 
# ks.test(PN_AllYears$fl[PN_AllYears$origin=="hatchery"],PN_AllYears$fl[PN_AllYears$origin=="wild"])
# 
# #Library needs to be called up to provide summary of data per origin
# library(plyr)
# library(dplyr)
# #Summary Data for FL
# PN_AllYears %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(fl,na.rm=TRUE),1),
#                    sd=round(sd(fl,na.rm=TRUE),1),
#                    min=min(fl,na.rm=TRUE),
#                    fQuant=quantile(fl,0.25,na.rm=TRUE),
#                    median=quantile(fl,0.5,na.rm=TRUE), 
#                    tQuant=quantile(fl,0.75,na.rm=TRUE),
#                    max=max(fl,na.rm=TRUE))
# 
# #Summary Data for tagb
# PN_AllYears %>%
#   dplyr::group_by(origin)%>%
#   dplyr::summarize(n=n(),
#                    mean=round(mean(tagb,na.rm=TRUE),3),
#                    sd=round(sd(tagb,na.rm=TRUE),3),
#                    min=min(tagb,na.rm=TRUE),
#                    fQuant=quantile(tagb,0.25,na.rm=TRUE),
#                    median=quantile(tagb,0.5,na.rm=TRUE), 
#                    tQuant=quantile(tagb,0.75,na.rm=TRUE),
#                    max=max(tagb,na.rm=TRUE))
# 
# 
# 
# ###########CJS Model Work############################
# 
# library(RMark)
# #STEP 2: set working directory
# setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
# #STEP 3: read INP (if INP was created using excel and notepad)
# 
# #For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length, tagb=tag burden, and year = year
# PNAYears=convert.inp('AllYears_H_W', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb", "year"))
# ##### code to get km intervals###
# ##Allsites#RKMs<-c(44.00, 43.23, 42.23, 40.75, 38.85, 34.60, 31.94, 30.12, 28.33, 26.68, 24.54, 23.61, 22.98, 22.35, 19.64, 17.61, 15.91, 14.72, 13.25, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# #10 Sites- Like G.Goulette Archive
# RKMs<-c(44.00, 38.85, 31.94, 23.61, 15.91, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
# 
# RKM2<-rep(NA,length(RKMs)-1)
# for (i in 1:length(RKM2)){
#   RKM2[i]<-RKMs[i]-RKMs[i+1]
# }
# RKM2
# ##############
# #STEP 4 create process object: include: model, groups name, time intervals
# PNAYears.process=process.data(PNAYears, time.intervals= RKM2, model= "CJS", groups="origin")
# #STEP 5: create design matrix, remove.unused = T (as all fish were released in the same site)
# PNAYears.ddl=make.design.data(PNAYears.process,remove.unused = T)
# 
# #STEP 6: create models 
# #time means reach (i.e which interval in being evaluated), a model that includes time explores if there are differences in Phi or P, between different reaches (intervals)
# # a model that includes time, gives a value for every interval
# 
# multi.models= function(){
#   Phi.dot=list(formula=~1) #this is a single survival value (one estimate for the whole river)
#   Phi.time.year=list(formula=~time+year) #difference in surv between reaches (time) and year 
#   Phi.time.origin.fl.year=list(formula=~time+origin+fl+year)  #diff between reaches (time),origin, fl by year  
#   #Phi.time.origin.tagb.year=list(formula=~time+origin+tagb+year)  #diff between reaches (time),origin, fl by year  
#   Phi.time.origin.fl.year=list(formula=~time+origin+fl+year)
#   #Phi.time.origin.tagb.year=list(formula=~time+origin+tagb+year)
#   Phi.time.origin.year=list(formula=~time+origin+year) #diff between reaches (time), origin, by year
#   #Phi.time.origin.tagb.fl.year=list(formula=~time+origin+tagb+fl+year) #diff between reaches (time), origin, tagb, fl, by year
#   Phi.time.fl.year=list(formula=~time+fl+year)  #diff between reaches (time), fl by year  
#   #Phi.time.tagb.year=list(formula=~time+tagb+year)
#   #Phi.time.tagb.fl.year=list(formula=~time+tagb+fl+year)
#   Phi.time.origin.fl=list(formula=~time+origin+fl)  #diff between reaches (time),origin, fl by year  
#   #Phi.time.origin.tagb=list(formula=~time+origin+tagb)  #diff between reaches (time),origin, fl by year  
#   Phi.time.origin=list(formula=~time+origin)
# 
#   
#   p.dot=list(formula~1) #this is a single detection value (one estimate for the whole river)
#   p.time.origin=list(formula=~time+origin)
#   p.time.year=list(formula=~time+year)
#   p.time.origin.year=list(formula=~time+origin+year)
#   
#   cml=create.model.list("CJS")
#   
#   results=mark.wrapper(cml,data=PNAYears.process,ddl=PNAYears.ddl,adjust=F)
# }   
# 
# PNAYears.results=multi.models() 
# 
# #run this to delete junk files
# 
# cleanup(ask=F)
# #This is to get the AIC table
# PNAYears.table<-as.data.frame(model.table(PNAYears.results))
# 
# #This is to get a specific model:
# #to get the model NAG2017.results$Phi$name.model$p$name.model$results$real
# ###in general model with lowest AIC is explored
# modeltoexplore<-PNAYears.results$Phi.time.origin.year.p.time.origin$results$real#everytime a $ is written, a menu comes out with all the models
# 
# #to save survival results into an excel file:
# write.csv(modeltoexplore,"resultsPNAYears.csv")
# 
# ###to explore a covariate
# minfl=min(PNAYears$fl)
# maxfl=max(PNAYears$fl)
# fl.values=minfl+(0:30)*(maxfl-minfl)/30
# fl.values
# 
# Phibyfl=covariate.predictions(PNAYears.results$Phi.time.origin.tagb.fl.p.time.origin,data = data.frame(fl=fl.values),indices=c(1:66))

#####################################################################################################
###########>>>>>>>>>>>> Group of year and Origin (6groups) w/ 6 covariates <<<<<<<<##################
#####################################################################################################

library(RMark)
#STEP 2: set working directory
#setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
#STEP 3: read INP (if INP was created using excel and notepad)

#For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length and tagb=tag burden
PN_AllYears7=convert.inp('ALLDATAnew.inp',
                         group.df = data.frame(origin=rep(c("wild","hatchery"),3),year=c(rep("2010",2),rep("2011",2),rep("2012",2))),
                         covariates = c("fl", "tagb", "doy", "diff", "dis", "ss"))
PN_AllYears7<-PN_AllYears7[(which(PN_AllYears7$ch!="10000000000")),]

PN_AllYears7$fl<-(PN_AllYears7$fl-mean(PN_AllYears7$fl))/sd(PN_AllYears7$fl)
PN_AllYears7$tagb<-(PN_AllYears7$tagb-mean(PN_AllYears7$tagb))/sd(PN_AllYears7$tagb)
PN_AllYears7$doy<-(PN_AllYears7$doy-mean(PN_AllYears7$doy))/sd(PN_AllYears7$doy)
PN_AllYears7$diff<-(PN_AllYears7$diff-mean(PN_AllYears7$diff))/sd(PN_AllYears7$diff)
PN_AllYears7$dis<-(PN_AllYears7$dis-mean(PN_AllYears7$dis))/sd(PN_AllYears7$dis)

PN_AllYears7$ss<-(PN_AllYears7$ss-mean(PN_AllYears7$ss))/sd(PN_AllYears7$ss)

PN_AllYears7<-rbind(PN_AllYears7,PN_AllYears7,PN_AllYears7,PN_AllYears7)


str(PN_AllYears7)
# REMOVE BANGOR HP
nrow(PN_AllYears7)

#Provides a summary of data within the dataframe
summary(PN_AllYears7)
PN_AllYears7[1:604,]

#Privides a count of the varibles within the group (Hatch vs Wild)
summary(PN_AllYears7$origin)


###########CJS Model Work############################

library(RMark)
#STEP 2: set working directory
#setwd("C:/Users/james.hawkes/Desktop/Rwork/PN_2010_12")
#STEP 3: read INP (if INP was created using excel and notepad)

#For this: 1 grouping variable (hatchery, wild), 2 covariates, FL= fork length, tagb=tag burden, and year = year
#PN_AllYears7=convert.inp('AllYears_H_WNew', group.df = data.frame(origin=c("wild","hatchery")),covariates = c("fl", "tagb", "year"))
##### code to get km intervals###
##Allsites#RKMs<-c(44.00, 43.23, 42.23, 40.75, 38.85, 34.60, 31.94, 30.12, 28.33, 26.68, 24.54, 23.61, 22.98, 22.35, 19.64, 17.61, 15.91, 14.72, 13.25, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)
#10 Sites- Like G.Goulette Archive
RKMs<-c(40.75, 38.85, 31.94, 23.61, 15.91, 10.12, 6.25, -0.18, -3.88, -19.88, -42.87)

RKM2<-rep(NA,length(RKMs)-1)
for (i in 1:length(RKM2)){
  RKM2[i]<-RKMs[i]-RKMs[i+1]
}
RKM2

write.table(RKM2, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
#Last RKM is 86.87
##############
#STEP 4 create process object: include: model, groups name, time intervals
PN_AllYears7.process=process.data(PN_AllYears7, time.intervals= RKM2, model= "CJS", groups=c("origin","year"))

#STEP 5: create design matrix, remove.unused = T (as all fish were released in the same site)
PN_AllYears7.ddl=make.design.data(PN_AllYears7.process)
#fix to .99
PN_AllYears7.ddl$p$fix[PN_AllYears7.ddl$p$Time==81.72]=0.99

#STEP 6: create models 
#time means reach (i.e which interval in being evaluated), a model that includes time explores if there are differences in Phi or P, between different reaches (intervals)
# a model that includes time, gives a value for every interval


multi.models= function(){
  Phi.dot=list(formula=~1) #this is a single survival value (one estimate for the whole river)
  Phi.time.year=list(formula=~time+year) #difference in surv between reaches (time) and year 
  Phi.time.origin.fl.year.2=list(formula=~time*origin+year)  #diff between reaches (time),origin, and year. Origin interactive effect  
  Phi.time.origin.fl.year=list(formula=~time+origin+year)  #diff between reaches (time),origin, fl by year  
  
 # Phi.time.origin.year=list(formula=~time+origin+year) #diff between reaches (time), origin, by year
  #Phi.time.fl.year=list(formula=~time+fl+year)  #diff between reaches (time), fl by year  
  #Phi.time.origin.fl=list(formula=~time+origin+fl)  #diff between reaches (time),origin, fl by year  
  Phi.time.origin=list(formula=~time+origin)
  Phi.time.origin.2=list(formula=~time*origin)
  
  Phi.time.origin.year.diff=list(formula=~time+origin+year+diff)
  Phi.time.origin.year.dis=list(formula=~time+origin+year+dis)
  Phi.time.origin.year.ss=list(formula=~time+origin+year+ss)
  
  Phi.time.origin.year.diff.2=list(formula=~time*origin*diff+year)
  Phi.time.origin.year.dis.2=list(formula=~time*origin*dis+year)
  Phi.time.origin.year.ss.2=list(formula=~time*origin*ss+year)
  
  Phi.time.origin.year.diff.4=list(formula=~time*origin+origin*diff+year)
  Phi.time.origin.year.dis.4=list(formula=~time*origin+origin*dis+year)
  Phi.time.origin.year.ss.4=list(formula=~time*origin+year+origin*ss)
  
  Phi.time.origin.year.diff.3=list(formula=~time+origin*diff+year)
  Phi.time.origin.year.dis.3=list(formula=~time+origin*dis+year)
  Phi.time.origin.year.ss.3=list(formula=~time+origin*ss+year)
  
  Phi.time.origin.diff=list(formula=~time*origin*diff)
  Phi.time.origin.dis=list(formula=~time*origin*dis)
  Phi.time.origin.ss=list(formula=~time*origin*ss)
  
  Phi.time.origin.diff2=list(formula=~time+origin*diff)
  Phi.time.origin.dis2=list(formula=~time+origin*dis)
  Phi.time.origin.ss2=list(formula=~time+origin*ss)
  
  Phi.time.dis=list(formula=~time+dis) #difference in surv between reaches (time) and year 
  Phi.time.dis.year=list(formula=~time+dis+year) #difference in surv between reaches (time) and year 
  
  Phi.time.ss=list(formula=~time+ss)
  Phi.time.year.ss=list(formula=~time+year+ss)#difference in surv between reaches (time) and year 
  
  Phi.time.diff=list(formula=~time+ss)
  Phi.time.year.diff=list(formula=~time+year+ss)
  
#  Phi.time.origin.fl.dis=list(formula=~time+origin+fl+dis)  #diff between reaches (time),origin, fl by year  
#  Phi.time.origin.fl.year.dis=list(formula=~time+origin+fl+year+dis)
  
  Phi.origin.year.dis=list(formula=~origin+year+dis)
  
  #p.dot=list(formula~1) #this is a single detection value (one estimate for the whole river)
  #p.time.origin=list(formula=~time+origin)
#  p.time.year=list(formula=~time+year)
  p.time.fixed=list(formula=~time+year)
  #p.time.origin.year=list(formula=~time+origin+year)
  
  cml=create.model.list("CJS")
  
  results=mark.wrapper(cml,data=PN_AllYears7.process,ddl=PN_AllYears7.ddl,adjust=F)
}    

models=multi.models() 
saveRDS(models,"newnewmodels.RDS")
#run this to delete junk files

cleanup(ask=F)
#This is to get the AIC table
models.table<-as.data.frame(model.table(models))
models.table
write.csv(models.table, "modelR2.csv")
#This is to get a specific model:
#to get the model NAG2017.results$Phi$name.model$p$name.model$results$real
###in general model with lowest AIC is explored
modeltoexplore<-models[[19]]$results$real#everytime a $ is written, a menu comes out with all the models
modeltoexplore<-PN_AllYears7.results$Phi.time.origin.fl.year.dis.p.time.fixed$results$real
#to save survival results into an excel file:
write.csv(modeltoexplore, file = "resultsnew.csv")

###to explore a covariate
mindis=min(PN_AllYears7$dis)
maxdis=max(PN_AllYears7$dis)
dis.values=mindis+(0:30)*(maxdis-mindis)/30
dis.values
Phibydis=covariate.predictions(models[[19]], data= data.frame(dis=dis.values),indices=c(10))

plot(x=c(min(Phibydis$estimates$covdata),min(Phibydis$estimates$covdata),max(Phibydis$estimates$covdata),max(Phibydis$estimates$covdata)),y=c(min(Phibydis$estimates$lcl),max(Phibydis$estimates$ucl),min(Phibydis$estimates$lcl),max(Phibydis$estimates$ucl)),type="n",
     xlab = "Discharge",ylab="Surv")

lines(x=Phibydis$estimates$covdata,y=Phibydis$estimates$estimate,type="l")
polygon(c(Phibydis$estimates$covdata,rev(Phibydis$estimates$covdata)),c(Phibydis$estimates$lcl,rev(Phibydis$estimates$ucl)),col=rgb(0.4, 0.4, 0.4,0.1),border=NA)

?covariate.predictions






mindis=min(PN_AllYears7$fl)
maxdis=max(PN_AllYears7$fl)
fl.values=mindis+(0:30)*(maxdis-mindis)/30
fl.values
Phibydis=covariate.predictions(PN_AllYears7.results$Phi.time.origin.fl.year.dis.p.time.fixed, data= data.frame(fl=fl.values),indices=c(10))
Phibydis2=covariate.predictions(PN_AllYears7.results$Phi.time.origin.fl.year.dis.p.time.fixed, data= data.frame(fl=fl.values),indices=c(65))

plot(x=c(min(Phibydis$estimates$covdata),min(Phibydis$estimates$covdata),max(Phibydis$estimates$covdata),max(Phibydis$estimates$covdata)),y=c(min(Phibydis2$estimates$lcl),max(Phibydis$estimates$ucl),min(Phibydis2$estimates$lcl),max(Phibydis$estimates$ucl)),type="n",
     xlab = "FL",ylab="Surv")

lines(x=Phibydis$estimates$covdata,y=Phibydis$estimates$estimate,type="l")
polygon(c(Phibydis$estimates$covdata,rev(Phibydis$estimates$covdata)),c(Phibydis$estimates$lcl,rev(Phibydis$estimates$ucl)),col=rgb(0.4, 0.4, 0.4,0.1),border=NA)


lines(x=Phibydis2$estimates$covdata,y=Phibydis2$estimates$estimate,type="l")
polygon(c(Phibydis2$estimates$covdata,rev(Phibydis2$estimates$covdata)),c(Phibydis2$estimates$lcl,rev(Phibydis2$estimates$ucl)),col=rgb(0.4, 0.2, 0.1,0.1),border=NA)

?covariate.predictions