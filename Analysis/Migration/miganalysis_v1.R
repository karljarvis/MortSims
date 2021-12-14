##########
# Code to extract migration results from CDPOP to test differences between dispersal functions
# Karl Jarvis August 2014
##########

options(scipen=999) # prevent scientific notation
library(ggplot2)
library(reshape2)
library(plyr)
library(MASS)

source("~/GoogleDriveNAU/MortSims/DataPrep.R")

Nruns = read.csv(invarPaths[1])$mcruns

################
# Extract migration numbers and population sizes from output files
### ONLY RUN ONCE - CREATES A SINGLE FILE CONTAINING ALL THE OUTPUT THAT CAN BE EASILY ACCESSED
# mm = matrix(NA,0,7)
# for (i in 1:length(dataDirs))
# {
#   m = matrix(NA, 0, 7)
#   runDirs = dir(dataPaths[i], pattern="batchrun", full.names = T)
#   for(j in 1:Nruns)
#   {
#     raw = read.csv(file.path(runDirs[j], "output.csv"), fill=T)
#     N = colsplit(raw$Population, "\\|", c("total","subA","subB"))[,1]
#     ms = colsplit(raw$SubpopImmigration, "\\|", c("subA","subB"))
#     ms$subB = gsub("\\|", "", ms$subB)
#     sums = as.numeric(ms$subA) + as.numeric(ms$subB)
#     mmean = mean(sums/N) # mean migration rate
#     ngen = length(N) # generation when sims stop
#     m = rbind(m, c(land=land[i], barr=barr[i], mort=mort[i], mcrun=j-1, 
#                    meann=round(mean(N),1), ngen=ngen, m=mmean))
#   }
#   mm = rbind(mm, m)
#   print(pars[i])
# }
# 
# mdf = data.frame(mm)
# write.csv(mdf, file.path(migPath, "mdf.csv"))

####################
# read stats in
mdf = read.csv(file.path(migPath, "mdf.csv"), row.names=1)

# means of all m data, long format
mmeans = ddply(.data = mdf, .variables = .(land, barr, mort), 
               summarize, mean=round(mean(m, na.rm=T),4), 
               sd=round(sd(m, na.rm=T),4), 
               meann=round(mean(meann, na.rm=T),0), 
               ngen=round(mean(ngen, na.rm=T),0))

# write mmeans to file
write.csv(mmeans, file.path(migPath, "mmeans.csv"))

# subset and scale means and write to file for use in sGD
mb0m0 = mmeans[mmeans$barr == 0 & mmeans$mort == 0,]
mb0mean = data.frame(mmeans[mmeans$barr == 0,], maxmig=rep(mb0m0$mean, each=5))
mscale = data.frame(mb0mean[,1:3], scalemig=mb0mean$mean/mb0mean$maxmig)
write.csv(mscale, mscalePath)

# ##########
# # plotting
# mmeans$mort = factor(mmeans$mort)
# col = scale_colour_brewer(palette="Set1")
# 
# 
# # All parameter sets
# plotAllMeans = ggplot(mmeans, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   facet_wrap(~ land, ncol=4) +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration \n All Sims")
# plotAllMeans
# pdf(paste0(migResults, "plotAllMeans.pdf"), width=10, height=12)
# print(plotAllMeans); dev.off()
# 
# # Panmictic 
# mpan = mmeans[mmeans$land == "pan0",]
# plotPANmeans = ggplot(mpan, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration \n PAN")
# pdf(paste0(migResults, "plotPANMeans.pdf"), width=10, height=10)
# print(plotPANmeans); dev.off()
# 
# # IBD 
# mibd = mmeans[mmeans$land == "ibd1",]
# plotIBDmeans = ggplot(mibd, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration \n IBD")
# pdf(paste0(migResults, "plotIBDmeans.pdf"), width=10, height=10)
# print(plotIBDmeans); dev.off()
# 
# # IBR2 
# mibr2 = mmeans[grep("ibr2", mmeans$land),]
# plotIBR2means = ggplot(mibr2, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   facet_wrap(~ land) +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration \n IBR2")
# pdf(paste0(migResults, "plotIBR2means.pdf"), width=12, height=10)
# print(plotIBR2means); dev.off()
# 
# # IBR4
# mibr4 = mmeans[grep("ibr4", mmeans$land),]
# plotIBR4means = ggplot(mibr4, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   facet_wrap(~ land) +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration \n IBR4")
# pdf(paste0(migResults, "plotIBR4means.pdf"), width=12, height=10)
# print(plotIBR4means); dev.off()
# 
# # PAN, IBD, high and low
# mmeanpan0 = grep("pan0", mmeans$land)
# mmeanibd1 = grep("ibd1", mmeans$land)
# mmeanibr27 = grep("ibr27", mmeans$land)
# mmeanibr47 = grep("ibr47", mmeans$land)
# mmeanibr23 = grep("ibr23", mmeans$land)
# mmeanibr43 = grep("ibr43", mmeans$land)
# mmeansamp = mmeans[c(mmeanpan0, mmeanibd1, mmeanibr27, mmeanibr47, mmeanibr23, mmeanibr43),]
# mmeansamp$land = factor(mmeansamp$land, levels=unique(mmeansamp$land))
# 
# plotSampMeans = ggplot(mmeansamp, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   facet_wrap(~ land, ncol=2) +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration")
# pdf(paste0(migResults, "plotSampMeans.pdf"), width=8, height=10)
# print(plotSampMeans); dev.off()
# 
# # IBD and IBR samps
# mmeanibd_r = mmeans[c(mmeanibd1, mmeanibr27, mmeanibr47, mmeanibr23, mmeanibr43),]
# mmeanibd_r$land = factor(mmeanibd_r$land, levels=unique(mmeanibd_r$land))
# 
# plotIBD_RMeans = ggplot(mmeanibd_r, aes(barr, mean, group=mort, color=mort)) + 
#   geom_line() + geom_point() + col +
#   facet_wrap(~ land, ncol=1) +
#   ylab("Mean Migration Rate") +
#   xlab("Road Resistance to Crossing (%)") +
#   ggtitle("Migration")
# pdf(paste0(migResults, "plotIBD_RMeans.pdf"), width=4, height=10)
# print(plotIBD_RMeans); dev.off()
# 
# 
# ###########
# # boxplots IBD
# # mibd = mlong[grep("ibd", mlong$land),]
# # mibd$mort = factor(mibd$mort, levels=unique(mibd$mort))
# # ggplot(mibd, aes(pars,m)) + 
# #   geom_boxplot(aes(fill=mort)) +
# #   scale_fill_manual(name = "Mortality", values = mort, labels = mort)
# # 
# 
# ##########
# # # model of migration rates
# # 
# # # log-log plot - slope is power 
# # # semi-log plot if linear, it's exponential
# # # plot most extereme landsacpes next to linear
# # # glm with landscape as random effect
# # 
# # # log of m
# # # Box-Cox transformation
# # mlm = lm(m+1e-10 ~ barr, data=mlong)
# # b = boxcox(mlm, data=mlong)
# # bdf = data.frame(b$x, b$y)
# # lambda = bdf[which.max(bdf[,2]),1] # pull out transformation with highest likelihood
# # 
# # # Do log transformations
# # mlog = ddply(.data=mlong, .variables=.(pars,land,barr,mort), mutate, mlog=log(m+1e-5), 
# #             barrlog=log(1/(barr+1e-5)), mbox=(m^lambda-1)/lambda)
# # mlogMean = ddply(.data=mlog, .variables=.(pars,land,barr,mort,barrlog), summarize, 
# #             meanm=mean(m), meanlog=mean(mlog), meanbox = mean(mbox))
# # mlogMeanIBD = mlogMean[grep("ibd", mlogMean$land),]
# # meanmelt = melt(data = mlogMeanIBD, id.vars = .(pars,land,barr,mort,barrlog), variable.name = "m")
# # 
# # # semi-log plots
# # ggplot(meanmelt, aes(barr, value, group=mort, color=mort)) +
# #   geom_line() + geom_point() + facet_wrap(~ m, scales="free")
# # 
# # # semi-log plots: barrier levels as factors
# # meanmelt$barr = as.factor(meanmelt$barr)
# # ggplot(meanmelt, aes(barr, value, group=mort, color=mort)) +
# #   geom_line() + geom_point() + facet_wrap(~ m, scales="free")
# # # ??? why is the boxcox worse than the others???
# # 
# # # Log-Log plot
# # ggplot(mlogMeanIBD, aes(barrlog, meanlog, group=mort, color=mort)) + 
# #   geom_line() + geom_point()  # Definitely not this one.
# # 
# ##########
# # modeling
# # library(package = MuMIn)
# # options(na.action = "na.fail")   #  prevent fitting models to different datasets
# # 
# # # IBD only
# # mibd = mlog[grep("ibd", mlog$land),]
# # mibdlong = melt(mibd, id=.(pars,land,barr,mort,mcrun,barrlog), variable.name="trans")
# # mibdlong$pars = factor(mibdlong$pars, levels=unique(mibdlong$pars))
# # mibdlong$land = factor(mibdlong$land, levels=unique(mibdlong$land))
# # 
# # # subsets
# # mibd = mibdlong[mibdlong$trans == "m",]
# # mibdlog = mibdlong[mibdlong$trans == "mlog",]
# # mibdbox = mibdlong[mibdlong$trans == "mbox",]
# # 
# # # lm
# # mibdmod = lm(value ~ barr + mort, mibd)
# # mibdlogmod = lm(value ~ barr + mort, mibdlog)
# # mibdboxmod = lm(value ~ barr + mort, mibdbox)
# # 
# # # AIC
# # dibdmod = dredge(mibdmod, subset = (barr&&mort) || !(barr||mort))
# # dibdlogmod = dredge(mibdlogmod, subset = (barr&&mort) || !(barr||mort))
# # dibdboxmod = dredge(mibdboxmod, subset = (barr&&mort) || !(barr||mort))
# # 
# # # Delta AIC full vs intercept only
# # dibdmod$delta[2]
# # dibdlogmod$delta[2]
# # dibdboxmod$delta[2]
# # 
# # ###########
# # # IBR2
# # mibr2sub = mlog[grep("ibr2", mlog$land),]
# # mibr2long = melt(mibr2sub, id=.(pars,land,barr,mort,mcrun,barrlog), variable.name="trans")
# # mibr2long$pars = factor(mibr2long$pars, levels=unique(mibr2long$pars))
# # mibr2long$land = as.numeric(substr(mibr2long$land, 5, 5))
# # 
# # # subsets
# # mibr2 = mibr2long[mibr2long$trans == "m",]
# # mibr2log = mibr2long[mibr2long$trans == "mlog",]
# # mibr2box = mibr2long[mibr2long$trans == "mbox",]
# # 
# # # lm
# # mibr2lm = lm(value ~ barr + mort + land, mibr2)
# # mibr2loglm = lm(value ~ barr + mort + land, mibr2log)
# # mibr2boxlm = lm(value ~ barr + mort + land, mibr2box)
# # 
# # # AIC
# # dibr2lm = dredge(mibr2lm, subset=(barr&&mort) || !(barr||mort||land))
# # dibr2loglm = dredge(mibr2loglm, subset=(barr&&mort) || !(barr||mort||land))
# # dibr2boxlm = dredge(mibr2boxlm, subset=(barr&&mort) || !(barr||mort||land))
# # 
# # # Delta AIC full vs intercept only
# # dibr2lm$delta[3]
# # dibr2loglm$delta[3]
# # dibr2boxlm$delta[3]

# sims that ended early
mmm = mdf[mdf$ngen < 500,]
median(mmm$ngen)

# general migration stats
panpan = mdf[grep("pan",mdf$land),]
ibdibd = mdf[mdf$land == "ibd1",]
ibr2 = mdf[grep("ibr2",mdf$land),]
ibr4 = mdf[grep("ibr4",mdf$land),]
ibr20 = mdf[grep("ibr20",mdf$land),]
ibr4nb = ibr4[ibr4$barr == 0 & ibr4$mort == 0,]

# relationship of road resistance to migration in pan
panrr = panpan[panpan$mort == 0,]
panrr.lm = lm(m ~ barr, data=panrr)
summary(panrr.lm)
plot(panrr$barr, panrr$m)
abline(panrr.lm)

# square of m to make it linear
panrr2.lm = lm(m^2 ~ barr, data=panrr)
summary(panrr2.lm)
plot(panrr$barr, panrr$m^2)
abline(panrr2.lm)

# relationship of roadkill to migration in pan
panrk = panpan[panpan$barr == 0,]
panrk.lm = lm(m ~ mort, data=panrk)
summary(panrk.lm)
plot(panrk$mort, panrk$m)
abline(panrk.lm)

# relationship of road resistance to migration in ibd
ibdrr = ibdibd[ibdibd$mort == 0,]
ibdrr.lm = lm(m ~ barr, data=ibdrr)
summary(ibdrr.lm)
plot(ibdrr$barr, ibdrr$m)
abline(ibdrr.lm)

ibdrrsq.lm = lm(sqrt(m) ~ barr, data=ibdrr)
summary(ibdrrsq.lm)
plot(ibdrr$barr, sqrt(ibdrr$m))
abline(ibdrrsq.lm)

# relationship of roadkill to migration in ibd
ibdrk = ibdibd[ibdibd$barr == 0,]
ibdrk.lm = lm(m ~ mort, data=ibdrk)
summary(ibdrk.lm)
plot(ibdrk$mort, ibdrk$m)
abline(ibdrk.lm)




ibrm = mdf[grep("ibr",mdf$land),]
ibrrr = ibrm[ibrm$mort == 0,]
ibrrr.lm = lm(m ~ barr, data=ibrrr)
summary(ibrrr.lm)
plot(ibrrr$barr, ibrrr$m)
abline(ibrrr.lm)

ibrrrsq.lm = lm(sqrt(m) ~ barr, data=ibrrr)
summary(ibrrrsq.lm)
plot(ibrrr$barr, sqrt(ibrrr$m))
abline(ibrrrsq.lm)

# relationship of roadkill to migration in ibd
ibrrk = ibrm[ibrm$barr == 0,]
ibrrk.lm = lm(m ~ mort, data=ibrrk)
summary(ibrrk.lm)
plot(ibrrk$mort, ibrrk$m)
abline(ibrrk.lm)

mmort = mdf[mdf$barr == 0,]
mbarr = mdf[mdf$mort == 0,]
mnc = rbind(mmort, mbarr[mbarr$barr != 0,])

# mixed models of migration
library(nlme)
lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc,
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("ib",mnc$land),],
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("ibr",mnc$land),],
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("ibr2",mnc$land),],
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("ibr4",mnc$land),],
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("ibd",mnc$land),],
    method = "ML")

lme(fixed = m ~ barr + mort, 
    random = ~ 1 | land, 
    data = mnc[grep("pan",mnc$land),],
    method = "ML")


# combine into table
mmodtab = rbind(c("PAN", "roadkill", summary(panrk.lm)$adj.r.squared, AICc(panrk.lm)),
                c("PAN", "roadresist", summary(panrr.lm)$adj.r.squared, AICc(panrr.lm)),
                c("PAN", "roadresist_2", summary(panrr2.lm)$adj.r.squared, AICc(panrr2.lm)),
                c("IBD", "roadkill", summary(ibdrk.lm)$adj.r.squared, AICc(ibdrk.lm)),
                c("IBD", "roadresist", summary(ibdrr.lm)$adj.r.squared, AICc(ibdrr.lm)),
                c("IBD", "roadresist_sqrt", summary(ibdrrsq.lm)$adj.r.squared, AICc(ibdrrsq.lm)))
write.csv(mmodtab, file.path(migPath, "mmodtable.csv"))




# both rk and rr
pan.lm = lm(m ~ mort + barr, data=panpan)
summary(pan.lm)
library(rgl)
plot3d(panpan$barr, panpan$mort, panpan$m)

# Landguth et al 2010 relationship of m, T and B
plot((3000-panrr$barr)/(6000-panrr$barr), panrr$m)
plot((100-panrk$mort)/(200-panrk$mort), panrk$m)
plot(panrk$mort, panrk$m)


#####
# tests to distinguish rk and rr
pan.aov = aov(m ~ barr*mort, data=panpan)
summary(pan.aov)


#####
# IBD


#####
# tests to distinguish rk and rr
ibd.aov = aov(m ~ barr*mort, data=ibdibd)
summary(ibd.aov)

ibr20.aov = aov(m ~ barr * mort, data=ibr20)
summary(ibr20.aov)

landlevs = levels(mdf$land)
maov = vector("list", length(landlevs))
names(maov) = landlevs
for(i in 1:length(landlevs))
{
  mdata = mdf[mdf$land == landlevs[i],]
  maov[[i]] = summary(aov(m ~ barr * mort, data=mdata))
}

# differences between no-barrier and road resistance
b0m0 = mmeans[mmeans$barr == 0 & mmeans$mort == 0,]
b50m0 = mmeans[mmeans$barr == 1500 & mmeans$mort == 0,]
b0m50 = mmeans[mmeans$barr == 0 & mmeans$mort == 50,]
b0m0$b50 = 1-b50m0$mean/b0m0$mean
b0m0$m50 = 1-b0m50$mean/b0m0$mean


