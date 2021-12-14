#################
# Genetic Diversity and Genetic Differentiation

library(gstudio)
library(reshape2)
library(plyr)
library(MASS)
options(scipen=3)

source("~/GoogleDriveNAU/MortSims/DataPrep.R")

#####
# # Load GD data
# GDmpath = file.path(anapath, "GD", "monsoon", "GD")
# GDpaths = dir(GDmpath, pattern = "GD_out", full.names = T)
# GDfilenames = dir(GDmpath, pattern = "GD_out", full.names = F)
# GDnames = substr(GDfilenames, 8, nchar(GDfilenames)-14)
# GDdata = data.frame() 
# for(i in 1:length(GDpaths))
# {
#   GD = read.csv(GDpaths[i])[,-1:-4]
#   
#   if(nrow(GD) > 0)
#   { GDdata = rbind(GDdata, cbind(land=substr(GDnames[i],0,nchar(GDnames[i])-10), 
#                                  barr=substr(GDnames[i],nchar(GDnames[i])-7,nchar(GDnames[i])-4), 
#                                  mort=substr(GDnames[i],nchar(GDnames[i])-2,nchar(GDnames[i])), 
#                                  mcrun=0:(nrow(GD)-1), GD)) 
#   } 
#   print(i)
# }
# GDdata$land = factor(GDdata$land, levels=landlevs)
# GDdata = arrange(GDdata, land)
# write.csv(GDdata, file.path(anapath, "GD", "GDdata.csv"))

######################################################################
# Analyze
######################################################################
# Import
sGDmeanxy = read.csv(file.path(sGDpath, "sGDmeanxy.csv"), row.names=1)
names(sGDmeanxy)[4:8] = c("mcrun",paste0("sGD",names(sGDmeanxy)[5:8]))
GDdata = read.csv(file.path(GDpath, "GDdata.csv"), row.names=1)
mdf = read.csv(file.path(migPath, "mdf.csv"), row.names=1)
names(mdf)[5] = "N"
GDsGD = merge(GDdata, sGDmeanxy, all=T)
mGD = merge(mdf, GDsGD, all=T)

######################################################################

head(GDdata)
tail(GDdata)

cor(GDdata$barr, GDdata$Gstp, use = "complete.obs")
cor(GDdata$mort, GDdata$Gstp, use = "complete.obs")

# Melt
GDmelt = melt(data = GDdata, id=c("land","barr","mort","mcrun"))

# Mean of each metric for each SCENario
GDscen = ddply(.data = GDmelt, 
               .variables = .(land,barr,mort,variable), 
               summarize, value=mean(value, na.rm=T))
scencast = dcast(GDscen, land + barr + mort ~ variable)
head(arrange(scencast, barr, mort), 22)
write.csv(scencast, file.path(anapath, "GD", "GDmeans.csv"))

# Correlations among means
GDcor = cor(scencast[,4:9], use="complete.obs")
write.csv(GDcor, file.path(anapath, "GD", "GDcor.csv"))
pdf(file.path(anapath, "GD", "GDpairs.pdf"), width = 10, height = 10)
pairs(scencast)
dev.off()

# Median of each metric by LANDscape
GDland = ddply(.data = GDmelt, .variables = .(land,variable), summarize, value=median(value, na.rm=T))
landcast = dcast(GDland, land ~ variable)
write.csv(landcast, file.path(anapath, "GD", "GDlandmedians.csv"))

# Median by landscape TYPE and scenario
levels(GDmelt$land)
GDmelt$type = GDmelt$land
levels(GDmelt$type) = c("ibd1",rep("ibr2",10),rep("ibr4",10),"pan0")
GDtypescen = ddply(.data = GDmelt, .variables = .(type,barr,mort,variable), summarize, value=median(value, na.rm=T))
typescencast = dcast(GDtypescen, type + barr + mort ~ variable)
arrange(typescencast, type, barr, mort)
write.csv(typescencast, file.path(anapath, "GD", "GDtypescenmedians.csv"))

# SD by landscape type and scenario
GDtypescensd = ddply(.data = GDmelt, .variables = .(type,barr,mort,variable), summarize, value=sd(value, na.rm=T))
typescencastsd = dcast(GDtypescensd, type + barr + mort ~ variable)
arrange(typescencastsd, barr, mort)

# Median by landscape type only
GDtype = ddply(.data = GDmelt, .variables = .(type,variable), summarize, value=median(value, na.rm=T))
typecast = dcast(GDtype, type ~ variable)

# SD by landscape type only
GDtypesd = ddply(.data = GDmelt, .variables = .(type,variable), summarize, value=sd(value, na.rm=T))
typecastsd = dcast(GDtypesd, type ~ variable)

######################################################################
# models

# general GD stats
GDpan = GDdata[grep("pan",GDdata$land),]
GDibd = GDdata[GDdata$land == "ibd1",]
GDibr2 = GDdata[grep("ibr2",GDdata$land),]
GDibr4 = GDdata[grep("ibr4",GDdata$land),]
GDibr20 = GDdata[grep("ibr20",GDdata$land),]
GDibr4nb = ibr4[ibr4$barr == 0 & ibr4$mort == 0,]

# relationship of road resistance to A in pan
GDpanrr = GDpan[GDpan$mort == 0,]
GDpanrr.lm = lm(A ~ barr, data=GDpanrr)
summary(GDpanrr.lm)
plot(GDpanrr$barr, GDpanrr$A)
abline(GDpanrr.lm)

##########
# Contrasts
# test if gen div results differ among landscapes
levels(GDdata$land) = landlevs
GD00 = GDdata[GDdata$barr == 0 & GDdata$mort ==0,]
GD00$landsamp = factor(GD00$land, levels=sample(levels(GD00$land), 
                                                size = 22, replace = F))
summary(lm(Gstp~land,data=GD00))
landlm = rbind(cbind("A",summary(lm(A~land,data=GD00))$coefficients),
               cbind("He",summary(lm(He~land,data=GD00))$coefficients),
               cbind("Ho",summary(lm(Ho~land,data=GD00))$coefficients),
               cbind("Fis",summary(lm(Fis~land,data=GD00))$coefficients))
write.csv(landlm, file.path(figPath, "Tab5_landlmbetas.csv"))

GDlowNpan = GDpan[GDpan$barr == 1500 & GDpan$mort == 50,]
GDlowNibr4 = GDibr4[GDibr4$barr == 1500 & GDibr4$mort == 50,]
GDlowN = rbind(GDlowNpan, GDlowNibr4)
GDlowN = GDdata[GDdata$barr == 1500 & GDdata$mort == 50,]
GDlowN$land = factor(GDlowN$land, levels=unique(GDlowN$land))

summary(lm(A~land,data=GDlowN))$coefficients
summary(lm(He~land,data=GDlowN))$coefficients
summary(lm(Ho~land,data=GDlowN))
summary(lm(Fis~land,data=GDlowN))



##########
# Models of influence of road resistance and roadkill
GDnames = names(mGD)[5:12]
GDmort = GDdata[GDdata$barr == 0,]
GDbarr = GDdata[GDdata$mort == 0,]
GDnc = rbind(GDmort, GDbarr[GDbarr$barr != 0,])
GDnc$barr = scale(GDnc$barr)
GDnc$mort = scale(GDnc$mort)

#####
# Run models on each landscape: m
GDlmm = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(m ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmm = rbind(GDlmm, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
}
row.names(GDlmm) = 1:nrow(GDlmm)
names(GDlmm) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmm, file.path(GDpath, "GDlmm.csv"))

GDlmmBetas = dcast(GDlmm, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmmBetas[-1] = sapply(GDlmmBetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmmBetas, file.path(GDpath, "GDlmmBetas.csv"))

GDlmmPvals = dcast(GDlmm, land ~ predictor, value.var = "p_val")
GDlmmPvals[-1] = sapply(GDlmmPvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmmPvals, file.path(GDpath, "GDlmmPvals.csv"))

median(GDlmmBetas$barr)
median(GDlmmBetas$mort)
median(GDlmmBetas$mort - GDlmmBetas$barr)

# Run models on each landscape: A
GDnames = names(GDnc)[5:12]
GDlmA = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(A ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmA = rbind(GDlmA, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
}
row.names(GDlmA) = 1:nrow(GDlmA)
names(GDlmA) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmA, file.path(GDpath, "GDlmA.csv"))

GDlmABetas = dcast(GDlmA, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmABetas[-1] = sapply(GDlmABetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmABetas, file.path(GDpath, "GDlmABetas.csv"))

GDlmAPvals = dcast(GDlmA, land ~ predictor, value.var = "p_val")
GDlmAPvals[-1] = sapply(GDlmAPvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmAPvals, file.path(GDpath, "GDlmAPvals.csv"))

median(GDlmABetas$barr)
median(GDlmABetas$mort)
median(GDlmABetas$mort - GDlmABetas$barr)

#####
# Run models on each landscape: He
GDnames = names(GDnc)[5:12]
GDlmHe = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(He ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmHe = rbind(GDlmHe, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
  print(i)
}
row.names(GDlmHe) = 1:nrow(GDlmHe)
names(GDlmHe) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmHe, file.path(GDpath, "GDlmHe.csv"))

GDlmHeBetas = dcast(GDlmHe, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmHeBetas[-1] = sapply(GDlmHeBetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmHeBetas, file.path(GDpath, "GDlmHeBetas.csv"))

GDlmHePvals = dcast(GDlmHe, land ~ predictor, value.var = "p_val")
GDlmHePvals[-1] = sapply(GDlmHePvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmHePvals, file.path(GDpath, "GDlmHePvals.csv"))

median(GDlmHeBetas$barr)
median(GDlmHeBetas$mort)
median(GDlmHeBetas$mort - GDlmHeBetas$barr)

#####
# Run models on each landscape: Ho
GDnames = names(GDnc)[5:12]
GDlmHo = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(Ho ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmHo = rbind(GDlmHo, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
  print(i)
}
row.names(GDlmHo) = 1:nrow(GDlmHo)
names(GDlmHo) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmHo, file.path(GDpath, "GDlmHo.csv"))

GDlmHoBetas = dcast(GDlmHo, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmHoBetas[-1] = sapply(GDlmHoBetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmHoBetas, file.path(GDpath, "GDlmHoBetas.csv"))

GDlmHoPvals = dcast(GDlmHo, land ~ predictor, value.var = "p_val")
GDlmHoPvals[-1] = sapply(GDlmHoPvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmHoPvals, file.path(GDpath, "GDlmHoPvals.csv"))

median(GDlmHoBetas$barr)
median(GDlmHoBetas$mort)
median(GDlmHoBetas$mort - GDlmHoBetas$barr)


#####
# Run models on each landscape: Fis
GDnames = names(GDnc)[5:12]
GDlmFis = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(Fis ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmFis = rbind(GDlmFis, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
  print(i)
}
row.names(GDlmFis) = 1:nrow(GDlmFis)
names(GDlmFis) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmFis, file.path(GDpath, "GDlmFis.csv"))

GDlmFisBetas = dcast(GDlmFis, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmFisBetas[-1] = sapply(GDlmFisBetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmFisBetas, file.path(GDpath, "GDlmFisBetas.csv"))

GDlmFisPvals = dcast(GDlmFis, land ~ predictor, value.var = "p_val")
GDlmFisPvals[-1] = sapply(GDlmFisPvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmFisPvals, file.path(GDpath, "GDlmFisPvals.csv"))

median(GDlmFisBetas$barr)
median(GDlmFisBetas$mort)
median(GDlmFisBetas$mort - GDlmFisBetas$barr)

#####
# Run models on each landscape: Gstp
GDnames = names(GDnc)[5:12]
GDlmGstp = data.frame()
for(i in 1:length(landlevs))
{
  data = GDnc[GDnc$land == landlevs[i],]
  lmsum = summary(lm(Gstp ~ barr + mort, data=data))
  coefs = lmsum$coefficients
  adjr2 = lmsum$adj.r.squared
  GDlmGstp = rbind(GDlmGstp, cbind(landlevs[i], row.names(coefs), coefs, adjr2))
  print(i)
}
row.names(GDlmGstp) = 1:nrow(GDlmGstp)
names(GDlmGstp) = c("land", "predictor","coefficient","SE","t_val","p_val","adj_r2")
write.csv(GDlmGstp, file.path(GDpath, "GDlmGstp.csv"))

GDlmGstpBetas = dcast(GDlmGstp, land + adj_r2 ~ predictor, value.var = "coefficient")
GDlmGstpBetas[-1] = sapply(GDlmGstpBetas[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmGstpBetas, file.path(GDpath, "GDlmGstpBetas.csv"))

GDlmGstpPvals = dcast(GDlmGstp, land ~ predictor, value.var = "p_val")
GDlmGstpPvals[-1] = sapply(GDlmGstpPvals[-1], function(x) as.numeric(as.character(x)))
write.csv(GDlmGstpPvals, file.path(GDpath, "GDlmGstpPvals.csv"))

median(GDlmGstpBetas$barr)
median(GDlmGstpBetas$mort)
median(GDlmGstpBetas$mort - GDlmGstpBetas$barr)








# Gstp differences between PAN and others on partial barriers
Gstptype = typescencast[c("type","barr","mort","Gstp")]
Gstpmelt = melt(Gstptype, id.vars = c("type","barr","mort"))
Gstpcast = dcast(Gstpmelt, barr + mort ~ type)
Gstpb = Gstpcast[Gstpcast$barr != 3000 & Gstpcast$barr != 0,]
Gstpm = Gstpb[Gstpb$mort != 100,]





median(Gstpm$ibd1 - Gstpm$pan0, na.rm = T)-Gstpcast$ibd1[1] - Gstpcast$pan0[1]
median(Gstpm$ibr2 - Gstpm$pan0, na.rm = T)-Gstpcast$ibr2[1] - Gstpcast$pan0[1]
median(Gstpm$ibr4 - Gstpm$pan0, na.rm = T)-Gstpcast$ibr4[1] - Gstpcast$pan0[1]


tsccast$ibr2 - tsccast$pan0
tsccast$ibr4 - tsccast$pan0

facts = [,1:4] 
fpan = merge(facts, GDpan[,5:12])
GDibd[,5:12]-GDpan[,5:12])


t=typescencast
b0m0 = t[t$barr == 0 & t$mort == 0,]
b50m0 = t[t$barr == 1500 & t$mort == 0,]
b0m50 = t[t$barr == 0 & t$mort == 50,]
b50m0$Gstp - b0m0$Gstp
b0m50$Gstp - b0m0[b0m0$type != "pan0","Gstp"]

# compare IBD with complete roadkill and complete avoidance
GDibdrk = GDibd[GDibd$mort == 100,]
GDibdrr = GDibd[GDibd$barr == 3000,]
t.test(GDibdrk$Fis, GDibdrr$Fis)

GDibr2rk = GDibr2[GDibr2$mort == 100,]
GDibr2rr = GDibr2[GDibr2$barr == 3000,]
t.test(GDibr2rk$Fis, GDibr2rr$Fis)

GDibr4rk = GDibr4[GDibr4$mort == 100,]
GDibr4rr = GDibr4[GDibr4$barr == 3000,]
t.test(GDibr4rk$Fis, GDibr4rr$Fis)

# Test for differences in G'st
GDibdrk = GDibd[GDibd$mort == 100,]
GDibdrr = GDibd[GDibd$barr == 3000,]
t.test(GDibdrk$Gstp, GDibdrr$Gstp)

GDibr2rk = GDibr2[GDibr2$mort == 100,]
GDibr2rr = GDibr2[GDibr2$barr == 3000,]
t.test(GDibr2rk$Gstp, GDibr2rr$Gstp)

GDibr4rk = GDibr4[GDibr4$mort == 100,]
GDibr4rr = GDibr4[GDibr4$barr == 3000,]
t.test(GDibr4rk$Gstp, GDibr4rr$Gstp)



# directory of results to analyze
# source("/Users/kjj/GoogleDriveNAU/MortSims/DataPrep.R")

# INex = read.csv(invarPaths[1])
# Ngen = INex$looptime
# Nmcruns = INex$mcruns
# Nloci = INex$loci

# load in results, calculate genetic divergence
# Gstats = matrix(NA, 0, 9)
# for(i in 1:length(dataPaths))
# {
#   genpaths = list.files(dataPaths[i], 
#                         recursive = T, 
#                         pattern = paste0("genepopgrid", Ngen, ".gen"), 
#                         full.names = T)
#   if(length(genpaths) > 0)
#   {
#     t0 = proc.time()
#     for(j in 1:length(genpaths))
#     {
#       gdata = read_population(genpaths[j], type = "genepop")
#       A = mean(A(gdata)$A)
#       He = mean(He(gdata)$He)
#       Ho = mean(Ho(gdata)$Ho)
#       Dps = suppressMessages(genetic_distance(gdata, "Population", mode="Dps")[1,2])
#       Gstp = Gst_prime(gdata)[Nloci+1,2]
#       gdf = c(land=land[i], barr=barr[i], mort=mort[i], mcrun=j-1, A=A, He=He, Ho=Ho, Dps=Dps, Gstp=Gstp)
#       Gstats = rbind(Gstats, gdf) 
#       print(dataPaths[i])
#     }
#     print(proc.time() - t0)
#   }
# } 
# 
# row.names(Gstats) = 1:nrow(Gstats)
# head(Gstats)

# test whether pan differs from ibr4 in any metrics
head(GDnc)
rbind(GDnc[grep("ibr4",GDnc$land),], GDnc[grep("pan",GDnc$land)]
