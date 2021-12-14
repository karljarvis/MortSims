# Extract genetic results from sGD results

# Stats: 
# N: number of individuals per neighborhood
# A: Mean number of alleles. How many alleles per locus on average in the neighborhood. Basically how many alleles divided by how many loci. Doesn"t take into account how common alleles are. If they"re there, they"re counted, whether 1% or 50%.
# He: Expected Heterozygosity/Nei"s Gene Diversity. Probability that 2 copies of a gene will have different alleles, average expected heterozygosity. Assumes H-W. If there were random mating, otherwise idealized populations, then we"d expect the alleles present to spread out in this way. I don"t see why we should trust He in a situation that explicitly rejects such kinds of assumptions.
# Ho: Observed Heterozygosity. Proportion of loci that are heterozygous. I like it better than He. It"s messier and less idealized, but it gives het without unrealistic assumptions. I don"t think we should expect the genes to spread out as extensively as He says.
# FIS: Fixation Index/Inbreeding Coefficient. Proportion of variance in subpopulation contained in an individual. High means more inbreeding. Or at least it means that Ho is less than He, and that could be because of fixation or because of rare alleles not having spread out and being localized. Basically more isolation due to landscape resistance. Formula is 1 - Ho/Hs. So it basically shows how different Ho and Hs are. Useful in showing how far the neighborhood is from being in HW. That"s useful, because we"d expect the more complex the landscape, the more different they are. Selkoe says "measure of heterozygote deficit"

#########################################
library(ggplot2)
library(plyr)
library(scales)
library(reshape2)

# Data and objects used by other scripts
source("/Users/kjj/GoogleDriveNAU/MortSims/DataPrep.R")

# Directories
dir.create(sGDanapath)
sGDvars = c("data","plots_A","plots_Ho","plots_He")

# All data
allPath = file.path(sGDanapath, "all")
allPathData = file.path(allPath, "data")
allPlotPaths = file.path(allPath, sGDvars)
dir.create(allPath)
dir.create(allPathData)
sapply(allPlotPaths, dir.create)

# Means of sGD values
meanPath = file.path(sGDanapath, "means")
meanDataPath = file.path(meanPath, "data")
meanPlotPath = file.path(meanPath, "plots")
dir.create(meanPath)
dir.create(meanDataPath)
dir.create(meanPlotPath)

# Percent difference in sGD values
diffPath = file.path(sGDanapath, "diff")
diffDataPath = file.path(diffPath, "data")
diffPlotPath = file.path(diffPath, "plots")
tranPlotPath = file.path(diffPath, "transect_plots")
dir.create(diffPath)
dir.create(diffDataPath)
dir.create(diffPlotPath)
dir.create(tranPlotPath)

# Find numbers of successful sGD runs for each parameter set
sGDn = numeric()
for(i in 1:length(sGDoutPaths))
{ sGDn = c(sGDn, length(list.files(sGDoutPaths[i], pattern=".csv"))) }

# indices of which folders have results
res = which(sGDn != 0)

# factors
sGDdirs = list.files(sGDrespath, full.names=F)[res]
sGDpars = substr(sGDdirs, 5, nchar(sGDdirs)-10)
sGDland = substr(sGDpars, 0, nchar(sGDpars)-10)
sGDbarr = substr(sGDpars, nchar(sGDpars)-7, nchar(sGDpars)-4)
sGDmort = substr(sGDpars, nchar(sGDpars)-2, nchar(sGDpars))

nullindres = grep("b0000m000", sGDpars) # "null" simulations (no barrier)
landlevs = sGDland[nullindres] # landscape types with results

cbind(sGDpars,sGDn[res])

# coordinates
xygridpts = data.frame(x=xygrid$XCOORD, y=xygrid$YCOORD)

#############################################
# Full dataset in one data frame
sGDdata = read.csv(file.path(sGDpath, "analysis/all/sGDdata.csv"), row.names=1)[,-1]
# Remove NA
sGDna = na.omit(sGDdata)
# Take means across simulations at same parameters, keeping points (collapsing runs)
sGDmeanruns = ddply(sGDna, .(land,barr,mort,x,y), summarize, 
                    N=mean(N), A=mean(A), He=mean(He), Ho=mean(Ho))
write.csv(sGDmeanruns, file.path(sGDpath, "sGDmeanruns.csv"))
# Take means by scenarios (collapsing x,y and runs)
sGDmeanxyruns = ddply(sGDna, .(land,barr,mort), summarize, 
                      N=mean(N), A=mean(A), He=mean(He), Ho=mean(Ho))
write.csv(sGDmeanxyruns, file.path(sGDpath, "sGDmeanxyruns.csv"))
# Take means by scenarios and landscape type (collapsing x, y, runs, and IBR2 and IBR4)
sGDna$type = sGDna$land
levels(sGDna$type)[grep("ibr2",levels(sGDna$type))] = "ibr2"
levels(sGDna$type)[grep("ibr4",levels(sGDna$type))] = "ibr4"
sGDtypes = ddply(sGDna, .(type,barr,mort), summarize, 
                 N=mean(N), A=mean(A), He=mean(He), Ho=mean(Ho))
sGDtypes$type = factor(sGDtypes$type,levels=typelevs)
sGDtypes = arrange(sGDtypes, type, barr, mort)
write.csv(sGDtypes, file.path(sGDpath, "sGDtypes.csv"))

sGDmeanruns = read.csv(file.path(sGDpath, "sGDmeanruns.csv"))

# mean across xy points in landscapes, keeping land, barr, mort, run
sGDmeanxy = ddply(sGDna, .(land,barr,mort,run), summarize, 
                  N=mean(N), A=mean(A), He=mean(He), Ho=mean(Ho))
write.csv(sGDmeanxy, file.path(sGDpath, "sGDmeanxy.csv"))

ibr4means = ddply(.data = sGDdata, 
                  .variables = .(land,barr,mort,metrics), 
                  summarize, 
                  value=median(value, na.rm=T))

# #########
# # Combine full data into list by landscape
# reslist = vector("list",length(landlevs))
# names(reslist) = landlevs
# for(i in 1:length(landlevs))
# {
#   reslist[[i]] = do.call(rbind, results[grep(landlevs[i], names(results))])
#   write.csv(reslist[[i]], file.path(allPath, paste0(landlevs[i], ".csv")))
#   print(i)
# }

# load in the best format
reslist = vector("list",length(landlevs))
names(reslist) = landlevs
for(i in 1:length(landlevs))
{
  reslist[[i]] = read.csv(file.path(allPath, paste0(landlevs[i], ".csv")), row.names=1)
  print(i)
}


# Run models on each landscape: A
Amodlist = vector("list",length(landlevs))
Amodres = matrix(NA, 0, 6)
names(Amodlist) = landlevs
for(i in 1:length(landlevs))
{
  Amodlist[[i]] = summary(lm(A ~ abs(x) + abs(y) + barr * mort, data=reslist[[i]]))
  coefs = Amodlist[[i]]$coefficients
  Amodres = rbind(Amodres, cbind(landlevs[i], row.names(coefs), coefs))
  print(i)
}
dimnames(Amodres) = list(1:nrow(Amodres), c("land", "predictor","coefficient","SE","t_val","p_val"))
Amoddf = as.data.frame(Amodres)
write.csv(Amoddf, file.path(sGDpath, "sGDlandmodsA.csv"))

Amodcast = dcast(Amoddf[c("land","predictor","coefficient")], land ~ predictor, value.var = "coefficient")
Amodcast[-1] = sapply(Amodcast[-1], as.numeric)
adjr2 = numeric()
for(i in 1:length(landlevs))
{
  adjr2 = c(adjr2, Amodlist[[i]]$adj.r.squared)
}
Amodcast$adjr2 = adjr2
median(Amodcast$barr)
median(Amodcast$mort)
median(Amodcast$mort - Amodcast$barr)
write.csv(Amodcast, file.path(sGDpath, "sGDlandmodsA.csv"))

# Run models on each landscape: Ho
Homodlist = vector("list",length(landlevs))
Homodres = matrix(NA, 0, 6)
names(Homodlist) = landlevs
for(i in 1:length(landlevs))
{
  Homodlist[[i]] = summary(lm(Ho ~ abs(x) + abs(y) + barr * mort, data=reslist[[i]]))
  coefs = Homodlist[[i]]$coefficients
  Homodres = rbind(Homodres, cbind(landlevs[i], row.names(coefs), coefs))
  print(i)
}
dimnames(Homodres) = list(1:nrow(Homodres), c("land", "predictor","coefficient","SE","t_val","p_val"))
Homoddf = as.data.frame(Homodres)
write.csv(Homoddf, file.path(sGDpath, "sGDlandmodsHo.csv"))

Homodcast = dcast(Homoddf[c("land","predictor","coefficient")], land ~ predictor, value.var = "coefficient")

Homodcast[-1] = sapply(Homodcast[-1], as.numeric)
Homodcast$barrdiff = Homodcast$mort - Homodcast$barr 
adjr2 = numeric()
for(i in 1:length(landlevs))
{
  adjr2 = c(adjr2, Homodlist[[i]]$adj.r.squared)
}
Homodcast$adjr2 = adjr2
median(Homodcast$mort)
median(Homodcast$barr)
median(Homodcast$mort - Homodcast$barr)

##########
# Find means of sGD results
resmeans = vector("list", length(res))
names(resmeans) = sGDpars

# calculate means for all sGD results by parameter set
for(i in 1:length(res))
{
  resmeans[[i]] = ddply(.data = results[[i]], 
    .variables=.(x,y,pars,land,barr,mort), summarize, 
    meanA=mean(A, na.rm=T), sdA=sd(A, na.rm=T),
    meanAe=mean(Ae, na.rm=T), sdAe=sd(Ae, na.rm=T), 
    meanHe=mean(He, na.rm=T), sdHe=sd(He, na.rm=T), 
    meanHo=mean(Ho, na.rm=T), sdHo=sd(Ho, na.rm=T), 
    meanN=mean(N, na.rm=T))
  print(paste("calculate means", sGDpars[i]))
}

# group by landscape type
meanlist = vector("list", length(nullindres))
names(meanlist) = landlevs
for(i in 1:length(meanlist))
{
  # collapses list into data frame if elements of list are identically sized dataframes
  meanlist[[i]] = do.call(rbind.data.frame, resmeans[grep(landlevs[i], names(resmeans))])
  write.csv(meanlist[[i]], file.path(meanDataPath, paste0(landlevs[i], ".csv")))
  print(paste("group by landscape type", landlevs[i]))
}

#############################################
# read in means
meanlist = vector("list", length(landlevs))
names(meanlist) = landlevs
for(i in 1:length(meanlist))
{
  meanlist[[i]] = read.csv(file.path(meanDataPath, paste0(landlevs[i], ".csv")), row.names=1)
  print(paste("read meanlist", landlevs[i]))
}

##########
# plot means
# A
(meanAcolmin = min(unlist(lapply(meanlist, function(x) {min(x$meanA, na.rm=T)}))))
(meanAcolmax = max(unlist(lapply(meanlist, function(x) {max(x$meanA, na.rm=T)}))))
meanAcols = scale_colour_gradientn(limits=c(meanAcolmin, meanAcolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(meanlist))
{
  meanAplot = ggplot(meanlist[[i]], aes(x, y, group = pars)) + meanAcols +
    geom_point(aes(colour = meanA, size = meanN)) + facet_wrap( ~ pars, ncol=5) +
    ggtitle(paste0("Mean Neighborhood A Values\n", landlevs[i])) + coord_equal()
  pdf(file.path(meanPlotPath, paste0(landlevs[i], "_A.pdf")), width=16, height=12)
  print(meanAplot); dev.off()
  print(names(meanlist)[i])
}

# He
(meanHecolmin = min(unlist(lapply(meanlist, function(x) {min(x$meanHe, na.rm=T)}))))
(meanHecolmax = max(unlist(lapply(meanlist, function(x) {max(x$meanHe, na.rm=T)}))))
meanHecols = scale_colour_gradientn(limits=c(meanHecolmin, meanHecolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(meanlist))
{
  meanHeplot = ggplot(meanlist[[i]], aes(x, y, group = pars)) + meanHecols +
    geom_point(aes(colour = meanHe, size = meanN)) + facet_wrap( ~ pars, ncol=5) +
    ggtitle(paste("Mean Neighborhood He Values\n", landlevs[i])) + coord_equal()
  pdf(file.path(meanPlotPath, paste0(landlevs[i], "_He.pdf")), width=16, height=12)
  print(meanHeplot); dev.off()
  print(names(meanlist)[i])
}

# Ho
(meanHocolmin = min(unlist(lapply(meanlist, function(x) {min(x$meanHo, na.rm=T)}))))
(meanHocolmax = max(unlist(lapply(meanlist, function(x) {max(x$meanHo, na.rm=T)}))))
meanHocols = scale_colour_gradientn(limits=c(meanHocolmin, meanHocolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(meanlist))
{
  meanHoplot = ggplot(meanlist[[i]], aes(x, y, group = pars)) + meanHocols +
    geom_point(aes(colour = meanHo, size = meanN)) + facet_wrap( ~ pars, ncol=5) +
    ggtitle(paste("Mean Neighborhood Ho Values\n", landlevs[i])) + coord_equal()
  pdf(file.path(meanPlotPath, paste0(landlevs[i], "_Ho.pdf")), width=16, height=12)
  print(meanHoplot); dev.off()
  print(names(meanlist)[i])
}
# 
# # FIS
# (meanFIScolmin = min(unlist(lapply(meanlist, function(x) {min(x$meanFIS, na.rm=T)}))))
# (meanFIScolmax = max(unlist(lapply(meanlist, function(x) {max(x$meanFIS, na.rm=T)}))))
# meanFIScols = scale_colour_gradientn(limits=c(meanFIScolmin, meanFIScolmax), colours=rainbow(n=100, start=0.15))
# 
# for(i in 1:length(meanlist))
# {
#   meanFISplot = ggplot(meanlist[[i]], aes(X, Y, group = pars)) + meanFIScols +
#     geom_point(aes(colour = meanFIS, size = meanN)) + facet_wrap( ~ pars, ncol=5) +
#     ggtitle(paste("Mean Neighborhood FIS Values\n", landlevs[i])) + coord_equal()
#   pdf(paste0(meanPlotPaths[6], "/", landlevs[i], ".pdf"), width=16, height=12)
#   print(meanFISplot); dev.off()
#   print(names(meanlist)[i])
# }
 
##########
# Normalized genetic diversity based on null landscapes
difflist = vector("list", length(landlevs))
names(difflist) = landlevs
for(i in 1:length(landlevs))
{
  nb = meanlist[[i]][grep("b0000m000", meanlist[[i]]$pars),]
  diff = mutate(meanlist[[i]], 
    diffA=(meanlist[[i]]$meanA-nb$meanA)/nb$meanA, 
    diffAe=(meanlist[[i]]$meanAe-nb$meanAe)/nb$meanAe, 
    diffHe=(meanlist[[i]]$meanHe-nb$meanHe)/nb$meanHe, 
    diffHo=(meanlist[[i]]$meanHo-nb$meanHo)/nb$meanHo)
    difflist[[i]] = diff[complete.cases(diff[,12:13]),] # remove NAs
  write.csv(difflist[[i]], file.path(diffDataPath, paste0(landlevs[i], ".csv")))

  # Show stats on difflist results
  print(paste(nrow(diff), "with NAs"))
  print(paste(nrow(difflist[[i]]), "without NAs"))
  print(difflist[[i]][c(1,nrow(difflist[[i]])),]) # show first and last rows
}

##################################################
# read difflist in
difflist = vector("list", length(landlevs))
names(difflist) = landlevs
for(i in 1:length(landlevs))
{
  difflist[[i]] = read.csv(file.path(diffDataPath, paste0(landlevs[i], ".csv")))
  print(paste("read difflist", landlevs[i])) # show first and last rows
}

diffsum = vector("list", length(landlevs))
diffdiff = diffmeanrr = diffmeanrk = meddiffrr = meddiffrk = numeric()
names(diffsum) = landlevs
for(i in 1:length(landlevs))
{
  diffsum[[i]] = dcast(data = difflist[[i]], barr ~ mort, value.var = "diffA", mean)
  diffdiff = c(diffvec, diffsum[[i]][1,6] - diffsum[[i]][5,2])
  diffmeanrr = c(diffmeanrr, diffsum[[i]][5,2])
  diffmeanrk = c(diffmeanrk, diffsum[[i]][1,6])
  meddiffrr = c(meddiffrr, diffsum[[i]][3,2]/diffsum[[i]][5,2])
  meddiffrk = c(meddiffrk, diffsum[[i]][1,4]/diffsum[[i]][1,6])
}

median(diffmeanrr, na.rm=T)
median(diffmeanrk, na.rm=T)
median(meddiffrr, na.rm=T)
median(meddiffrk, na.rm=T)


rdnondiff = vector("list", length(landlevs))
names(rdnondiff) = landlevs
rdnonrr = rdnonrk = numeric()
for(i in 16:length(landlevs))
{
  rd = difflist[[i]][abs(difflist[[i]]$x) < 500,]
  non = difflist[[i]][abs(difflist[[i]]$x) > 500,]
  rddiff = dcast(data = rd, barr ~ mort, value.var = "diffA", mean)
  nondiff = dcast(data = non, barr ~ mort, value.var = "diffA", mean)
  rdnondiff[[i]] = rddiff-nondiff
  rdnonrr = c(rdnonrr, rdnondiff[[i]][5,2])
  rdnonrk = c(rdnonrk, rdnondiff[[i]][1,6])
}

median(rdnonrr, na.rm=T)
median(rdnonrk, na.rm=T)
###########
# plot normalized means

# A
(diffAcolmin = min(unlist(lapply(difflist, function(x) {min(x$diffA, na.rm=T)}))))
(diffAcolmax = max(unlist(lapply(difflist, function(x) {max(x$diffA, na.rm=T)}))))
diffAcols = scale_colour_gradientn(limits=c(diffAcolmin, diffAcolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(difflist))
{
  diffplot = ggplot(difflist[[i]], aes(x, y, group = pars)) + 
    geom_point(aes(colour = diffA, size=10)) + diffAcols + 
    facet_wrap( ~ pars, ncol=5) + coord_equal() +
    ggtitle(paste0("Normalized sGD Mean Values\n", landlevs[i]))
  pdf(file.path(diffPlotPath, paste0(landlevs[i], "_A.pdf")), width=16, height=12)
    print(diffplot); dev.off()
  print(landlevs[i])
}


# He
(diffHecolmin = min(unlist(lapply(difflist, function(x) {min(x$diffHe, na.rm=T)}))))
(diffHecolmax = max(unlist(lapply(difflist[-17], function(x) {max(x$diffHe, na.rm=T)}))))
diffHecols = scale_colour_gradientn(limits=c(diffHecolmin, diffHecolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(difflist))
{
  diffplot = ggplot(difflist[[i]], aes(x, y, group = pars)) + 
    geom_point(aes(colour = diffHe, size=10)) + diffHecols + 
    facet_wrap( ~ pars, ncol=5) + coord_equal() +
    ggtitle(paste0("Normalized sGD Mean Values\n", landlevs[i]))
  pdf(file.path(diffPlotPath, paste0(landlevs[i], "_He.pdf")), width=16, height=12)
    print(diffplot); dev.off()
  print(landlevs[i])
}

# Ho
(diffHocolmin = min(unlist(lapply(difflist, function(x) {min(x$diffHo, na.rm=T)}))))
(diffHocolmax = max(unlist(lapply(difflist[-17], function(x) {max(x$diffHo, na.rm=T)}))))
diffHocols = scale_colour_gradientn(limits=c(diffHocolmin, diffHocolmax), colours=rainbow(n=100, start=0.15))

for(i in 1:length(difflist))
{
  diffplot = ggplot(difflist[[i]], aes(x, y, group = pars)) + 
    geom_point(aes(colour = diffHo, size=10)) + diffHocols + 
    facet_wrap( ~ pars, ncol=5) + coord_equal() +
    ggtitle(paste0("Normalized sGD Mean Values\n", landlevs[i]))
  pdf(file.path(diffPlotPath, paste0(landlevs[i], "_Ho.pdf")), width=16, height=12)
    print(diffplot); dev.off()
  print(landlevs[i])
}

  
###########
# Moustache plots

#####
# Plot overview of all the data
diffdf = do.call(rbind.data.frame, difflist)

# A
allTranAPlot = ggplot(diffdf, aes(x, diffA, color=barr)) + 
  geom_point(size=1) + facet_wrap(~ land, ncol=5)
pdf(file.path(tranPlotPath, "allTranA.pdf"), width=16, height=12)
  print(allTranAPlot); dev.off()

# He
allTranHePlot = ggplot(diffdf, aes(X, diffHe)) + 
  geom_point() + facet_wrap(~ land, ncol=5)
pdf(file.path(tranPlotPath, "allTranHe.pdf"), width=16, height=12)
  print(allTranHePlot); dev.off()

# N
allTranNPlot = ggplot(diffdf, aes(X, meanN)) + geom_point() +
  facet_wrap(~ land, ncol=5)
pdf(paste0(tranPlotPath, "allTranN.pdf"), width=16, height=12)
print(allTranNPlot); dev.off()


#####
# Plot all in separate figures, color by transect
for(i in 1:length(difflist))
{
  difflist[[i]]$y = factor(difflist[[i]]$y)

  # Normalized A by transect
  tranAplot = ggplot(difflist[[i]], aes(x, diffA, group=y)) + 
    geom_line(aes(colour = y)) + facet_wrap(~ pars)
  pdf(file.path(tranPlotPath, paste0(landlevs[i],"TranA.pdf")), width=16, height=12)
  print(tranAplot); dev.off()
  
  # Normalized He by transect
  tranHeplot = ggplot(difflist[[i]], aes(x, diffHe, group=y)) + 
    geom_line(aes(colour = y)) + facet_wrap(~ pars)  
  pdf(file.path(tranPlotPath, paste0(landlevs[i],"TranHe.pdf")), width=16, height=12)
  print(tranHeplot); dev.off()

  # Normalized Ho by transect
  tranHoplot = ggplot(difflist[[i]], aes(x, diffHo, group=y)) + 
    geom_line(aes(colour = y)) + facet_wrap(~ pars)  
  pdf(file.path(tranPlotPath, paste0(landlevs[i],"TranHo.pdf")), width=16, height=12)
  print(tranHoplot); dev.off()

  # Mean N by transect
  tranNplot = ggplot(difflist[[i]], aes(x, meanN, group=y)) + 
    geom_line(aes(colour = y)) + facet_wrap(~ pars)
  pdf(file.path(tranPlotPath, paste0(landlevs[i],"meanN.pdf")), width=16, height=12)
  print(tranNplot); dev.off()
  
  print(landlevs[i])
}

