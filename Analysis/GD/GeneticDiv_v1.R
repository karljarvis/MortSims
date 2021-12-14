#################
# gstudio approach
library(gstudio)

# directory of results to analyze
source("/Users/kjj/GoogleDriveNAU/MortSims/DataPrep.R")

# load in results, calculate genetic divergence
Gstats = matrix(NA,0,8)
for(i in 1:length(dataPaths))
{
  batchmcrunnames = list.dirs(dataPaths[i])[-1]
#   print(pars[i])
  for(j in 1:Nbatchruns)
  {
    Ngen = INex$looptime[j]
    Nmcruns = INex$mcruns[j]
    Nloci = INex$loci[j]
    for(k in 1:Nmcruns)
    {
      mcrun = batchmcrunnames[grep(paste0("batchrun", j-1), batchmcrunnames)]
      genfile = paste0(mcrun, "/genepopgrid", Ngen, ".gen")
      if (file.exists(genfile) == TRUE) {
      gdata = read_population(genfile, type = "genepop")
      Dps = suppressMessages(genetic_distance(gdata, "Population", mode="Dps"))[1,2]
      Gst = Gst(gdata)[Nloci+1,2]
      Gstp = Gst_prime(gdata)[Nloci+1,2]
      gdf = c(land=land[i], barr=barr[i], mort=mort[i], batchrun=j-1, mcrun=k-1, Dps=Dps, Gst=Gst, Gstp=Gstp)
      Gstats = rbind(Gstats, gdf) 
      }
      print(mcrun)
    }
  }
} 

Gstats
# write.csv(Gstats, file.path(outDir, "GeneticDivergenceResults.csv"), row.names=F)


GDpath = "/Users/kjj/GoogleDriveNAU/MortSims/Analysis/GD/GD_20150323"
GDpaths = dir(GDpath, pattern = "GD_out", full.names = T)
GDfilenames = dir(GDpath, pattern = "GD_out", full.names = F)
GDnames = substr(GDfilenames, 8, nchar(GDfilenames)-14)
GDdata = data.frame()
for(i in 1:length(GDpaths))
{
  GD = read.csv(GDpaths[i])
  if(nrow(GD) > 0)
  { GDdata = rbind(GDdata, cbind(pars=GDnames[i], GD)) }
  print(i)
}
GDdata$pars = as.character(GDdata$pars)
GDdata$land = substr(GDdata$pars, 0, nchar(GDdata$pars)-10)
GDdata$barr = as.numeric(substr(GDdata$pars, nchar(GDdata$pars)-7, nchar(GDdata$pars)-4))
GDdata$mort = as.numeric(substr(GDdata$pars, nchar(GDdata$pars)-2, nchar(GDdata$pars)))
write.csv(GDdata, file.path(GDpath, "GDall.csv"))

GDdata = read.csv(file.path(GDpath, "GDall.csv"), row.names=1)
head(GDdata)

GDmean = read.csv(file.path(GDpath, "GDmean.csv"), row.names=1)

library(dplyr)
GD0 = filter(GDmean, barr==0, mort==0)


dim(panpan)
dim(pangd)
pangd = GDdata[grep("pan",GDdata$pars),]
pangd
