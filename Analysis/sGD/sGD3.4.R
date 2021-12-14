##########
# sGD 3.4 
# Based on sGD from Andrew Shirk, University of Washington, Climate Impacts Group, ashirk@uw.edu

##########
# Modified by Karl Jarvis
# June-Aug 2014 changes:
# loop through output of CDPOP simulations
# analyze only cropped subset of CDPOP output
# deal with missing data in CDPOP output
# run via shell scripting using Rscript
# run on NAU High Performance Computing Cluster "monsoon"

# Nov 2014 changes:
# remove individuals from neighborhoods next to roadkill barriers

# 3.1 changes Apr 2015:
# limit radius of neighborhoods to max dispersal (3000)
# rarify neighborhood to standard level (30)
# use gstudio to calculate genetic diversity instead of diveRsity

# 3.2 changes
# don't crop. edge effect is both due to edge of pop and due to edge of sGD analysis
# no extra function to remove roadkill from neighborhoods, because individuals 
# killed have already been removed from populations.

# 3.3 changes
# no sampling within neighborhoods. size already limited because it's the same as 
#   dispersal distance. also, no cropping means edges not dropping off so quickly

# 3.4 changes (2015-06-04)
# change CD mats to take into account effective barrier strength - use CD matrix 
#   that simulates resistance of both avoidance and mortality

##########

# load required packages
library(gstudio)

# # input files
# src = "~/sGD/DataPrep.R"
# dataPath = "/Volumes/JarvisBackup/mmdata_20150323"
# sGDpath = "~/sGD"
# dataDir = "out_ibd1_b0000m0001427556389"
# dataDir = "out_ibr43_b0000m0501427556372"
# dataDir = "out_ibr43_b1500m0501427556371"

args = commandArgs(TRUE)  # Run in Rscript through shell
src = args[1]
dataPath = args[2]
sGDpath = args[3]
dataDir = args[4]

source(src)

# directory of results to analyze
sGDdirPath = file.path(sGDpath, dataDir)
dir.create(sGDdirPath)
dataPathDir = file.path(dataPath, dataDir) # full path to data
dataDirVec = dir(dataPathDir, pattern="batchrun", full.names = T)
invarFile = substr(dataDir, 5, nchar(dataDir)-10) # input variables file name
invarPath = file.path(dataPath, paste0("invar_",invarFile,".csv")) # input variables full path
print(paste("CDPOP Directory", dataPathDir))
print(paste("input variables file", invarPath))

# input variables
invar = read.csv(invarPath) # read input variables file
Ngen = invar$looptime
Nmcruns = invar$mcruns
Nloci = invar$loci
dispdist = invar$matemovethresh

# set CD matrix based on cumulative effect of barr and mort
print("reading cost distance matrix")
levs$netresist = as.character(levs$barr*30)
levs[levs$netresist == 0,3] = "0000"
levs[levs$netresist == 750,3] = "0750"
levs$totresist = as.character(migprob*3000)
levs[levs$totresist == 0,4] = "0000"
levs[levs$totresist == 750,4] = "0750"
barrlev = substr(dataDir, nchar(dataDir)-17, nchar(dataDir)-14)
mortlev = as.numeric(substr(dataDir, nchar(dataDir)-12, nchar(dataDir)-10))
landlev = substr(dataDir, 5, nchar(dataDir)-20)
totbarr = levs[levs$netresist == barrlev & levs$mort == mortlev,4]
CDfilePath = file.path(dataPath, paste0("CD_", landlev, "_b", totbarr, ".csv"))
CDgrid = read.csv(CDfilePath,header=F) # read CD matrix file
      
# effective mort rates derived from migration
mscale = read.csv(mscalePath, row.names=1)
mortmig = mscale[mscale$land == landlev & mscale$mort == mortlev,]$scalemig

# sGD inputs
minhoodsize = 15

# sGD function
sGD = function(sGDpath,sGDfilePath,nhoodFilePath,Nloci,genepopFilePath,xygrid,CDgrid,minhoodsize)
{
# read in CDPOP file and rename loci
  print(paste("reading genetic data", dataFilePath))
  CDPOPgrid = read_population(dataFilePath, type="cdpop") 
  locind = (ncol(CDPOPgrid)-Nloci+1):ncol(CDPOPgrid) 
  occind = which(CDPOPgrid$ID != "OPEN") # indices of occupied points
  names(CDPOPgrid)[locind] = gsub("-","",names(CDPOPgrid)[locind])
  CDPOPgridocc = CDPOPgrid[occind,] # remove missing points
  NgridOcc = nrow(CDPOPgridocc)
  
# define neighborhoods
  print("defining neighborhoods")
  CDgridocc = CDgrid[occind,occind] # remove missing points in CD matrix
  nhoods = matrix(NA,0,NgridOcc)  # create an empty list to fill with the indices of neighborhood members
  nhoodsize = numeric()
  for (i in 1:NgridOcc)
  {
    nbrs = CDgridocc[i,] < dispdist
    nhoodsize = c(nhoodsize,sum(nbrs))
    if (sum(nbrs) < minhoodsize)
    {
      nhoods = rbind(nhoods, FALSE)
#      print(paste("neighborhood < minhoodsize at pt",i))
    } else
    {
      nhoods = rbind(nhoods, nbrs)
#      print(paste("neighborhood for pt",i))
    }
  }

# Calculate genetic diversity
  print("calculating genetic diversity")
  sGDdata = data.frame()
  for(i in 1:nrow(CDPOPgridocc))
  {
    nhood = CDPOPgridocc[nhoods[i,],]
    xy = data.frame(subpop=CDPOPgridocc[i,1],
                    x=CDPOPgridocc[i,2],
                    y=CDPOPgridocc[i,3])
    if (nrow(nhood) > 0)
    {
      data = data.frame(xy,
                        N=nhoodsize[i],
                        A=mean(A(nhood)$A),
                        Ae=mean(Ae(nhood)$Ae),
                        Ho=mean(Ho(nhood)$Ho),
                        He=mean(He(nhood)$He))
    } else
    {
      data = data.frame(xy, N=nhoodsize[i], 
                        A=NA, Ae=NA, Ho=NA, He=NA)
    }
    sGDdata = rbind(sGDdata,data)
#   print(paste("diversity for point",row.names(CDPOPgridocc)[i]))
  }
  
  # write genetic diversity results to csv file
  write.table(sGDdata, sGDfilePath, row.names=F, sep=",", na="")
}

# Loop through CDPOP runs and run sGD function
for(i in 1:Nmcruns)
{
  dataFilePath = list.files(dataDirVec[i], pattern=paste0("grid",Ngen,".csv"), full.names = T)
  sGDfilePath = file.path(sGDdirPath, paste0("sGD_","b",0,"m",i-1,"g",Ngen,".csv"))
    
  sGD(sGDpath,sGDfilePath,nhoodFilePath,Nloci,genepopFilePath,xygrid,CDgrid,minhoodsize)
}


