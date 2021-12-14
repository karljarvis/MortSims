#################
# gstudio approach
library(gstudio)

args = commandArgs(TRUE)  # Run in Rscript through shell
# dataPath = "/scratch/kj375/CDPOP/mmdata"
# GDdir = "/scratch/kj375/GD"
# CDPOPoutDir = "out_ibd1_b0000m0001427556389"

dataPath = args[1]
GDdir = args[2]
CDPOPoutDir = args[3]

# directory of results to analyze
CDPOPfullDir = file.path(dataPath, CDPOPoutDir)
invar = paste0("invar_", substr(CDPOPoutDir, 5, nchar(CDPOPoutDir)-10), ".csv")
invarFile = file.path(dataPath, invar)
print(paste("CDPOP Directory", CDPOPfullDir))
print(paste("input variables file", invarFile))

INex = read.csv(invarFile)
Nbatchruns = nrow(INex)

dataDirs = list.dirs(dataPath, recursive=F, full.names=F)
pars = substr(dataDirs, 5, nchar(dataDirs)-10)
land = substr(pars, 0, nchar(pars)-10)
barr = substr(pars, nchar(pars)-7, nchar(pars)-4)
mort = substr(pars, nchar(pars)-2, nchar(pars))

# load in results, calculate genetic divergence
GDstats = matrix(NA,0,12)
for(i in 1:length(CDPOPfullDir))
{
  allrunnames = list.dirs(CDPOPfullDir[i])[-1]
  #   print(pars[i])
  for(j in 1:Nbatchruns)
  {
    Ngen = INex$looptime[j]
    Nmcruns = INex$mcruns[j]
    Nloci = INex$loci[j]
    batchrun = allrunnames[grep(paste0("batchrun", j-1), allrunnames)]
    for(k in 1:Nmcruns)
    {
      mcrun = batchrun[grep(paste0("mcrun",k-1), batchrun)]
      genfile = paste0(mcrun, "/genepopgrid", Ngen, ".gen")
      if (file.exists(genfile) == TRUE) {
        gdata = read_population(genfile, type = "genepop")
        
        A = mean(A(gdata)$A)
        He = mean(He(gdata)$He)
        Ho = mean(Ho(gdata)$Ho)
        Fis = mean(Fis(gdata)$Fis)
        
        Dps = suppressMessages(genetic_distance(gdata, "Population", mode="Dps"))[1,2]
        Gst = Gst(gdata)[Nloci+1,2]
        Gstp = Gst_prime(gdata)[Nloci+1,2]
        Djost = Dest(gdata)[Nloci+1,2]
        
        gdf = c(land=land[i], 
                barr=barr[i], 
                mort=mort[i], 
                mcrun=j-1, 
                A=A,  
                He=He, 
                Ho=Ho, 
                Fis=Fis,
                Dps=Dps, 
                Gst=Gst,
                Gstp=Gstp,
                Djost=Djost)
        GDstats = rbind(GDstats, gdf)
      }
      print(mcrun)
    }
  }
}
write.csv(GDstats, file.path(GDdir, paste0("GD_", CDPOPoutDir,".csv")), row.names=F)

