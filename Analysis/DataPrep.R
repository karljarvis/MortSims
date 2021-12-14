# Folder locations
mortsimpath = "/Users/kjj/GoogleDriveNAU/MortSims"
anapath = file.path(mortsimpath, "Analysis")
simpath = file.path(mortsimpath, "Sims")
CDPOPpath = file.path(simpath, "CDPOP_v1.2.26_20150321")
date = "20150323"
dataPath = paste0("/Volumes/JarvisBackup/mmdata_", date)
dataPaths = list.dirs(dataPath, recursive=F)        # results files
dataDirs = list.dirs(dataPath, recursive=F, full.names=F)

# sGD folders
sGDpath = paste0(anapath, "/sGD/sGD_", date)
sGDrespath = "/Volumes/JarvisBackup/sGD_20150323"
sGDoutPaths = list.dirs(sGDrespath)[-1]
sGDanapath = file.path(sGDpath, "analysis")
sGDrastPath = file.path(anapath, "sGD", "Rasters")
sGDcdpath = file.path(anapath, "sGD", "CDmats")
sGDmodpath = file.path(anapath, "sGD", "Models")

# Rasters
rastPath = file.path(CDPOPpath, paste0("sims_", date), "Rasters")
rasterFiles <- list.files(rastPath, full.names=TRUE)

# Coordinates for individuals
xygridfile = "xygrid800.csv"
# xycropfile = "xycrop476.csv"
xygrid <- read.csv(file.path(CDPOPpath, "data_20150323", xygridfile))
# xycrop <- read.csv(file.path(dataPath, xycropfile))
Ngrid <- nrow(xygrid)
# Ncrop <- nrow(xycrop)
library(sp)
xy.sp = SpatialPoints(xygrid[,c("XCOORD","YCOORD")])
colnames(xy.sp@coords) = c("x","y")

# Input variables
invarPaths = list.files(dataPath, pattern="invar", full.names=TRUE ) 

# Migration rates
migPath = file.path(anapath, "Migration", "Migration_20150323")
mscalePath = file.path(migPath, "mscale.csv")

# Genetic Differentiation and Diversity analysis
GDpath = file.path(anapath, "GD")

# Figures
figPath = file.path(mortsimpath, "Manuscript/figures")

# Factors
pars = substr(dataDirs, 5, nchar(dataDirs)-10)
land = substr(pars, 0, nchar(pars)-10)
barr = substr(pars, nchar(pars)-7, nchar(pars)-4)
mort = substr(pars, nchar(pars)-2, nchar(pars))
landlevs = c("pan0", "ibd1", paste0("ibr2",0:9), paste0("ibr4",0:9))
typelevs = c("pan0", "ibd1", "ibr2", "ibr4")
