# Folder locations
CDPOPpath = "/scratch/kj375/CDPOP"
date = "20150323"
dataPath = "/scratch/kj375/CDPOP/mmdata"
dataPaths = list.dirs(dataPath, recursive=F)        # results files
dataDirs = list.dirs(dataPath, recursive=F, full.names=F)

# rastPath = file.path(CDPOPpath, paste0("sims_", date), "Rasters")
# rasterFiles <- list.files(rastPath, full.names=TRUE) 

# figPath = "~/Projects/MortSims/Manuscript/figures"

invarPaths = list.files(dataPath, pattern="invar", full.names=TRUE ) 

migPath = "/scratch/kj375/Migration"
mscalePath = file.path(migPath, "mscale.csv")

GDpath = "/scratch/kj375/GD"

xygridfile = "xygrid800.csv"
xycropfile = "xycrop476.csv"
xygrid <- read.csv(file.path(dataPath, xygridfile))
xycrop <- read.csv(file.path(dataPath, xycropfile))

Ngrid <- nrow(xygrid)
Ncrop <- nrow(xycrop)

# Factors
pars = substr(dataDirs, 5, nchar(dataDirs)-10)
land = substr(pars, 0, nchar(pars)-10)
barr = substr(pars, nchar(pars)-7, nchar(pars)-4)
mort = substr(pars, nchar(pars)-2, nchar(pars))