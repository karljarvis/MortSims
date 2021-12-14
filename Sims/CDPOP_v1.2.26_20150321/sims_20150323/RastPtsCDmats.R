#####
# Importing and processing fractal patterned rasters from QRULE
# Creating all PAN, IBD, and IBR landscapes for roadkill simulations
# KJJ Jan 29, 2014
#####

library(gdistance)
library(rgdal)
library(ggplot2)

# Folder locations
source("/Users/kjj/GoogleDriveNAU/MortSims/DataPrep.R")

# Import QRULE output and use to set values for rasters
sims = as.matrix(read.table(file.path(simpath, "QRULE/qrule_out10.txt")))
ndim = ncol(sims) # Number of cells per side in the square QRULE outputs
nibr = 10   # Number of landscapes from QRULE
ext = 12800 # Extent in each dimension from center

#########
# Convert to raster format
ibrList = vector(mode="list",length=nibr)
for(i in 1:nibr)
{
  ibrList[[i]] = raster(sims[(ndim*(i-1)+1):(ndim*(i)), 1:ndim])
  extent(ibrList[[i]]) = c(-ext,ext,-ext,ext)		
  proj4string(ibrList[[i]]) = "+proj=utm +datum=WGS84"
}

# Clip extent to be half as tall as wide
ibr = ibrList
extc = c(-ext, ext, -ext/2, ext/2)
for(i in 1:nibr) 
{
  ibr[[i]] = crop(ibr[[i]],extc)
}

##########
# Create rasters with all combinations of parameters
pan0 = ibd1 = ibr[[1]]
values(pan0) = 1e-10  # PAN
values(ibd1) = 1  # IBD

# IBR (patch resistance of 2 and 4)
ibr2 = ibr4 = ibr
for(i in 1:nibr)
{ 
  values(ibr2[[i]]) = as.numeric(gsub(pattern=0, replacement=2, x=values(ibr2[[i]])))
  values(ibr4[[i]]) = as.numeric(gsub(pattern=0, replacement=4, x=values(ibr4[[i]])))
}

# Combine all 22 landscapes into one list
lands = c(pan0, ibd1, ibr2, ibr4)
names(lands) = c("pan0", "ibd1", paste0("ibr2",0:(nibr-1)), paste0("ibr4",0:(nibr-1)))

# barrier location
barrL = seq(ndim/2, ndim/2*ndim, by=ndim)
barrLoc = sort(c(barrL, barrL+1))

# barrier strengths
Tval = 15 # 15*2 (2 pixel wide road) * 100 (100m resolution)
Bval = seq(0, Tval, length.out = 5)
mortvals = seq(0,100,length.out = 5)

BvalsGDall = numeric()
for(i in 1:length(Bval))
{
  for(j in 1:length(mortvals))
  {
    BvalsGDall = c(BvalsGDall, Bval[i] + (max(Bval) - Bval[i]) * mortvals[j]/100)
  }
}

Bvalmat = matrix(BvalsGDall, 5, 5, dimnames = list(Bval, mortvals))
write.csv(Bvalmat, file.path(sGDpath, "effective_Bval_sGD.csv"))

BvalsGD = unique(BvalsGD)[6:10]

(barrVals = BvalsGD*2*res(ibd1)[1])

# Create variable barrier resistances, including
barrNames = paste0("b", round(barrVals,0))
names(barrNames)[1:2] = c("b0000","b0750")

# insert barrier values
rastList = vector("list", length(BvalsGD))
rastList[[1]] = lands
for(i in 2:length(BvalsGD))
{
  l = lands
  for(j in 1:length(lands))
  {
    values(l[[j]])[barrLoc] = Bval[i] 
  }
  rastList[[i]] = l
}
names(rastList) = barrNames

# Write rasters to file
for(i in 1:length(rastList))
{
  for(j in 1:length(lands))
  {
    writeRaster(rastList[[i]][[j]], file.path(sGDrastPath, paste0(names(rastList[[i]][j]), 
                                                                  "_", names(rastList[i]), ".asc")), format="ascii", overwrite=T)
  }
  print(paste0(barrNames[i]))
}

# proportion of landscape in cropped version
ibrs = stack(rastList$b0000[grep("ibr2",names(rastList$b0000))])
patchprop = colSums(ibrs[] == 1)/nrow(ibrs[])
names(patchprop) = gsub("ibr2","ibr_",names(patchprop))
patchpropdf = data.frame(patchprop)
names(patchpropdf) = "prop_in_patch"
write.csv(patchprop, file.path(CDPOPpath, "sims_20150323", "raster_stats_prop_patch.csv"))

##################################
# Generate coordinates for sites #
##################################

# Import sample raster
rast = raster(file.path(rastPath,list.files(rastPath)[6]))
proj4string(rast) = "+proj=utm +datum=WGS84"

###############################
# Simulate points on a grid
# Suggested by Andrew Shirk as a way to ensure that we can distinguish 
# between random variation in sgd due to clustering of points and variation 
# due to actual differences in diversity.

# Full grid
nindside = 20 # number of points per side of square gride on each side
space = 600 # meters between points
locs = seq(space/2, ext, by=space)[1:nindside] # spacing of points
xvals = sort(c(locs, -locs)) 
yvals = sort(c(locs[1:(nindside/2)], -locs[1:(nindside/2)]))
xy = expand.grid(XCOORD=xvals, YCOORD=yvals)
xy = xy[order(xy$XCOORD),]

# CDPOP
Subpopulation = (xy$XCOORD > 0) + 1
ID = paste0("initial", seq(0, nrow(xy)-1))
sex = rep(0:1, length.out=nrow(xy))
xygrid = data.frame(Subpopulation, xy, ID, sex)
write.table(xygrid, file.path(dataDir, paste0("xygrid", nrow(xy), ".csv")), 
            sep=",", eol="\r\n", row.names=FALSE, quote=FALSE)

# Crop outside layer by 3
xycrop = xygrid[abs(xygrid$XCOORD) < 10500 & abs(xygrid$YCOORD) < 4500,]
write.table(xycrop, file.path(dataDir, paste0("xycrop", nrow(xycrop), ".csv")), 
            sep=",", eol="\r\n", row.names=FALSE, quote=FALSE)

# Load in points
xygrid = read.csv(list.files(dataDir, pattern="xygrid", full.names = T))
xycrop = read.csv(list.files(dataDir, pattern="xycrop", full.names = T))

# Plot
plot(rast)
points(xygrid$XCOORD, xygrid$YCOORD, col=xygrid$Subpopulation+3)
points(xycrop$XCOORD, xycrop$YCOORD)


################################################################
# Create cost distance matrices based on simulated landscapes 
################################################################

# Set names of rasters and corresponding CD matrices
rastFiles = list.files(rastPath, pattern=".asc")
rastNames = gsub(".asc", "", rastFiles)

# Set empty list to contain loaded rasters and other steps toward CD matrices
rastList = vector("list",length=length(rastNames))
names(rastList) = rastNames
transList = gcList = CDlist = rastList

# Loop through raster names ane load them into the list
for(i in 1:length(rastNames)) 
{ 
  rastList[[i]] = raster(file.path(rastPath, paste0(rastNames[[i]], ".asc")))
  proj4string(rastList[[i]]) = "+proj=utm +datum=WGS84"
  print(rastNames[i])
}

# Load coordinates for sites
xygrid = read.csv(list.files(dataDir, pattern="xygrid", full.names = T))
sites = SpatialPointsDataFrame(coords=xygrid[,2:3], data=data.frame(xygrid[,4]))

############################
# Create Transition Objects 

# Transition all the rasters
for(i in 1:length(transList))
{ 
  # use 1/mean as transition function to convert resistance to conductance
  transList[[i]] = transition(rastList[[i]], transitionFunction=function(x) 1/mean(x), directions=8) 
  print(rastNames[i])
}

##########################
# Geo Correction 

# Determine geocorrection factor to use on all the transition objects	
gcFactor = geoCorrection(transList[[1]], type="c", multpl=TRUE)

# Geocorrect the transition objects
for(i in 1:length(gcList)) 
{ 
  gcList[[i]] = gcFactor*transList[[i]]
  print(rastNames[i])
}

############################
# Make CD matrices 

# Calculate Cost Distance matrices by least cost paths and write to file
for(i in 1:length(CDlist))
{ 
#   CDlist[[i]] = round(as.matrix(costDistance(gcList[[i]], sites)),1)
  write.table(CDlist[[i]], file.path(dataDir, paste0("CD_", rastNames[i],".csv")), 
              sep=",", eol="\n", row.names=FALSE, col.names=FALSE) 
  print(rastNames[i])
}

###############################
##### Analyze CD matrices #####
###############################
# Read in CDmats as matrix objects

for(i in 1:length(CDlist)) 
{ 
  CDlist[[i]] = as.matrix(read.csv(file.path(dataDir, paste0("CD_", rastNames[i],".csv")), header=F))
  print(rastNames[i])
}

# Give CD matrices the proper names
names(CDmats) = CDnames

# Make matrix of sample of CD matrices 
# bring CD dist objects into one matrix
mat =  matrix(NA, N*(N-1)/2, length(CDnames), dimnames=list(NULL, CDnames))
for (i in 1:length(CDnames)) { mat[,i] = as.dist(CDmats[[i]], diag=F, upper=F) }

# Index values for extremes of each landscape type  
maxno = c(1,5,6,10,56,60,106,110)

# sample 10000 of the 124750 elements of each CDmat to make them manageable
samp = matrix(NA, 10000, length(maxno), dimnames=list(NULL, CDnames[maxno]))
for(i in 1:10000) 
{ 
  samp[i,] = mat[sample(1:nrow(mat),1),maxno] 
}

### Plots	
# Density plots
dev.new()
par(mfrow=c(2,4))
for(i in 1:length(maxno)) 
{	
  plot(density(samp[,i], adjust=0.1), main=dimnames(samp)[[2]][i])	
}

# Pairwise comparisons of CD matrices 
dev.new(); pairs(samp,cex=0.3)

# find max 
max(CDmats$ibdm5)
lapply(CDmats, max)