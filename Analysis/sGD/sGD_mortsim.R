library(sGD)
library(adegenet)
library(gstudio)

# genetic data
panpaths = dataPaths[grep("pan",dataPaths)]
reppaths = dir(panpaths[1], full.names = T, pattern = "batchrun")
gppaths = dir(reppaths, full.names = T, pattern = "genepopgrid500.gen")
gi = read.genepop(gppaths[1], ncode=3L)
gs = read_population(gppaths[1],type = "genepop")
mean(He(gs)$He)
mean(Ho(gs)$Ho)
mean(A(gs)$A)

# coordinates
xy = xygrid[,c(ID="ID",x="XCOORD",y="YCOORD")]

# cost distance matrix
cdpath = dir(dataPath, full.names = T, pattern = "CD_pan0")[1]
cdmat = read.csv(cdpath, header = F)

# radius
rad = 3000

# min_N
min_N = 12

# output file
sGDout = file.path(sGDpath, "test", "test.csv")

# run sGD
source(file.path(sGDpath, "sGD_verbose.R"))
library(hierfstat)

sGD_verbose(genind_obj = gi, xy = xy, dist.mat = cdmat, radius = rad, min_N = min_N, GD_ans = T, file_name = )

NH_hierfstat = hierfstat:::.genind2hierfstat(gi)
nb.alleles(NH_hierfstat)
allelic.richness(NH_hierfstat)
basic.stats(NH_hierfstat)
