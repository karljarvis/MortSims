##########
# make batch files to use on Monsoon cluster

library(reshape2)
CDPOPdir = "/Users/kjj/Projects/MortSims/CDPOP_v1.2.26_20150321"
date = "20150323"
dataDir = file.path(CDPOPdir, paste0("data_", date, "_matemort"))

# cost distance matrices
CDmats = gsub(".csv", "", list.files(dataDir, pattern = "CD_"))
landbarr = data.frame(colsplit(CDmats, "_", c("land","barr")), CDmats)
lands = levels(as.factor(landbarr$land))
barrs = levels(as.factor(landbarr$barr))

# mortality levels
subpopmortpercs = c("0|0","25|25","50|50","75|75","100|100")
morts = paste0("m", c("000","025","050","075","100"))

disp = 3000

# static levels
xyfilename = "xygrid800"
agefilename = "Agevars_nonOverlap.csv"
mcruns = 10
looptime = 500
nthfile_choice = "list"
nthfile_list = "0|500"
nthfile_seq = 0
gridformat = "genepop"
gendmatans = "Dps"
cdclimgentime = 0

matemoveno = 1
matemoveparA = 0
matemoveparB = 0
matemoveparC = 0
matemovethresh = disp
sexans = "Y"
Freplace = "N"
Mreplace = "N"
philopatry = "N"
multiple_paternity = 0
selfans = "N"
Fdispmoveno = 1
FdispmoveparA = 0
FdispmoveparB = 0
FdispmoveparC = 0
Fdispmovethresh = disp
Mdispmoveno = 1
MdispmoveparA = 0
MdispmoveparB = 0
MdispmoveparC = 0
Mdispmovethresh = disp
offno = 2
Femalepercent = 50
EqualsexratioBirth = "N"
popModel = "exp"
r = 0
K_env = 0

muterate = 0.001
mutationtype = "random"
loci = 10
intgenesans = "file"
allefreqfilename = "allelefreqs_20150323.csv"
alleles = 10
mtdna = "N"
startGenes = 50
cdevolveans = "N"
startSelection = 0
Fitness_AA = NA
Fitness_Aa = NA
Fitness_aa = NA
Fitness_AABB = NA
Fitness_AaBB = NA
Fitness_aaBB = NA
Fitness_AABb = NA
Fitness_AaBb = NA
Fitness_aaBb = NA
Fitness_AAbb = NA
Fitness_Aabb = NA
Fitness_aabb = NA
cdinfect = "N"
transmissionprob = 0

for(i in 1:length(CDmats))
{
  matecdmat = dispcdmat = CDmats[i]
  CDname = gsub("CD_","",CDmats[i])
  for(j in 1:length(morts))
  {
    subpopmortperc = subpopmortpercs[j]
    invars = data.frame(xyfilename,agefilename,mcruns,looptime,nthfile_choice,nthfile_list,nthfile_seq,
                        gridformat,gendmatans,cdclimgentime,matecdmat,dispcdmat,matemoveno,matemoveparA,
                        matemoveparB,matemoveparC,matemovethresh,sexans,Freplace,Mreplace,philopatry,
                        multiple_paternity,selfans,Fdispmoveno,FdispmoveparA,FdispmoveparB,FdispmoveparC,
                        Fdispmovethresh,Mdispmoveno,MdispmoveparA,MdispmoveparB,MdispmoveparC,
                        Mdispmovethresh,offno,Femalepercent,EqualsexratioBirth,popModel,r,K_env,
                        subpopmortperc,muterate,mutationtype,loci,intgenesans,allefreqfilename,alleles,
                        mtdna,startGenes,cdevolveans,startSelection,Fitness_AA,Fitness_Aa,Fitness_aa,
                        Fitness_AABB,Fitness_AaBB,Fitness_aaBB,Fitness_AABb,Fitness_AaBb,Fitness_aaBb,
                        Fitness_AAbb,Fitness_Aabb,Fitness_aabb,cdinfect,transmissionprob)
    write.table(invars, file.path(dataDir, paste0("invar_", CDname, morts[j], ".csv")), 
                sep=",", eol="\n", quote=F, row.names=F)          
  }
  print(CDname)
}
