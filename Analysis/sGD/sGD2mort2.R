##########
# sGD 2.0.1 30 May 2014 - Andrew Shirk, University of Washington, Climate Impacts Group, ashirk@uw.edu

##########
# Modified by Karl Jarvis, June-Aug 2014
# Changes:
# loop through output of CDPOP simulations
# analyze only cropped subset of CDPOP output
# deal with missing data in CDPOP output
# run via shell scripting using Rscript
# run on NAU High Performance Computing Cluster "monsoon"
##########

rm(list=ls())

# load required packages
library(diveRsity)
library(adegenet)
library(ecodist)

# How to run
# args = commandArgs(TRUE)  # Run in Rscript through shell
# args = c("ibd1b0000m100d30001419985604", "invar_ibd1b0000m100d3000.csv") # Use with regular R (first one in list, ibd with complete barrier)
args = c("ibd1b3000m000d30001419985605", "invar_ibd1b3000m100d3000.csv") # Use with regular R (first one in list, ibd with complete barrier)
# args = c("ibd1b3000m100d30001419985605", "invar_ibd1b3000m100d3000.csv") # Use with regular R (first one in list, ibd with complete barrier)
# args = c("ibr40b0000m050d30001414623460", "invar_ibr40b0000m050d3000.csv") # Use with regular R (first one in list, ibd with complete barrier)


# Run on monsoon
# CDdir = "/scratch/kj375/CDPOP/mdata/"  # path of folder containing CDPOP output
# sGDdir = "/scratch/kj375/sGD/results/"  # path of folder to contain sGD output

# Run on laptop
CDdir = "/Volumes/JarvisBackup/mdata20141230/" # path of folder containing CDPOP output
sGDdir = "~/Projects/MortSims/sGD/sGD20141230/testresults/" # path of folder to contain sGD output

# Input files
CDoutResults = args[1] # CDPOP output folders like "ibdb0m0d24001405192897"
CDoutDir = paste0(CDdir,CDoutResults) # path of CDPOP results
invarFile = paste0(CDdir, args[2]) # path of CDPOP input variable files like "invar_ibdb0m0d2400.csv"
output_dir = paste0(sGDdir, CDoutResults) # path of output from each sGD run
xy_file = paste0(CDdir, "xy722crop450.csv") # path of file with points of cropped landscape
xy_all = paste0(CDdir, "xygrid722.csv") # path of file with all points in full landscape
fullPop = 722 # maximum population size in full landscape

# input variables
popMax = nrow(read.csv(xy_file)) # find maximum population size in cropped landscape
IN = read.csv(invarFile) # read input variables file
Ngen = IN$looptime
costdistance_file = paste0(CDdir,IN$matecdmat,".csv")
Nseq = IN$nthfile_seq
Nrun = IN$mcruns
outputNames = list.dirs(CDoutDir)[-1]
dir.create(output_dir)
numloci = IN$loci

# crossing mortality rate
mort = as.numeric(substr(args[2], nchar(args[2])-11, nchar(args[2])-9))
land = substr(args[2], 7, nchar(args[2])-18)
# mortlev = as.numeric(substr(mort, 1, (nchar(mort)+1)/2-1))
# mortperc = (1-mortlev)/100

# effective mort rates derived from migration
# migs = read.csv("/home/kj375/sGD/mscale.csv")[c(3:5,9)]
migs = read.csv("/Users/kjj/Projects/MortSims/Migration/results20141230/mscale.csv")[c(3:5,9)]
mortmig = migs[migs$land == land & migs$mort == mort,4]

# Log file
# sink(paste0(output_dir, "/", args[1], ".log"))

# sGD inputs
minhoodsize = 12 
CRSproj = "+proj=utm +zone=12 +datum=WGS84"
GD_ans = T

# sGD function
sGD = function(output_dir,outfilename,neighborhood_file,numloci,genepop_file,xy_file,xy_all,costdistance_file,minhoodsize,CRSproj,GD_ans,corrfile)
{
  # read in points
  xyfile = read.csv(xy_file)
  xyall = read.csv(xy_all)
  Npops = length(unique(xyfile$Subpopulation))
  
  # core area that we want to keep so we can crop the edge effect off
  core = abs(xyall$YCOORD) < 4700 & abs(xyall$XCOORD) < 9200 # 15x15x2=450
  corePop = sum(core) # maximum pop size in core of landscape

  # read in genepop file
  # core reformatted to apply to the right indices of genepop file
  genepopCore = c(rep(TRUE, numloci+2), core[1:(length(core)/2)], TRUE, core[(length(core)/2+1):length(core)])
  genepopfile_in = readLines(genepop_file_N)          # read in genepop file
  genepopfileCoreN  = genepopfile_in[genepopCore]   # subset of genepop in core area
  writeLines(gsub("NANA","000000",genepopfileCoreN), genepop_file) # rewrite genepop, cropped and subbed 0 for NA
  genepopfile = read.genepop(genepop_file)        # read in cropped genepop file with 0 instead of NA
  # the read.genepop function will give a warning about deleting an individual but it does not affect the outcome.

  # begin writing genepop file for sGD output
  genepopfile_lines = readLines(genepop_file)[genepopfileCoreN != "POP"] # read in genepop, removing 2 POP lines
  header_length = 1+numloci                   
  genepop_header = genepopfile_lines[1:header_length]
  write.table(genepop_header,neighborhood_file,sep="\t",quote=F,col.names=F,row.names=F)

  # read in cost distance matrix
  genepop_data_core = genepopfile_lines[(header_length+1):(header_length+corePop)]
  genepop_data_all = genepopfile_in[genepopfile_in != "POP"][-1:-(1+numloci)]
  numindivs = length(genepopfile@ind.names) # number of occupied positions in core area
  openAll = grep("OPEN", genepop_data_all)  # indices of open positions in whole landscape
  openCore = grep("OPEN", genepop_data_core) # indices of open positions in core

  if (numindivs == corePop) # if there are no open positions in the core, it's simpler:
  { 
    occCore = core # indices of all positions in core (all are occupied)
    costdistances = read.csv(costdistance_file,header=F)[core,core]
    genepop_data = genepop_data_core
  } else 
  { 
    # but if there are open positions in the core, we need to crop them out
    occCore = core[-openAll]
    costdistances = read.csv(costdistance_file,header=F)[core,core][-openCore,-openCore]
    genepop_data = genepop_data_core[-openCore]
  }
  cd_dist = as.dist(as.matrix(costdistances))

  # read in genetic distance matrix
  GDist = read.csv(gendistance_file,header=F)     # read gen dist file
  gendistances = GDist[,1:nrow(GDist)][occCore,occCore] # chop off extra blank column, crop to core area
  gd_dist = as.dist(as.matrix(gendistances))      # format to dist object

  # find the limit of mantel correlation and divide in half for neighborhood size
  print("Calculating mantel correlation")
  if (max(costdistances) > 20000) 
  {
    corrbreaks = c(seq(400, 19600, by=400), seq(20000, max(costdistances), by=5000))
  } else
  {
    corrbreaks = c(seq(400, max(costdistances), by=400))
  }
  correlogram = mgram(gd_dist, cd_dist, breaks=corrbreaks, nperm = 99, mrank = F)
  pdf(corrfile); plot(correlogram); dev.off() # plot to file
  mgramAct = correlogram$mgram[correlogram$mgram[,2] > 0 ,] # save only points with observations that support them
  ndist = which(mgramAct[,3] <= 0)[1] - 1 # index of point where correlation hits zero or below
  maxcostdist = as.numeric(mgramAct[ndist,1]) # distance where correlation reaches zero
  halfmax = maxcostdist/2 # half of maxcostdist is suggested by Shirk, keeps corr very high
  write.table(halfmax, file.path(output_dir, paste0(outfilename,".txt")), row.names=F, col.names=F)
  
  # define neighborhoods and write to genepop file
  print("Defining neighborhoods")
  neighborhoods = list()  # create an empty list to fill with the indices of neigbhorhood members
  neighborhood_N = NULL   # create an empty vector to fill with neighborhood sizes  
  counter = 1 # counter needed to track neighborhoods with minhoodsize individuals (i.e. npops in gen file) 
  for (i in 1:numindivs)
  {
    # find indices of neighbors without roadkill
    nbrs = which(costdistances[i,] < halfmax)
    
    # subset neighbors from same/other side, kill off from other
    nsame = nbrs[which(xyfile[nbrs,1] == xyfile[i,1])] # same side
    nother = nbrs[which(xyfile[nbrs,1] != xyfile[i,1])] # other side
    nkill = sample(nother, round(length(nother)*mortmig, 0)) # other side - roadkill
    print(nother)
    print(nkill)
    nact = sort(c(nsame,nkill))   # combine
    
    # add neighborhood size of neighborhood minus roadkill
    neighborhood_N = append(neighborhood_N, length(nact))
    
    # if neighborhoods are big enough 
    if (length(nact)>=minhoodsize)
    {
      # add on roadkill-affected neighborhood
      neighborhoods[[counter]] = nact
      counter = counter+1
    }
  }

  # npops tracks the number of neighborhoods with minhoodsize individuals
  npops = length(neighborhoods)
  
  for (pop in 1:npops)
  {
    write.table("POP", neighborhood_file,sep="\t",quote=F,col.names=F,row.names=F,append=T)
    write.table(genepop_data[neighborhoods[[pop]]], neighborhood_file,sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
  
  # Beginning of genetic diversity output file. Writes names, coordinates, and N of neighborhood. 
  if (length(openCore) == 0) 
  { 
    sGD_output = data.frame(genepopfile@ind.names,xyfile[,c(2,3)],neighborhood_N,stringsAsFactors=F) 
  } else
  { 
    sGD_output = data.frame(genepopfile@ind.names,xyfile[-openCore,c(2,3)],neighborhood_N,stringsAsFactors=F) 
  }
  names(sGD_output) = c("ID","X","Y","N")
  
  # Calculate genetic diversity
  print("Calculating genetic diversity indices")
  if (GD_ans == T)
  {
    GD = divBasic(neighborhood_file, outfile = NULL, gp = 3)
    A = as.numeric(tail(GD[[2]],1)/numloci)
    Ap = as.numeric(tail(GD[[3]],1))
    Ar = as.numeric(tail(GD[[4]],1))
    Ho = as.numeric(tail(GD[[5]],1))
    He = as.numeric(tail(GD[[6]],1))
    FIS = round((He-Ho)/He,3)
    HWE = as.numeric(tail(GD[[7]],1))
    
    A.final = NULL
    Ap.final = NULL
    Ar.final = NULL
    He.final = NULL
    Ho.final = NULL
    FIS.final = NULL
    HWE.final = NULL
    
    counter=1
    
    for (i in 1:numindivs)
    {
      if(neighborhood_N[i]<minhoodsize)
      {
        A.final = append(A.final,NA)
        Ap.final = append(Ap.final,NA)
        Ar.final = append(Ar.final,NA)
        He.final = append(He.final,NA)
        Ho.final = append(Ho.final,NA)
        FIS.final = append(FIS.final,NA)
        HWE.final = append(HWE.final,NA)  
      }
      
      else
      {
        A.final = append(A.final,A[counter])
        Ap.final = append(Ap.final,Ap[counter])
        Ar.final = append(Ar.final,Ar[counter])
        He.final = append(He.final,He[counter])
        Ho.final = append(Ho.final,Ho[counter])
        FIS.final = append(FIS.final,FIS[counter])
        HWE.final = append(HWE.final,HWE[counter])
        counter = counter + 1  
      }
    } 
    sGD_output = data.frame(sGD_output,data.frame(A.final,Ap.final,Ar.final,He.final,Ho.final,FIS.final,HWE.final))
    names(sGD_output)[5:11] = c("A","Ap","Ar","He","Ho","FIS","HWE_p")
  }
  
  # write genetic diversity results to csv file
  write.table(sGD_output,file.path(output_dir,paste0(outfilename,".csv")),row.names=F,sep=",",na="")
}

# Loop through CDPOP runs and run sGD function
for(i in 1:Nrun)
{
  end = paste0("b", 0, "m", i-1, "g", Ngen)
  outfilename = paste0("sGD_", end)
  neighborhood_file = paste0(output_dir, "/", outfilename,"nhood.gen")
  genepop_file_N = paste0(CDoutDir, "/", "batchrun", 0, "mcrun", i-1, "/", "genepopgrid", Ngen, ".gen")
  genepop_file = gsub("[.]gen", "_0.gen", genepop_file_N)
  gendistance_file = paste0(CDoutDir, "/", "batchrun", 0, "mcrun", i-1, "/", "Gdmatrix", Ngen, ".csv")
  corrfile = paste0(output_dir, "/corr_", end, ".pdf")
  
  sGD(output_dir,outfilename,neighborhood_file,numloci,genepop_file,xy_file,xy_all,costdistance_file,minhoodsize,CRSproj,GD_ans,corrfile)
}

