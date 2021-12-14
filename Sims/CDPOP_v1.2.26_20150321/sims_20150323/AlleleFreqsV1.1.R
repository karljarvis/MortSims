#########################################
# Determine ideal starting allele frequencies for simulations

library(reshape2)
library(ggplot2)

# CDPOP output folder name
CDPOPdir = "/Users/kjj/Projects/MortSims/CDPOP_v1.2.26_20150321"
date = "20150323"
run = "ibd1_b0000m000"
dataDir = file.path(CDPOPdir, paste0("data_", date, "_matemort"))
CDoutDir = dir(dataDir, pattern=run, full.names=T)[1]
invarFile = dir(dataDir, pattern=run, full.names=T)[2]

# pull in input variables
IN = read.csv(invarFile)
Ngen = IN$looptime[1]
# Nseq = IN$nthfile_seq[1]
Nseq = 100
popMax = as.numeric(gsub("xygrid","",IN$xyfilename))
Nalleles = IN$alleles[1]
Nloci = IN$loci[1]
# nthlist = seq(1, Ngen, by=Nseq)
# nthlist = seq(0, Ngen, by=Nseq)
nthlist = c(0, 500)

# pull in data, get means
allFreqs = matrix(NA, Nalleles * Nloci, length(nthlist))
freqMeans = matrix(NA, Nalleles, length(nthlist))
for (i in 1:length(nthlist))
{
  genout = read.csv(file.path(dir(CDoutDir, pattern="batchrun0mcrun0", full.names=T), paste0("grid",nthlist[i],".csv")))[,-1:-8]
  simEndN = popMax - sum(is.na(genout[1]))  # finds ending population size
  allFreqs[,i] = colSums(genout, na.rm=T)/(simEndN*2)

  # order the alleles, find average values 
  freqMat = matrix(allFreqs[,i], Nalleles, Nloci, byrow=F)
  freqSort = apply(freqMat, 2, sort)
  freqMeans[,i] = rowMeans(freqSort)
}

# plot changes in allele frequencies over generations
freqDF = melt(freqMeans)
names(freqDF) = c("allele", "gen", "freq")
freqDF$gen = rep(nthlist, each=Nalleles)
ggplot(freqDF, aes(x=gen, y=freq, group=allele, colour=allele)) + 
  geom_line() + 
  geom_point()

# Make file with means of alleles across all 30 loci at last generation
# freqs = round(freqMeans[,ncol(freqMeans)], 5)
# alleleNames = paste0("L", rep(0:9, each=10), "A", rep(0:9, times=10))
# df = data.frame("Allele List"=alleleNames, "Frequency"=freqs)
# write.table(df, file.path(dataDir, paste0("allelefreqs_", date, ".csv")), sep=",", eol="\n", row.names=F)

output = read.csv(file.path(dir(CDoutDir, pattern="batchrun0mcrun0", full.names=T), "output.csv"))
A = colsplit(output$Alleles, "\\|", c("total","sub1","sub2"))[1]
He = colsplit(output$He, "\\|", c("total","sub1","sub2"))[1]
Ho = colsplit(output$Ho, "\\|", c("total","sub1","sub2"))[1]

plot(A, type="l")
plot(1:500, A$total, type="l")
plot(1:500, He$total, type="l")
points(1:500, Ho$total, type="l", col=2)

