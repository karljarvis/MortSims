# Extract Fragstats results

library(Rfrag)
fsin = frag.combine(file.path(anapath, "Fragstats"), 
                    inland="Sim_Metrics.land", 
                    inclass="Sim_Metrics.class")
names(fsin)[1] = "land"
fsin$land = as.character(fsin$land)
fsin$land = substr(fsin$land, nchar(fsin$land)-14, nchar(fsin$land)-10)
write.csv(fsin, file.path(anapath, "Fragstats", "Sim_Metrics.csv"))

# Find correlations among fragstats metrics
fsc2 = fsin[,c(11:15)]
write.csv(cor(fsc2), file.path(anapath, "Fragstats", "Sim_Corr_Class2.csv"))

# Keep top 3 least correlated metrics
fs2 = fsin[c("land","GYRATE_AM","CLUMPY.cls_2","COHESION.cls_2")]
names(fs2) = gsub(".cls_2","",names(fs2))
fs = rbind(fs2, cbind(land=gsub("ibr2","ibr4",fs2$land),fs2[,2:4]))

# Add in patch density in road zone
ibrnames = rasterFiles[grep("ibr2", rasterFiles)]
ibrnobar = ibrnames[grep("b0000", ibrnames)]
ibrstack = stack(ibrnobar)
e = extent(-3000, 3000, -6400, 6400)
pdensvals = extract(ibrstack, e)
pdens = apply(pdensvals, 2, function(x) {length(grep(1, x))})/nrow(pdensvals)

# Write to file
fs$pdens = c(pdens,pdens)
write.csv(fs, file.path(anapath, "Fragstats", "Sim_Metrics_for_Model.csv"))
