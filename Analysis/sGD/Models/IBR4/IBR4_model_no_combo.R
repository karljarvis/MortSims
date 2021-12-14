############################################################
# Model the effect of road resistance and roadkill on genetic diversity
# taking into account spatial heterogeneity using mixed models
############################################################
# Sep 8, 2015: Omit cases where road resistance and roadkill are both
# greater than zero.

# Load and subset data
library(raster)
library(reshape2)
library(car)
library(nlme)
library(plyr)
library(AICcmodavg)
source("/Users/kjj/GoogleDriveNAU/MortSims/DataPrep.R")

############################################################
# Read in means
sGDmeanruns = read.csv(file.path(sGDpath, "sGDmeanruns.csv"))
sGDmeanruns["X"] = factor(paste(sGDmeanruns$barr, sGDmeanruns$mort, sep="_"))
names(sGDmeanruns)[1] = "road"

# Omit instances where both road resistance and roadkill are > 0 
sGDmort = sGDmeanruns[sGDmeanruns$barr == 0,]
sGDbarr = sGDmeanruns[sGDmeanruns$mort == 0,]
sGDnc = rbind(sGDmort, sGDbarr[sGDbarr$barr != 0,])

# Extract resistances from rasters at points
rb0000 = rasterFiles[grep("b0000",rasterFiles)]
rstack = stack(rb0000)
names(rstack) = gsub("_b0000","",names(rstack))
xy = SpatialPoints(cbind(xygrid$XCOORD, xygrid$YCOORD))
xyvals = extract(rstack, xy)
xydf = data.frame(x=xygrid$XCOORD, y=xygrid$YCOORD, xyvals)
xymelt = melt(xydf, id=c("x","y"), variable.name = "land", value.name = "cost")
df = merge(sGDnc, xymelt, by=c("x","y","land"), all=T)
df = na.omit(df)

# E side of landscapes
df = df[df$x > 0,]
scalevars = c("x","y","barr","mort","N","cost")
# Subset and scale IBR4
ibr4 = df[grep("ibr4",df$land),]
ibr4[,scalevars] = data.frame(scale(ibr4[,scalevars]))
ibr4$cost = factor(ibr4$cost)

# Multicollinearity
ibr4.lm = lm(Ho ~ barr + mort + x + y + N + cost, data=ibr4)
vif(ibr4.lm)
cor(ibr4[,c("x","y","barr","mort","N")])

############################################################
# explore distributions
# A is pretty normal
hist(ibr4$A, 1000, main="Histogram IBR4 A")
qqnorm(ibr4$A)

# Ho
hist(ibr4$Ho, 1000, main="Histogram IBR4 Ho")
qqnorm(ibr4$Ho)

ibr4t = ibr4
ibr4t$Ho = 100^ibr4$Ho

hist(ibr4t$Ho, 1000, main="Histogram IBR4 Ho transformed")
qqnorm(ibr4t$Ho)

############################################################
# Fit model of spatial autocorrelation 
# Help from these pages:
# http://www.ats.ucla.edu/stat/r/faq/variogram_lme.htm
# http://www.ats.ucla.edu/stat/r/faq/spatial_regression.htm

######
# Dealing with neighborhood size
l = lme(fixed = He ~ 1, random = ~ 1 | land, data = ibr4, 
        method = "ML", control=lmeControl(opt = "optim"))
l_xy = update(l, correlation = corExp(form = ~ x + y | land/road))
rd_l = update(l, fixed = He ~ barr + mort)
rd_l_xy = update(l_xy, fixed = He ~ barr + mort)

# Incorporate neighborhood size
Nl = update(l, random = ~ N | land)
Nl_xy = update(l_xy, random = ~ N | land)
rd_Nl = update(rd_l, random = ~ N | land)
rd_Nl_xy = update(rd_l_xy, random = ~ N | land)

# Incorporate cost
costl = update(l, random = ~ cost | land)
costl_xy = update(l_xy, random = ~ cost | land)
rd_costl = update(rd_l, random = ~ cost | land)
rd_costl_xy = update(rd_l_xy, random = ~ cost | land)

# Incorporate cost and N
costNl = update(l, random = ~ cost + N | land)
costNl_xy = update(l_xy, random = ~ cost + N | land)
rd_costNl = update(rd_l, random = ~ cost + N | land)
rd_costNl_xy = update(rd_l_xy, random = ~ cost + N | land)

# lme1 = lme(fixed = A ~ barr + mort, 
#            random = ~ cost + N | land, 
#            correlation = corExp(form = ~ x + y | land/road),
#            data = ibr4, 
#            method = "ML", 
#            control=lmeControl(opt = "optim"))
# 
# lme2 = lme(fixed = A ~ barr + mort, 
#            random = ~ cost + N | land, 
#            correlation = corExp(form = ~ x + y | land:road),
#            data = ibr4, 
#            method = "ML", 
#            control=lmeControl(opt = "optim"))
# 
# lme3 = lme(fixed = He ~ 1, 
#            random = ~ cost + N | land:road, 
#            correlation = corExp(form = ~ x + y | land:road),
#            data = ibr4, 
#            method = "ML", 
#            control=lmeControl(opt = "optim"))

##########
# save models
modnames = c("l","l_xy","rd_l","rd_l_xy",
             "Nl","Nl_xy","rd_Nl","rd_Nl_xy",
             "costl","costl_xy","rd_costl","rd_costl_xy",
             "costNl","costNl_xy","rd_costNl","rd_costNl_xy")

# A models
# save(list=modnames, file=file.path(sGDmodpath,"IBR4","Amods.Rdata"))
load(file.path(sGDmodpath,"IBR4","Amods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Amods_resids.pdf"))
plot(rd_costNl_xy); dev.off()
mods = list(l,l_xy,rd_l,rd_l_xy,
            Nl,Nl_xy,rd_Nl,rd_Nl_xy,
            costl,costl_xy,rd_costl,rd_costl_xy,
            costNl,costNl_xy,rd_costNl,rd_costNl_xy)
modtab = aictab(mods, modnames=modnames)
modtab = arrange(modtab, Delta_AICc)
sink(file.path(sGDmodpath, "IBR4", "A_full_model.txt"))
summary(rd_costNl_xy); sink()
write.csv(modtab, file.path(sGDmodpath,"IBR4","A_model_selection.csv"))
write.csv(summary(rd_costNl_xy)$tTable, file.path(sGDmodpath,"IBR4","A_full_model.csv"))

# He models
# save(list=modnames, file=file.path(sGDmodpath,"IBR4","Hemods.Rdata"))
load(file.path(sGDmodpath,"IBR4","Hemods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Hemods_resids.pdf"))
plot(rd_costNl_xy); dev.off()
mods = list(l,l_xy,rd_l,rd_l_xy,
            Nl,Nl_xy,rd_Nl,rd_Nl_xy,
            costl,costl_xy,rd_costl,rd_costl_xy,
            costNl,costNl_xy,rd_costNl,rd_costNl_xy)
modtab = aictab(mods, modnames=modnames)
modtab = arrange(modtab, Delta_AICc)
write.csv(modtab, file.path(sGDmodpath,"IBR4","He_model_selection.csv"))
write.csv(summary(rd_costNl_xy)$tTable, file.path(sGDmodpath,"IBR4","He_full_model.csv"))

# Ho models
# save(list=modnames, file=file.path(sGDmodpath,"IBR4","Homods.Rdata"))
load(file.path(sGDmodpath,"IBR4","Homods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Homods_resids.pdf"))
plot(rd_costNl_xy); dev.off()
mods = list(l,l_xy,rd_l,rd_l_xy,
            Nl,Nl_xy,rd_Nl,rd_Nl_xy,
            costl,costl_xy,rd_costl,rd_costl_xy,
            costNl,costNl_xy,rd_costNl,rd_costNl_xy)
modtab = aictab(mods, modnames=modnames)
modtab = arrange(modtab, Delta_AICc)
write.csv(modtab, file.path(sGDmodpath,"IBR4","Ho_model_selection.csv"))
write.csv(summary(rd_costNl_xy)$tTable, file.path(sGDmodpath,"IBR4","Ho_full_model.csv"))

# Ho models, data transformed by 100^Ho (residuals not looking any better)
# save(list=modnames, file=file.path(sGDmodpath,"IBR4","Homods_transformed.Rdata"))
load(file.path(sGDmodpath,"IBR4","Homods_transformed.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Homods_transformed_resids.pdf"))
plot(rd_costNl_xy); dev.off()
mods = list(l,l_xy,rd_l,rd_l_xy,
            Nl,Nl_xy,rd_Nl,rd_Nl_xy,
            costl,costl_xy,rd_costl,rd_costl_xy,
            costNl,costNl_xy,rd_costNl,rd_costNl_xy)
modtab = aictab(mods, modnames=modnames)
modtab = arrange(modtab, Delta_AICc)
write.csv(modtab, file.path(sGDmodpath,"IBR4","Ho_transformed_model_selection.csv"))
write.csv(summary(rd_costNl_xy)$tTable, file.path(sGDmodpath,"IBR4","Ho_transformed_full_model.csv"))

