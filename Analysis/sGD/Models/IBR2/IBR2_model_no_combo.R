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
scalevars = c("x","y","barr","mort","N")
# Subset and scale IBR2
ibr2 = df[grep("ibr2",df$land),]
ibr2[,scalevars] = data.frame(scale(ibr2[,scalevars]))
ibr2$cost = factor(ibr2$cost)

# Multicollinearity
ibr2.lm = lm(A ~ barr + mort + x + y + N + cost, data=ibr2)
vif(ibr2.lm)
cor(ibr2[,c("x","y","barr","mort","N")])

############################################################
# explore distributions
# A is pretty normal
pdf(file.path(sGDmodpath, "IBR2_Anorm.pdf"), width=8)
par(mfrow=c(1,2))
hist(ibr2$A, 1000, main="Histogram IBR2 A")
qqnorm(ibr2$A)
dev.off()

# Ho
pdf(file.path(sGDmodpath, "IBR2_Honorm.pdf"), width=8)
par(mfrow=c(1,2))
hist(ibr2$Ho, 1000, main="Histogram IBR2 Ho")
qqnorm(ibr2$Ho)
dev.off()

ibr2t = ibr2
ibr2t$Ho = 100^ibr2$Ho

pdf(file.path(sGDmodpath, "IBR2_Ho_transformed_norm.pdf"), width=8)
par(mfrow=c(1,2))
hist(ibr2t$Ho, 1000, main="Histogram IBR2 Ho transformed")
qqnorm(ibr2t$Ho)
dev.off()

############################################################
# Fit model of spatial autocorrelation 
# Help from these pages:
# http://www.ats.ucla.edu/stat/r/faq/variogram_lme.htm
# http://www.ats.ucla.edu/stat/r/faq/spatial_regression.htm

######
# Dealing with neighborhood size
l = lme(fixed = Ho ~ 1, 
        random = ~ 1 | land, 
        data = ibr2t, 
        method = "ML", 
        control=lmeControl(opt = "optim"))
l_xy = update(l, correlation = corExp(form = ~ x + y | land/road))
rd_l = update(l, fixed = Ho ~ barr + mort)
rd_l_xy = update(l_xy, fixed = Ho ~ barr + mort)

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

##########
# save models
save(l, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "l.Rdata"))
save(l_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "l_xy.Rdata"))
save(rd_l, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_l.Rdata"))
save(rd_l_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_l_xy.Rdata"))
save(Nl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "Nl.Rdata"))
save(Nl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "Nl_xy.Rdata"))
save(rd_Nl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_Nl.Rdata"))
save(rd_Nl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_Nl_xy.Rdata"))
save(costl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "costl.Rdata"))
save(costl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "costl_xy.Rdata"))
save(rd_costl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costl.Rdata"))
save(rd_costl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costl_xy.Rdata"))
save(costNl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "costNl.Rdata"))
save(costNl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "costNl_xy.Rdata"))
save(rd_costNl, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costNl.Rdata"))
save(rd_costNl_xy, file=file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costNl_xy.Rdata"))

# load in models
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "l.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "l_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_l.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_l_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "Nl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "Nl_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_Nl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_Nl_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "costl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "costl_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costl_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "costNl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "costNl_xy.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costNl.Rdata"))
load(file.path(sGDmodpath, "IBR2", "Ho_transformed", "rd_costNl_xy.Rdata"))

# Model selection table
mods = list(l,l_xy,rd_l,rd_l_xy,
            Nl,Nl_xy,rd_Nl,rd_Nl_xy,
            costl,costl_xy,rd_costl,rd_costl_xy,
            costNl,costNl_xy,rd_costNl,rd_costNl_xy)
modnames = c("l","l_xy","rd_l","rd_l_xy",
             "Nl","Nl_xy","rd_Nl","rd_Nl_xy",
             "costl","costl_xy","rd_costl","rd_costl_xy",
             "costNl","costNl_xy","rd_costNl","rd_costNl_xy")
modtab = aictab(mods, modnames=modnames)
modtab = arrange(modtab, Delta_AICc)

summary(rd_costNl_xy)

write.csv(modtab, file.path(sGDmodpath, "IBR2", "Ho_transformed_model_selection.csv"))
write.csv(summary(rd_costNl_xy)$tTable, file.path(sGDmodpath, "IBR2", "Ho_transformed_full_model.csv"))
