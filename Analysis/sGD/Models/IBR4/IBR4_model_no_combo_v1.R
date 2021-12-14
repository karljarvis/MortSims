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
ribr4 = rb0000[grep("ibr4",rb0000)]
rstack = stack(ribr4)
names(rstack) = gsub("_b0000","",names(rstack))
xy = SpatialPoints(cbind(xygrid$XCOORD, xygrid$YCOORD))
xyvals = extract(rstack, xy)
xydfvals = data.frame(x=xygrid$XCOORD, y=xygrid$YCOORD, xyvals)
xymeltvals = melt(xydfvals, id=c("x","y"), variable.name = "land", value.name = "cost")
df = merge(sGDnc, xymeltvals, by=c("x","y","land"), all=T)
df = na.omit(df)
df$lr = factor(paste(df$land, df$barr, df$mort, sep="_"))
df$land = droplevels(df$land)

# E side of landscapes
dfs = df[df$x %in% c(300,1500,2700,3900,5100,6300,7500,8700,9900,11100) & 
     df$y %in% c(300,1500,2700,3900,5100),]

# Subset and scale IBR4
scalevars = c("x","y","barr","mort","N","cost")
ibr4 = dfs[grep("ibr4",dfs$land),]
ibr4[,scalevars] = data.frame(scale(ibr4[,scalevars]))

# Multicollinearity
ibr4.lm = lm(Ho ~ barr + mort + x + y + N + cost, data=ibr4)
vif(ibr4.lm)
cor(ibr4[, c("x","y","barr","mort","N","cost")])

############################################################

######
# A 
Am1.1 = lme(fixed = A ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 11 pars
Am1.2 = lme(fixed = A ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 8 pars
Am1.3 = lme(fixed = A ~ barr + mort, random = ~ N | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 10 pars
Am1.4 = lme(fixed = A ~ barr + mort, random = ~ N | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 7 pars
Am1.5 = lme(fixed = A ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 16 pars
Am1.6 = lme(fixed = A ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim")) # 10 pars

# Fit model of spatial autocorrelation 
Am1.1exp = update(Am1.1, correlation = corExp(1, form = ~ x + y | land/road))
Am1.5exp = update(Am1.5, correlation = corExp(1, form = ~ x + y | land/road))
Am1.1gau = update(Am1.1, correlation = corGaus(1, form = ~ x + y | land/road))
Am1.5gau = update(Am1.5, correlation = corGaus(1, form = ~ x + y | land/road))
Am1.1sph = update(Am1.1, correlation = corSpher(1, form = ~ x + y | land/road))
Am1.5sph = update(Am1.5, correlation = corSpher(1, form = ~ x + y | land/road))
Am1.1lin = update(Am1.1, correlation = corLin(1, form = ~ x + y | land/road))
Am1.5lin = update(Am1.5, correlation = corLin(1, form = ~ x + y | land/road))
Am1.1rat = update(Am1.1, correlation = corRatio(1, form = ~ x + y | land/road))
Am1.5rat = update(Am1.5, correlation = corRatio(1, form = ~ x + y | land/road))

aictab(list(Am1.1,Am1.5,
            Am1.1exp,Am1.5exp,
            Am1.1gau,Am1.5gau,
            Am1.1sph,Am1.5sph,
            Am1.1lin,Am1.5lin,
            Am1.1rat,Am1.5rat),
       c("Am1.1","Am1.5",
         "Am1.1exp","Am1.5exp",
         "Am1.1gau","Am1.5gau",
         "Am1.1sph","Am1.5sph",
         "Am1.1lin","Am1.5lin",
         "Am1.1rat","Am1.5rat"))

Am2.1 = lme(fixed = A ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | land/road)) # 12 pars
Am2.2 = lme(fixed = A ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | lr)) # 9 pars
Am3.1 = lme(fixed = A ~ barr + mort, random = ~ N | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | land/road)) # 11 pars
Am3.2 = lme(fixed = A ~ barr + mort, random = ~ N | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | lr)) # 8 pars
Am5.1 = lme(fixed = A ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | land/road)) # 17 pars
Am5.2 = lme(fixed = A ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | lr)) # 11 pars
Am5.3 = lme(fixed = A ~ 1, random = ~ N + cost | land/road, data = ibr4, 
            method = "ML", control=lmeControl(opt = "optim"), 
            correlation = corExp(1e20, form = ~ x + y | land/road)) # 15 pars

# He
Hem1.1 = lme(fixed = He ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 11 pars
Hem1.2 = lme(fixed = He ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 8 pars
Hem1.3 = lme(fixed = He ~ barr + mort, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 10 pars
Hem1.4 = lme(fixed = He ~ barr + mort, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 7 pars
Hem1.5 = lme(fixed = He ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 16 pars
Hem1.6 = lme(fixed = He ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 10 pars

# Fit model of spatial autocorrelation 
Hem1.1exp = update(Hem1.1, correlation = corExp(10, form = ~ x + y | land/road))
Hem1.5exp = update(Hem1.5, correlation = corExp(10, form = ~ x + y | land/road))
Hem1.1gau = update(Hem1.1, correlation = corGaus(10, form = ~ x + y | land/road))
Hem1.5gau = update(Hem1.5, correlation = corGaus(10, form = ~ x + y | land/road))
Hem1.1sph = update(Hem1.1, correlation = corSpher(10, form = ~ x + y | land/road))
Hem1.5sph = update(Hem1.5, correlation = corSpher(10, form = ~ x + y | land/road))
Hem1.1lin = update(Hem1.1, correlation = corLin(10, form = ~ x + y | land/road))
Hem1.5lin = update(Hem1.5, correlation = corLin(10, form = ~ x + y | land/road))
Hem1.1rat = update(Hem1.1, correlation = corRatio(10, form = ~ x + y | land/road))
Hem1.5rat = update(Hem1.5, correlation = corRatio(10, form = ~ x + y | land/road))

aictab(list(Hem1.1,Hem1.5,
            Hem1.1exp,Hem1.5exp,
            Hem1.1gau,Hem1.5gau,
            Hem1.1sph,#Hem1.5sph,
            Hem1.1lin,#Hem1.5lin,
            Hem1.1rat,Hem1.5rat),
       c("Hem1.1","Hem1.5",
         "Hem1.1exp","Hem1.5exp",
         "Hem1.1gau","Hem1.5gau",
         "Hem1.1sph",#"Hem1.5sph",
         "Hem1.1lin",#"Hem1.5lin",
         "Hem1.1rat","Hem1.5rat"))

Hem2.1 = lme(fixed = He ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 12 pars
Hem2.2 = lme(fixed = He ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 9 pars
Hem3.1 = lme(fixed = He ~ barr + mort, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 11 pars
Hem3.2 = lme(fixed = He ~ barr + mort, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 8 pars
Hem5.1 = lme(fixed = He ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 17 pars
Hem5.2 = lme(fixed = He ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 11 pars
Hem5.3 = lme(fixed = He ~ 1, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 15 pars

# Ho
Hom1.1 = lme(fixed = Ho ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 11 pars
Hom1.2 = lme(fixed = Ho ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 8 pars
Hom1.3 = lme(fixed = Ho ~ barr + mort, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 10 pars
Hom1.4 = lme(fixed = Ho ~ barr + mort, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 7 pars
Hom1.5 = lme(fixed = Ho ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 16 pars
Hom1.6 = lme(fixed = Ho ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim")) # 10 pars

# Fit model of spatial autocorrelation 
Hom1.1exp = update(Hom1.1, correlation = corExp(1, form = ~ x + y | land/road))
Hom1.5exp = update(Hom1.5, correlation = corExp(1, form = ~ x + y | land/road))
Hom1.1gau = update(Hom1.1, correlation = corGaus(1, form = ~ x + y | land/road))
Hom1.5gau = update(Hom1.5, correlation = corGaus(1, form = ~ x + y | land/road))
Hom1.1sph = update(Hom1.1, correlation = corSpher(1, form = ~ x + y | land/road))
Hom1.5sph = update(Hom1.5, correlation = corSpher(1, form = ~ x + y | land/road))
Hom1.1lin = update(Hom1.1, correlation = corLin(1, form = ~ x + y | land/road))
Hom1.5lin = update(Hom1.5, correlation = corLin(1, form = ~ x + y | land/road))
Hom1.1rat = update(Hom1.1, correlation = corRatio(1, form = ~ x + y | land/road))
Hom1.5rat = update(Hom1.5, correlation = corRatio(1, form = ~ x + y | land/road))

aictab(list(Hom1.1,Hom1.5,
            Hom1.1exp,Hom1.5exp,
            #Hom1.1gau,#Hom1.5gau,
            Hom1.1sph,#Hom1.5sph,
            #Hom1.1lin,#Hom1.5lin,
            Hom1.1rat#,Hom1.5rat
            ),
       c("Hom1.1","Hom1.5",
         "Hom1.1exp","Hom1.5exp",
         #"Hom1.1gau",#"Hom1.5gau",
         "Hom1.1sph",#"Hom1.5sph",
         #"Hom1.1lin",#"Hom1.5lin",
         "Hom1.1rat"#,"Hom1.5rat"
         ))

Hom2.1 = lme(fixed = Ho ~ barr + mort + cost, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 12 pars
Hom2.2 = lme(fixed = Ho ~ barr + mort + cost, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 9 pars
Hom3.1 = lme(fixed = Ho ~ barr + mort, random = ~ N | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 11 pars
Hom3.2 = lme(fixed = Ho ~ barr + mort, random = ~ N | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 8 pars
Hom5.1 = lme(fixed = Ho ~ barr + mort, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 17 pars
Hom5.2 = lme(fixed = Ho ~ barr + mort, random = ~ N + cost | lr, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | lr)) # 11 pars
Hom5.3 = lme(fixed = Ho ~ 1, random = ~ N + cost | land/road, data = ibr4, 
             method = "ML", control=lmeControl(opt = "optim"), 
             correlation = corExp(1e20, form = ~ x + y | land/road)) # 15 pars

##########
# Model names
Anames = c("Am1.1","Am1.2","Am1.3","Am1.4","Am1.5","Am1.6",
            "Am2.1","Am2.2","Am3.1","Am3.2","Am5.1","Am5.2","Am5.3")
Henames = c("Hem1.1","Hem1.2","Hem1.3","Hem1.4","Hem1.5","Hem1.6",
            "Hem2.1","Hem2.2","Hem3.1","Hem3.2","Hem5.1","Hem5.2","Hem5.3")
Honames = c("Hom1.1","Hom1.2","Hom1.3","Hom1.4","Hom1.5","Hom1.6",
            "Hom2.1","Hom2.2","Hom3.1","Hom3.2","Hom5.1","Hom5.2","Hom5.3")

##########
# Save models
save(list=Anames, file=file.path(sGDmodpath,"IBR4","Amods.Rdata"))
save(list=Henames, file=file.path(sGDmodpath,"IBR4","Hemods.Rdata"))
save(list=Honames, file=file.path(sGDmodpath,"IBR4","Homods.Rdata"))


##########
# Load models
# A models
load(file.path(sGDmodpath,"IBR4","Amods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Amods_resids.pdf"))
plot(Am5.1); dev.off()
Atab = aictab(list(Am1.1,Am1.2,Am1.3,Am1.4,Am1.5,Am1.6,
            Am2.1,Am2.2,Am3.1,Am3.2,Am5.1,Am5.2,Am5.3), Anames)
sink(file.path(sGDmodpath, "IBR4", "A_full_model.txt"))
summary(Am5.1); sink()
write.csv(Atab, file.path(sGDmodpath,"IBR4","A_model_selection.csv"))
write.csv(summary(Am5.1)$tTable, file.path(sGDmodpath,"IBR4","A_full_model.csv"))

# He models
load(file.path(sGDmodpath,"IBR4","Hemods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Hemods_resids.pdf"))
plot(Hem5.1); dev.off()
Hetab = aictab(list(Hem1.1,Hem1.2,Hem1.3,Hem1.4,Hem1.5,Hem1.6,
                   Hem2.1,Hem2.2,Hem3.1,Hem3.2,Hem5.1,Hem5.2,Hem5.3), Henames)
sink(file.path(sGDmodpath, "IBR4", "He_full_model.txt"))
summary(Hem5.1); sink()
write.csv(Hetab, file.path(sGDmodpath,"IBR4","He_model_selection.csv"))
write.csv(summary(Hem5.1)$tTable, file.path(sGDmodpath,"IBR4","He_full_model.csv"))

# Ho models
load(file.path(sGDmodpath,"IBR4","Homods.Rdata"))
pdf(file.path(sGDmodpath, "IBR4", "Homods_resids.pdf"))
plot(Hom5.1); dev.off()
Hotab = aictab(list(Hom1.1,Hom1.2,Hom1.3,Hom1.4,Hom1.5,Hom1.6,
                   Hom2.1,Hom2.2,Hom3.1,Hom3.2,Hom5.1,Hom5.2,Hom5.3), Honames)
sink(file.path(sGDmodpath, "IBR4", "Ho_full_model.txt"))
summary(Hom5.1); sink()
write.csv(Hotab, file.path(sGDmodpath,"IBR4","Ho_model_selection.csv"))
write.csv(summary(Hom5.1)$tTable, file.path(sGDmodpath,"IBR4","Ho_full_model.csv"))

