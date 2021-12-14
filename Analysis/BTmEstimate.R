m = matrix(1:100,10,10)
d = distance(m)
dm = as.matrix(d)
l = lower(d)
dm[,1]
cg = data.frame(names=DY14names, side=DY14sides, cd=DY14cdmat[,1], gd=DY14gdmat[,1])
v = as.vector(dm)

B=1500
R=0.5
Tmax=3000

Befun = function(B,R,Tmax) { B + R * (Tmax - B)}
Befun(B,R,Tmax)



Tmax=3000
m=0
m=0.5

# Estimate B based on migration rate and dispersal distance for PAN
# Landguth et al 2010
Bl = function(Tmax,m) {Tmax*(2*m-1)/(m-1)}
Tmax-Bl(Tmax,m)

# How to convert to IBD. To develop:
# L (landscape) factor which effectively increases barrier strength
# proportional to landscape-relevant factors
# Probably dependent on the extent of the landscape relative to dispersal distance
m=0.0387

12000
L=2879/Tmax

m=0
Tmax*(2*m-1)/(m-1)

m=0.5
Tmax*(2*m-1)/(m-1)


B=1500
(Tmax-B)/((Tmax^2)-B)
