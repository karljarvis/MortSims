# MortSims
This is a set of simulations I put together as part of my dissertation. I used [CDPOP](https://github.com/ComputationalEcologyLab), a Python program written by Erin Landguth for individual-based, spatially-explicit population genetic simulations. The goal was to determine the effects of roads on genetic diversity. I modeled roads in two major ways: as avoidance barriers (wildlife simply bouncing off of the road barrier and continuing on with life), as mortality barriers (wildlife die when in contact with the road). Then I modeled many levels of avoidance and mortality, and many permutations of combined avoidance and mortality. I also modeled these effects across many types of complex, realistic simulated landscapes with a variety of resistances to movement.

The quick version of my results: in terms of genetic differentiation, mortality and avoidance are not hugely different, but they are definitely different in terms of genetic diversity, especially in genetic neighborhoods near roads, and even more so in complex realistic landscapes.

## Analysis
Analyses of the simulations

### Fragstats
I used Kevin McGarigal's [FRAGSTATS](https://www.fs.usda.gov/treesearch/pubs/3064) program to calculate spatial statistics of the complex landscapes I generated.

### GD
Genetic diversity statistics of the simulations. Those with a _slurm_ prefix were analyzed on a parallel computing cluster.

### Migration
Migration statistics: how often the simulated wildlife made it across the road

### sGD
The most advanced approach to genetic diversity -- spatial genetic diversity ([sGD](https://github.com/Andrew-Shirk/sGD)), written by Andrew Shirk to calculate genetic diversity in small local neighborhoods around each individual in a landscape. Requires cost distance matrices (CDmats), models (isolation by distance and isolation by resistance, essentially the rules for how fast an individual can travel across any given part of the landscape), and rasters (rasters of the landscapes, that the models are based on).

## Sims
The spatially explicit wildlife genetic simulations were run here.

### CDPOP_v1.2.26_20150321
The CDPOP program, documentation, the code for setting up the simulations, and some data output from it. The bulk of the results are not on Github due to size constraints. I ran them on a parallel processor to cut ~1 month computing time down to a few hours.

### QRULE
Robert Gardner's [QRULE](https://www.umces.edu/qrule) program for generating rasters for landscape simulations, which I used to generate the complex landscapes.
