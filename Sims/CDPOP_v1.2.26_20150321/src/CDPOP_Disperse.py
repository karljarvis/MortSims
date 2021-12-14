# -------------------------------------------------------------------------------------------------
# CDPOP_Disperse.py
# Author: Erin L Landguth
# Created: December 2010
# Description: This is the function/module file for dispersal processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# Python specific functions
import pdb, random, copy, os, sys

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# --------------------------------------------------------------------------
def logMsg(outf,msg):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print("%s"%(msg))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	# The case where all of the values in lst are the same
	if len(lst) == count:
		count = count-1
	return item,count
	
	#End::w_choice_general()	

# ---------------------------------------------------------------------------------------------------	 
def w_choice_item(lst):
	'''
	w_choice_item()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(lst)
	n=random.uniform(0,wtotal)
	for i in xrange(len(lst)):
		if n < lst[i]:
			break
		n = n-lst[i]
	return i
	
	#End::w_choice_item()
	
# ---------------------------------------------------------------------------------------------------	
def GetProbArray(Fxycdmatrix,Mxycdmatrix,tempoffspring,index,freegrid,\
philopatry,females,females_nomate,males,males_nomate):
	'''
	GetProbArray()
	This function gets indices for F and M specific cdmatrix values
	'''
	
	# Index into offspring array
	currentoff = tempoffspring[index]
	
	# If not philopatry
	if philopatry == 'N':	
		# Append the freegrid probabilities for the offspring choice
		if tempoffspring[index][4] == 0: # Female offspring
			probarray = Fxycdmatrix[currentoff[0]][freegrid]
		elif tempoffspring[index][4] == 1: # Male offspring
			probarray = Mxycdmatrix[currentoff[0]][freegrid]
		else:
			print('Invalid offspring list.')
			sys.exit(-1)		
	
	# If female philopatry
	elif philopatry == 'F' or philopatry == 'f' or philopatry == 'Female' or philopatry == 'female' or philopatry == '0':	
		'''
		Females are phylopatric and therefore stay at their birthplace
		(mother place), while a male can disperse to a dead male location.
		If by chance a female did not mate (no female offspring to take
		over), then male from nearest can disperse to fill.
		'''
		# Keep track of females that did not mate, first fill up those grid spots with nearest
		# Careful here, birth rate must be greater than 2.
		probarray = []
		# Loop through free grids appending sex prob to storage list
		for j in xrange(len(freegrid)):	
			# Append the freegrid probabilities for the offspring choice
			if currentoff[4] == 0:	# Offspring is a female			
				# Get the freegrid spots for just females
				find = False
				for fspot in xrange(len(females)):
					if freegrid[j] == females[fspot]:
						# Append probability female cdmatrix
						probarray.append(Fxycdmatrix[currentoff[0]][freegrid[j]])
						find = True
						break
				if find:
					continue
				else:
					probarray.append(0.0)
			elif currentoff[4] == 1: # Offspring is a male
				# Get the freegrid spots minus the females
				find = False
				for fspot in xrange(len(females)):
					if freegrid[j] == females[fspot]:
						probarray.append(0.0)
						find = True
						break
				if find:
					continue
				else:
					# Append the males can disperse with male cdmatrix
					probarray.append(Mxycdmatrix[currentoff[0]][freegrid[j]])
			else:
				print('Invalid offspring list.')
				sys.exit(-1)
				
		# If there are females that did not mate and female offspring probarray is zero
		if len(females_nomate) != 0 and currentoff[4] == 0 and sum(probarray) == 0.0:
			# Then change probability array so that this female is forced to disperse to nomate locations use male dispersal
			for nmspot in xrange(len(females_nomate)):
				fg_nomatespot = np.where(np.asarray(freegrid)==females_nomate[nmspot])[0]
				if len(fg_nomatespot) != 0:
					probarray[fg_nomatespot[0]] = Mxycdmatrix[currentoff[0]][females_nomate[nmspot]]
	
	# If male philopatry
	elif philopatry == 'M' or philopatry == 'm' or philopatry == 'male' or philopatry == 'Male' or philopatry == '1':	
		'''
		Males are phylopatric and therefore stay at their birthplace
		(mother place), while a female can disperse to a dead male location.
		If by chance a female did not mate (no male offspring to take
		over), then male from nearest can disperse to fill.
		'''
		# Keep track of females that did not mate, first fill up those grid spots with nearest male
		# Careful here, birth rate must be greater than 2.
		probarray = []
		
		# Loop through free grids appending sex prob to storage list
		for j in xrange(len(freegrid)):	
			# Append the freegrid probabilities for the offspring choice
			if currentoff[4] == 1: # Offspring is a male				
				# Get the freegrid spots for just females
				find = False
				for fspot in xrange(len(females)):
					if freegrid[j] == females[fspot]:
						# Append probability male cdmatrix
						probarray.append(Mxycdmatrix[currentoff[0]][freegrid[j]])
						find = True
						break
				if find:
					continue
				else:
					probarray.append(0.0)
			elif currentoff[4] == 0: # Offspring is a female				
				# Get the freegrid spots minus the females
				find = False
				for fspot in xrange(len(females)):
					if freegrid[j] == females[fspot]:
						probarray.append(0.0)
						find = True
						break
				if find:
					continue
				else:
					# Append probabilily of the female cdmatrix
					probarray.append(Fxycdmatrix[currentoff[0]][freegrid[j]])
			else:
				print('Invalid offspring list.')
				sys.exit(-1)
				
		# If there are females that did not mate and male offspring probarray is zero
		if len(females_nomate) != 0 and currentoff[4] == 1 and sum(probarray) == 0.0:
			# Then change probability array so that this male is forced to disperse to nomate locations use female dispersal
			for nmspot in xrange(len(females_nomate)):
				fg_nomatespot = np.where(np.asarray(freegrid)==females_nomate[nmspot])[0]
				if len(fg_nomatespot) != 0:
					probarray[fg_nomatespot[0]] = Fxycdmatrix[currentoff[0]][females_nomate[nmspot]]
	
	return probarray
	
	# End::GetProbArray()
	
# ---------------------------------------------------------------------------------------------------	
def Do1LocusSelection(offspring,fitvals1,tempfreegrid,iteminlist):
	'''
	Do1LocusSelection()
	This function calculates offsprings differential mortality, ie,
	offspring viability selection, for the 1-locus selection model.
	'''
	
	# If L0A0|L0A0 -- loci under selection:
	if offspring[6][0] == 2:

		# The grab it's fitness values
		differentialmortality = fitvals1[tempfreegrid[iteminlist]][0]/100
																
	# If L0A0|L0A1 -- loci under selection:
	elif offspring[6][0] == 1:

		# The grab it's fitness values
		differentialmortality = fitvals1[tempfreegrid[iteminlist]][1]/100
																															
	# If L0A1|L0A1 -- loci under selection
	else:
		
		# The grab it's fitness values
		differentialmortality = fitvals1[tempfreegrid[iteminlist]][2]/100
	
	return differentialmortality
	
	# End::Do1LocusSelection()
	
# ---------------------------------------------------------------------------------------------------	
def Do2LocusSelection(offspring,fitvals2,tempfreegrid,iteminlist):
	'''
	Do2LocusSelection()
	This function calculates offsprings differential mortality, ie,
	offspring viability selection, for the 2-locus selection model.
	'''
	# If L0A0|L0A0|L1A0|L1A0 - AABB -- loci under selection:
	if offspring[6][0] == 2 and offspring[6][2] == 2:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][0]/100
									
	# If L0A0|L0A1|L1A0|L1A0 - AaBB -- loci under selection:
	elif offspring[6][0] == 1 and offspring[6][2] == 2:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][1]/100
																															
	# If L0A1|L0A1|L1A0|L1A0 - aaBB -- loci under selection
	elif offspring[6][0] == 0 and offspring[6][2] == 2:
		
		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][2]/100
									
	# If L0A0|L0A0|L1A0|L1A1 - AABb -- loci under selection:
	elif offspring[6][0] == 2 and offspring[6][2] == 1:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][3]/100
									
	# If L0A0|L0A1|L1A0|L1A1 - AaBb -- loci under selection:
	elif offspring[6][0] == 1 and offspring[6][2] == 1:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][4]/100
																	
	# If L0A1|L0A1|L1A0|L1A1 - aaBb -- loci under selection
	elif offspring[6][0] == 0 and offspring[6][2] == 1:
		
		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][5]/100
	
	# If L0A0|L0A0|L1A1|L1A1 - AAbb -- loci under selection:
	elif offspring[6][0] == 2 and offspring[6][2] == 0:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][6]/100
									
	# If L0A0|L0A1|L1A1|L1A1 - Aabb -- loci under selection:
	elif offspring[6][0] == 1 and offspring[6][2] == 0:

		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][7]/100
																															
	# If L0A1|L0A1|L1A1|L1A1 - aabb -- loci under selection
	elif offspring[6][0] == 0 and offspring[6][2] == 0:
		
		# The grab it's fitness values
		differentialmortality = fitvals2[tempfreegrid[iteminlist]][8]/100
		
	return differentialmortality
	
	# End::Do2LocusSelection()

# ---------------------------------------------------------------------------------------------------	
def DoEmigration(offspring,freegrid,Migrants,Open,loci,alleles,\
Fxycdmatrix,Mxycdmatrix,gen,\
offspringno,cdevolveans,fitvals1,\
fitvals2,subpop,subpopmigration,DisperseDeaths,CouldNotDisperse,\
subpopmortperc,philopatry,females,subpopemigration,females_nomate,\
males,males_nomate,burningen):
	'''
	DoEmigration()
	This function enforces emigration when there are
	more offspring than open grid spots.
	'''	
	
	# Create variable to store offspring that disperse inside grid
	OffDisperseIN=[]
	
	# Do a shuffle on offpspring
	shuffle(offspring)
		
	# Initialize the while loop
	dispcount = 0
	offcount = 0
	
	# Add spot to track dispersing deaths for cdevolve
	DisperseDeaths.append([])
	CouldNotDisperse.append([])
		
	# Deep copy the freegrid spots to delete from
	tempfreegrid = copy.deepcopy(freegrid)
	
	# This while loop makes sure loop stops at carrying capacity (ie, total number of freegrids)
	#	or stops at end of offpsring list
	while dispcount < len(freegrid) and offcount < len(offspring):
		
		# Loop through offspring that are shuffled
		for i in xrange(len(offspring)):
			
			# Create a function here that gets indices for male and female
			probarray = GetProbArray(Fxycdmatrix,Mxycdmatrix,offspring,i,\
			tempfreegrid,philopatry,females,females_nomate,males,males_nomate)
									
			# If statement to check if there are spots for offpsring to disperse to
			if sum(probarray) != 0.0:
				
				# CDEVOLVE
				if cdevolveans == '1' and gen >= burningen:
											
					# Select the w_choice item
					iteminlist = w_choice_item(probarray)
					
					# Call 1-locus selection model
					differentialmortality = Do1LocusSelection(offspring[i],fitvals1,tempfreegrid,iteminlist)
											
					# Then flip the coin to see if offspring survives its location
					randcheck = rand()
					
					# If offspring did not survive: break from loop, move to next offspring
					if randcheck < differentialmortality:
						offcount = offcount + 1
						DisperseDeaths[gen].append(1)
						CouldNotDisperse[gen].append(0)
						continue
												
				# CDEVOLVE - 2 loci
				elif cdevolveans == '2' and gen >= burningen:
					
					# Select the w_choice item
					iteminlist = w_choice_item(probarray)

					# Call 2-locus selection model
					differentialmortality = Do2LocusSelection(offspring[i],fitvals2,tempfreegrid,iteminlist)
											
					# Then flip the coin to see if offspring survives its location
					randcheck = rand()
					
					# If offspring did not survive: break from loop, move to next offspring
					if randcheck < differentialmortality:
						offcount = offcount + 1
						DisperseDeaths[gen].append(1)
						CouldNotDisperse[gen].append(0)
						continue
						
				# If subpopulation differential mortality is on
				elif sum(subpopmortperc) != 0.0:
											
					# Select the w_choice item
					iteminlist = w_choice_item(probarray)
											
					# What subpopulation is offspring coming from
					fromsubpop = subpop[offspring[i][0]]
					# Where is subpopulation proposing to go
					tosubpop = subpop[tempfreegrid[iteminlist]]
					
					# If it is dispersing to another subpopulation
					if fromsubpop != tosubpop:
						
						# Grab its mortality percentage
						differentialmortality = subpopmortperc[int(tosubpop)-1]
						
						# Then flip the coin to see if offspring survives its location
						randcheck = rand()
						
						# If offspring did not survive: break from loop, move to next offspring
						if randcheck < differentialmortality:
							offcount = offcount + 1
							DisperseDeaths[gen].append(1)
							CouldNotDisperse[gen].append(0)
							continue
						
				# If not cdevolve or if cdevolve but it is in burn in gen
				else:
					
					# Select the w_choice item
					iteminlist = w_choice_item(probarray)				
				
				# Append information to variable [offspring, grid it dispersed to, and name]
				recd = [offspring[i],tempfreegrid[iteminlist],'T'+str(gen)+\
				'M'+str(offspring[i][0])+'F'+str(offspring[i][1])+\
				'Pop'+str(subpop[offspring[i][0]])]
							
				# Record offspring disperse information	
				OffDisperseIN.append(recd)
									
				# Update count for freegrid filling up
				dispcount = dispcount + 1
				offcount = offcount + 1
				DisperseDeaths[gen].append(0)
				CouldNotDisperse[gen].append(0)
				
				# Store the subpop dispersing to another subpop...what if 1, 3, 5 labeled subpop. either error intially or fix here....
				dispersingfrom = subpop[recd[0][0]]
				dispersingto = subpop[tempfreegrid[iteminlist]]
				if dispersingto != dispersingfrom:
					subpopmigration[gen][int(dispersingto)-1].append(1)
					subpopemigration[gen][int(dispersingfrom)-1].append(1)
					
				# And then delete freegrid spot from temp variable
				del(tempfreegrid[iteminlist])
										
			# If statement to check if there were not spots to disperse to
			elif sum(probarray) == 0.0:
				
				# Then Break from the loop and move to next offspring
				offcount = offcount + 1
				CouldNotDisperse[gen].append(1)
				continue
										
			
	# The store the number of disperses to separate subpops
	subpopmigration.append([]) # This adds a spot for next generation
	subpopemigration.append([])
	for i in xrange(len(subpopmigration[0])):
		subpopmigration[gen][i]=sum(subpopmigration[gen][i])
		subpopmigration[gen+1].append([0]) # This adds spots for subpops in next generation
		subpopemigration[gen][i]=sum(subpopemigration[gen][i])
		subpopemigration[gen+1].append([0]) # This adds spots for subpops in next generation
	
	# Store numbers
	DisperseDeaths[gen] = sum(DisperseDeaths[gen])
	CouldNotDisperse[gen] = sum(CouldNotDisperse[gen])
	Migrants.append(len(OffDisperseIN))
	Open.append(len(tempfreegrid))
	 
	# Variables returned
	tupDoEmi = OffDisperseIN,tempfreegrid
	return tupDoEmi
	
	# End::DoEmigration()
	
# ---------------------------------------------------------------------------------------------------	
def CalculateDispersalMetrics(OffDisperseIN,xgridcopy,ygridcopy,\
Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistED,\
MDispDistED,FDispDistEDstd,MDispDistEDstd,FDispDistCD,MDispDistCD,\
FDispDistCDstd,MDispDistCDstd,Fthreshold,Mthreshold):
	'''
	CalculateDispersalMetrics()
	This function calculates how far disperses are moving.
	'''		
	# Store the average dispersal distance offspring went
	# temp variable to store offspring dispersal distance
	FtempAvgDispDistED = []
	MtempAvgDispDistED = []
	FtempAvgDispDistCD = []
	MtempAvgDispDistCD = []
	OffDispDistCD = [] # Combined dispersal distance for storage
	Fcount = 0
	Mcount = 0
	
	# Loop through each OffDisperseIN
	for ioffspring in xrange(len(OffDisperseIN)):
		# Store the ED/CD distance offspring went - split up into sex
		if OffDisperseIN[ioffspring][0][4] == 0:
			FtempAvgDispDistED.append(np.sqrt((xgridcopy[OffDisperseIN[ioffspring][0][2]]-xgridcopy[OffDisperseIN[ioffspring][1]])**2+(ygridcopy[OffDisperseIN[ioffspring][0][2]]-ygridcopy[OffDisperseIN[ioffspring][1]])**2))
			Fcount = Fcount + 1
			
			# If panmictic
			if Fdispmoveno == '4' or Fdispmoveno == '6':
				cdval = 0.0
				FtempAvgDispDistCD.append(cdval)
			# If linear
			elif Fdispmoveno == '1':
				probval = Fxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				cdval = (probval - 1.) * (-Fthreshold)
				FtempAvgDispDistCD.append(cdval)
			# If inverse square
			elif Fdispmoveno == '2':
				probval = Fxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				if probval == 1.0:
					cdval = 0.0
				else:	
					cdval = np.sqrt(1. / (probval * (- 1./(pow(Fthreshold,2))) + (1./(pow(Fthreshold,2)))))
				FtempAvgDispDistCD.append(cdval)
			# These function not converted yet....####
			elif Fdispmoveno == '5' or Fdispmoveno == '7' or Fdispmoveno == '8':
				probval = Fxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				cdval = probval
				FtempAvgDispDistCD.append(cdval)
				
			# Then add to combined storage
			OffDispDistCD.append(cdval)
			
		# Else a male
		else:
			MtempAvgDispDistED.append(np.sqrt((xgridcopy[OffDisperseIN[ioffspring][0][2]]-xgridcopy[OffDisperseIN[ioffspring][1]])**2+(ygridcopy[OffDisperseIN[ioffspring][0][2]]-ygridcopy[OffDisperseIN[ioffspring][1]])**2))
			Mcount = Mcount + 1
			# If panmictic
			if Mdispmoveno == '4' or Mdispmoveno == '6':
				cdval = 0.0
				MtempAvgDispDistCD.append(cdval)
			# If linear
			elif Mdispmoveno == '1':
				probval = Mxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				cdval = (probval - 1.) * (-Mthreshold)
				MtempAvgDispDistCD.append(cdval)
			# If inverse square
			elif Mdispmoveno == '2':
				probval = Mxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				if probval == 1.0:
					cdval = 0.0
				else:	
					cdval = np.sqrt(1. / (probval * (- 1./(pow(Mthreshold,2))) + (1./(pow(Mthreshold,2)))))
				MtempAvgDispDistCD.append(cdval)
			# These function not converted yet....####
			elif Mdispmoveno == '5' or Mdispmoveno == '7' or Mdispmoveno == '8':
				probval = Mxycdmatrix[OffDisperseIN[ioffspring][0][2]][OffDisperseIN[ioffspring][1]]
				cdval = probval
				MtempAvgDispDistCD.append(cdval)
			
			# Then add to combined storage
			OffDispDistCD.append(cdval)
			
	# If at least 1 Female offspring dispersed
	if Fcount > 0:		
		# And append to DispDistED
		FDispDistED.append(sum(FtempAvgDispDistED)/Fcount)
		FDispDistEDstd.append(np.std(FtempAvgDispDistED))
		# And append to DispDistCD
		FDispDistCD.append(sum(FtempAvgDispDistCD)/Fcount)
		FDispDistCDstd.append(np.std(FtempAvgDispDistCD))
	else:
		# And append to DispDistED
		FDispDistED.append(0)
		FDispDistEDstd.append(0)
		# And append to DispDistCD
		FDispDistCD.append(0)
		FDispDistCDstd.append(0)
	
	# If at least 1 Male offspring dispersed
	if Mcount > 0:		
		# And append to DispDistED
		MDispDistED.append(sum(MtempAvgDispDistED)/Mcount)
		MDispDistEDstd.append(np.std(MtempAvgDispDistED))
		# And append to DispDistCD
		MDispDistCD.append(sum(MtempAvgDispDistCD)/Mcount)
		MDispDistCDstd.append(np.std(MtempAvgDispDistCD))
	else:
		# And append to DispDistED
		MDispDistED.append(0)
		MDispDistEDstd.append(0)
		# And append to DispDistCD
		MDispDistCD.append(0)
		MDispDistCDstd.append(0)

	return OffDispDistCD
	
	# End::CalculateDispersalMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def DoDisperse(offspringno,freegrid,offspring,Fdispmoveno,Mdispmoveno,\
Fxycdmatrix,Mxycdmatrix,gen,Migrants,Open,loci,alleles,\
xgridcopy,ygridcopy,FDispDistED,MDispDistED,FDispDistCD,MDispDistCD,
logfHndl,cdevolveans,fitvals1,fitvals2,FDispDistEDstd,MDispDistEDstd,\
FDispDistCDstd,MDispDistCDstd,subpop,subpopmigration,DisperseDeaths,CouldNotDisperse,\
subpopmortperc,philopatry,females,subpopemigration,females_nomate,males,males_nomate,\
burningen,Fthreshold,Mthreshold):
	'''
	DoDisperse()
	Disperse the new offspring to empty spots on grid
	Input: Units of dipsersal, movement function,
	offspring, freegrid, cdmatrix 
	Output: OffDisperseIN = [offspring,freegrid,name,[offspringgenes]] 
	'''		
		
	tupDoEmi = DoEmigration(offspring,freegrid,Migrants,Open,\
	loci,alleles,Fxycdmatrix,Mxycdmatrix,gen,offspringno,\
	cdevolveans,fitvals1,fitvals2,subpop,subpopmigration,\
	DisperseDeaths,CouldNotDisperse,subpopmortperc,philopatry,\
	females,subpopemigration,females_nomate,males,males_nomate,\
	burningen)
	
	OffDisperseIN = tupDoEmi[0]
	opengrids = tupDoEmi[1]
	
	# Calculate Dispersal Metrics
	OffDispDistCD = CalculateDispersalMetrics(OffDisperseIN,xgridcopy,ygridcopy,\
	Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistED,\
	MDispDistED,FDispDistEDstd,MDispDistEDstd,\
	FDispDistCD,MDispDistCD,FDispDistCDstd,MDispDistCDstd,Fthreshold,Mthreshold)
		
	# Return variables from this argument
	tupDoDisp = OffDisperseIN,opengrids,OffDispDistCD
	return tupDoDisp
	
	# End::DoDisperse()