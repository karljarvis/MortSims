# -------------------------------------------------------------------------------------------------
# CDPOP_Mate.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for mate processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Python specific functions
import pdb, random, os, sys, copy
from sets import Set

# --------------------------------------------------------------------------
def countDuplicatesInList(dupedList):
	'''
	countDuplicatesInList() - Counts dupicates in lists
	'''
	uniqueSet = Set(item for item in dupedList)
	return [dupedList.count(item) for item in uniqueSet]
	
	# End::countDuplicatesInList()

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()
	
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
def DoSexualNY(nomales,xycdmatrix,females,count,\
males,matemovethresh,Bearpairs,subpopmatemort,subpopmortperc,subpop,gen):
	'''
	DoSexualNY()
	This function is the mating function for:
	sexual reproduction
	females	without replacement
	males with replacement
	'''	
	# Create an empty probability array to be appended to
	probarray = []
						
	# Extract each male grid probability
	probarray.append(xycdmatrix[females[count]][males])
				
	# If statement to check if there were individuals in probarray:
	if sum(probarray[0]) != 0.0:
	
		# If subpopulation differential mortality is on
		if sum(subpopmortperc) != 0.0:
		
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
			
			# What subpopulation are each mate pair
			fromsubpop = subpop[females[count]] # female
			tosubpop = subpop[males[itemselect]] # male
			
			# If there is a mate selection to another subpopulation
			if fromsubpop != tosubpop:
			
				# Grab its mortality percentage
				differentialmortality = subpopmortperc[int(tosubpop)-1]
				
				# Then flip the coin to see if mating occurs
				randcheck = rand()
				
				# If mate did not occur: 
				if randcheck < differentialmortality:
					# Track this 'mortality'
					subpopmatemort[gen].append(1)
					# Store the female that did not mate with -9999 designation
					Bearpairs.append([females[count],-9999])
					# Then delete that male from the male list
					males = np.delete(males,itemselect)
					
				# If mate did occur, record
				else:
					Bearpairs.append([females[count],males[itemselect]])
					
			# If they are in the same population: mate occured
			else:
				Bearpairs.append([females[count],males[itemselect]])

		# If no subpopmortperc
		else:		
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
	
			# And store the mated pair information.						
			Bearpairs.append([females[count],males[itemselect]])
							
	# If statement to check if there were not males < matemovethresh:
	elif sum(probarray[0]) == 0.0:
		
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[count],-9999])	
	
	return Bearpairs,males
	
	# End::DoSexualNY()
	
# ---------------------------------------------------------------------------------------------------	
def DoSexualYY(nomales,xycdmatrix,females,males,\
matemovethresh,Bearpairs,nofemales,subpopmatemort,subpopmortperc,subpop,gen):
	'''
	DoSexualYY()
	This function is the mating function for: 
	sexual reproduction
	females	with replacement
	males with replacement.
	'''	
	# Create an empty probability array to be appended to
	probarray = []					

	# Randomly grab a female
	intfemale = int((nofemales)*rand())
	
	# Extract each male grid probability
	probarray.append(xycdmatrix[females[intfemale]][males])
		
	# If statement to check if there were individuals in probarray:
	if sum(probarray[0]) != 0.0:

		# If subpopulation differential mortality is on
		if sum(subpopmortperc) != 0.0:
		
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
			
			# What subpopulation are each mate pair
			fromsubpop = subpop[females[count]] # female
			tosubpop = subpop[males[itemselect]] # male
			
			# If there is a mate selection to another subpopulation
			if fromsubpop != tosubpop:
			
				# Grab its mortality percentage
				differentialmortality = subpopmortperc[int(tosubpop)-1]
				
				# Then flip the coin to see if mating occurs
				randcheck = rand()
				
				# If mate did not occur: 
				if randcheck < differentialmortality:
					# Track this 'mortality'
					subpopmatemort[gen].append(1)
					# Store the female that did not mate with -9999 designation
					Bearpairs.append([females[count],-9999])
					# Then delete that male from the male list
					males = np.delete(males,itemselect)
					
				# If mate did occur, record
				else:
					Bearpairs.append([females[count],males[itemselect]])
					
			# If they are in the same population: mate occured
			else:
				Bearpairs.append([females[count],males[itemselect]])
		
		# If no subpopmortperc
		else:
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
			
			# And store the mated pair information.						
			Bearpairs.append([females[intfemale],males[itemselect]])
						
	# If statement to check if there were not males < matemovethresh:
	elif sum(probarray[0]) == 0.0:
	
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[intfemale],-9999])				
	
	# Return Variables from this function
	return Bearpairs,males
	
	# End::DoSexualYY()

# ---------------------------------------------------------------------------------------------------	
def DoSexualNN(nomales,xycdmatrix,females,count,\
males,matemovethresh,Bearpairs,subpopmatemort,subpopmortperc,subpop,gen):
	'''
	DoSexualNN()
	This function is the mating function for
	sexual reproduction
	females	with replacement
	males with replacement
	'''	
	# Create an empty probability array to be appended to
	probarray = []
					
	# Extract each male grid probability
	probarray.append(xycdmatrix[females[count]][males])
				
	# If statement to check if there were individuals in probarray:
	if sum(probarray[0]) != 0.0:

		# If subpopulation differential mortality is on
		if sum(subpopmortperc) != 0.0:
		
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
			
			# What subpopulation are each mate pair
			fromsubpop = subpop[females[count]] # female
			tosubpop = subpop[males[itemselect]] # male
			
			# If there is a mate selection to another subpopulation
			if fromsubpop != tosubpop:
			
				# Grab its mortality percentage
				differentialmortality = subpopmortperc[int(tosubpop)-1]
				
				# Then flip the coin to see if mating occurs
				randcheck = rand()
				
				# If mate did not occur: 
				if randcheck < differentialmortality:
					# Track this 'mortality'
					subpopmatemort[gen].append(1)
					# Store the female that did not mate with -9999 designation
					Bearpairs.append([females[count],-9999])
					# Then delete that male from the male list
					males = np.delete(males,itemselect)
					
				# If mate did occur, record
				else:
					Bearpairs.append([females[count],males[itemselect]])
					# Then delete that male from the male list
					males = np.delete(males,itemselect)
					
			# If they are in the same population: mate occured
			else:
				Bearpairs.append([females[count],males[itemselect]])
				# Then delete that male from the male list
				males = np.delete(males,itemselect)
		
		# If no subpopmortperc
		else:
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
		
			# And store the mated pair information.						
			Bearpairs.append([females[count],males[itemselect]])
			
			# Then delete that male from the male list
			males = np.delete(males,itemselect)
							
	# If statement to check if there were not males < matemovethresh:
	elif sum(probarray[0]) == 0.0:
	
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[count],-9999])	
	
	return Bearpairs,males
	
	# End::DoSexualNN()	

# ---------------------------------------------------------------------------------------------------	
def DoAsexualNY(xycdmatrix,matemovethresh,Bearpairs,selfans,count,breedgrids):
	'''
	DoAsexualNY()
	This function is the mating function for: 
	asexual reproduction
	females	without replacement
	males with replacement.
	'''		
	# Create an empty probability array to be appended to
	probarray = []
	
	# First mate
	firstmate = breedgrids[count]
	
	# Make array of individuals, removing itself unless selfing on
	if selfans == 'Y':
		indspots = np.asarray(copy.deepcopy(breedgrids))
	else:
		indspots = np.asarray(copy.deepcopy(breedgrids))
		delspot = np.where(indspots==firstmate)[0]
		indspots = np.delete(indspots,delspot)
	shuffle(indspots)
	
	# Extract each male grid probability
	probarray.append(xycdmatrix[firstmate][indspots])
					
	# If statement to check if there were individuals in probarray:
	if sum(probarray[0]) != 0.0:

		# Select the w_choice item
		itemselect = w_choice_item(probarray[0])
	
		# And store the mated pair information.						
		Bearpairs.append([firstmate,indspots[itemselect]])
							
	# If statement to check if there were not males < matemovethresh:
	elif sum(probarray[0]) == 0.0:
		
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([firstmate,-9999])
					
	# Return Variables from this function
	return Bearpairs
	
	# End::DoAsexualNY()	
	
# ---------------------------------------------------------------------------------------------------	
def DoAsexualYY(xycdmatrix,matemovethresh,Bearpairs,\
selfans,breedgrids):
	'''
	DoAsexualYY()
	This function is the mating function for: 
	asexual reproduction
	females	with replacement
	males with replacement.
	'''		
	# Create an empty probability array to be appended to
	probarray = []
	
	# Random number for grab of first grid - first mate
	intgrid = int((len(breedgrids))*rand())
	firstmate = breedgrids[intgrid]
	# Select second mate, unless selfing on
	if selfans == 'Y':
		indspots = np.asarray(copy.deepcopy(breedgrids))
	else:
		indspots = np.asarray(copy.deepcopy(breedgrids))
		delspot = np.where(indspots==firstmate)[0]
		indspots = np.delete(indspots,delspot)
	shuffle(indspots)
	# Extract each male grid probability
	probarray.append(xycdmatrix[firstmate][indspots])
					
	# If statement to check if there were individuals in probarray:
	if sum(probarray[0]) != 0.0:

		# Select the w_choice item
		itemselect = w_choice_item(probarray[0])
	
		# And store the mated pair information.						
		Bearpairs.append([firstmate,indspots[itemselect]])
							
	# If statement to check if there were not males < matemovethresh:
	elif sum(probarray[0]) == 0.0:
	
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([firstmate,-9999])
					
	# Return Variables from this function
	return Bearpairs
	
	# End::DoAsexualYY()	
		
# ---------------------------------------------------------------------------------------------------	 
def DoMate(nogrids,sex,lmbda,age,\
freplace,mreplace,matemoveno,matemovethresh,\
xycdmatrix,MateDistED,\
MateDistCD,xgridcopy,ygridcopy,ToTMales,\
ToTFemales,BreedMales,BreedFemales,\
sexans,selfans,FID,matemoveparA,matemoveparB,\
MateDistEDstd,MateDistCDstd,FAvgMate,MAvgMate,\
FSDMate,MSDMate,filledgrids,Female_BreedEvents,\
gen,subpop,BreedFemales_age,agemort,subpopmatemort,subpopmortperc):

	'''
	DoMate()
	This is the mating function for choosing
	individual mate pairs. Two calls here:
	sexual and asexual mating.	
	'''
	
	# --------------------------------------------------------
	# Preliminary: Needed for both sexual and asexual routines	
	# --------------------------------------------------------
	
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	
	# Get reproduction age from lmbda vals
	for i in xrange(len(lmbda)):
		if lmbda[i] != 0.0:
			reproage = i+1 # Add one to age class lmbda 1, 2, ...
			break
		
	# ---------------------------------------
	# Step10a: Call DoMateSexual()
	# ---------------------------------------
	if (sexans == 'Y'):
		'''
		DoMateSexual()
		This function is the mating function for 
		sexual reproduction.
		'''
		# ---------------------------------------------------
		# Select males and females for mating
		# ---------------------------------------------------

		# Storage variables for breeding age males and females
		females = []		# These are the zeros
		males = []			# These are the ones

		# Storage variables for total numbers of males and females
		allfemales = []
		allmales = []
		
		# Loop through and grab each index for each catagory for reproage and sex
		for i in xrange(nogrids):
			if sex[i] != 'NA':
				if sex[i] == '0': 
					allfemales.append(i)
					if str(age[i]) >= str(reproage):
						females.append(i)				
				elif sex[i] == '1':
					allmales.append(i)
					if str(age[i]) >= str(reproage):
						males.append(i)
		
		# Then get the length of each sex that are reproducing
		nomales = len(males)
		nofemales = len(females)
		
		# And then sum them up - and store numbers
		subpopmatemort.append([]) # For generation
		ToTMales.append([])
		ToTFemales.append([])
		BreedMales.append([])
		BreedFemales.append([])
		BreedFemales_age.append([])
		ToTMales[gen].append(len(allmales))
		ToTFemales[gen].append(len(allfemales))
		BreedMales[gen].append(nomales)
		BreedFemales[gen].append(nofemales)
		for i in xrange(nosubpops):
			# Get all subpop sex
			tempallpopsex = np.asarray(sex)[np.where(np.asarray(subpop) == str(i+1))[0]]
			ToTMales[gen].append(len(np.where(tempallpopsex == '1')[0]))
			ToTFemales[gen].append(len(np.where(tempallpopsex == '0')[0]))
			# Get all breeding subpop sex
			tempallpopage = np.asarray(age)[np.where(np.asarray(subpop) == str(i+1))[0]]
			tempmalepopage = np.asarray(tempallpopage[np.where(tempallpopsex == '1')[0]],dtype='int')
			tempfemalepopage = np.asarray(tempallpopage[np.where(tempallpopsex == '0')[0]],dtype='int')			
			BreedMales[gen].append(len(np.where(tempmalepopage>=reproage)[0]))
			BreedFemales[gen].append(len(np.where(tempfemalepopage>=reproage)[0]))
		# Get count of female age breeders
		Fsex = np.where(np.asarray(sex) == '0')[0]
		Fsexage = np.asarray(age)[Fsex]
		Fsexage = np.asarray(Fsexage,dtype = int)
		Fsexreproage = Fsexage[np.where(Fsexage >= int(reproage))[0]]
		for i in xrange(len(agemort)):
			BreedFemales_age[gen].append([])
		countFage = count_unique(Fsexreproage)
		
		count = 0
		for item in countFage[0]:
			if item != 'NA':
				if int(item) > len(agemort):
					count = len(agemort)-1				
					BreedFemales_age[gen][count].append(countFage[1][count])
				else:
					BreedFemales_age[gen][count].append(countFage[1][count])
					count = count+1
		# Sum up population age tracker
		for i in xrange(len(agemort)):
			BreedFemales_age[gen][i] = sum(BreedFemales_age[gen][i])
		
		# Choose mate for each female or individual
		Bearpairs = []	# Empty matrix: xy indexes
		
		# Shuffle the females
		shuffle(females)
		
		# If there were no reproducing males or females
		if nomales == 0 or nofemales == 0:
			Bearpairs.append([-9999,-9999])
		
		# If there were reproducing males and females
		if nomales != 0 and nofemales != 0:
		
			# For the case of a Female without replacement and a male with replacement
			if freplace == 'N' and mreplace == 'Y':
							
				# Loop through while loop until all females paired up.		
				count = 0		# Initialize the while loop
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)
				while count < nofemales:
								
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNY(nomales,xycdmatrix,females,count,\
					tempmales,matemovethresh,Bearpairs,subpopmatemort,subpopmortperc,subpop,gen)
						
					# Update count
					count = count + 1
					
			# For the case of a Female with replacement and a male with replacement
			elif freplace == 'Y' and mreplace == 'Y':
					
				# Loop through while loop until all females paired up, but do this nogrid times.		
				count = 0		# Initialize the while loop
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)
				while count < filledgrids:
							
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualYY(nomales,xycdmatrix,females,\
					tempmales,matemovethresh,Bearpairs,nofemales,subpopmatemort,subpopmortperc,subpop,gen)
														
					# Update count
					count = count + 1
									
			# For the case of Female with replacement and male without replacement
			elif freplace == 'Y' and mreplace == 'N':
			
				print('Female with replacement and Male without replacement not coded yet.')
				print('Email Erin.')
				sys.exit(-1)
				
			# For the case of Female without replacement and male without replacement
			elif freplace == 'N' and mreplace == 'N':
				
				# Loop through while loop until all male female pairs occur		
				count = 0		# Initialize the while loop for females
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)				
				while count < nofemales:
								
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNN(nomales,xycdmatrix,females,count,tempmales,matemovethresh,Bearpairs,subpopmatemort,subpopmortperc,subpop,gen)
											
					# Update count
					count = count + 1
			
			# Error check
			else:
				print('This Female/Male mating structure does not exist. Must be Y/N combinations.')
				sys.exit(-1)		
		
		# Sort the females for later functions
		females = list(np.sort(females))
		
	# ---------------------------------------
	# Step10b: Call DoMateAsexual()
	# ---------------------------------------
	if (sexans=='N'):	
		'''
		DoMateAsexual()
		This function is the mating function for 
		asexual reproduction.
		'''
		
		# If subpopulation differential mortality is on
		if sum(subpopmortperc) != 0.0:
			print('Subpopulation mortlaity specified and not currently implemented with sexans==N.')
			sys.exit(-1)
		
		# Storage variable for total number of grids
		allgrids = []
		breedgrids = []
		
		# Loop through and grab each index for each catagory for reproage
		for i in xrange(nogrids):
			if sex[i] == '0' or sex[i] == '1':
				allgrids.append(i)
				if age[i] >= reproage:
					breedgrids.append(i)				
		nobreedgrids = len(breedgrids)
		
		# And then sum them up - and store numbers
		ToTMales.append([])
		ToTFemales.append([])
		BreedMales.append([])
		BreedFemales.append([])
		ToTMales[gen].append(len(allgrids))
		ToTFemales[gen].append(len(allgrids))
		BreedMales[gen].append(nobreedgrids)
		BreedFemales[gen].append(nobreedgrids)
		for i in xrange(nosubpops):
			# Get all subpop sex
			tempallpopsex = np.asarray(sex)[np.where(np.asarray(subpop) == str(i+1))[0]]
			ToTMales[gen].append(len(tempallpopsex))
			ToTFemales[gen].append(len(tempallpopsex))
			# Get all breeding subpop sex
			tempallpopage = np.asarray(age)[np.where(np.asarray(subpop) == str(i+1))[0]]
			BreedMales[gen].append(len(np.where(tempallpopage>=reproage)[0]))
			BreedFemales[gen].append(len(np.where(tempallpopage>=reproage)[0]))
				
		# Choose mate for each female or individual
		Bearpairs = []	# Empty matrix: xy indexes
		
		# Shuffle the grids
		shuffle(breedgrids)
		
		# If there were no reproducing grids
		if nobreedgrids == 0:
			Bearpairs.append([-9999,-9999])	

		# If there were reproducing grids
		else:
					
			# For the case of a first mate without replacement and second mate with replacement
			if freplace == 'N' and mreplace == 'Y':
			
				# Loop through while loop until all grids are paired up according to function		
				count = 0		# Initialize the while loop
				while count < filledgrids:
													
					# Get probability function of user defined input number
					Bearpairs = DoAsexualNY(xycdmatrix,matemovethresh,Bearpairs,\
					selfans,count,breedgrids)		
															
					# Update count
					count = count + 1
					
			# For the case of a first mate with replacement and second mate with replacement
			elif freplace == 'Y' and mreplace == 'Y':
			
				# Loop through filledgrids times.		
				count = 0		# Initialize the while loop
				while count < filledgrids:
				
					# Get probability function of user defined input number
					Bearpairs = DoAsexualYY(xycdmatrix,matemovethresh,Bearpairs,\
					selfans,breedgrids)											
					
					# Update count
					count = count + 1					
				
			# For the case of first mate without replacement and second mate without replacement
			elif freplace == 'N' and mreplace == 'N':
			
				print('This mate function is not coded up yet, email Erin.')
				sys.exit(-1)
				
			# For the case of the first mate with replacement and second mate without
			elif freplace == 'Y' and mreplace == 'N':				
				print('This is not coded up yet, but is the same case as NY replacement.')
				sys.exit(-1)
				
			# Error check
			else:
				print('This asexual mating structure does not exist. Must be Y/N combinations.')
				sys.exit(-1)
		
		# For passing on only! Careful here....
		females = breedgrids
		males = breedgrids
		
	# ----------------------------------------
	# Summary Stats on Mate functions
	# ----------------------------------------
	# Store the average distance mates were choosen 
	# temp variable to store the number of -9999 numbers
	tempImmStorage = []
	# temp variable to store mate distance
	tempAvgMateED = []
	tempAvgMateCD = []
	
	# Loop through each CDpair
	for ipair in xrange(len(Bearpairs)):
		
		# if -9999 store
		if Bearpairs[ipair][1]==-9999:
			tempImmStorage.append(1)
		
		# else calculate average mate distance
		else:
			tempAvgMateED.append(np.sqrt((xgridcopy[Bearpairs[ipair][0]]-xgridcopy[Bearpairs[ipair][1]])**2+(ygridcopy[Bearpairs[ipair][0]]-ygridcopy[Bearpairs[ipair][1]])**2))
			if matemoveno == '4':
				tempAvgMateCD.append(0)
			else:
				tempAvgMateCD.append(xycdmatrix[Bearpairs[ipair][0]][Bearpairs[ipair][1]])
			
	# If at least some individuals mated
	if len(Bearpairs) > sum(tempImmStorage):
		
		# And append to MateDistED
		MateDistED.append(sum(tempAvgMateED)/(len(Bearpairs)-sum(tempImmStorage)))
		MateDistEDstd.append(np.std(tempAvgMateED))

		# And append to MateDistCD
		MateDistCD.append(sum(tempAvgMateCD)/(len(Bearpairs)-sum(tempImmStorage)))
		MateDistCDstd.append(np.std(tempAvgMateCD))
		
	# If all immigrants
	else:
	
		# And append to MateDistED
		MateDistED.append(0)
		MateDistEDstd.append(0)

		# And append to MateDistCD
		MateDistCD.append(0)
		MateDistCDstd.append(0) 
	
	# Track actual number of breeding events of females.
	Female_BreedEvents.append(len(Bearpairs)-sum(tempImmStorage))
		
	# Delete temp storage
	del(tempImmStorage)
	del(tempAvgMateED)
	del(tempAvgMateCD)
	
	# Get the avgerage and sdev for female and male mating times
	#	For variance in reproductive succes.
	Bearpairs_array = np.asarray(Bearpairs)
	femalesmated = Bearpairs_array[:,0][np.where(Bearpairs_array[:,1]>=0)]
	femalesmated = list(femalesmated)
	malesmated = Bearpairs_array[:,1][np.where(Bearpairs_array[:,1]>=0)]
	malesmated = list(malesmated)
	females_nomate = list(Bearpairs_array[:,0][np.where(Bearpairs_array[:,1]==-9999)])
	males_nomate = list(Bearpairs_array[:,1][np.where(Bearpairs_array[:,0]==-9999)])
	# Append to tracker
	FAvgMate.append(np.mean(countDuplicatesInList(femalesmated)))
	MAvgMate.append(np.mean(countDuplicatesInList(malesmated)))
	FSDMate.append(np.std(countDuplicatesInList(femalesmated)))
	MSDMate.append(np.std(countDuplicatesInList(malesmated)))
	# No mate due to subpopmortperc
	subpopmatemort[gen] = sum(subpopmatemort[gen])
		
	# Return variables from this function
	tupMate = Bearpairs, females, females_nomate, males, males_nomate
	return tupMate
	
	#End::DoMate()