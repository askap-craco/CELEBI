import os,sys
import numpy as np
sys.path.insert(0,'/fred/oz002/askap/craft/craco/ankflag_craft')
from convertfits import *

ANTS		=	36
ANT0		=	1
SHOWTF		=	1

infilename	=	sys.argv[1]
badantfilename	=	sys.argv[2]

infile	=	fits.open(infilename)
data	=	infile[0].data

blid	=	[]
for a in range (ANT0,ANTS):
	for b in range (a+1,ANTS+1):
		blid.append([a,b,256*a+b])
blid	=	np.array(blid)
nbase	=	len(blid)

print('\nIdeally total baselines =	%d'%nbase)

flaggingstatus	=	[]
for i in range (0,nbase):	
	flaggingstatus.append(returnbadbaselines(blid[i],data))

flaggingstatus	=	np.array(flaggingstatus)

listbadants		=	[]

for a in range (ANT0,ANTS):
	antinds	=	np.concatenate((np.where(blid[:,0]==a)[0],np.where(blid[:,1]==a)[0]), axis=0)
	#print (a,blid[antinds])
	antstatus	=	flaggingstatus[antinds]
	#print (a, antstatus)
	if(len(antstatus[antstatus==0])==0):
		listbadants.append(a)
		#print(a)

listbadants		=	np.array(listbadants)
print("Bad antennas")
print(listbadants)

np.savetxt(badantfilename,listbadants.T,fmt='%d')

infile.close()

	
	
































