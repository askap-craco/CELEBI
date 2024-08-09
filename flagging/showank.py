import os,sys
import numpy as np
from convertfits import *

ANTS		=	25
ANT0		=	1
SHOWTF		=	1

infilename	=	sys.argv[1]
infile2name	=	sys.argv[2]

infile2	=	fits.open(infile2name)		
data2	=	infile2[0].data	

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
	flaggingstatus.append(showbasecomparison(blid[i],data,data2,SHOWTF))

flaggingstatus	=	np.array(flaggingstatus)
avgflag			=	np.mean(flaggingstatus,axis=0)
print('\n\nAverage flagging fraction	%.3f %.3f	%.3f %.3f'%(avgflag[0],avgflag[1],avgflag[2],avgflag[3]))
print('\nFlagged data			%.3f %.3f\n'%((avgflag[2]-avgflag[0])/(1.0-avgflag[0]),(avgflag[3]-avgflag[1])/(1.0-avgflag[1])))

infile.close()
infile2.close()
	
	
































