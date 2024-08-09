import os,sys
import numpy as np
import time as tm
from convertfits import *
from inputs import *

#	----------------------------------------------

if(len(sys.argv)<6):
    print("Inputs are -- infile, outfile, flagstage, badchanfile, flagautocorr")
    sys.exit()

infilename		=	sys.argv[1]					#	Name of the input FITS file
outfilename		=	sys.argv[2]					#	Name of the output file
flagstage		=	int(sys.argv[3])			#	stage - 1/2/3/4
allbadchanfile	=	sys.argv[4]					#	list of bad channels
FLAGAUTOCORR	=	int(sys.argv[5])
ANKDIR          =   os.path.split(os.path.abspath(__file__))[0]
scratchfitsname	=	'scratchfits.fits'

print("runank running with -- "+infilename+" "+outfilename+" "+str(flagstage)+" "+allbadchanfile+" "+str(flagautocorr))

#	--------------		Convert inputs to numbers

if (CLEARSCRATCH):
    print('\nClearing scratch directory....\n')
    os.system('rm -rf '+scratchdir)
    os.system('mkdir '+scratchdir)


exmode		=	['baseline', 'uvbin']
flagwhat	=	['vis_ind', 'chan_ind', 'rec_ind', 'vis_block', 'chan_block', 'rec_block']
flagon		=	['mean', 'rms', 'mean_rms']
statused	=	['median', 'mean']
datype		=	['re', 'im', 'am', 'ph']
blkorder	=	['ascending', '', 'descending']

if(flagstage==1):
    FLAGMODE	=	'baseline'
    FLAGPARS	=	FLAGPARS1
elif(flagstage==2):
    FLAGMODE	=	'uvbin'
    FLAGPARS	=	FLAGPARS2
elif(flagstage==3):
    FLAGMODE	=	'baseline'
    FLAGPARS	=	FLAGPARS3
elif(flagstage==4):
    FLAGMODE	=	'uvbin'
    FLAGPARS	=	FLAGPARS4
else:
    print ('Unknown flagstage !!')

print('Total flagging rounds	=	%d'%len(FLAGPARS))

flagparfile	=	open(scratchdir+'flagpars.pars','w')
flagparfile.write('%d	%d	%d	%d	%d	'%(ANTS, exmode.index(FLAGMODE), ugrids, flagon.index(SCANFLAGMEAN[0])+1, vgrids))
flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(SCANFLAGMEAN[1], SCANFLAGMEAN[2], len(FLAGPARS), THREADS, WRITEOUT, BLOCKPOW, SCANFLAGMEAN[3]))

for flpar in FLAGPARS:
    flagparfile.write('%d	%d	%d	%d	%d	'%(flagwhat.index(flpar[0]), flagon.index(flpar[1])+1, statused.index(flpar[2]), datype.index(flpar[3]), -(blkorder.index(flpar[7])-1)))
    flagparfile.write('%f	%f	%d	%d	%d	%f	%f\n'%(flpar[4], flpar[5], flpar[6]+1, flpar[8], flpar[9], flpar[10], flpar[11]))

flagparfile.close()

#	----------------------------	Convert FITS to binary files
start0	=	tm.time()	

if(FLAGAUTOCORR):
    infile		=	fits.open(infilename)
    data		=	infile[0].data	
    nauto 		= 	flagautocorr(ANTS,data)
    infile.writeto(scratchdir+scratchfitsname,output_verify='warn',overwrite=True)
	
if (CONVERTFITS):
    if(FLAGAUTOCORR):
        infile		=	fits.open(scratchdir+scratchfitsname)
    else:
        infile		=	fits.open(infilename)
    data		=	infile[0].data
		
    #	-----------	Flagging hopeless channels	---------------------------
    if(os.path.exists(allbadchanfile)):
        hopelesschans	=	[]
        allbadchans		=	np.loadtxt(allbadchanfile).astype('int')

        if(len(allbadchans.shape)==2):
            if(allbadchans.shape[1]==2):
                print("Bad channel list")
                for lk in range(0,allbadchans.shape[0]):
                    print("%d	%d"%(allbadchans[lk,0],allbadchans[lk,1]+1))
                    for bc in range(allbadchans[lk,0],allbadchans[lk,1]+1):
                        hopelesschans.append(bc)
                hopelesschans	=	np.array(hopelesschans)	
                nhopelesschans	=	flagchanlist(ANTS,data,hopelesschans)
                print("Flagged %d bad channels everywhere..."%nhopelesschans)
            else:
                print("Useless bad channel list... ignoring")
        else:
            print("Useless bad channel list... ignoring")
    else:
        print("No bad channel list found... continuing")
	
    #	-------------------------------------------------------------------
	
    if (FLAGMODE==exmode[1]):		
        uvfitstobinary(data,scratchdir,ugrids,vgrids,plotuv)

    elif (FLAGMODE==exmode[0]):		
        baselinestobinary(ANTS,data,scratchdir,BREAKTIMESEC)
			
    else:
        print("Unknown flagging mode !!!!			Please tell me how to execute it ........")
	
    infile.close()	
    print("\nConvertion done in 		%d seconds\n"%(tm.time()-start0))
	
#	------------------------------		Flag data
start1	=	tm.time()

if (DOFLAG):
    print("runank - running core flagging")
    status	=	os.system('module load gcc/12.2.0 && module load gsl/2.7 && '+ANKDIR+'/ankflag')	
    print("\nFlagging done in 		%d seconds\n"%(tm.time()-start1))
	
#	------------------------------		Convert back binary files to FITS	

if (READBACK):
    print("runank - reading back flagged data")
    if(FLAGAUTOCORR):
        infile2		=	fits.open(scratchdir+scratchfitsname)
    else:
        infile2		=	fits.open(infilename)		
    data2		=	infile2[0].data	

    if (FLAGMODE==exmode[1]):
        bintofits	=	uvfitsfrombinary(data2,scratchdir,ugrids,vgrids)

    elif (FLAGMODE==exmode[0]):
        baselinesfrombinary(ANTS,data2,scratchdir)

    #	---------------------------------------------------------------------
    #					Plot Baselines 
    #	---------------------------------------------------------------------

    if (SHOWBASE):

        infile	=	fits.open(infilename)
        data	=	infile[0].data

        blid	=	[]
        for a in range (1,ANTS):
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
		
        print('\n\nAverage flagging fraction	'),
        for p in range (0,npols):
            print('%.3f %.3f	'%(avgflag[p],avgflag[p+npols])),
        print('\n')
		
        print('\nFlagged data			'),
        for p in range (0,npols):
            print('%.3f	'%((avgflag[p+npols]-avgflag[p]))),
        print('\n')

        infile.close()

    if (WRITEOUT):
        infile2.writeto(outfilename,output_verify='warn',overwrite=True)

    infile2.close()
    print("\nEverything done in 		%d seconds\n"%(tm.time()-start0))
	#	----------------------------------------------------------------------------

if (CLEARSCRATCH):
    print('\nClearing scratch directory....\n')
    os.system('rm -rf '+scratchdir)






























