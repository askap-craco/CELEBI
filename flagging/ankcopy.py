import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import multiprocessing as mp
import time as tm

ANTS		=	30							#	Total number of antennas
parallel	=	False						#	Process baselines parallelly		#	SET TO "FALSE"
recflag		=	False						#	Copy flags for whole records
chanflag	=	False						#	Copy flags for whole channels	
visflag		=	True						#	Copy flags for individual visibilities
writeout	=	True						#	Write the output file

srcname		=	'/data/apurba/ptf10hgi/gmrt/band_3_uvsub_out.fits'			#	Input FITS file (FG TABLE IS NOT READ!)
targetname	=	'/data/apurba/ptf10hgi/gmrt/band_3_cal.fits'				#	Target FITS file
outname		=	'/data/apurba/ptf10hgi/gmrt/band_3_cal_out.fits'	#	Output FITS file


#	--------------------------------------------------------------------	Functions	----------------------------------------------------
#	Function to copy flags from one dataset to another for a baseline
def copybase(bl,datasrc,datatar):
	print('Baseline	%d-%d'%(bl[0],bl[1]))	
	bindxsrc	=	np.where(datasrc.par('BASELINE')==bl[2])
	bdatasrc	=	datasrc[bindxsrc[0]]
	bindxtar	=	np.where(datatar.par('BASELINE')==bl[2])
	bdatatar	=	datatar[bindxtar[0]]
	
	if(len(bdatasrc)==0):
		print('No data in source file .... flagging entire baseline in target file')
		if(len(bdatatar)==0):
			print('No data in target file either.... nothing to do')
			return [2.0,[[1.0,1.0,1.0],[1.0,1.0,1.0]],0]
		for p in range (0,2):		
			tfdatatar	=	bldatatar[:,0,0,0,:,p,:]		#	Time frequency data target
			tfretar		=	tfdatatar[:,:,0]				#	Real part of target
			tfimtar		=	tfdatatar[:,:,1]				#	Imaginary Part of target
			tfwtar		=	tfdatatar[:,:,2]				#	Weight of target
			tfwtar		=	0.0
			tfwre		=	0.0
			tfwim		=	0.0
		return [1.0,[[1.0,1.0,1.0],[1.0,1.0,1.0]],0]
	bldatasrc	=	bdatasrc.data
	bldatatar	=	bdatatar.data
	flgfrac		=	[]
	for p in range (0,2):		
		tfdatasrc	=	bldatasrc[:,0,0,0,:,p,:]		#	Time frequency data source
		tfdatatar	=	bldatatar[:,0,0,0,:,p,:]		#	Time frequency data target
		tfresrc		=	tfdatasrc[:,:,0]				#	Real part of source
		tfretar		=	tfdatatar[:,:,0]				#	Real part of target
		tfimsrc		=	tfdatasrc[:,:,1]				#	Imaginary Part of source
		tfimtar		=	tfdatatar[:,:,1]				#	Imaginary Part of target
		tfwsrc		=	tfdatasrc[:,:,2]				#	Weight of source
		tfwtar		=	tfdatatar[:,:,2]				#	Weight of target
		
		if(np.size(tfwsrc)!=np.size(tfwtar)):
			print('Caution !	Sizes do not match	src = %d	target = %d'%(np.size(tfwsrc),np.size(tfwtar)))
		
		srcfrac		=	1.0-float(np.count_nonzero(tfwsrc>0.0))/np.size(tfwsrc)
		tarfrac		=	1.0-float(np.count_nonzero(tfwtar>0.0))/np.size(tfwtar)
		print('Pol	%d	flagged fractions	src 	= %.2f	target = %.2f'%(p,srcfrac,tarfrac))
		
		if(recflag):
			#print('Copyting record flags ...........')
			nrecs	=	len(tfwsrc)
			nrect	=	len(tfwtar)
			if(nrecs!=nrect):
				print('Error !	Records do not natch !!!!')
				return 3
			for r in range(0,nrecs):
				if(np.count_nonzero(tfwsrc[r])==0):
					tfwtar[r]	=	0.0
					tfretar[r]	=	0.0
					tfimtar[r]	=	0.0		
		
		if(chanflag):
			#print('Copyting channel flags ...........')
			nchans	=	len(tfwsrc[0])
			nchant	=	len(tfwtar[0])
			if(nchans!=nchant):
				print('Error !	Channels do not natch !!!!')
				return 3
			for c in range(0,nchans):
				if(np.count_nonzero(tfwsrc[:,c])==0):
					tfwtar[:,c]		=	0.0
					tfretar[:,c]	=	0.0
					tfimtar[:,c]	=	0.0	
		
		if(visflag):
			#print('Copyting visibilities flags ...........')
			nrecs	=	len(tfwsrc)
			nrect	=	len(tfwtar)
			nchans	=	len(tfwsrc[0])
			nchant	=	len(tfwtar[0])
			if((nrecs!=nrect) or (nchans!=nchant)):
				print('Error !	Vsibilities do not natch !!!!')
				return 3
			for r in range(0,nrecs):
				for c in range(0,nchans):
					if(tfwsrc[r,c]<=0.0):
						tfwtar[r,c]		=	0.0
						tfretar[r,c]	=	0.0
						tfimtar[r,c]	=	0.0		
		
		tarfracf		=	1.0-float(np.count_nonzero(tfwtar>0.0))/np.size(tfwtar)
		print('					target	= %.2f'%(tarfracf))
		flgfrac.append([srcfrac,tarfrac,tarfracf])
	return [0.0,flgfrac,bldatatar,bindxsrc[0],bindxtar[0],bl]
#	---------------------------------------------------------------------------------------------




#	---------------------------------------------------	Main programme	-----------------------------------------------------------------	

b_start	=	tm.time()

srcfile		=	fits.open(srcname)		
targetfile	=	fits.open(targetname)	

datasrc		=	srcfile[0].data	
nrecsrc		=	len(datasrc)
datatar		=	targetfile[0].data	
nrectar		=	len(datatar)

blid	=	[]
for a in range (1,ANTS):
	for b in range (a+1,ANTS+1):
		blid.append([a,b,256*a+b])
blid	=	np.array(blid)
nbase	=	len(blid)

print('\nIdeally total baselines = %d'%nbase)
print('Source records	=	%d	\nTarget records	=	%d\n'%(nrecsrc,nrectar))
print('\nCopying flags for each baseline\n')

if(parallel):
	#	parallel
	print('Sorry !	Parallel execution is not yet supported ...........')
else:
	#	Serial
	blstat	=	[]
	for i in range (0,nbase):
		bls	=	copybase(blid[i],datasrc,datatar)
		blstat.append(bls)

#	Copying baseline
		
print('\nCopying baselines\n')
for i in range (0,nbase):
	if(blstat[i][0]<1.5):
		bldata	=	blstat[i][2]
		lenbase	=	len(bldata)
		if(lenbase!=len(blstat[i][4])):
			print('Lengths do not match !!!!')
		for k in range (0,lenbase):
			if(datatar[blstat[i][4][k]].par('BASELINE')!=blstat[i][5][2]):
				print('Baselines do not match !!!!	%d	%d	%d'%(k,datatar[blstat[i][4][k]].par('BASELINE'),blstat[i][5][2]))
			np.copyto(datatar[blstat[i][4][k]].data,bldata[k])
		print('Baseline	%d-%d	Copied'%(blid[i][0],blid[i][1]))
	else:
		print('Baseline	%d-%d	No data'%(blid[i][0],blid[i][1]))		

flstat	=	[]
for i in range (0,nbase):
	flstat.append(blstat[i][1])
flstat	=	np.array(flstat)
srcf	=	[np.mean(flstat[:,0,0]),np.mean(flstat[:,1,0])]
tarfb	=	[np.mean(flstat[:,0,1]),np.mean(flstat[:,1,1])]
tarfa	=	[np.mean(flstat[:,0,2]),np.mean(flstat[:,1,2])]

print('\nSource		%.2f	%.2f'%(srcf[0],srcf[1]))
print('\nTarget before	%.2f	%.2f'%(tarfb[0],tarfb[1]))
print('\nTarget after	%.2f	%.2f'%(tarfa[0],tarfa[1]))

if(writeout):		
	targetfile.writeto(outname,output_verify='warn',overwrite=True)

srcfile.close()
targetfile.close()

b_end	=	tm.time()
print('\nDone in %d seconds\n'%(int(b_end-b_start)))





















































































































