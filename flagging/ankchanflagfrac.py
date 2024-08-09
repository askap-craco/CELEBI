import numpy as np
from astropy.io import fits

ANTS			=	30
breaktimesec	=	20000.0
intimesec		=	15.0
npols			=	2

bands	=	['y','a','b','c','d','e','f','g','h','i','j']

for band in bands:

	infilename		=	'/data/apurba/egs_all/egs_q/fits/egs_102_'+band+'_wt.fits'

	infile	=	fits.open(infilename)
	data	=	infile[0].data
	nchan	=	infile[0].header['NAXIS4']
	freqghz	=	((np.arange(0,nchan,1.0,dtype=float)-infile[0].header['CRPIX4'])*infile[0].header['CDELT4'] + infile[0].header['CRVAL4'])/1.0e6

	blid	=	[]
	for a in range (1,ANTS):
		for b in range (a+1,ANTS+1):
			blid.append([a,b,256*a+b])
	blid	=	np.array(blid)
	nbase	=	len(blid)
	nvis	=	len(data)

	print('\nTotal channels	=	%d'%nchan)
	print('\nWriting Total baselines = %d	Visibilities = %ld'%(nbase,nvis))

	datepars	=	[]
	for ipar in range (0,len(data.parnames)):
		if (data.parnames[ipar]=='DATE'):
			datepars.append(ipar)

	print ("Date parameters "+str(datepars))
	if (len(datepars)>2):
		print ("More than two DATE parameters !!! EXITING ...")

	mjdates	=	data.par(datepars[1]) #+ data.par(datepars[0]) 

	missed	=	0
	for dt in range (0,len(mjdates)-1):
		if (((mjdates[dt+1]-mjdates[dt])*86400.0>intimesec) and ((mjdates[dt+1]-mjdates[dt])*86400.0<breaktimesec)):			
			missed	=	missed + int((mjdates[dt+1]-mjdates[dt])*86400.0/8.05)

	print('Missed %d records'%missed)

	breaks	=	[]
	for dt in range (0,len(mjdates)-1):
		if ((mjdates[dt+1]-mjdates[dt])*86400.0>breaktimesec):			
			breaks.append(dt)
		elif ((mjdates[dt+1]-mjdates[dt])<0.0):# and (data.par(datepars[0])[dt+1] - data.par(datepars[0])[dt+1] > 0.0)):
			breaks.append(dt)
	breaks	=	np.array(breaks)
			
	if(len(breaks)==0):
		breaks	=	np.array([len(data)])

	scanlens	=	breaks[1:]-breaks[:-1]	
	print ("Found total		%d scans"%(len(breaks)+1))
	print np.concatenate((np.array([breaks[0]]),scanlens,np.array([nvis-breaks[-1]-1])),axis=0)

	chanstat	=	np.zeros(nchan,dtype=int)
	chantot		=	0
	scstart		=	0
	scend		=	0

	for scan in range (0,len(breaks)+1):
		
		if (scan < len(breaks)):
			scend	=	breaks[scan] + 1
		else :
			scend	=	nvis
			
		print ("\n\nScan %d / %d		Visibilities	%ld - %ld"%(scan,len(breaks)+1,scstart,scend))
		scdata	=	data[scstart:scend]	
		print ("Total visibilities	= %ld"%len(scdata))	
		
		basemaxrec	=	0
		for bline in range (0,nbase):

			bl		=	blid[bline]		
			bindx	=	np.where(scdata.par('BASELINE')==bl[2])[0]
			bdata	=	scdata[bindx]
		
			if(len(bdata)==0):
				print('Scan	%d / %d		Baseline	%d-%d		No data'%(scan,len(breaks)+1,bl[0],bl[1]))
				continue
				
			bldata	=	bdata.data
			
			for p in range (0,npols):			
				tfdata	=	bldata[:,0,0,0,:,p,:]		#	Time frequency data			
				basemaxrec	=	max(basemaxrec,len(tfdata))
				for c in range (0,nchan):
					tffdata		=	tfdata[:,c]
					#print len(tffdata[tffdata[:,2]>0.0])	
					chanstat[c]	=	chanstat[c] + len(tffdata[tffdata[:,2]>0.0])			
		
		chantot		=	chantot + nbase*basemaxrec*2
		scstart		=	scend			

	chantot		=	chantot + (missed - basemaxrec/10)*nbase*2

	flfrac		=	np.array([freqghz,chanstat/float(chantot)]).T

	infile.close()	
	
	np.savetxt('/data/apurba/egs_all/egs_q/fits/egs_102_'+band+'_wt.txt',flfrac,fmt='%.6f	%.6f')
	
































