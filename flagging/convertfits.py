import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from inputs import *




#	----------------------------------------------------------------------------------------
#
#	Function to convert UVFITS file to binary files and write in scratch direcotory
#
#	----------------------------------------------------------------------------------------

def uvfitstobinary(data,scratchdir,ugrids,vgrids,plotuv):
	
	ngroups	=	len(data)
	recids	=	np.arange(ngroups)
	
	#print('Parameters	'+str(data.parnames))
	print('Total no of records	=	%d'%ngroups)
	
	uvwarr		=	np.transpose(np.array([data.par(ukey),data.par(vkey),data.par(wkey)]))

	usorted		=	np.argsort(uvwarr[:,0])
	uvwarrus	=	uvwarr[usorted]
	recdataus	=	data[usorted]
	recidus		=	recids[usorted]
	
	uvbins		=	[]

	ustart		=	0
	uglen	=	1 + int(ngroups/ugrids)
	
	for ug in range (0,ugrids):
	
		ustop		=	ustart + uglen
		
		if(ustop>(ngroups+1)):
			ustop	=	ngroups+1
			
		ugridarr	=	uvwarrus[ustart:ustop]
		ugsize		=	len(ugridarr)
		#print('%d		%d		%d	%d'%(ug,ugsize,ustart,ustop))
		vsorted		=	np.argsort(ugridarr[:,1])
		ugarrvs		=	ugridarr[vsorted]
	
		urecdata	=	recdataus[ustart:ustop]
		urecvs		=	urecdata[vsorted]
		
		urecid		=	recidus[ustart:ustop]
		urecidvs	=	urecid[vsorted]
	
		vstart		=	0
		vglen		=	1 + int(ugsize/vgrids)
	
		for vg in range (0,vgrids):
			vstop		=	vstart + vglen
			if(vstop>(ugsize+1)):
				vstop	=	ugsize+1
			vgridarr	=	ugarrvs[vstart:vstop]	
			vgsize		=	len(vgridarr)
			#print('%d	%d	%d		%d	%d'%(ug,vg,vgsize,vstart,vstop))
			print('%d	%d		writing...'%(ug,vg))
			
			vrecdata	=	urecvs[vstart:vstop]
			vrecid		=	urecidvs[vstart:vstop]
			uvbins.append([ug,vg,vrecdata,vrecid])
		
			#print(vgridarr[100],vrecdata.par(ukey)[100],vrecdata.par(vkey)[100],vrecdata.par(wkey)[100],vrecdata.par('BASELINE')[100])
			if(plotuv):
				colsel		=	'b.'
				if((ug+vg)%2):
					colsel	=	'r.'
				plt.plot(vgridarr[:,0],vgridarr[:,1],colsel,markersize=0.5)
		
			vstart		=	vstop
		
		ustart		=	ustop
	
	if(plotuv):
		plt.xlim([-0.0001,0.0001])
		plt.ylim([-0.0001,0.0001])
		plt.show()
	
	uvbinstats	=	[]
		
	for ug in range (0,ugrids):
		for vg in range (0,vgrids):
			
			datauvg		=	uvbins[ug*vgrids+vg]
			
			datarray	=	datauvg[2].data
			#print(np.shape(datarray))
			
			uvbinstats.append([ug,vg,1,len(datarray),len(datarray[0,0,0,0])])
			
			bindetails	=	[datauvg[3],datauvg[2].par(ukey),datauvg[2].par(vkey),datauvg[2].par(wkey),datauvg[2].par('BASELINE')]
			bindetails	=	np.array(bindetails)
			bindetails	=	np.transpose(bindetails)
			
			np.savetxt(scratchdir+'uvbindetails_%d_%d_0.txt'%(ug,vg),bindetails,fmt='%d	%e	%e	%e	%d')
			
			fname		=	open(scratchdir+'uvbin_%d_%d_0.array'%(ug,vg),'wb')		
			for p in range (0,npols):
				for d in range (0,3):					
					fname.write(datarray[:,0,0,0,:,p,d].astype('float32').tobytes())
			fname.close()	
	
	uvbinstats		=	np.array(uvbinstats)
	np.savetxt(scratchdir+'uvbin_status.txt', uvbinstats, fmt='%4d	%4d	%4d	%4d	%4d')
		
	return

#	------------------------------------------------------------------------------------------------------







#	----------------------------------------------------------------------------------------
#
#	Function to convert binary files from the scratch directory to UVFITS file
#
#	----------------------------------------------------------------------------------------

def uvfitsfrombinary(data,scratchdir,ugrids,vgrids):
		
	for ug in range (0,ugrids):
		for vg in range (0,vgrids):
		
			bindetails	=	np.loadtxt(scratchdir+'uvbindetails_%d_%d_0.txt'%(ug,vg))
			recids		=	bindetails[:,0].astype('int32')
			baselines	=	bindetails[:,4]
			
			temparr	=	np.fromfile(scratchdir+'uvbin_%d_%d_0_f.array'%(ug,vg),dtype='float32',count=-1,sep="")
			
			if((len(temparr)%len(recids)) or (len(temparr)%npols) or (len(temparr)%3)):
				print('Length Mismatch !!!!!!! ....... exiting...')
				return 1
			
			nchan	=	int(len(temparr)/(npols*3*len(recids)))
			
			temparr	=	np.reshape(temparr,(npols,3,len(recids),nchan))
			
			for p in range (0,npols):
				for d in range (0,3):						
					#print("Records,channels		"+str(np.shape(temparr)))
				
					for l in range (0,len(recids)):
						#print(data[recids[l]].par(ukey),data[recids[l]].par(vkey),data[recids[l]].par(wkey),bindetails[l,1],bindetails[l,2],bindetails[l,3])
						np.copyto( data.data[ recids[l] ,0,0,0,:,p,d] , temparr[ p, d, l ] )	
			
			print('Copied	%d %d	of	%d %d'%(ug,vg,ugrids,vgrids))

	return 0

#	----------------------------------------------------------------------------------------








#	----------------------------------------------------------------------------------------
#
#	Function to convert baselines to binary files and write in scratch direcotory
#
#	----------------------------------------------------------------------------------------

def baselinestobinary(ANTS,data,scratchdir,breaktimesec):
	
	#print (data.par('BASELINE'))
	#sys.exit()
	
	blid	=	[]
	for a in range (1,ANTS):
		for b in range (a+1,ANTS+1):
			blid.append([a,b,256*a+b])
	blid	=	np.array(blid)
	nbase	=	len(blid)
	nvis	=	len(data)

	print('\nWriting Total baselines = %d	Visibilities = %ld'%(nbase,nvis))
	
	datepars	=	[]
	for ipar in range (0,len(data.parnames)):
		if (data.parnames[ipar]=='DATE'):
			datepars.append(ipar)
	
	print ("Date parameters "+str(datepars))
	if (len(datepars)>2):
		print ("More than two DATE parameters !!! EXITING ...")
		return 1
	
	mjdates	=	data.par(datepars[1]) + data.par(datepars[0]) 
	
	breaks	=	[]
	scanlens=	[]
	prebrk	=	-1
	
	for dt in range (0,len(mjdates)-2):
		if ((mjdates[dt+1]-mjdates[dt])*86400.0>breaktimesec):			
			breaks.append(dt)
			scanlens.append(dt-prebrk)
			prebrk	=	dt
	
	breaks.append(len(mjdates)-1)
	scanlens.append(len(mjdates)-1-prebrk)
	breaks		=	np.array(breaks)
	scanlens	=	np.array(scanlens)	
	totscans	=	len(breaks)
	
	print ("Found total		%d scans"%totscans)
	print (breaks)
	print (scanlens)
	
	scandetails	=	np.transpose(np.concatenate((np.array([totscans]),breaks),axis=0))
	np.savetxt(scratchdir+'scandetails.txt',scandetails, fmt = "%ld")
		
	blstatus	=	[]
	scstart		=	0
	scend		=	0
	
	for scan in range (0,totscans):
		
		scend	=	breaks[scan] + 1
			
		print ("\n\nScan %d / %d		Visibilities	%ld - %ld"%(scan,totscans,scstart,scend))
		scdata	=	data[scstart:scend]	
		print ("Total visibilities	= %ld"%len(scdata))	
		
		for bline in range (0,nbase):
	
			bl		=	blid[bline]		
			bindx	=	np.where(scdata.par('BASELINE')==bl[2])[0]
			bdata	=	scdata[bindx]
		
			if(len(bdata)==0):
				print('Writing scan	%d / %d		Baseline	%d-%d		No data'%(scan,totscans,bl[0],bl[1]))
				blstatus.append([bl[0],bl[1],0,0,0,scan])
				continue
				
			print('Writing scan	%d / %d		Baseline	%d-%d		writing...'%(scan,totscans,bl[0],bl[1]))
		
			bindetails	=	[bindx,bdata.par(ukey),bdata.par(vkey),bdata.par(wkey),bdata.par('BASELINE'),scan*np.ones(bindx.shape,dtype=int)]
			bindetails	=	np.array(bindetails)
			bindetails	=	np.transpose(bindetails)
			
			np.savetxt(scratchdir+'baselinedetails_%d_%d_%d.txt'%(bl[0],bl[1],scan),bindetails,fmt='%d	%e	%e	%e	%d	%d')
		
			bldata	=	bdata.data
			blstatus.append([bl[0],bl[1],1,len(bldata),len(bldata[0,0,0,0]),scan])
			
			fname		=	open(scratchdir+'baseline_%d_%d_%d.array'%(bl[0],bl[1],scan),'wb')
			for p in range (0,npols):
				
				tfdata	=	bldata[:,0,0,0,:,p,:]		#	Time frequency data			
				lent	=	len(tfdata)
				
				for d in range (0,3):
					for r in range (0,lent):
						fname.write(tfdata[r,:,d].astype('float32').tobytes())
			fname.close()
		
		scstart		=	scend			
	
	blstatus	=	np.array(blstatus)		
	np.savetxt(scratchdir+'baseline_status.txt', blstatus, fmt='%4d	%4d	%4d	%4d	%4d	%4d')	
		
	return 0

#	------------------------------------------------------------------------------------------------------







#	----------------------------------------------------------------------------------------
#
#	Function to convert baselines from binary files in scratch direcotory to UVFITS
#
#	----------------------------------------------------------------------------------------

def baselinesfrombinary(ANTS,data,scratchdir):
	
	blid	=	[]
	for a in range (1,ANTS):
		for b in range (a+1,ANTS+1):
			blid.append([a,b,256*a+b])
	blid	=	np.array(blid)
	nbase	=	len(blid)
	nvis	=	len(data)
	
	scanbreaks	=	np.loadtxt(scratchdir+'scandetails.txt',dtype=int,ndmin=1)
	
	nscan		=	scanbreaks[0]
	scanbreaks	=	scanbreaks[1:]
		
	blstatus	=	np.loadtxt(scratchdir+'baseline_status.txt', dtype=int)
	baselineflag=	np.loadtxt(scratchdir+'badbase.list',dtype=int)
		
	print('\nReading Total baselines = %d	Scans = %d\n'%(nbase,nscan))
	
	scanstart	=	0
	for scan in range (0,nscan):	
		
		if (scan > 0) :	
			scanstart	=	scanbreaks[scan-1] + 1
		if (scan < (nscan-1)) :
			scanend		=	scanbreaks[scan] + 1
		else :
			scanend		=	nvis + 1
			
		scdata	=	data[scanstart:scanend]
		
		for bline in range (0,nbase):
			
			bl		=	blid[bline]			
			if(blstatus[scan*nbase + bline,2]==0):
				#print('\nReading scan	%d / %d		Baseline	%d-%d		No data'%(scan,nscan,bl[0],bl[1]))
				continue
			
			bindetails	=	np.loadtxt(scratchdir+'baselinedetails_%d_%d_%d.txt'%(bl[0],bl[1],scan))
			if (len(bindetails.shape)==1):
				bindetails	=	np.array([bindetails])
			if (len(bindetails.shape)<1):
				print ("Blanck details for baseline	%d %d	scan %d ... CONTINUING !!!"%(bl[0],bl[1],scan))
				continue
			
			recids		=	bindetails[:,0].astype('int32')
			baselines	=	bindetails[:,4]
		
			temparr		=	np.fromfile(scratchdir+'baseline_%d_%d_%d_f.array'%(bl[0],bl[1],scan),dtype='float32',count=-1,sep="")
			
			if((len(temparr)%len(recids)) or (len(temparr)%npols) or (len(temparr)%3)):
				print('Length Mismatch !!!!!!! ....... exiting...')
				return 1
			
			nchan	=	int(len(temparr)/(npols*3*len(recids)))
			
			temparr	=	np.reshape(temparr,(npols,3,len(recids),nchan))	
			temp00	=	np.zeros(np.shape(temparr), dtype='float32')
										
			for p in range (0,npols):
				#for d in range (0,3):
				for d in range (2,3):				
					#print("Records,channels		"+str(np.shape(temparr)))
					if (baselineflag[scan*nbase + bline,3+p]):
						for l in range (0,len(recids)):
							#print(data[recids[l]].par('UU'),data[recids[l]].par('VV'),data[recids[l]].par('WW'),bindetails[l,1],bindetails[l,2],bindetails[l,3])
							np.copyto( scdata.data[ recids[l] ,0,0,0,:,p,d] , temparr[ p, d, l ] )
					else:
						for l in range (0,len(recids)):
							#print(data[recids[l]].par('UU'),data[recids[l]].par('VV'),data[recids[l]].par('WW'),bindetails[l,1],bindetails[l,2],bindetails[l,3])
							np.copyto( scdata.data[ recids[l] ,0,0,0,:,p,d] , temp00[ p, d, l ] )	
		
			del (temp00)
			del (temparr)
			
			print('\nCopied	Scan	%d / %d		baseline	%d - %d'%(scan, nscan,bl[0],bl[1]))		
				
	return 0

#	------------------------------------------------------------------------------------------------------





#	----------------------------------------------------------------------------------------
#
#	Function auto-correlations if present
#
#	----------------------------------------------------------------------------------------

def flagautocorr(ANTS,data):
	
	blid	=	[]
	for a in range (1,ANTS+1):
		blid.append([a,a,256*a+a])
		
	blid	=	np.array(blid)
	nbase	=	len(blid)
	nvis	=	len(data)
	
	print("Flagging %d auto-correlations..."%nbase)
	
	for bline in range (0,len(data.par('BASELINE'))):	
		if(np.isin(data.par('BASELINE')[bline],blid)):
			data[bline].data[:,:,:,:,:,2]	=	0.0
	
	return(nbase)

#	----------------------------------------------------------------------------------------



#	----------------------------------------------------------------------------------------
#
#	Function to flag a list of channels in all baselines
#
#	----------------------------------------------------------------------------------------

def flagchanlist(ANTS,data,hopelesschans):
	
	nhopeless		=	len(hopelesschans)
							
	if(nhopeless>0):		
		print('\nFlagging %d hopeless channels')
		
		nchan	=	data.data.shape[4]
		
		for hc in hopelesschans:
			if(hc<nchan):
				data.data[:,:,:,:,hc,:,2] = 0.0
				
	else:
		print('\nNo bad channels to flag')
		
	return(nhopeless)

#	------------------------------------------------------------------------------------------------------





#	-----------------------------------------------------------------------------------------------------
#
#	Function to show time-frequency plot for one single baseline
#
#	-----------------------------------------------------------------------------------------------------

def showbasecomparison(bl,data,data2,dopl):
	
	flfrac	=	[]
	
	bdata	=	data[data.par('BASELINE')==bl[2]]
	
	if(len(bdata)==0):
		print('Baseline	%d-%d		No data'%(bl[0],bl[1]))
		return (np.ones(2*npols))
	
	if(dopl):
		fig=plt.figure(figsize=(12,8))
		
	bldata	=	bdata.data
	
	for p in range (0,npols):
			
		tfdata	=	bldata[:,0,0,0,:,p,:]			#	Time frequency data
		tfre	=	tfdata[:,:,0]					#	Real part
		tfim	=	tfdata[:,:,1]					#	Imaginary Part
		tfw		=	tfdata[:,:,2]					#	Weight
		tfamp	=	np.sqrt(tfre*tfre+tfim*tfim)	#	Amplitude
		
		tfamp[tfw<=0.0]	=	0.0
		
		flfrac.append(1.0-float(np.count_nonzero(tfw>0.0))/np.size(tfw))
		
		if(dopl):
			ax=fig.add_subplot(npols,2,npols*p+1)
			plt.imshow(tfamp,interpolation='none',aspect='auto')
			plt.colorbar()
	
	bdata	=	data2[data.par('BASELINE')==bl[2]]
		
	bldata	=	bdata.data
	
	for p in range (0,npols):
			
		tfdata	=	bldata[:,0,0,0,:,p,:]			#	Time frequency data
		tfre	=	tfdata[:,:,0]					#	Real part
		tfim	=	tfdata[:,:,1]					#	Imaginary Part
		tfw		=	tfdata[:,:,2]					#	Weight
		tfamp	=	np.sqrt(tfre*tfre+tfim*tfim)	#	Amplitude
		
		tfamp[tfw<=0.0]	=	0.0
		
		flfrac.append(1.0-float(np.count_nonzero(tfw>0.0))/np.size(tfw))
		
		if(dopl):
			ax=fig.add_subplot(npols,2,npols*p+2)
			plt.imshow(tfamp,interpolation='none',aspect='auto')
			plt.colorbar()
			
	print('Baseline	%d-%d	flagged fractions	'%(bl[0],bl[1])),
	for p in range (0,npols):
		print('%.2f %.2f	'%(flfrac[p],flfrac[p+npols])),
	print('\n')
	
	if(dopl):
		plt.show()
	
	return (np.array(flfrac))

#	----------------------------------------------------------------------------------------------------


#	-----------------------------------------------------------------------------------------------------
#
#	Function to find out completely flagged baselines
#
#	-----------------------------------------------------------------------------------------------------

def returnbadbaselines(bl,data):
	
	baseflagged	=	1
	
	bdata	=	data[data.par('BASELINE')==bl[2]]
	
	if(len(bdata)==0):
		baseflagged	=	1
		return (baseflagged)
			
	bldata	=	bdata.data
	
	for p in range (0,npols):
			
		tfdata	=	bldata[:,0,0,0,:,p,:]			#	Time frequency data
		tfw		=	tfdata[:,:,2]					#	Weight
		if(len(tfw[tfw>0.0])>0):
			baseflagged	=	0
		
	return (baseflagged)

#	----------------------------------------------------------------------------------------------------


	

































