#
#	Script for estimating DM of polcal by maximizing S/N
#
#								AB, July 2024

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt

print("\n******************************************************\n")
print("This is a DM estimation by maximizing S/N")
print("\n******************************************************\n")

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FRB name string> <current DM> <search range> <DM step> <t_avg> <fcen_mhz> <bw mhz>\n")
		
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<8):
	print_instructions()
	sys.exit()

frbname		=	sys.argv[1]					#	FRB name string (YYMMDDx)
dm0s		=	sys.argv[2]					#	Current DM (string)
edm			=	float(sys.argv[3])			#	DM serach range
ddm			=	float(sys.argv[4])			#	DM step
tavg		=	int(sys.argv[5])			#	Time averaging factor
f0mhz		=	float(sys.argv[6])			#	Central frequency in MHz
bwmhz		=	float(sys.argv[7])			#	Bandwidth in MHz

#	-------------------------	Initialize thing	-----------------------

dcon		=	0.0							#	Dispersion constant
dfmhz		=	1.0							#	Channel width in MHz
dtns		=	tavg*1.0e3/bwmhz			#	Time resolution in ns
dm0			=	float(dm0s)

dmarr		=	np.arange(-edm,edm,ddm)
pkarr		=	np.zeros(dmarr.shape,dtype=float)

dspec0		=	np.load(frbname+'_polcal_I_dynspec_'+dm0s+'.npy')	
qdspec0		=	np.load(frbname+'_polcal_Q_dynspec_'+dm0s+'.npy')	
udspec0		=	np.load(frbname+'_polcal_U_dynspec_'+dm0s+'.npy')	
vdspec0		=	np.load(frbname+'_polcal_V_dynspec_'+dm0s+'.npy')	
dspec		=	np.zeros(dspec0.shape, dtype=float)	
qdspec		=	np.zeros(qdspec0.shape, dtype=float)	
udspec		=	np.zeros(udspec0.shape, dtype=float)	
vdspec		=	np.zeros(vdspec0.shape, dtype=float)		
tser		=	np.zeros(dspec0.shape[1], dtype=float)	

nc0			=	dspec0.shape[0]
nt0			=	dspec0.shape[1]
nt			=	tavg*int(nt0/tavg)

print("Originally - %d channels - %d times - DM = %.2f"%(nc0, nt0, dm0))
print("Frequency resolution = %.2f MHz"%(dfmhz))
print("Time resolution = %.2f ns"%(dtns))

for i in range(0,len(dmarr)):
	for c in range(0,nc0):
		cfmhz	=	f0mhz + (c-(float(nc0)/2) + 0.5)*dfmhz
		dlns	=	4.15 * dmarr[i] * 1.0e6 * ((1.0e3/cfmhz)**2 - (1.0e3/f0mhz)**2)
		dldt	=	int(np.rint(dlns/dtns))
		#print(c,dmarr[i],dldt)	
		dspec[c]=	np.roll(dspec0[c], dldt)
	tser	=	np.nanmean(dspec, axis=0)
	pkarr[i]=	np.nanmax(tser)

plt.plot(pkarr)
plt.savefig('polcal_dm_snmax.png')
plt.close()

optdm		=	dmarr[np.argmax(pkarr)]

print("Optimum -- deltaDM = %.3f - DM = %.3f"%(optdm, dm0+optdm))

for c in range(0,nc0):
	cfmhz	=	f0mhz + (c-(float(nc0)/2) + 0.5)*dfmhz
	dlns	=	4.15 * optdm * 1.0e6 * ((1.0e3/cfmhz)**2 - (1.0e3/f0mhz)**2)
	dldt	=	int(np.rint(dlns/dtns))
	#print(c,dmarr[i],dldt)	
	dspec[c]=	np.roll(dspec0[c], dldt)
	qdspec[c]=	np.roll(qdspec0[c], dldt)
	udspec[c]=	np.roll(udspec0[c], dldt)
	vdspec[c]=	np.roll(vdspec0[c], dldt)

np.save(frbname+'_polcal_I_dynspec_snmax.npy',dspec)	
np.save(frbname+'_polcal_Q_dynspec_snmax.npy',qdspec)	
np.save(frbname+'_polcal_U_dynspec_snmax.npy',udspec)	
np.save(frbname+'_polcal_V_dynspec_snmax.npy',vdspec)	

plt.imshow(dspec, aspect='auto')
plt.savefig('polcal_I_dynspec_snmax.png')
plt.close()










































