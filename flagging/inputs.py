infilename		=	''						#	Name of the input FITS file
outfilename		=	''						#	Name of the output file
allbadchanfile	=	''
scratchdir		=	'scratch/'				#	Scratch directory with at least twice the size of the inputs FITS file of space

CLEARSCRATCH	=	True											#	Clear the scratch directory
ukey			=	'UU---SIN'
vkey			=	'VV---SIN'
wkey			=	'WW---SIN'
npols			=	2												#	Number of polarizations

ugrids			=	1
vgrids			=	1
plotuv			=	False

FLAGAUTOCORR	=	0											#	Flag auto-correlations?

CONVERTFITS		=	1											#	Convert FITS to binary ?	
DOFLAG			=	1											#	Do flagging ?
ANTS			=	36											#	Total number of antennas
THREADS			=	4											#	Number of parallel threads
FLAGMODE		=	'baseline'		 									#	'baseline' / 'uvbin' 
SCANFLAGMEAN	=	['mean_rms',	10.1,	10.1,	0.01]		#	[FLAGON,	tolerance_mean,	tolearnce_rms,	min fraction (should be > 0.01)]]			ONLY for 'baseline'
BREAKTIMESEC	=	300.0
READBACK		=	1											#	Read back baselines ?
SHOWBASE		=	0											#	SHOW baseline stats ?
SHOWTF			=	0											#	Show time-frequency plots ?
WRITEOUT		=	1											#	Write output ?
BLOCKPOW		=	0.8											#	Power law for Block non-Gaussianity (DON'T CHANGE UNLESS YOU KNOW WHAT IT IS !)

#	For bandpass calibrated data

FLAGPARS1		=	[	[ 'chan_ind',	'rms',		'median',	'am',	2.0,	0.1,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'median',	'am',	2.0,	0.1,	0,	'',				0,	0,	0.0,	0.0],
						[ 'vis_ind',	'mean',		'mean',		'am',	3.0,	0.1,	0,	'',				0,	0,	0.0,	0.0]	]

FLAGPARS2		=	[	[ 'rec_ind',	'rms',		'median',	'am',	5.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'median',	'am',	5.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'median',	'am',	2.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'median',	'am',	2.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'mean',		'mean',		'am',	3.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'mean',		'mean',		'am',	3.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'rec_ind',	'rms',		'mean',		'am',	2.0,	0.01,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'mean',		'am',	2.0,	0.01,	0,	'',				0,	0,	0.0,	0.0]	]
						
FLAGPARS3		=	[	[ 'chan_ind',	'mean',		'mean',		'am',	1.5,	0.1,	0,	'',				0,	0,	0.0,	0.0],
						[ 'chan_ind',	'rms',		'median',	'am',	1.2,	0.1,	0,	'',				0,	0,	0.0,	0.0]	]

FLAGPARS4		=	[	[ 'chan_ind',	'rms',		'mean',		'am',	1.5,	0.1,	0,	'',				0,	0,	0.0,	0.0]	]

#
#	FLAGON			=	'mean'		/	'rms'		/	'mean_rms'
#	STATYPE			=	'mean'		/	'median'
#	DATATYPE		=	're'		/	'im'		/	'am'		(/	'ph'	---- NOT YET SUPPORTED)
#	ORDER			=	'ascending'			(/	'descending' --- NOT YET SUPPORTED)

#	Vis individual	=	[ 'vis_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	0, 				fit order,  '', 	0, 			0, 			0, 				0			 	]

#	Chan individual	=	[ 'chan_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	'',		0,			0,			0,				0			 	]

#	Rec individual	=	[ 'rec_ind',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,			0,	'',		0,			0,			0,				0			 	]

#	Vis block		=	[ 'vis_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	chan_block,	rec_block,	chan_max_frac,	rec_max_frac 	]

#	Chan block		=	[ 'chan_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	chan_block,	0,			chan_max_frac,	0 				]

#	Rec block		=	[ 'rec_block',	FLAGON,	STATYPE,	DATATYPE,	tolerance,	min fraction,	fit order,	ORDER,	0,			rec_block,	0,				rec_max_frac 	]







































































