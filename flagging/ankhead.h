typedef	struct	{			//	Structure containing baseline parameters
	
	int		exists;				//	Whether the baseline data exists
	int		anta;				//	Antenna A
	int		antb;				//	Antenna B
	int		recsize;			//	Record length
	int		chansize;			//	Channel number
	int		polsize;			//	Polarizations
	int		datsize;			//	Number of data arrays
	float	****data;			//	Data arrays
	float	flfrac[10];			//	Initial flag fraction
	float	flfrac_after[10];	//	Final flag fraction

}	BaseParType;



typedef struct {				//	Structure containing flagging parameters
	
	int		whatflag;			//	Do flagging? ------------  	0 = No, 1 = meanonly, 2 = rmsonly, 3 = mean+rms
	int		doflag;				//	Do flagging? ------------  	0 = No, 1 = meanonly, 2 = rmsonly, 3 = mean+rms
	int		domean;				//	Do flagging based on mean?	0 = median, 1 = mean 
	int		qtype;				//	Data array to fo flagging based on ---- 0 = Real, 1 = Imaginary, 2 = Amplitude, 3 = Phase
	float	tolerance;			//	Tolerance factor
	float	min_fraction;		//	Minimum fraction to be kept
	int		fitorder;			//	Fit order for channels
	int 	ascending;			//	Ascending (1) or descending (0)
	int		chanblockfac;		//	Channel block factor
	int		recblockfac;		//	Record block factor
	float	chanmaxfrac;		//	Channel maximum fraction
	float 	recmaxfrac;			//	Record maximum fraction

}	FlagParType;


//	-----------------------------------------------------------------------------------
//					Function prototypes
//	-----------------------------------------------------------------------------------

//													exc functions

int processbaselines(int ants, char *statusfile, int *nbaselines, char *fnames, char *flagparfile, FILE *badbasefile, char *scanfile);

//													flag functions

int	flagbaseline (BaseParType *basepar, float bl_minfrac, float blockpow, int flagrounds, FlagParType *flagparams, float **datarray, float *goutarr,
					float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d, float **statsarray);

//													Init functions

int init_gout (float *goutarr);
int init_data (BaseParType *basepar, float **datarray, int i, int qtype);
int	init_baseflag(FILE *fpflagpar, int *bl_doflag, float *mean_tolerance, float *rms_tolerance, float *bl_minfrac, float *blockpow, int *flagrounds, int *threads, int uvgrid[2]);
int	init_flagpars (FILE *fpflagpar, int  curound, FlagParType *flagpar);

//													basic functions

int	flag_chan_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, 
							float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d);
int	flag_rec_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, float **buffarrc, 
												float **buffarrr, float **buffdev1d);
int	flag_vis_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, 
							float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d);
int	flag_chan_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int cw, float blockpow, 
										float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d);
int	flag_rec_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int rw, float blockpow,
									float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d);
int	flag_vis_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int rw, int cw, float blockpow,
									float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d);
float findbasestat (BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *statsarray, 
										float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d);
int findbadbase (BaseParType **baseparams, int nscan, int nbase, float ****basestatsarray, int bl_doflag, float mean_tolerance, float rms_tolerance, 
														float bl_minfrac, int maxize, float *goutarr, FILE *badbasefile);

//													auxiliary functions

int		count_nonzero (float *data, int nels);
int		count_nonzero_2d (float **data, int n1, int n2);
float	find_mean (float *data, int nels);
float	find_rms (float *data, int nels, float mean);
float	find_median (float *data, int nels, float *temparr);
float	find_mad (float *data, int nels, float median, float *temparr);
int	 	dopolyfit (float *xarr, float *yarr, float *yerr, double *coeffs, int lendata, int order);

//	-----------------------------------------------------------------------------------

#define	GOUTFILE	"glogout.dat"
#define	GOUTLEN		10000
#define	DATARRS		3
#define	POLARRS		2















