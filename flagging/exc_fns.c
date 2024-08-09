#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>
#	include <omp.h>


//	---------------------------------------------------------------------------------------
//				Function to read a single baseline
//															AB	26 July 2018															
//	---------------------------------------------------------------------------------------

int	read_a_baseline (BaseParType basepar,	char *namestring,	int iscan) {
	
	int		p,d,r,c;
	char	fname[100];
	
	FILE	*fp;
	sprintf(fname,"%s_%d_%d_%d.array",namestring,basepar.anta,basepar.antb,iscan);
	
	if((fp	=	fopen(fname,"r"))==NULL){
		printf("File not found	%s\n",fname);
		return 1;
	}
		
	for (p = 0; p < basepar.polsize; p++)
		for (d = 0; d < basepar.datsize; d++)			
			for (r=0;r<basepar.recsize;r++)
				fread(basepar.data[p][d][r], sizeof(float), basepar.chansize, fp);
	//printf("Test value	%f\n",basepar.data[0][2][100][12]);	
	fclose(fp);
			
	return 0;
}

//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to write a single baseline
//															AB	26 July 2018
//	---------------------------------------------------------------------------------------

int	write_a_baseline (BaseParType basepar,	char *namestring, int iscan){
	
	int		p,d,r,c;
	char	fname[100];
	float	test = 0.0;
	
	FILE	*fp;
	sprintf(fname,"%s_%d_%d_%d_f.array",namestring,basepar.anta,basepar.antb,iscan);
	
	if((fp	=	fopen(fname,"w"))==NULL){
		printf("Cannot write to	%s\n",fname);
		return 1;
	}
	
	for (p = 0; p < basepar.polsize; p++)
		for (d = 0; d < basepar.datsize; d++)			
			for (r=0;r<basepar.recsize;r++)
				fwrite(basepar.data[p][d][r], sizeof(float), basepar.chansize, fp);
	//printf("Test value	%f\n",basepar.data[0][2][100][12]);
	fclose(fp);
				
	return 0;
}

//	-----------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to read baseline status file
//															AB	26 July 2018
//															AB	26 November 2019
//	---------------------------------------------------------------------------------------

int read_base_status (int ants, char *statusfile, int nbase, int nscans, BaseParType **baseparams)	{
	
	int			i, j, k;
				
	FILE		*fp;
	if((fp		=	fopen(statusfile,"r"))==NULL) {
		printf("Error reading status file .........!!!");
		return (1);
	}
	
	for (j=0; j<nscans; j++) {
		for (i = 0; i < nbase; i++) {
		
			fscanf(fp,"%d",&baseparams[j][i].anta);
			fscanf(fp,"%d",&baseparams[j][i].antb);
			fscanf(fp,"%d",&baseparams[j][i].exists);
			fscanf(fp,"%d",&baseparams[j][i].recsize);
			fscanf(fp,"%d",&baseparams[j][i].chansize);
			fscanf(fp,"%d",&k);
			if (j!=k) {
				printf("Unrecognized scan structure ... Contact AB ...\n");
				return (1);
			}
			baseparams[j][i].polsize			=	POLARRS;
			baseparams[j][i].datsize			=	DATARRS;
			baseparams[j][i].flfrac[0]			=	baseparams[j][i].flfrac[1] = 0.0;
			baseparams[j][i].flfrac_after[0]	=	baseparams[j][i].flfrac_after[1] = 0.0;
			if(baseparams[j][i].exists==0){
				baseparams[j][i].flfrac[0]			=	baseparams[j][i].flfrac[1] = 1.0;
				baseparams[j][i].flfrac_after[0]	=	baseparams[j][i].flfrac_after[1] = 1.0;
			}
		}
	}
	
	fclose(fp);
	
	return (0);
}

//	-----------------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to read uvbin status file
//															AB	6 August 2018
//	---------------------------------------------------------------------------------------

int read_uvbin_status (int uvgrids[2], char *statusfile, int nbase, int nscans, BaseParType **baseparams) {
	
	int			i, j, k;
				
	FILE		*fp;
	if((fp		=	fopen(statusfile,"r"))==NULL) {
		printf("Error reading status file .........!!!");
		return (1);
	}
	
	for (j=0; j<nscans; j++) {
		for (i = 0; i < nbase; i++){
		
			fscanf(fp,"%d",&baseparams[j][i].anta);
			fscanf(fp,"%d",&baseparams[j][i].antb);
			fscanf(fp,"%d",&baseparams[j][i].exists);
			fscanf(fp,"%d",&baseparams[j][i].recsize);
			fscanf(fp,"%d",&baseparams[j][i].chansize);
			baseparams[j][i].polsize			=	POLARRS;
			baseparams[j][i].datsize			=	DATARRS;
			baseparams[j][i].flfrac[0]			=	baseparams[j][i].flfrac[1] = 0.0;
			baseparams[j][i].flfrac_after[0]	=	baseparams[j][i].flfrac_after[1] = 0.0;
			if(baseparams[j][i].exists==0){
				baseparams[j][i].flfrac[0]			=	baseparams[j][i].flfrac[1] = 1.0;
				baseparams[j][i].flfrac_after[0]	=	baseparams[j][i].flfrac_after[1] = 1.0;
			}
		}
	}
	
	fclose(fp);
	
	return (0);
}

//	-----------------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to process baselines 
//															AB	5 August 2018
//															AB	26 November 2019
//	---------------------------------------------------------------------------------------

int processbaselines(int ants, char *statusfile, int *nbaselines, char *fnames, char *flagparfile, FILE *badbasefile, char *scanfile) {
	
	int			i,p,d,r,c,t,j,k,curound,donex;
	int			nbase, flagmode, nbadbase, nscans;
	int			maxrecsize=0,maxchansize=0;
	int			uvgrid[2];
	BaseParType	**baseparams;
	FlagParType	*flagparams;
	int			bl_doflag, flagrounds, threads;
	float		mean_tolerance, rms_tolerance, bl_minfrac, blockpow;
	FILE		*fpflagpar, *fpscans;
	
	fpflagpar	=	fopen (flagparfile,"r");
	flagmode	=	init_baseflag (fpflagpar, &bl_doflag, &mean_tolerance, &rms_tolerance, &bl_minfrac, &blockpow, &flagrounds, &threads, uvgrid);
	//printf("%d	%d	%d	%d	%.2f	%.2f	%.2f	%.2f	%d %d\n", flagmode,bl_doflag,threads,flagrounds,mean_tolerance,rms_tolerance,blockpow,bl_minfrac,uvgrid[0],uvgrid[1]);
	
	flagparams	=	(FlagParType *) malloc (flagrounds * sizeof(FlagParType));
	
	for (curound = 0; curound < flagrounds; curound++)	{
		init_flagpars (fpflagpar, curound, &flagparams[curound]);
		//	printf("\n%d %d %d %d %.2f %.2f %d %d %d %d %.2f %.2f\n",flagparams[curound].whatflag,flagparams[curound].doflag,flagparams[curound].domean,
			//	flagparams[curound].qtype,flagparams[curound].tolerance,flagparams[curound].min_fraction,flagparams[curound].fitorder,
			//	flagparams[curound].ascending,flagparams[curound].chanblockfac,flagparams[curound].recblockfac,flagparams[curound].chanmaxfrac,flagparams[curound].recmaxfrac);
	}
	
	fclose (fpflagpar);
	
	if (flagmode) {
		nscans		=	1;
		nbase		=	uvgrid[0]*uvgrid[1];
	}
	else {
		fpscans		=	fopen (scanfile, "r");
		fscanf (fpscans,"%d",&nscans);		
		fclose (fpscans);	
		nbase		=	ants*(ants-1)/2;
	}
	
	nbaselines	=	&nbase;		
	printf ("Total scans	=	%d	baselines / uvbins	= %d\n",nscans,*nbaselines);	
	
	baseparams	=	(BaseParType **) malloc ( nscans * sizeof(BaseParType *));
	
	for (i=0;i<nscans;i++)
		baseparams[i]	=	(BaseParType *) malloc ( nbase * sizeof(BaseParType));
	
	if (flagmode)
		donex	=	read_uvbin_status(uvgrid, statusfile, nbase, nscans, baseparams);
	else
		donex	=	read_base_status(ants, statusfile, nbase, nscans, baseparams);
	
	if (donex) {
		printf("Exiting ... \n");
		return (1);
	}	
	
	for (j=0; j< nscans; j++) {
		for (i = 0; i < nbase; i++) {
			if (baseparams[j][i].chansize > maxchansize) 
				maxchansize = baseparams[j][i].chansize; 
			if (baseparams[j][i].recsize > maxrecsize) 
				maxrecsize = baseparams[j][i].recsize;
		} 
	}
	
	printf("\nMaximum channels	= %d	Records	= %d\n\n",maxchansize,maxrecsize);
	
	//	---------------------------------------------------------------------------------------------------
	//								Allocate buffer arrays
	//	---------------------------------------------------------------------------------------------------
	
	float		***datarray, *goutarr, ***buffarrc, ***buffarrr, *****baserawdata, ***buffdev1d, ****buffdev2d, ****basestatsarray;
	double		**coeffs;	
	
	basestatsarray	=	(float ****) malloc (nscans * sizeof(float ***));
	for (j=0; j< nscans; j++) {
		basestatsarray[j]	=	(float ***) malloc (nbase * sizeof(float **));
		for (i = 0; i < nbase; i++) {
			basestatsarray[j][i]	=	(float **) malloc (POLARRS * sizeof(float *));
			for (p=0; p<POLARRS; p++) 
				basestatsarray[j][i][p]	=	(float *) malloc (5 * sizeof(float));
		}
	}
	
	goutarr		=	(float *) malloc (GOUTLEN * sizeof(float));
	if ( init_gout (goutarr))	{
		free (goutarr);
		printf("Hmm... May be you forgot it at the clinic !!!\n");
		return (1);
	}	
	
	baserawdata	=	(float *****) malloc ( threads * sizeof(float ****) );
	datarray	=	(float ***) malloc ( threads * sizeof(float **) );
	buffarrc	=	(float ***) malloc ( threads * sizeof(float **) );
	buffarrr	=	(float ***) malloc ( threads * sizeof(float **) );	
	coeffs		=	(double **) malloc ( threads * sizeof(double *) );
	buffdev1d	=	(float ***) malloc ( threads * sizeof(float **) );
	buffdev2d	=	(float ****) malloc ( threads * sizeof(float ***) );
	
	for (t = 0; t < threads; t++) {
	
		buffarrc[t]		=	(float **) malloc (12 * sizeof(float *));
		buffarrr[t]		=	(float **) malloc (5 * sizeof(float *));
		datarray[t]		=	(float **) malloc ( maxrecsize * sizeof(float *) );
		
		for (r = 0; r < maxrecsize; r++)
			datarray[t][r]	=	(float *) malloc ( maxchansize * sizeof(float) );
			
		for (j = 0; j < 12; j++)
			buffarrc[t][j]	=	(float *) malloc (maxchansize * sizeof(float));
			
		for (j = 0; j < 5; j++)
			buffarrr[t][j]	=	(float *) malloc (maxrecsize * sizeof(float));
			
		coeffs[t]			=	(double *) malloc (10 * sizeof(double));
		
		baserawdata[t]		=	(float ****) malloc ( POLARRS * sizeof(float ***) );
		for (p=0; p<POLARRS; p++) {
			baserawdata[t][p]	=	(float ***)malloc( DATARRS *sizeof(float **));
			
			for (d=0; d<DATARRS; d++) {
				baserawdata[t][p][d]	=	(float **)malloc( maxrecsize *sizeof(float *));
				
				for (r=0; r<maxrecsize; r++)
					baserawdata[t][p][d][r]	=	(float *)malloc( maxchansize *sizeof(float));
			}
		}
		
		buffdev1d[t]	=	(float **) malloc ( 4 * sizeof(float *) );
		buffdev2d[t]	=	(float ***) malloc ( 6 * sizeof(float **) );
		
		for (j = 0; j < 4; j++)
			buffdev1d[t][j]	=	(float *) malloc ( maxrecsize*maxchansize * sizeof(float) );
		
		for (j = 0; j < 6; j++) {
			buffdev2d[t][j]	=	(float **) malloc ( maxrecsize * sizeof(float *) );
			
			for (r = 0; r < maxrecsize; r++)
				buffdev2d[t][j][r]	=	(float *) malloc ( maxchansize * sizeof(float) );
		}
	}	
	
	
	for (j=0; j < nscans; j++) {	
		//	---------------------------------------------------------------------------------------------------
		//							Process baselines in parallel
		//	---------------------------------------------------------------------------------------------------
	
		#	pragma omp parallel for num_threads (threads) shared (baseparams,fnames,fpflagpar,bl_minfrac,blockpow,flagrounds,flagparams,j) private (i,p,d,r,t)
			for (i = 0; i < nbase; i++) {
	
				if (baseparams[j][i].exists==0) {
					printf("Baseline	%d %d		Scan %d		Flagged......\n",baseparams[j][i].anta,baseparams[j][i].antb,j);
					continue;
				}
			
				t		=	omp_get_thread_num();
			
				baseparams[j][i].data	=	baserawdata[t];		
		
				if(read_a_baseline (baseparams[j][i], fnames, j)) 
					printf("Baseline	%d %d		Scan %d		Reading Error...... !!!!!\n",baseparams[j][i].anta,baseparams[j][i].antb,j);
		
				//	----------------------	Actual flagging function	------------------------------------------------		
			
				printf("\nBaseline	%d %d		Scan %d		taken by thread	%d\n\n",baseparams[j][i].anta,baseparams[j][i].antb,j,t);
			
				flagbaseline (&baseparams[j][i], bl_minfrac, blockpow, flagrounds, flagparams, datarray[t], goutarr, buffarrc[t], buffarrr[t], coeffs[t], 
																													buffdev1d[t], buffdev2d[t], basestatsarray[j][i]);			
		
				//	----------------------	Flagging done	---------------------------------------------------
		
				if(write_a_baseline (baseparams[j][i], fnames, j)) 
					printf("Baseline	%d %d		Scan %d		Writing Error...... !!!!!\n",baseparams[j][i].anta,baseparams[j][i].antb,j);
		
				printf("\nBaseline	%d %d	Scan %d		Fractions	",baseparams[j][i].anta,baseparams[j][i].antb,j);					
				for (p=0; p<POLARRS; p++)
					printf("%.3f %.3f	",baseparams[j][i].flfrac[p],baseparams[j][i].flfrac_after[p]);
				printf("\n\n");
			}
	}
		
	//	--------------------------------------------------------------------------------------------------------
	//							Flag bad baselines
	//	-------------------------------------------------------------------------------------------------------
	
	if (flagmode==0) {
		for (j=0; j < nscans; j++) {
			for (i = 0; i < nbase; i++) {
				if (baseparams[j][i].exists==0)
					for (p=0; p<POLARRS; p++)
						basestatsarray[j][i][p][4]	=	1.0;
				//printf("%f	%f	%f	%f\n",basestatsarray[i][0][0],basestatsarray[i][1][0],basestatsarray[i][0][1],basestatsarray[i][1][1]);
			}
		}	
	
		nbadbase	=	findbadbase (baseparams, nscans, nbase, basestatsarray, bl_doflag, mean_tolerance, rms_tolerance, bl_minfrac, maxchansize*maxrecsize, goutarr, badbasefile);
	}
	
	//	-------------------------------------------------------------------------------------------------------
	
	free (basestatsarray);
	free (buffdev1d);
	free (buffdev2d);
	free (baserawdata);
	free (buffarrc);
	free (buffarrr);
	free (coeffs);
	free (datarray);
	free (goutarr);
	free (flagparams);
	
	free (baseparams);
	
	return (0);
}

//	-------------------------------------------------------------------------------------------------------
































































