#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>


//	---------------------------------------------------------------------------------------
//				Function to initialize the gout array
//															AB	3 August 2018															
//	---------------------------------------------------------------------------------------

int init_gout (float *goutarr)	{
	
	int i;
	FILE	*fp	=	fopen(GOUTFILE,"r");
	if (fp==NULL)	{
		printf("GOUTFILE not found !!!\n");
		return 1;
	}
	for (i=0; i<GOUTLEN; i++)
		fscanf(fp,"%f	",&goutarr[i]);
	fclose(fp);	
	
	return 0;
}

//	---------------------------------------------------------------------------------------






//	---------------------------------------------------------------------------------------
//				Function to initialize the data array
//															AB	3 August 2018															
//	---------------------------------------------------------------------------------------

int init_data (BaseParType *basepar, float **datarray, int i, int qtype)	{
	
	int		r,c;
	
	switch (qtype)	{
	
		case 0:
			for (r = 0; r < basepar->recsize; r++)
				for (c = 0; c < basepar->chansize; c++)	
					datarray[r][c]	=	basepar->data[i][0][r][c];
			//printf("\nWorking on real array \n\n");
			break;
		
		case 1:
			for (r = 0; r < basepar->recsize; r++)
				for (c = 0; c < basepar->chansize; c++)	
					datarray[r][c]	=	basepar->data[i][1][r][c];
			//printf("\nWorking on imaginary array \n\n");
			break;
			
		case 2:
			for (r = 0; r < basepar->recsize; r++)
				for (c = 0; c < basepar->chansize; c++)	
					datarray[r][c]	=	sqrt (basepar->data[i][1][r][c]*basepar->data[i][1][r][c] + basepar->data[i][0][r][c]*basepar->data[i][0][r][c]);
			//printf("\nWorking on amplitude array \n\n");
			break;	

		default:
			printf("\nUnknown recipe !!! Cannot prepare the requested item....\n");
			return 1;
	}
				
	return 0;
}

//	--------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to initialize baseline flagging parameters
//															AB	5 August 2018															
//	---------------------------------------------------------------------------------------

int	init_baseflag(FILE *fpflagpar, int *bl_doflag, float *mean_tolerance, float *rms_tolerance, float *bl_minfrac, float *blockpow, int *flagrounds, int *threads, int uvgrid[2]) {

	int		junkint, flagmode;
	
	fscanf(fpflagpar,"%d",&junkint);
	fscanf(fpflagpar,"%d",&flagmode);
	fscanf(fpflagpar,"%d",&uvgrid[0]);
	fscanf(fpflagpar,"%d",bl_doflag);
	fscanf(fpflagpar,"%d",&uvgrid[1]);
	fscanf(fpflagpar,"%f",mean_tolerance);
	fscanf(fpflagpar,"%f",rms_tolerance);
	fscanf(fpflagpar,"%d",flagrounds);
	fscanf(fpflagpar,"%d",threads);
	fscanf(fpflagpar,"%d",&junkint);
	fscanf(fpflagpar,"%f",blockpow);
	fscanf(fpflagpar,"%f",bl_minfrac);
	
	return flagmode;
}






//	---------------------------------------------------------------------------------------
//				Function to initialize flagging parameters
//															AB	5 August 2018															
//	---------------------------------------------------------------------------------------

int	init_flagpars (FILE *fpflagpar, int  curound, FlagParType *flagpar)	{
	
	fscanf(fpflagpar,"%d",&flagpar->whatflag);		
	fscanf(fpflagpar,"%d",&flagpar->doflag);
	fscanf(fpflagpar,"%d",&flagpar->domean);
	fscanf(fpflagpar,"%d",&flagpar->qtype);
	fscanf(fpflagpar,"%d",&flagpar->ascending);	
	fscanf(fpflagpar,"%f",&flagpar->tolerance);
	fscanf(fpflagpar,"%f",&flagpar->min_fraction);
	fscanf(fpflagpar,"%d",&flagpar->fitorder);	
	fscanf(fpflagpar,"%d",&flagpar->chanblockfac);
	fscanf(fpflagpar,"%d",&flagpar->recblockfac);
	fscanf(fpflagpar,"%f",&flagpar->chanmaxfrac);
	fscanf(fpflagpar,"%f",&flagpar->recmaxfrac);

	return 0;
}

























































































