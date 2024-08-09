#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>
#	include <gsl/gsl_poly.h>

//	---------------------------------------------------------------------------------------
//				Function to find bad channels and flag them
//															AB	28 July 2018															
//	---------------------------------------------------------------------------------------

int	flag_chan_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, 
							float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d) {

	int		i,j,r,c,templen,tempnon0,status;
	int		flagged	=	0;
	int		nonzero;
	float	*chmednon0[4];
	float	*temparrc, *meddev, *maddev, *medfit, *madfit;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
	
	for (i = 0; i < 4; i++)
		chmednon0[i]	=	buffarrc[i];			
	
	temparrc	=	buffarrr[0];
	meddev		=	buffarrc[5];
	maddev		=	buffarrc[6];
	medfit		=	buffarrc[7];
	madfit		=	buffarrc[8];	
	
	//	--------------------------------------------	Find channel statistics			
	
	tempnon0	=	0;
	
	for (c = 0; c < basepar.chansize; c++) {
		
		templen		=	0;
		
		for (r = 0; r < basepar.recsize; r++)			
			if ( flagarr[r][c] > 0.0 ) {				
				temparrc[ templen ] = datarray[r][c];
				templen++;
			}
		
		if ( templen > (int) (flagpar.min_fraction * basepar.recsize)) {	
			
			if (flagpar.domean) {
				chmednon0[0][tempnon0]	=	find_mean (temparrc, templen);
				chmednon0[1][tempnon0]	=	find_rms (temparrc, templen, chmednon0[0][tempnon0]);
			}
			else {
				chmednon0[0][tempnon0]	=	find_median (temparrc, templen, buffdev1d[2]);
				chmednon0[1][tempnon0]	=	find_mad (temparrc, templen, chmednon0[0][tempnon0], buffdev1d[2]);
			}
			chmednon0[2][tempnon0]	=	(float) c;
			chmednon0[3][tempnon0]	=	(float) templen;
			tempnon0++;
		
		}		
		else {
			for (r = 0; r < basepar.recsize; r++)
				flagarr[r][c] = -1.0 ;
		}	
	}	
	
	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
	//printf("Nonzero	=	%d	of	%d\n",nonzero,basepar.recsize * basepar.chansize);	
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		return (flagged);
		
	}	
	
	//	-----------------------------------------		Flag based on statistics
	
	outlim	=	goutarr[(int)round(log10(tempnon0)/0.001)];
	//printf("N = %d	outlim = %f\n",tempnon0,outlim);
	
	if (tempnon0) {
	
		if (flagpar.doflag % 2) {
		
			status		=		dopolyfit (chmednon0[2], chmednon0[0], chmednon0[1], coeffs, tempnon0, flagpar.fitorder);
			
			for (c = 0; c < basepar.chansize; c++)
				medfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
						
			for (i = 0; i < tempnon0; i++)
				meddev[i]	=	chmednon0[0][i] - medfit[(int) round(chmednon0[2][i])];					
			
			if (flagpar.domean)	{
				medval		=		find_mean (meddev, tempnon0);		
				madval		=		find_rms (meddev, tempnon0, medval);
			}
			else {
				medval		=		find_median (meddev, tempnon0, buffdev1d[2]);
				madval		=		find_mad (meddev, tempnon0, medval, buffdev1d[2]);
			}
		}
		//printf("Medval = %f Madval = %f	Medmad = %f Madmad = %f\n",medval,madval,medmad,madmad);		
		if (flagpar.doflag > 1){
		
			status		=		dopolyfit (chmednon0[2], chmednon0[1], chmednon0[1], coeffs, tempnon0, flagpar.fitorder);
			
			for (c = 0; c < basepar.chansize; c++)
				madfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
						
			for (i = 0; i < tempnon0; i++)
				maddev[i]	=	chmednon0[1][i] - madfit[(int) round(chmednon0[2][i])];
			
			if (flagpar.domean) {
				medmad		=		find_mean (maddev, tempnon0);
				madmad		=		find_rms (maddev, tempnon0, medmad);	
			}
			else {
				medmad		=		find_median (maddev, tempnon0, buffdev1d[2]);
				madmad		=		find_mad (maddev, tempnon0, medmad, buffdev1d[2]);	
			}		
		}
		//printf("Medval = %f Madval = %f	Medmad = %f Madmad = %f\n",medval,madval,medmad,madmad);				
		for (i = 0; i < tempnon0; i++) {
		
			c	=	(int) round(chmednon0[2][i]);
			
			if (chmednon0[3][i] > 0.0) {
				exfac	=	sqrt( (float) basepar.recsize / chmednon0[3][i] );

				if (flagpar.doflag % 2) 					
					if ( fabs(chmednon0[0][i]-medfit[c]-medval) > exfac*madval* outlim *flagpar.tolerance) {						
						for (r = 0; r < basepar.recsize; r++)
							flagarr[r][c] = -1.0 ;
						flagged++;
						continue;
					}
					
				if (flagpar.doflag > 1) 					
					if ( (chmednon0[1][i]-madfit[c]-medmad) > exfac*madmad* outlim *flagpar.tolerance) {						
						for (r = 0; r < basepar.recsize; r++)
							flagarr[r][c] = -1.0 ;
						flagged++;
					}				
			}
		}
	}
	//printf("Flagged	%d cahnnels\n",flagged);
	return (flagged);
}

//	---------------------------------------------------------------------------------------






//	---------------------------------------------------------------------------------------
//				Function to find bad records and flag them
//															AB	31 July 2018															
//	---------------------------------------------------------------------------------------

int	flag_rec_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, float **buffarrc, 
												float **buffarrr, float **buffdev1d) {

	int		i,j,r,c,templen,tempnon0;
	int		flagged	=	0;
	int		nonzero;
	float	*rmednon0[4];
	float	*temparr;
	float	medval, madval, medmad, madmad, exfac, outlim;

	//	------------------------------------------		Allocatig working arrays
	
	for (i = 0; i < 4; i++)
		rmednon0[i]	=	buffarrr[i];	
	
	temparr		=	buffarrc[0];	

	//	--------------------------------------------	Find record statistics			
	
	tempnon0	=	0;
	
	for (r = 0; r < basepar.recsize; r++){

		templen		=	0;
		
		for (c = 0; c < basepar.chansize; c++)			
			if ( flagarr[r][c] > 0.0 ){				
				temparr[ templen ] = datarray[r][c];
				templen++;
			}

		if( templen > (int) (flagpar.min_fraction * basepar.chansize)) {	

			if (flagpar.domean){
				rmednon0[0][tempnon0]	=	find_mean (temparr, templen);
				rmednon0[1][tempnon0]	=	find_rms (temparr, templen, rmednon0[0][tempnon0]);
			}
			else {
				rmednon0[0][tempnon0]	=	find_median (temparr, templen, buffdev1d[2]);
				rmednon0[1][tempnon0]	=	find_mad (temparr, templen, rmednon0[0][tempnon0], buffdev1d[2]);
			}
			
			rmednon0[2][tempnon0]	=	(float) r;
			rmednon0[3][tempnon0]	=	(float) templen;
			tempnon0++;
		}

		else{
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0 ;
		}	
	}

	//	---------------------------  	If the baseline is flagged more than critical fraction, flag it completely
	
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
	
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ) {

		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;

		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
		
		return (flagged);
	}
		
	//	-----------------------------------------		Flag based on statistics
	
	outlim	=	goutarr[(int)round(log10(tempnon0)/0.001)];
	//printf("N = %d	outlim = %f\n",tempnon0,outlim);
	
	if (tempnon0) {
		if (flagpar.doflag % 2) {
			if (flagpar.domean){
				medval		=		find_mean (rmednon0[0], tempnon0);
				madval		=		find_rms (rmednon0[0], tempnon0, medval);
			}
			else {
				medval		=		find_median (rmednon0[0], tempnon0, buffdev1d[2]);
				madval		=		find_mad (rmednon0[0], tempnon0, medval, buffdev1d[2]);
			}
		}	
		if (flagpar.doflag > 1){
			if (flagpar.domean){
				medmad		=		find_mean (rmednon0[1], tempnon0);
				madmad		=		find_rms (rmednon0[1], tempnon0, medmad);
			}
			else {
				medmad		=		find_median (rmednon0[1], tempnon0, buffdev1d[2]);
				madmad		=		find_mad (rmednon0[1], tempnon0, medmad, buffdev1d[2]);
			}
		}
	
		for (i = 0; i < tempnon0; i++) {
			
			r	=	(int) round(rmednon0[2][i]);
			
			if (rmednon0[3][i] > 0.0) {
				exfac	=	sqrt( (float) basepar.chansize / rmednon0[3][i] );	
		
				if (flagpar.doflag % 2) 					
					if ( fabs(rmednon0[0][i]-medval) > exfac*madval* outlim *flagpar.tolerance) {						
						for (c = 0; c < basepar.chansize; c++)
							flagarr[r][c] = -1.0 ;
						flagged++;
						continue;
					}		
				if (flagpar.doflag > 1) 					
					if ( (rmednon0[1][i]-medmad) > exfac*madmad* outlim *flagpar.tolerance) {						
						for (c = 0; c < basepar.chansize; c++)
							flagarr[r][c] = -1.0 ;
						flagged++;
					}				
			}
		}
	}

	return (flagged);
}

//	-------------------------------------------------------------------------------------------






//	---------------------------------------------------------------------------------------
//				Function to find bad visibilities and flag them
//															AB	31 July 2018															
//	---------------------------------------------------------------------------------------

int	flag_vis_individual (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, 
							float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d) {

	int		i,j,r,c,templen,tempnon0,vistotal,status;
	int		flagged	=	0;
	int		nonzero;
	float	*chmednon0[4];
	float	**devdata, *temparrc, *meddev, *medfit;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
	
	for (i = 0; i < 4; i++)
		chmednon0[i]	=	buffarrc[i];			
	
	temparrc	=	buffarrr[0];
	medfit		=	buffarrc[7];	
	meddev		=	buffdev1d[0];
	devdata		=	buffdev2d[0];
	
	//	--------------------------------------------	Find channel statistics			
	
	tempnon0	=	0;
	vistotal	=	0;
	
	for (c = 0; c < basepar.chansize; c++) {
		
		templen		=	0;
		
		for (r = 0; r < basepar.recsize; r++)			
			if ( flagarr[r][c] > 0.0 ) {				
				temparrc[ templen ] = datarray[r][c];
				templen++;
			}
		vistotal	+=	templen;
		
		if ( templen > (int) (flagpar.min_fraction * basepar.recsize)) {	
			
			if (flagpar.domean) {
				chmednon0[0][tempnon0]	=	find_mean (temparrc, templen);
				chmednon0[1][tempnon0]	=	find_rms (temparrc, templen, chmednon0[0][tempnon0]);
			}
			else {
				chmednon0[0][tempnon0]	=	find_median (temparrc, templen, buffdev1d[2]);
				chmednon0[1][tempnon0]	=	find_mad (temparrc, templen, chmednon0[0][tempnon0], buffdev1d[2]);
			}
			chmednon0[2][tempnon0]	=	(float) c;
			chmednon0[3][tempnon0]	=	(float) templen;
			tempnon0++;
		
		}		
		else {
			for (r = 0; r < basepar.recsize; r++)
				flagarr[r][c] = -1.0 ;
		}	
	}	
		
	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	//printf("Baseline	%d %d	Counting nonzero for %d\n",basepar.anta,basepar.antb,basepar.recsize * basepar.chansize);
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
		
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		return (flagged);
	}	
	
	//	-----------------------------------------		Flag based on statistics
	
	if (tempnon0) {
	
		if (flagpar.doflag % 2) {
		
			status		=		dopolyfit (chmednon0[2], chmednon0[0], chmednon0[1], coeffs, tempnon0, flagpar.fitorder);
			
			for (c = 0; c < basepar.chansize; c++)
				medfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );	
			
			j	=	0;
			for (r = 0; r < basepar.recsize; r++)
				for (c = 0; c < basepar.chansize; c++) {
					
					devdata[r][c]	=	datarray[r][c] - medfit[c];
					
					if (flagarr[r][c] > 0.0) {
						meddev[j]	=	devdata[r][c];
						j++;
					}
				}	
			
			outlim	=	goutarr[(int)round(log10(j)/0.001)];
			//printf("N = %d	outlim = %f\n",j,outlim);
				
			if (flagpar.domean)	{		
				medval		=		find_mean (meddev, j);				
				madval		=		find_rms (meddev, j, medval);
			}
			else {
				medval		=		find_median (meddev, j, buffdev1d[2]);			
				madval		=		find_mad (meddev, j, medval, buffdev1d[2]);
			}
		}			
			
		else
			printf("Hmm... You want RMS from a single point... Try again when the effect of alcohol wears off...\n");			
			
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)	
				if(flagarr[r][c] > 0.0)
					if( fabs(devdata[r][c]-medfit[c]-medval) > madval* outlim *flagpar.tolerance){
						flagarr[r][c] 	= 	-1.0;
						flagged++;	
					}	
	}
	
	return (flagged);
}

//	---------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to smooth a 2d array with a tophat kernel
//															AB	1 August 2018															
//	---------------------------------------------------------------------------------------

int	smooth_2d (float **inarr, float **flagarr, int a, int b, float **outarr, float **rmsarr, float **junkarr, BaseParType basepar, FlagParType flagpar,
						double *coeffs,	float *arr0, float *arr1, float *arr2, float *arr3, float *temparrc, float *medfit, float **copyarr, float *temparr2d, float *tempflag) {
	
	int		i,j,r,c,templen,tempnon0,status,nonzero;
	float	*chmednon0[4];	
	float	tempmean;
	
	//	------------------------------------------		Allocatig working arrays
	
	chmednon0[0]=	arr0;
	chmednon0[1]=	arr1;
	chmednon0[2]=	arr2;
	chmednon0[3]=	arr3;				
		
	//	--------------------------------------------	Find channel statistics	
	
	tempnon0	=	0;
	
	for (c = 0; c < basepar.chansize; c++) {
		
		templen		=	0;
		
		for (r = 0; r < basepar.recsize; r++)			
			if ( flagarr[r][c] > 0.0 ) {				
				temparrc[ templen ] = inarr[r][c];
				templen++;
			}	
	
		if ( templen > (int) (flagpar.min_fraction * basepar.recsize)) {	
			chmednon0[0][tempnon0]	=	find_mean (temparrc, templen);
			chmednon0[1][tempnon0]	=	find_rms (temparrc, templen, chmednon0[0][tempnon0]);	
			chmednon0[2][tempnon0]	=	(float) c;
			chmednon0[3][tempnon0]	=	(float) templen;
			tempnon0++;
		}		
		else {
			for (r = 0; r < basepar.recsize; r++)
				flagarr[r][c] = -1.0 ;
		}	
	}
	
	if ( tempnon0 <= 0)		
		return 0;	
	
	status		=		dopolyfit (chmednon0[2], chmednon0[0], chmednon0[1], coeffs, tempnon0, flagpar.fitorder);
			
	for (c = 0; c < basepar.chansize; c++)
		medfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
	
	for (r = 0; r < basepar.recsize; r++)
		for (c = 0; c < basepar.chansize; c++) {
			if (flagarr[r][c] > 0.0)
				copyarr[r][c]	=	inarr[r][c];
			else
				copyarr[r][c]	=	medfit[c];
		}
	
	for (r = 0; r < basepar.recsize - a + 1; r++) 	
		for (c = 0; c < basepar.chansize - b + 1; c++) {			
			for (i = 0; i < a; i++) 
				for (j = 0; j < b; j++) {
					temparr2d[i*b+j]	=	copyarr[r+i][c+j];
					tempflag[i*b+j]		=	flagarr[r+i][c+j];
				}
			
			tempmean		=	find_mean(temparr2d,a*b);	
			rmsarr[r][c]	=	find_rms (temparr2d,a*b,tempmean);
			outarr[r][c]	=	tempmean;
			
			nonzero		=	count_nonzero (tempflag, a*b);
			//nonzero			=	count_nonzero_2d (tempflag,a,b);
			if ( nonzero < flagpar.min_fraction*a*b)
				junkarr[r][c] = 0.0;	
			else
				junkarr[r][c] = (float) nonzero;
	
		}	
		
	return 1;
}

//	---------------------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to find bad blocks of channels and flag them
//															AB	2 August 2018															
//	---------------------------------------------------------------------------------------

int	flag_chan_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int cw, float blockpow, 
										float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d) {

	int		i,j,r,c,templen,tempnon0,status,status0;
	int		flagged	=	0;
	int		nonzero;
	float	*chnums, *meddev, *maddev, *medfit, *madfit, **smootharr, **smoothrms, **junkarr;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
	
	chnums		=	buffarrc[0];
	for (c = 0; c < basepar.chansize - cw + 1; c++) 			
		chnums[c]		=	(float) c;		
	
	meddev		=	buffarrc[1];
	maddev		=	buffarrc[2];
	medfit		=	buffarrc[3];
	madfit		=	buffarrc[4];	
	smootharr	=	buffdev2d[0];
	smoothrms	=	buffdev2d[1];
	junkarr		=	buffdev2d[2]; 
		
	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
		
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		return (flagged);
	}	
	
	int		localen;
	float	*localarr;
			
	localarr=	buffdev1d[2];	
	
	outlim	=	goutarr[(int)round(log10((float) (basepar.chansize)/pow(cw,blockpow))/0.001)];
	//printf("outlim = %f\n",outlim);
	
	status	=	smooth_2d (datarray, flagarr, basepar.recsize, cw, smootharr, smoothrms, junkarr, basepar, flagpar, coeffs,
											buffarrc[7], buffarrc[8], buffarrc[9], buffarrc[10], buffarrr[1], buffarrc[11], buffdev2d[3], buffdev1d[0], buffdev1d[1]);

	if (status) {
		
		if (flagpar.doflag % 2) {	
			
			status0		=		dopolyfit (chnums, smootharr[0], smoothrms[0], coeffs, (basepar.chansize - cw + 1), flagpar.fitorder);		
			
			for (c = 0; c < basepar.chansize; c++)
				medfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
			
			for (c = 0; c < basepar.chansize - cw + 1; c++)
				meddev[c]	=	smootharr[0][c] - medfit[c];
			
			/*FILE	*fptemp;
			fptemp	=	fopen("fittest.dat","w");
			for (c = 0; c < basepar.chansize - cw + 1; c++)
				fprintf(fptemp,"%f %f %f\n",(float) c, smootharr[0][c], medfit[c]);			
			fclose(fptemp);*/
			
			localen	=	0;	
			for (c = 0; c < basepar.chansize - cw + 1; c++)	{
				if (junkarr[0][c] > 0.0)	{
					localarr[localen]	=	meddev[c];
					localen++;
				}	
			}
			
			if (flagpar.domean)	{			
				if (localen){
					medval		=		find_mean (localarr, localen);	
					madval		=		find_rms (localarr, localen, medval);
				}
				else {	
					medval		=		find_mean (meddev, basepar.chansize - cw + 1);	
					madval		=		find_rms (meddev, basepar.chansize - cw + 1, medval);
				}
			}
			else	{							
				if (localen)	{
					medval		=		find_median (localarr, localen, buffdev1d[3]);	
					madval		=		find_mad (localarr, localen, medval, buffdev1d[3]);
				}
				else {	
					medval		=		find_median (meddev, basepar.chansize - cw + 1, buffdev1d[3]);	
					madval		=		find_mad (meddev, basepar.chansize - cw + 1, medval, buffdev1d[3]);
				}
			}
		}
		
		if (flagpar.doflag > 1) {	
		
			status0		=		dopolyfit (chnums, smoothrms[0], smoothrms[0], coeffs, (basepar.chansize - cw + 1), flagpar.fitorder);
			
			for (c = 0; c < basepar.chansize; c++)
				madfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
			
			for (c = 0; c < basepar.chansize - cw + 1; c++)
				maddev[c]	=	smoothrms[0][c] - madfit[c];
			
			localen	=	0;			
			for (c = 0; c < basepar.chansize - cw + 1; c++)	{
				if (junkarr[0][c] > 0.0)	{
					localarr[localen]	=	maddev[c];
					localen++;
				}	
			}
			
			if (flagpar.domean)	{							
				if (localen) {
					medmad		=		find_mean (localarr, localen);	
					madmad		=		find_rms (localarr, localen, medmad);
				}
				else {	
					medmad		=		find_mean (maddev, basepar.chansize - cw + 1);	
					madmad		=		find_rms (maddev, basepar.chansize - cw + 1, medmad);
				}
			}
			else	{			
				if (localen)	{
					medmad		=		find_median (localarr, localen, buffdev1d[3]);	
					madmad		=		find_mad (localarr, localen, medmad, buffdev1d[3]);
				}
				else {	
					medmad		=		find_median (maddev, basepar.chansize - cw + 1, buffdev1d[3]);	
					madmad		=		find_mad (maddev, basepar.chansize - cw + 1, medmad, buffdev1d[3]);
				}
			}	
		}
	}
	
	//printf("Mean	%f	RMS		%f\n",medval,madval);
	
	if (status)	{
		for (c = 0; c < basepar.chansize - cw + 1; c++) {
			if (junkarr[0][c]<=0.0)
				continue;	
			exfac	=	sqrt ( cw * basepar.recsize/junkarr[0][c] );			
			
			if (flagpar.doflag % 2) 					
				if ( (meddev[c]-medval) > exfac*madval* outlim *flagpar.tolerance) {						
					for (r = 0; r < basepar.recsize; r++)
						for (j = 0; j < cw; j++)
							flagarr[r][c+j] = -1.0 ;
					flagged++;
					continue;
				}
				
			if (flagpar.doflag > 1) 					
				if ( (maddev[c]-medmad) > exfac*madmad* outlim *flagpar.tolerance) {						
					for (r = 0; r < basepar.recsize; r++)
						for (j = 0; j < cw; j++)
							flagarr[r][c+j] = -1.0 ;
					flagged++;
				}
		}
	}
	
	return (flagged);
}

//	---------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to find bad blocks of records and flag them
//															AB	3 August 2018															
//	---------------------------------------------------------------------------------------

int	flag_rec_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int rw, float blockpow,
									float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d) {

	int		i,j,r,c,templen,tempnon0,status;
	int		flagged	=	0;
	int		nonzero;
	float	**smootharr, **smoothrms, **junkarr;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
		
	smootharr	=	buffdev2d[0];
	smoothrms	=	buffdev2d[1];
	junkarr		=	buffdev2d[2];
			
	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
		
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		return (flagged);
	}	
	
	int		localen;
	float	*localarr;
			
	localarr=	buffdev1d[2];
	
	outlim	=	goutarr[(int)round(log10((float) (basepar.recsize)/pow(rw,blockpow))/0.001)];
	//printf("outlim = %f\n",outlim);
	
	status	=	smooth_2d (datarray, flagarr, rw, basepar.chansize, smootharr, smoothrms, junkarr, basepar, flagpar, coeffs,
											buffarrc[7], buffarrc[8], buffarrc[9], buffarrc[10], buffarrr[1], buffarrc[11], buffdev2d[3], buffdev1d[0], buffdev1d[1]);
	
	/*
	FILE	*fptemp;
	fptemp	=	fopen("fittest.dat","w");
	for (r = 0; r < basepar.recsize - rw + 1; r++)	
		fprintf(fptemp,"%f %f\n",(float) r, smootharr[r][0]);			
	fclose(fptemp);*/
		
	if (status) {
		
		if (flagpar.doflag % 2) {	
			
			localen	=	0;	
			for (r = 0; r < basepar.recsize - rw + 1; r++)	{
				if (junkarr[r][0] > 0.0)	{
					localarr[localen]	=	smootharr[r][0];
					localen++;
				}	
			}
			
			if (flagpar.domean)	{			
				if (localen){
					medval		=		find_mean (localarr, localen);	
					madval		=		find_rms (localarr, localen, medval);
				}
				else {	
					medval		=		find_mean ((float *)smootharr, basepar.recsize - rw + 1);	
					madval		=		find_rms ((float *)smootharr, basepar.recsize - rw + 1, medval);
				}
			}
			else	{							
				if (localen)	{
					medval		=		find_median (localarr, localen, buffdev1d[3]);	
					madval		=		find_mad (localarr, localen, medval, buffdev1d[3]);
				}
				else {	
					medval		=		find_median ((float *)smootharr, basepar.recsize - rw + 1, buffdev1d[3]);	
					madval		=		find_mad ((float *)smootharr, basepar.recsize - rw + 1, medval, buffdev1d[3]);
				}
			}
		}
		
		if (flagpar.doflag > 1) {	
			
			localen	=	0;		
			for (r = 0; r < basepar.recsize - rw + 1; r++)	{
				if (junkarr[r][0] > 0.0)	{
					localarr[localen]	=	smoothrms[r][0];
					localen++;
				}	
			}
			
			if (flagpar.domean)	{							
				if (localen) {
					medmad		=		find_mean (localarr, localen);	
					madmad		=		find_rms (localarr, localen, medmad);
				}
				else {	
					medmad		=		find_mean ((float *)smoothrms, basepar.recsize - rw + 1);	
					madmad		=		find_rms ((float *)smoothrms, basepar.recsize - rw + 1, medmad);
				}
			}
			else	{			
				if (localen)	{
					medmad		=		find_median (localarr, localen, buffdev1d[3]);	
					madmad		=		find_mad (localarr, localen, medmad, buffdev1d[3]);
				}
				else {	
					medmad		=		find_median ((float *)smoothrms, basepar.recsize - rw + 1, buffdev1d[3]);	
					madmad		=		find_mad ((float *)smoothrms, basepar.recsize - rw + 1, medmad, buffdev1d[3]);
				}
			}	
		}
	}
	
	//printf("Mean	%.3f	RMS	%.3f	Lim	%.3f\n",medval,madval,madval* outlim *flagpar.tolerance);
	
	if (status)	{
		for (r = 0; r < basepar.recsize - rw + 1; r++) {
			if (junkarr[r][0]<=0.0)
				continue;	
			exfac	=	sqrt (((float) (rw * basepar.chansize))/junkarr[r][0] );			
			
			if (flagpar.doflag % 2) 					
				if ( (smootharr[r][0]-medval) > exfac*madval* outlim *flagpar.tolerance) {					
					for (c = 0; c < basepar.chansize; c++)
						for (j = 0; j < rw; j++)
							flagarr[r+j][c] = -1.0 ;
					flagged++;
					continue;
				}
				
			if (flagpar.doflag > 1) 					
				if ( (smoothrms[r][0]-medmad) > exfac*madmad* outlim *flagpar.tolerance) {						
					for (c = 0; c < basepar.chansize; c++)
						for (j = 0; j < rw; j++)
							flagarr[r+j][c] = -1.0 ;
					flagged++;
				}
		}
	}
	
	return (flagged);
}

//	---------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to find bad blocks of visibilities and flag them
//															AB	3 August 2018															
//	---------------------------------------------------------------------------------------

int	flag_vis_block (FlagParType flagpar, BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *goutarr, int rw, int cw, float blockpow,
									float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d) {

	int		i,j,r,c,templen,tempnon0,status,status0,localen;
	int		flagged	=	0;
	int		nonzero;
	float	*chmeds[2];
	float	*chnums, *temparr, *temprms, **meddev, **maddev, *medfit, *madfit, **smootharr, **smoothrms, **junkarr, *localarr;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
	
	chnums		=	buffarrc[0];
	for (c = 0; c < basepar.chansize - cw + 1; c++) 			
		chnums[c]		=	(float) c;		
	
	chmeds[0]	=	buffarrc[1];
	chmeds[1]	=	buffarrc[2];
	medfit		=	buffarrc[3];
	madfit		=	buffarrc[4];
	
	temparr		=	buffarrr[0];
	temprms		=	buffarrr[1];
	
	localarr	=	buffdev1d[2];	
		
	smootharr	=	buffdev2d[0];
	smoothrms	=	buffdev2d[1];
	junkarr		=	buffdev2d[2];	
	meddev		=	buffdev2d[3];
	maddev		=	buffdev2d[4];
	
	status	=	smooth_2d (datarray, flagarr, rw, cw, smootharr, smoothrms, junkarr, basepar, flagpar, coeffs,
											buffarrc[7], buffarrc[8], buffarrc[9], buffarrc[10], buffarrr[2], buffarrc[11], buffdev2d[5], buffdev1d[0], buffdev1d[1]);
		
	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	
	//nonzero	=	count_nonzero ((float *)flagarr, basepar.recsize * basepar.chansize);
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
		
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		return (flagged);
	}		
	
	outlim	=	goutarr[(int)round(log10((float) (basepar.chansize*basepar.recsize)/pow(rw*cw,blockpow))/0.001)];
	//printf("outlim = %f\n",outlim);
		
	for (c = 0; c < basepar.chansize - cw + 1; c++) {		
		for (r = 0; r < basepar.recsize - rw + 1; r++){							
			temparr[ r ] = smootharr[r][c];
			temprms[ r ] = smoothrms[r][c];
		}		
		if (flagpar.domean) {	
			chmeds[0][c]	=	find_mean (temparr, basepar.recsize - rw + 1);
			chmeds[1][c]	=	find_mean (temprms, basepar.recsize - rw + 1);
		}
		else {	
			chmeds[0][c]	=	find_median (temparr, basepar.recsize - rw + 1, buffdev1d[3]);
			chmeds[1][c]	=	find_median (temprms, basepar.recsize - rw + 1, buffdev1d[3]);
		}
	}		
	
	if (status) {
		
		if (flagpar.doflag % 2) {	
			
			status0		=		dopolyfit (chnums, chmeds[0], chmeds[1], coeffs, (basepar.chansize - cw + 1), flagpar.fitorder);		
			
			for (c = 0; c < basepar.chansize - cw + 1; c++)
				medfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
			
			for (r = 0; r < basepar.recsize - rw + 1; r++)
				for (c = 0; c < basepar.chansize - cw + 1; c++)
					meddev[r][c]	=	smootharr[r][c] - medfit[c];
			
			/*FILE	*fptemp;
			fptemp	=	fopen("fittest.dat","w");
			for (r = 0; r < basepar.recsize - rw + 1; r++){
				for (c = 0; c < basepar.chansize - cw + 1; c++)
					fprintf(fptemp,"%f ",meddev[r][c]);	
				fprintf(fptemp,"\n");
			}		
			fclose(fptemp);*/
			
			localen	=	0;	
			for (r = 0; r < basepar.recsize - rw + 1; r++)
				for (c = 0; c < basepar.chansize - cw + 1; c++)	{
					if (junkarr[r][c] > 0.0)	{
						localarr[localen]	=	meddev[r][c];
						localen++;
					}	
				}
			
			if (flagpar.domean)	{			
				if (localen){
					medval		=		find_mean (localarr, localen);	
					madval		=		find_rms (localarr, localen, medval);
				}
				else {	
					medval		=		find_mean ((float *) meddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1));	
					madval		=		find_rms ((float *) meddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), medval);
				}
			}
			else	{							
				if (localen)	{
					medval		=		find_median (localarr, localen, buffdev1d[3]);	
					madval		=		find_mad (localarr, localen, medval, buffdev1d[3]);
				}
				else {	
					medval		=		find_median ((float *) meddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), buffdev1d[3]);	
					madval		=		find_mad ((float *) meddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), medval, buffdev1d[3]);
				}
			}
		}
		
		if (flagpar.doflag > 1) {	
		
			status0		=		dopolyfit (chnums, chmeds[1], chmeds[1], coeffs, (basepar.chansize - cw + 1), flagpar.fitorder);		
			
			for (c = 0; c < basepar.chansize - cw + 1; c++)
				madfit[c]	=	(float) gsl_poly_eval ( coeffs, flagpar.fitorder, (double) c );
			
			for (r = 0; r < basepar.recsize - rw + 1; r++)
				for (c = 0; c < basepar.chansize - cw + 1; c++)
					maddev[r][c]	=	smoothrms[r][c] - madfit[c];
			
			localen	=	0;	
			for (r = 0; r < basepar.recsize - rw + 1; r++)
				for (c = 0; c < basepar.chansize - cw + 1; c++)	{
					if (junkarr[r][c] > 0.0)	{
						localarr[localen]	=	maddev[r][c];
						localen++;
					}	
				}
			
			if (flagpar.domean)	{			
				if (localen){
					medmad		=		find_mean (localarr, localen);	
					madmad		=		find_rms (localarr, localen, medmad);
				}
				else {	
					medmad		=		find_mean ((float *) maddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1));	
					madmad		=		find_rms ((float *) maddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), medmad);
				}
			}
			else	{							
				if (localen)	{
					medmad		=		find_median (localarr, localen, buffdev1d[3]);	
					madmad		=		find_mad (localarr, localen, medmad, buffdev1d[3]);
				}
				else {	
					medmad		=		find_median ((float *) maddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), buffdev1d[3]);	
					madmad		=		find_mad ((float *) maddev, (basepar.chansize - cw + 1)*(basepar.recsize - rw + 1), medmad, buffdev1d[3]);
				}
			}	
		}
	}
		
	//printf("Mean	%f	RMS		%f\n",medval,madval);
	
	if (status)	{
		for (r = 0; r < basepar.recsize - rw + 1; r++)
			for (c = 0; c < basepar.chansize - cw + 1; c++) {
				if (junkarr[r][c]<=0.0)
					continue;	
				exfac	=	sqrt ( cw * rw/junkarr[r][c] );			
				
				if (flagpar.doflag % 2) 					
					if ( (meddev[r][c]-medval) > exfac*madval* outlim *flagpar.tolerance) {						
						for (i = 0; i < rw; i++)
							for (j = 0; j < cw; j++)
								flagarr[r+i][c+j] = -1.0 ;
						flagged++;
						continue;
					}
				
				if (flagpar.doflag > 1) 					
					if ( (maddev[r][c]-medmad) > exfac*madmad* outlim *flagpar.tolerance) {						
						for (i = 0; i < rw; i++)
							for (j = 0; j < cw; j++)
								flagarr[r+i][c+j] = -1.0 ;
						flagged++;
					}
			}
	}
	
	return (flagged);
}

//	---------------------------------------------------------------------------------------






//	---------------------------------------------------------------------------------------
//				Function to find statistics for a given baseline
//															AB	8 August 2018															
//	---------------------------------------------------------------------------------------

float findbasestat (BaseParType basepar, float **datarray, float **flagarr, float bl_minfrac, float *statsarray, 
										float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d) {

	int		i,j,r,c,templen,tempnon0,vistotal,statusmean,statusmed;
	int		flagged	=	0;
	int		nonzero;
	float	*chmednon0[6];
	float	**devdata, *temparrc, *meddev, *medfit;
	float	medval, madval, medmad, madmad, exfac, outlim;	
	
	//	------------------------------------------		Allocatig working arrays
	
	for (i = 0; i < 6; i++)
		chmednon0[i]	=	buffarrc[i];			
	
	temparrc	=	buffarrr[0];
	medfit		=	buffarrc[6];	
	meddev		=	buffdev1d[0];
	devdata		=	buffdev2d[0];

	//	--------------------------------------------	Find channel statistics			
	
	tempnon0	=	0;
	vistotal	=	0;
	
	for (c = 0; c < basepar.chansize; c++) {
		
		templen		=	0;
		
		for (r = 0; r < basepar.recsize; r++)			
			if ( flagarr[r][c] > 0.0 ) {				
				temparrc[ templen ] = datarray[r][c];
				templen++;
			}
		vistotal	+=	templen;
		if (templen) {
			chmednon0[0][tempnon0]	=	find_mean (temparrc, templen);
			chmednon0[1][tempnon0]	=	find_rms (temparrc, templen, chmednon0[0][tempnon0]);
		
			chmednon0[4][tempnon0]	=	find_median (temparrc, templen, buffdev1d[2]);
			chmednon0[5][tempnon0]	=	find_mad (temparrc, templen, chmednon0[0][tempnon0], buffdev1d[2]);
			
			chmednon0[2][tempnon0]	=	(float) c;
			chmednon0[3][tempnon0]	=	(float) templen;
			tempnon0++;	
		}					
	}
	
	//	------------------------------------------		Find statistics

	if (tempnon0) {
		
		statusmean		=		dopolyfit (chmednon0[2], chmednon0[0], chmednon0[1], coeffs, tempnon0, 3);
		for (c = 0; c < basepar.chansize; c++)
			medfit[c]	=	(float) gsl_poly_eval ( coeffs, 3, (double) c );
		
		j	=	0;
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++) {
				
				devdata[r][c]	=	datarray[r][c] - medfit[c];
				
				if (flagarr[r][c] > 0.0) {
					meddev[j]	=	devdata[r][c];
					j++;
				}
			}
		
		medval		=		find_mean (meddev, j);				
		madval		=		find_rms (meddev, j, medval);	
			
		statsarray[0] 	=	medval;
		statsarray[2]	=	madval;		
		
		statusmed		=		dopolyfit (chmednon0[2], chmednon0[4], chmednon0[5], coeffs, tempnon0, 3);
		for (c = 0; c < basepar.chansize; c++)
			medfit[c]	=	(float) gsl_poly_eval ( coeffs, 3, (double) c );
		
		j	=	0;
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++) {
				
				devdata[r][c]	=	datarray[r][c] - medfit[c];
				
				if (flagarr[r][c] > 0.0) {
					meddev[j]	=	devdata[r][c];
					j++;
				}
			}
		
		medval		=		find_median (meddev, j, buffdev1d[2]);			
		madval		=		find_mad (meddev, j, medval, buffdev1d[2]);	
		
		statsarray[1] 	=	medval;
		statsarray[3]	=	madval;			
	}

	//	------------------------------------------   	If the baseline is flagged more than critical fraction, flag it completely
	
	nonzero	=	count_nonzero_2d (flagarr, basepar.recsize, basepar.chansize);
		
	if (nonzero < (int) (bl_minfrac * basepar.recsize * basepar.chansize) ){
		
		for (r = 0; r < basepar.recsize; r++)
			for (c = 0; c < basepar.chansize; c++)
				flagarr[r][c] = -1.0;	
		
		printf("Baseline	%d %d		Flagged off......\n",basepar.anta,basepar.antb);
				
		statsarray[4] 	= 	1.0;
	}
	else
		statsarray[4]	=	(1.0 - ((float) nonzero)/(basepar.recsize * basepar.chansize) );
	
	return (statsarray[4]);
}

//	---------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Function to find bad baselines
//															AB	8 August 2018	
//															AB	26	November 2019														
//	---------------------------------------------------------------------------------------

int findbadbase (BaseParType **baseparams, int nscan, int nbase, float ****basestatsarray, int bl_doflag, float mean_tolerance, float rms_tolerance, 
														float bl_minfrac, int maxize, float *goutarr, FILE *badbasefile) {

	int			nbad	=	0;
	int			p,i,j,k,l,m;
	float		meanval, medval, rmsval, madval, outlim;
	
	int			***baselines;
	float		**temparr;	
	temparr		=	(float **) malloc (5 * sizeof(float *));
	for (i = 0; i < 5; i++)
		temparr[i]		=	(float *) malloc (nscan*nbase * sizeof(float));	
	
	baselines	=	(int ***) malloc (nscan * sizeof(int **));
	for (l = 0; l < nscan; l++) {
		baselines[l]	=	(int **) malloc (nbase * sizeof(int *));
		for (i = 0; i < nbase; i++) {
			baselines[l][i]	=	(int *) malloc (POLARRS * sizeof(int));	
			for (p = 0; p < POLARRS; p++)
				baselines[l][i][p]	=	1;	
		}
	}
	
	for (p = 0; p < POLARRS; p++) {
	
		j	=	0;
		for (l = 0; l < nscan; l++) {
			for (i = 0; i < nbase; i++) {
				if (basestatsarray[l][i][p][4] < 1.0) {
					for (m=0; m < 4; m++)
						temparr[m][j]	=	basestatsarray[l][i][p][m];
					j++;
				}
				else
					baselines[l][i][p]	=	0;	
			}
		}
		
		if (j==0)
			j	=	1;
		
		outlim	=	goutarr[(int)round(log10((float) j)/0.001)];
		
		printf("Total baseline scans %d\n",j);
		
		if (bl_doflag%2) {
			//	------	Outliers based on Median
		
			meanval	=	find_median(temparr[1],j,temparr[4]);
			rmsval	=	find_mad(temparr[1],j,meanval,temparr[4]);
			//printf("Baseline median	= %f	MAD	= %f\n",meanval,rmsval);
			for (l = 0; l < nscan; l++) {
				for (i = 0; i < nbase; i++) {
			
					if (basestatsarray[l][i][p][1] > (meanval + outlim*rmsval*mean_tolerance*sqrt(maxize/((float) (baseparams[l][i].chansize*baseparams[l][i].recsize))))) {
						printf ("Baseline	%d %d	Scan %d	Pol	%d	median flagged...\n",baseparams[l][i].anta,baseparams[l][i].antb,l,p);
						baselines[l][i][p]	=	0;
					}			
				}
			}
		
			//	------	Outliers based on Mean
		
			meanval	=	find_mean(temparr[0],j);
			rmsval	=	find_rms(temparr[0],j,meanval);
			//printf("Baseline mean	= %f	RMS	= %f\n",meanval,rmsval);
			for (l = 0; l < nscan; l++) {
				for (i = 0; i < nbase; i++) {
			
					if (basestatsarray[l][i][p][0] > (meanval + outlim*rmsval*mean_tolerance*sqrt(maxize/((float) (baseparams[l][i].chansize*baseparams[l][i].recsize))))) {
						printf ("Baseline	%d %d	Scan %d	Pol	%d	mean flagged...\n",baseparams[l][i].anta,baseparams[l][i].antb,l,p);
						baselines[l][i][p]	=	0;
					}			
				}
			}
		}
		
		if (bl_doflag>1) {
			//	------	Outliers based on MAD
		
			meanval	=	find_median(temparr[3],j,temparr[4]);
			rmsval	=	find_mad(temparr[3],j,meanval,temparr[4]);
			//printf("Baseline median	= %f	MAD	= %f\n",meanval,rmsval);
			for (l = 0; l < nscan; l++) {
				for (i = 0; i < nbase; i++) {
			
					if (basestatsarray[l][i][p][3] > (meanval + outlim*rmsval*rms_tolerance*sqrt(maxize/((float) (baseparams[l][i].chansize*baseparams[l][i].recsize))))) {
						printf ("Baseline	%d %d	Scan %d	Pol	%d	MAD flagged...\n",baseparams[l][i].anta,baseparams[l][i].antb,l,p);
						baselines[l][i][p]	=	0;
					}			
				}
			}
		
			//	------	Outliers based on RMS
		
			meanval	=	find_mean(temparr[2],j);
			rmsval	=	find_rms(temparr[2],j,meanval);
			//printf("Baseline mean	= %f	RMS	= %f\n",meanval,rmsval);
			for (l = 0; l < nscan; l++) {
				for (i = 0; i < nbase; i++) {
			
					if (basestatsarray[l][i][p][2] > (meanval + outlim*rmsval*rms_tolerance*sqrt(maxize/((float) (baseparams[l][i].chansize*baseparams[l][i].recsize))))) {
						printf ("Baseline	%d %d	Scan %d	Pol	%d	RMS flagged...\n",baseparams[l][i].anta,baseparams[l][i].antb,l,p);
						baselines[l][i][p]	=	0;
					}			
				}
			}
		}			
	}
	
	for (l = 0; l < nscan; l++) {
		for (i = 0; i < nbase; i++)
				fprintf(badbasefile,"%d	%d	%d	%d	%d	%d\n",baseparams[l][i].anta,baseparams[l][i].antb,baseparams[l][i].exists,baselines[l][i][0],baselines[l][i][1],l);
	}
	
	free (baselines);
	free (temparr);
	return	(nbad);
}

























































