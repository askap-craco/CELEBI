#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>
#	include <gsl/gsl_multifit.h>


//	---------------------------------------------------------------------------------------
//				Function to count nonzero elements of an array
//															AB	27 July 2018															
//	---------------------------------------------------------------------------------------

int	count_nonzero (float *data, int nels) {

	int		i;
	int		nonzerocount = 0;

	for ( i = 0; i < nels; i++)
		if( data[i] > 0.0)
				nonzerocount++;
	
	return (nonzerocount);
}

//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to count nonzero elements of a 2D array
//															AB	6 August 2018															
//	---------------------------------------------------------------------------------------

int	count_nonzero_2d (float **data, int n1, int n2) {

	int		i,j;
	int		nonzerocount = 0;

	for ( i = 0; i < n1; i++)
		for ( j = 0; j < n2; j++)
			if( data[i][j] > 0.0)
					nonzerocount++;
	
	return (nonzerocount);
}

//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to find mean of an array
//															AB	28 July 2018															
//	---------------------------------------------------------------------------------------

float	find_mean (float *data, int nels){

	int		i;
	float	mean = 0.0;

	for ( i = 0; i < nels; i++)
		mean += data[i];
	
	mean	=	mean/nels;
	
	return (mean);
}

//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to find RMS of an array
//															AB	28 July 2018															
//	---------------------------------------------------------------------------------------

float	find_rms (float *data, int nels, float mean){

	int		i;
	float	rms = 0.0;

	for ( i = 0; i < nels; i++)
		rms += (data[i] - mean) * (data[i] - mean);
	
	rms	=	sqrt(rms/nels);
	
	return (rms);
}

//	---------------------------------------------------------------------------------------





//	---------------------------------------------------------------------------------------
//				Compare function
int cmpare (const void * a, const void * b) {
	float	fa	=	*(const float *)a;
	float	fb	=	*(const float *)b;
				
   	return ((int)copysign(1.0, (fa-fb)));
}
//	---------------------------------------------------------------------------------------


//	---------------------------------------------------------------------------------------
//				Function to find median of an array
//															AB	28 July 2018															
//	---------------------------------------------------------------------------------------

float	find_median (float *data, int nels, float *temparr){
	
	float	median;	
	
	memcpy (temparr, data, nels * sizeof(float));	
	qsort (temparr, nels, sizeof(float), cmpare);
	
	if (nels%2)
		median	=	temparr[nels/2];
	else
		median	=	0.5*(temparr[(nels/2)-1] + temparr [nels/2]);

	return (median);
}
//	---------------------------------------------------------------------------------------



//	---------------------------------------------------------------------------------------
//				Function to find MAD of an array
//															AB	28 July 2018															
//	---------------------------------------------------------------------------------------

float	find_mad (float *data, int nels, float median, float *temparr){
	
	int		i;
	float 	mad;
		
	memcpy (temparr, data, nels * sizeof(float));
	for (i=0;i<nels;i++)
		temparr[i] = fabs(temparr[i] - median);
		
	qsort (temparr, nels, sizeof(float), cmpare);
	
	if (nels%2)
		mad	=	temparr[nels/2];
	else
		mad	=	0.5*(temparr[(nels/2)-1] + temparr [nels/2]);

	return (1.48*mad);
}
//	---------------------------------------------------------------------------------------




//	---------------------------------------------------------------------------------------
//				Function to fit a polynomial to the data
//															AB	29 July 2018															
//	---------------------------------------------------------------------------------------

int dopolyfit (float *xarr, float *yarr, float *yerr, double *coeffs, int lendata, int order) {

	int 		i, j, status;
	double 		chisq;
	
	if((order > 0) && (lendata > 2*order)) {
	    gsl_multifit_linear_workspace 	*ws;
     	gsl_matrix 						*cov, *xx;
      	gsl_vector 						*y, *c, *wt;
      	
	    xx 		= gsl_matrix_alloc (lendata, order);
      	y 		= gsl_vector_alloc (lendata);
      	wt 		= gsl_vector_alloc (lendata);
      	c 		= gsl_vector_alloc (order);
      	cov 	= gsl_matrix_alloc (order, order);

	    for(i=0; i < lendata; i++) {
		    for(j=0; j < order; j++) 
		      gsl_matrix_set( xx, i, j, pow(xarr[i], j));
		
		    gsl_vector_set(y, i, yarr[i]);
		    if (yerr==NULL)
			    gsl_vector_set(wt, i, 1.0);
		    else
			    gsl_vector_set(wt, i, 1.0/pow(yerr[i],2.0));
	    }
	
	    ws		= gsl_multifit_linear_alloc (lendata, order);
	    status	= gsl_multifit_wlinear (xx, wt, y, c, cov, &chisq, ws);

	    for(i=0; i < order; i++)
        	coeffs[i] = gsl_vector_get(c, i);

	    gsl_multifit_linear_free(ws);
      	gsl_matrix_free(xx);
      	gsl_matrix_free(cov);
      	gsl_vector_free(y);
     	gsl_vector_free(c);
     	gsl_vector_free(wt);
    }
    else    {
        status      =   0;
        coeffs[0]   =   find_mean(yarr,lendata);
    }
	return status;
}


























































































