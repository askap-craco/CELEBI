#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>

int	flagbaseline (BaseParType *basepar, float bl_minfrac, float blockpow, int flagrounds, FlagParType *flagparams, float **datarray, float *goutarr,
					float **buffarrc, float **buffarrr, double *coeffs, float **buffdev1d, float ***buffdev2d, float **statsarray)	{
		
	int			i,j,k,l,pol,r,c,curound,cw,rw;
	int			flagstat;
	float		initflag, currflag;
	FlagParType	flagpar;
			
//	------------------------------------------		Flag iterating on flag commands	
	
	for (i = 0; i < basepar->polsize; i++) {
		currflag			=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
		basepar->flfrac[i]	=	currflag;
	}
	
	for (curound = 0; curound < flagrounds; curound++)	{
		
		//printf("\n%d %d %d %d %.2f %.2f %d %d %d %d %.2f %.2f\n",flagparams[curound].whatflag,flagparams[curound].doflag,flagparams[curound].domean,
					//flagparams[curound].qtype,flagparams[curound].tolerance,flagparams[curound].min_fraction,flagparams[curound].fitorder,
					//flagparams[curound].ascending,flagparams[curound].chanblockfac,flagparams[curound].recblockfac,flagparams[curound].chanmaxfrac,flagparams[curound].recmaxfrac);
		
		for (i = 0; i < basepar->polsize; i++) {
			
			if (init_data (basepar, datarray, i, flagparams[curound].qtype)) {
				printf("\nWait till I learn it .....\n");
				return 1;
			}
		
			currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);			
			initflag	=	currflag;
			
			switch (flagparams[curound].whatflag)	{
		
				case 0:
					flagstat	=	1;
					while (flagstat) {	
						flagstat	=	flag_vis_individual (flagparams[curound], *basepar, datarray, basepar->data[i][2], bl_minfrac, goutarr, buffarrc, buffarrr, coeffs,
																															buffdev1d, buffdev2d);	
					}	
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging vis 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);						
					break;
				
				case 1:
					flagstat	=	1;
					while (flagstat) {	
						flagstat	=	flag_chan_individual (flagparams[curound], *basepar, datarray, basepar->data[i][2], bl_minfrac, goutarr, buffarrc, buffarrr, coeffs, buffdev1d);
					}
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging chan 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);							
					break;
				
				case 2:
					flagstat	=	1;
					while (flagstat) {	
						flagstat	=	flag_rec_individual (flagparams[curound], *basepar, datarray, basepar->data[i][2], bl_minfrac, goutarr, buffarrc, buffarrr, buffdev1d);
					}
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging rec 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);							
					break;
				
				case 3:
					rw		=	1;
					while (rw < (int)(basepar->recsize*flagparams[curound].recmaxfrac)) {
						cw	=	1;
						while (cw < (int)(basepar->chansize*flagparams[curound].chanmaxfrac)) {
							if ((rw != 1) || (cw != 1)) {
								//printf("Baseline	%d-%d	pol %d flagging vblock 		rw %d	cw %d\n",basepar->anta,basepar->antb,i,rw,cw);
								flagstat	=	1;
								while (flagstat) {							
									flagstat = flag_vis_block(flagparams[curound],*basepar,datarray,basepar->data[i][2],bl_minfrac,goutarr,rw,cw,blockpow, buffarrc, buffarrr, coeffs, 
																																						buffdev1d, buffdev2d);
								}
							}
							cw	=	cw * flagparams[curound].chanblockfac;
						}	
						rw	=	rw * flagparams[curound].recblockfac;
					}
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging vblock 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);							
					break;
					
				case 4:
					cw			=	flagparams[curound].chanblockfac;
					while (cw < (int)(basepar->chansize*flagparams[curound].chanmaxfrac)) {
						flagstat	=	1;
						while (flagstat) {							
							flagstat = flag_chan_block(flagparams[curound],*basepar,datarray,basepar->data[i][2],bl_minfrac,goutarr,cw,blockpow, buffarrc, buffarrr, coeffs, 
																																					buffdev1d, buffdev2d);
						}
						cw	=	cw * flagparams[curound].chanblockfac;
					}
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging cblock 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);							
					break;
				
				case 5:
					rw			=	flagparams[curound].recblockfac;
					while (rw < (int)(basepar->recsize*flagparams[curound].recmaxfrac)) {
						flagstat	=	1;
						while (flagstat) {							
							flagstat = flag_rec_block(flagparams[curound],*basepar,datarray,basepar->data[i][2],bl_minfrac,goutarr,rw,blockpow, buffarrc, buffarrr, coeffs, 
																																					buffdev1d, buffdev2d);
						}
						cw	=	cw * flagparams[curound].chanblockfac;
						rw	=	rw * flagparams[curound].recblockfac;
					}
					currflag	=	1.0 - ((float) count_nonzero_2d(basepar->data[i][2], basepar->recsize, basepar->chansize)) / (basepar->recsize*basepar->chansize);
					//printf("Baseline	%d-%d	pol %d flagging rblock 	%d	%d =	%.3f\n",
												//basepar->anta,basepar->antb,i,flagparams[curound].domean,flagparams[curound].qtype,currflag-initflag);							
					break;
					
				default:
					printf("\nUnknown language !!  Are you speaking Hebrew? \n");
			}			
		}
	}
	
	/*for (i = 0; i < basepar->polsize; i++)
		for (r = 0; r < basepar->recsize; r++)
			for (c = 0; c < basepar->chansize; c++)
				basepar->data[i][2][r][c] = 0.0;*/
		
	for (i = 0; i < basepar->polsize; i++) {
		if (init_data (basepar, datarray, i, 2)) {
			printf("\nWait till I learn it .....\n");
			return 1;
		}
		basepar->flfrac_after[i]	=	findbasestat(*basepar,datarray,basepar->data[i][2],bl_minfrac,statsarray[i],buffarrc,buffarrr,coeffs,buffdev1d,buffdev2d);
	}
	
	/*for (i = 0; i < basepar->polsize; i++)
		for (r = 0; r < basepar->recsize; r++)
			for (c = 0; c < basepar->chansize; c++)
				//if (basepar->data[i][2][r][c] <= 0.0)
					basepar->data[i][0][r][c] = basepar->data[i][1][r][c] = basepar->data[i][2][r][c] = 0.0;*/		
						
	return 0;
}






































































































































































