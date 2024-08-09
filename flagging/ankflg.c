#	include <stdio.h>
#	include <stdlib.h>	
#	include	<ankhead.h>
#	include <math.h>
#	include <string.h>


//	-----------------------------------------------------------------------------------

int main()
{
	int				i,j,k;
	
	int				ants, flagmode, done;
	char			statusfile[100],binfiles[100],flagparfile[100],scanfile[100],badbasename[100];
	int				nbase	=	0;
	
	sprintf(flagparfile,"scratch/flagpars.pars");
	sprintf(badbasename,"scratch/badbase.list");
	sprintf(scanfile,"scratch/scandetails.txt");
	
	FILE		*fp,*badbasefile;
	fp		=	fopen(flagparfile,"r");
	badbasefile	=	fopen(badbasename,"w");
	
	fscanf (fp,"%d	%d",&ants,&flagmode);
	
	fclose(fp);
	
	if (flagmode) {
		sprintf(statusfile,"scratch/uvbin_status.txt");
		sprintf(binfiles,"scratch/uvbin");	
	}
	else {
		sprintf(statusfile,"scratch/baseline_status.txt");
		sprintf(binfiles,"scratch/baseline");
	}
	
	printf("\nFlagmode	= %d	Antennas	= %d\n\n",flagmode,ants);
	
	done	=	processbaselines (ants, statusfile, &nbase, binfiles, flagparfile, badbasefile, scanfile);	
	
	printf("Total baselines / uvbins		%d\n",nbase);
		
	fclose (badbasefile);
	
	return 0;
}

























