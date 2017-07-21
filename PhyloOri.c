/* PhyloOri.c */

# include "libhdr"
# include <string.h>

/* ************************************************************************* */

int GEN, TGEN, g, gen, NSAM, r, i, j, k, x, y, z, pop, REP, oo[100];
int ii, is, it;
int FtreeLoc[100], FtreeEco[100], FtreeOther[100];
int HtreeLoc[100], HtreeEco[100], HtreeOther[100];
char ch;
double w, d[10];

FILE *fdatF, *fdatH, *ftre;

/* ************************************************************************** */

main()
{
	char outnameF[30]="DISTANCES_FREQ_";
	char outnameH[30]="DISTANCES_HAPLO_";
	char id[15];

	tracestart();

	ftre = fopen("TREEFILE","w");

	// getinputs

	getintandskip("Number of subpopulations (max 100):",&pop,1,100);
	getintandskip("Number of generations:",&GEN,1,infinity);
	getintandskip("Number of generations between trees:",&TGEN,1,infinity);
	getintandskip("Number of sampled individuals:",&NSAM,1,500);
	getintandskip("Number of replicates:",&REP,1,infinity);
	
	for (gen=0; gen<=GEN; gen+=TGEN)
	{
		if (tracelevel != 0)	fprintf (ftre, "***** gen = %d *****\n\n", gen);

		g = gen / TGEN;

		// DISTANCES FROM FREQUENCIES

		sprintf (id, "%d", gen);
		strcat(outnameF,id);
		fdatF = fopen( outnameF, "r" );
		memset(outnameF,'\0',strlen(outnameF));
		strcat(outnameF,"DISTANCES_FREQ_"); 

		for (r=1; r<=REP; r++)
		{
			if (tracelevel != 0)	fprintf (ftre, "*** FREQ rep = %d ***\n\n", r);

			fscanf(fdatF,"%d", &x);

			for (j=1; j<=13; j++)	fscanf(fdatF,"%c", &ch);
			for (j=1; j<=2; j++)		fscanf(fdatF,"%lf", &w);
			d[1] = w;
			fscanf(fdatF,"%lf", &w);
			d[2] = w;
			fscanf(fdatF,"%lf", &w);
			d[3] = w;
			for (j=1; j<=13; j++)	fscanf(fdatF,"%c", &ch);
			for (j=1; j<=3; j++)		fscanf(fdatF,"%lf", &w);
			d[4] = w;
			fscanf(fdatF,"%lf", &w);
			d[5] = w;
			for (j=1; j<=13; j++)	fscanf(fdatF,"%c", &ch);
			for (j=1; j<=3; j++)		fscanf(fdatF,"%lf", &w);
			fscanf(fdatF,"%lf", &w);
			d[6] = w;
			for (j=1; j<=13; j++)	fscanf(fdatF,"%c", &ch);
			for (j=1; j<=4; j++)		fscanf(fdatF,"%lf", &w);
			fscanf(fdatF,"%c", &ch);

			for (i=1; i<=6; i++)    if (d[i] == (-1))	d[i] = 10.0;

			if (tracelevel != 0)	fprintf (ftre, "FREQ d[1]=%f d[2]=%f d[3]=%f d[4]=%f d[5]=%f d[6]=%f\n\n", d[1], d[2], d[3], d[4], d[5], d[6]);

			// order distances

			for (i=1; i<=6; i++)    oo[i]=i;

			for (is=1; is<=6; is++)
			{
				ii=is;
				for (it=is+1; it<=6; it++)
					if (d[oo[ii]] > d[oo[it]])     ii=it;
			
				x=oo[is]; oo[is]=oo[ii]; oo[ii]=x;
			}

			//if (tracelevel != 0)	fprintf (ftre, "oo[1]=%d oo[2]=%d \n", oo[1], oo[2]);

			// frequency tree type

			if ( ( (oo[1]==1) && (oo[2]==6) ) || ( (oo[1]==6) && (oo[2]==1) ) )		FtreeLoc[g] ++;
			else if ( ( (oo[1]==2) && (oo[2]==5) ) || ( (oo[1]==5) && (oo[2]==2) ) )	FtreeEco[g] ++;
			else											FtreeOther[g] ++;
		}

		fclose(fdatF);									

		// DISTANCES FROM HAPLOTYPES

		sprintf (id, "%d", gen);
		strcat(outnameH,id);
		fdatH = fopen( outnameH , "r" );
		memset(outnameH,'\0',strlen(outnameH));
		strcat(outnameH,"DISTANCES_HAPLO_"); 

		for (r=1; r<=REP; r++)
		{
			if (tracelevel != 0)	fprintf (ftre, "*** HAPLO rep = %d ***\n\n", r);

			// TAKE ALL DISTANCES FOR ALL GROUPS OF SAMPLED INDIVIDUALS (NSAM)

			for (i=1; i<=6; i++)	d[i]=0.0;

			for (i=1; i<=NSAM; i++)
			{
				fscanf(fdatH,"%d", &x);

				for (j=1; j<=13; j++)	fscanf(fdatH,"%c", &ch);
				for (j=1; j<=2; j++)		fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[1] += w;
				fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[2] += w;
				fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[3] += w;
				for (j=1; j<=13; j++)	fscanf(fdatH,"%c", &ch);
				for (j=1; j<=3; j++)		fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[4] += w;
				fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[5] += w;
				for (j=1; j<=13; j++)	fscanf(fdatH,"%c", &ch);
				for (j=1; j<=3; j++)		fscanf(fdatH,"%lf", &w);
				fscanf(fdatH,"%lf", &w);
				if (w == (-1)) w=10.0;
				d[6] += w;
				for (j=1; j<=13; j++)	fscanf(fdatH,"%c", &ch);
				for (j=1; j<=4; j++)		fscanf(fdatH,"%lf", &w);
				fscanf(fdatH,"%c", &ch);
			}

			// AVERAGE VALUES AMONG GROUPS OF SAMPLED INDIVIDUALS (NSAM)

			for (i=1; i<=6; i++)	 d[i] = d[i] / (double)NSAM;
			
			if (tracelevel != 0)	 fprintf (ftre, "\nHAPLO d[1]=%f d[2]=%f d[3]=%f d[4]=%f d[5]=%f d[6]=%f\n\n", d[1], d[2], d[3], d[4], d[5], d[6]);

			// order distances
	
			for (i=1; i<=6; i++)    oo[i]=i;

			for (is=1; is<=6; is++)
			{
				ii=is;
				for (it=is+1; it<=6; it++)
					if (d[oo[ii]] > d[oo[it]])     ii=it;
				
				x=oo[is]; oo[is]=oo[ii]; oo[ii]=x;
			}

			//if (tracelevel != 0)	fprintf (ftre, "oo[1]=%d oo[2]=%d \n", oo[1], oo[2]);

			// haplotype tree type

			if ( ( (oo[1]==1) && (oo[2]==6) ) || ( (oo[1]==6) && (oo[2]==1) ) )		HtreeLoc[g] ++;
			else if ( ( (oo[1]==2) && (oo[2]==5) ) || ( (oo[1]==5) && (oo[2]==2) ) )	HtreeEco[g] ++;
			else											HtreeOther[g] ++;
		}

		fclose(fdatH);									

	}

	// printoutput
	
//	for (i=1; i<=6; i++)
//	fprintf (ftre, "d[%d]=%f\n", oo[i], d[oo[i]]);

	fprintf (ftre, "GEN 	L  E  O\n\n");

	fprintf (ftre, "FREQ\n");
	g=0;					fprintf (ftre, "%d     %d  %d  %d\n", g*10000, FtreeLoc[g], FtreeEco[g], FtreeOther[g]);
	for (g=1; g<=(GEN/TGEN); g++)	fprintf (ftre, "%d  %d  %d  %d\n", g*10000, FtreeLoc[g], FtreeEco[g], FtreeOther[g]);

	fprintf (ftre, "\nHAPLO\n");
	g=0;					fprintf (ftre, "%d     %d  %d  %d\n", g*10000, HtreeLoc[g], HtreeEco[g], HtreeOther[g]);
	for (g=1; g<=(GEN/TGEN); g++)	fprintf (ftre, "%d  %d  %d  %d\n", g*10000, HtreeLoc[g], HtreeEco[g], HtreeOther[g]);

	return(0);
}

/* ************************************************************************** */

