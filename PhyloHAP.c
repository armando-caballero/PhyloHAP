
/* PhyloHAP.c (n haploid subpopulations RB and SU; 1 selected locus s=0.99;
1 tag original haplotype; and x neutral; study only extreme populations; NSAM samples for frequency and haplotype tree) */

# include "libhdr"
# include <time.h>
# include <string.h>

# define NmaxS 101 // max pop
# define NmaxI 10001 // max nind and haplotypes
# define NmaxL 3001 // loci
# define maxmpt 5001
# define normalthreshold 30

/* ************************************************************************* */

int INITYPE, nind, nhap, NIND, NSAM, pop, z, h, i, j, k, l, m, n, a;
int gen, GEN, TREEGEN, l1, l2;
int muts, lastinmutantspoissontable, lastinrecombinantpoissontable;
int migs, lastinmigrantsECOpoissontable, lastinmigrantsLOCpoissontable;
int NLOCI, NCRO, TOTLOCI, RM[NmaxL];
int gm[NmaxS][NmaxI][201], sm[NmaxS][NmaxI][201];
int ran_i, ran_l, ran_k, ran_m;
int EE[201], FF[201], p1, p2, marker;
int haplotype[NmaxI][201], type, typehaplotype[NmaxI];
int numberrecs, nr, pointrec[201][31], ncrorec[201], rndk, rndl;

double u, mutantspoissontable[maxmpt];
double freqRB1[201][31], freqRB2[201][31], freqSU1[201][31], freqSU2[201][31];
double mrECO, migrantsECOpoissontable[maxmpt];
double mrLOC, migrantsLOCpoissontable[maxmpt];
double L, recombinantpoissontable[maxmpt];
double s, pm_s[NmaxS][NmaxI], numRec;

void recombination_masks();

FILE *fsal, *fdat, *fdat0, *fphylF, *fphylH, *fichero;

/* ************************************************************************** */

main()
{
	tracestart();
	getseed();	
	getinputs();
	recombination_masks(RM);
	headings();
	initialize();

	for (gen=1; gen<=GEN; gen++)
	{
		mutation();
		migration();
		fitnesses();
		mating();
		if ((gen == 1) || (gen%TREEGEN == 0))
		{
			frequencies();
			writing_phylfile();
			writing_phylfreqfile();
		}
	}
	finding_haplotypes();
	output_file();
	printf("numRec=%f\n", numRec);

	writeseed();
	return(0);
}

/* ************************************************************************** */

getinputs()
{
	fdat = fopen("datafile","w");
	fdat0 = fopen("datafile0","w");
	if (tracelevel!=0)	fsal = fopen("salfile","w");

	getintandskip("Number of subpopulations (max 100):",&pop,1,100);
	getintandskip("Number of individuals (max 10000) for subpopulation :",&nind,1,10000);
	NIND = nind * pop;
	getrealandskip("Migration rate (ECO):",&mrECO,0.0,1.0);
	getrealandskip("Migration rate (LOC):",&mrLOC,0.0,1.0);
	getintandskip("Number of pieces of chromosome (min 1, max 200) :",&NCRO,1,200);
	getintandskip("Number of loci per piece (min 1, max 30) :",&NLOCI,1,30);
	TOTLOCI = NCRO * NLOCI;
	getrealandskip("Mutation rate:",&u,0.0,1.0);
	getrealandskip("Length of sequence in Morgans :",&L,0.0,99.0);
	getrealandskip("Selection coefficient:",&s,0.0,1.0);
	getintandskip("Number of generations:",&GEN,1,infinity);
	getintandskip("Number of generations between trees :",&TREEGEN,1,infinity);
	getintandskip("Number of sampled individuals (max 100):",&NSAM,1,100);
	getintandskip("Random(0), habitat occupancy (1), locality occupancy (2):",&INITYPE,0,2);

  	/* MIGRATION RATE WITH POISSON (mrECO) NEW MIGRANTS PER GENERATION */
	if ( (exp(-(double)nind*(double)mrECO) != 0.0)&&((double)nind*(double)mrECO < normalthreshold) )
	generatepoissontable((double)nind*(double)mrECO, &lastinmigrantsECOpoissontable, migrantsECOpoissontable, maxmpt-1);

  	/* MIGRATION RATE WITH POISSON (mrLOC) NEW MIGRANTS PER GENERATION */
	if ( (exp(-(double)nind*(double)mrLOC) != 0.0)&&((double)nind*(double)mrLOC < normalthreshold) )
	generatepoissontable((double)nind*(double)mrLOC, &lastinmigrantsLOCpoissontable, migrantsLOCpoissontable, maxmpt-1);

  	/* BIALLELIC LOCI WITH POISSON (Nu) NEW MUTATIONS PER LOCUS */
	if ( (exp(-(double)NIND*((double)TOTLOCI-1)*u) != 0.0)&&((double)NIND*((double)TOTLOCI-1)*u < normalthreshold) )
	generatepoissontable((double)NIND*((double)TOTLOCI-1)*u, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

	/* NUMBER OF RECOMBINATIONS POISSON WITH MEAN L*/
	if ( (exp(-L) != 0.0) && (L < normalthreshold) )
	generatepoissontable(L, &lastinrecombinantpoissontable, recombinantpoissontable, maxmpt-1);
}

/* ************************************************************************** */

void recombination_masks (RM)
int RM[];
{
	for (m=0; m<NLOCI; m++)   RM[m]=pow(2.0,(double)m);
}

/* ************************************************************************** */

headings()
{
	printf("\n\n pop=%d    nind=%d     TOTLOCI=%d    gens=%d    INITYPE=%d\n", pop, nind, TOTLOCI, GEN, INITYPE);
	printf(" u=%7.6f    L=%f    mECO=%f    mLOC=%f    \n", u, L, mrECO, mrLOC);
}

/* ************************************************************************** */

initialize()
{
	if (INITYPE == 0)
	{
		for (i=1; i<=pop; i++)
		for (l=1; l<=nind; l++)
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	if (uniform() < 0.5) gm[i][l][k] = (gm[i][l][k] | RM[m]);
	}
	else if (INITYPE == 1)
	{
		// ********* Allele 0 in RB populations (i=1,3,5,...), allele 1 in SU populations (habitat)
		for (i=2; i<=pop; i+=2)
		for (l=1; l<=nind; l++)
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	gm[i][l][k] = (gm[i][l][k] | RM[m]);
	}
	else
	{
		// ********* Allele 1 in left part (i=1,2,3,...,pop/2), allele 0 in SU populations (locality)
		for (i=1; i<=pop/2; i++)
		for (l=1; l<=nind; l++)
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	gm[i][l][k] = (gm[i][l][k] | RM[m]);
	}

	if (tracelevel!=0)
	{
		fprintf(fsal,"\n Genotypes of initial generation\n");
		printing_for_check();
	}

	output_file_gen0();
}

/* ************************************************************************** */

mutation()
{
	muts = mutationnumber();

	if ((tracelevel!=0)&&(gen == GEN))   fprintf(fsal,"\n New mutants\n");

	if (muts > 0)
	for (m=1; m<=muts; m++)
	{
		ran_i = (int)(uniform()*pop)+1;
		ran_l = (int)(uniform()*nind)+1;

		// There is no mutation in the first locus
		do
		{
			ran_k = (int)(uniform()*NCRO);
			ran_m = (int)(uniform()*NLOCI);
		}
		while ( (ran_k==0) && (ran_m==0) );

		if ( (gm[ran_i][ran_l][ran_k] & RM[ran_m])==RM[ran_m] )
			gm[ran_i][ran_l][ran_k]=(gm[ran_i][ran_l][ran_k] & (~RM[ran_m]));
		else	gm[ran_i][ran_l][ran_k]=(gm[ran_i][ran_l][ran_k] | RM[ran_m]);

		if ((tracelevel!=0)&&(gen == GEN))    fprintf(fsal,"ran_i=%d  ran_l=%d  ran_k=%d  ran_m=%d\n", ran_i, ran_l, ran_k, ran_l);
	}
}

/* ************************************************************************** */

int mutationnumber()
{
	int r;
	if (((double)NIND*((double)TOTLOCI-1)*u < normalthreshold) && (exp(-(double)NIND*((double)TOTLOCI-1)*u) != 0.0) )
	{
		r= poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal((double)NIND*((double)TOTLOCI-1)*u, sqrt((double)NIND*((double)TOTLOCI-1)*u)) );
	return(r);
}

/* ************************************************************************** */

migration()
{
	int k1, k2;

	// ********* Migration (mrLOC) between localities (i=1-3,2-4,...)

	if (mrLOC > 0.0)	migs = migrationLOCnumber();
	else			migs = 0;

	for (i=1; i<=(pop-2); i+=1)
	{
		j = i+2;

		if ((tracelevel!=0)&&(gen == GEN))		fprintf(fsal,"\ni=%d j=%d RBmigs=%d\n", i, j, migs);

		for (z=1; z<=migs; z++)
		{
			l1 = (int)(uniform()*nind)+1;
			l2 = (int)(uniform()*nind)+1;

			for (k=0; k<NCRO; k++)
			{
				k1 = gm[i][l1][k];
				gm[i][l1][k] = gm[j][l2][k];
				gm[j][l2][k] = k1;
			}
		}
	}

	// ********* Migration (mrECO) between RB and SU populations (i=1-2,3-4,5-6,...)

	if (mrECO > 0.0)	migs = migrationECOnumber();
	else			migs = 0;

	for (i=1; i<pop; i+=2)
	{
		j = i+1;

		if ((tracelevel!=0)&&(gen == GEN))		fprintf(fsal,"\ni=%d j=%d RB-SUmigs=%d\n", i, j, migs);

		for (z=1; z<=migs; z++)
		{
			l1 = (int)(uniform()*nind)+1;
			l2 = (int)(uniform()*nind)+1;

			for (k=0; k<NCRO; k++)
			{
				k1 = gm[i][l1][k];
				gm[i][l1][k] = gm[j][l2][k];
				gm[j][l2][k] = k1;
			}
		}
	}

	if ((tracelevel!=0)&&(gen == GEN))
	{
		fprintf(fsal,"\n Genotypes of progeny after migration\n");
		printing_for_check();
	}
}

/* ************************************************************************** */

int migrationLOCnumber()
{
	int r;
	if (((double)nind*(double)mrLOC < normalthreshold) && (exp(-(double)nind*(double)mrLOC) != 0.0) )
	{
		r= poisson(lastinmigrantsLOCpoissontable, migrantsLOCpoissontable);
	}
	else r = (int)( normal((double)nind*(double)mrLOC, sqrt((double)nind*(double)mrLOC)) );
	return(r);
}

/* ************************************************************************** */

int migrationECOnumber()
{
	int r;
	if (((double)nind*(double)mrECO < normalthreshold) && (exp(-(double)nind*(double)mrECO) != 0.0) )
	{
		r= poisson(lastinmigrantsECOpoissontable, migrantsECOpoissontable);
	}
	else r = (int)( normal((double)nind*(double)mrECO, sqrt((double)nind*(double)mrECO)) );
	return(r);
}

/* ************************************************************************** */

fitnesses()
{
	for (i=1; i<=pop; i++)
	for (l=1; l<=nind; l++)	pm_s[i][l] = 1.0;

	for (i=1; i<=pop; i++)
	{
		if (i%2 != 0)
		{
			// RB populations, where allele 1 is deleterious

			for (l=1; l<=nind; l++)
			{
				if ((gm[i][l][0] & RM[0]) != RM[0])	pm_s[i][l] *= 1.0;		/*0*/	
				else    					pm_s[i][l] *= (1.0 - s);	/*1*/
			}
		}
		else
		{
			// SU populations, where allele 0 is deleterious
			for (l=1; l<=nind; l++)
			{
				if ((gm[i][l][0] & RM[0]) != RM[0])	pm_s[i][l] *= (1.0 - s);	/*0*/	
				else    					pm_s[i][l] *= 1.0;		/*1*/
			}
		}
	}

	if ((tracelevel!=0) && (gen == GEN))
	{
		fprintf(fsal,"\n Fitnesses\n");
		for (i=1; i<=pop; i++)
		for (l=1; l<=nind; l++)
		if (l<10)	fprintf(fsal,"pm_s[%d][%d]=%f\n", i, l, pm_s[i][l]);
	}

	return(0);
}

/* ************************************************************************** */

mating()
{
	double cum, cum_pm_s[NmaxI];

	for (i=1; i<=pop; i++)
	for (l=1; l<=nind; l++)
	for (k=0; k<NCRO; k++)	sm[i][l][k] = gm[i][l][k];

	if ((tracelevel!=0)&&(gen == GEN))
	{
		fprintf(fsal,"\n Genotypes of parents\n");
		for (i=1; i<=pop; i++)
		for (l=1; l<=nind; l++)
		{
			if (l<10)
			{
				fprintf(fsal,"\nsm[%d][%d][0]=%d\n", i, l, sm[i][l][0]);

				fprintf (fsal,"i=%d l=%d   ", i, l);

				for (k=0; k<NCRO; k++)
				for (m=0; m<NLOCI; m++)
				{
					if ((sm[i][l][k] & RM[m]) == RM[m])	fprintf (fsal,"1 ");
					else						fprintf (fsal,"0 ");
				}
			}
		}
	}

	for (i=1; i<=pop; i++)
	{
		// Relative fitnesses

		cum = 0.0;
		for (l=1; l<=nind; l++)
		{
			cum += pm_s[i][l];
			cum_pm_s[l] = cum;
		}
		for (l=1; l<=nind; l++)	cum_pm_s[l] = cum_pm_s[l] / cum;

		if ((tracelevel!=0)&&(gen == GEN))
		{
			fprintf(fsal,"\ncum = %f\n", cum);
			fprintf(fsal,"\n Scaled relative fitnesses subpop=%d\n", i);
			for (l=1; l<=nind; l++)   fprintf(fsal,"cum_pm_s = %f\n", cum_pm_s[l]);
		}

		if ((tracelevel!=0)&&(gen == GEN))	fprintf(fsal,"\n Parents of subpopulation %d\n", i);

		for (l=1; l<=nind; l++)
		{
			p1 = binarysearchcumulative(cum_pm_s, nind);
			p2 = binarysearchcumulative(cum_pm_s, nind);

			if ((tracelevel!=0)&&(gen == GEN))		if (l<10)	fprintf (fsal,"p1=%d\tp2=%d\n", p1, p2);

			/* ******************* Neutral linked ******************* */

			RECOMBINANT_CHROMOSOME();

			if (uniform() < 0.5)		for (k=0; k<NCRO; k++)	gm[i][l][k]=((EE[k]&sm[i][p1][k])|(FF[k]&sm[i][p2][k]));
			else				for (k=0; k<NCRO; k++)	gm[i][l][k]=((FF[k]&sm[i][p1][k])|(EE[k]&sm[i][p2][k]));


/**************************************************/
//CHECK: COMMENT IN OUT FILE

if (numberrecs==4)
{
	printf ("\nINDIVIDUAL i=%d l=%d   \n", i, l);

	printf ("\n\nPARENTS\n\n");

	printf ("\ni=%d l=%d   \n", i, p1);

	for (k=0; k<NCRO; k++)
	for (m=0; m<NLOCI; m++)
	{
		if ((sm[i][p1][k] & RM[m]) == RM[m])	printf ("1 ");
		else						printf ("0 ");
	}

	printf ("\ni=%d l=%d   \n", i, p2);

	for (k=0; k<NCRO; k++)
	for (m=0; m<NLOCI; m++)
	{
		if ((sm[i][p2][k] & RM[m]) == RM[m])	printf ("1 ");
		else						printf ("0 ");
	}

	//EE and FF

	printf ("\n\nnumberrecs=%d\n", numberrecs);
	printf ("EE\n");
	for (k=0; k<NCRO; k++)
	{
		for (m=0; m<NLOCI; m++)
		{
			if((EE[k]&RM[m])==RM[m])  	printf ("1 ");
			else			   	printf ("0 ");
		}
	}
	printf ("\nFF\n");
	for (k=0; k<NCRO; k++)
	{
		for (m=0; m<NLOCI; m++)
		{
			if((FF[k]&RM[m])==RM[m])  	printf ("1 ");
			else			   	printf ("0 ");
		}
	}

	printf ("\n\nOFFSPRING\n\n");

	printf ("\ni=%d l=%d   \n", i, l);

	for (k=0; k<NCRO; k++)
	for (m=0; m<NLOCI; m++)
	{
		if ((gm[i][l][k] & RM[m]) == RM[m])	printf ("1 ");
		else						printf ("0 ");
	}
}
/****************************************************/

		}
	}

	/* *************************** Printing ******************************* */

	if ((tracelevel!=0)&&(gen == GEN))
	{
		fprintf(fsal,"\n Genotypes of progeny\n");
		printing_for_check();
	}
}

/* ************************************************************************** */

int recombinationnumber()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecombinantpoissontable, recombinantpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ************************************************************************** */

RECOMBINANT_CHROMOSOME()
{
	/* ******************* Neutral sequence with recombination c between SNPs ******************* */

	if (L == 0.0)	 	for (k=0; k<NCRO; k++)	EE[k]=0;
	else if (L == 99)	for (k=0; k<NCRO; k++)	EE[k]=(int)(uniform()*(pow(2.0,(double)NLOCI)));
	else
	{
		numberrecs = recombinationnumber();
		numRec += numberrecs/(double)(GEN*pop*nind);

//		if (numberrecs != 0)	printf("numberrecs=%d\n", numberrecs);

		if (tracelevel!=0)   fprintf (fsal,"numberrecs=%d\n", numberrecs);

		if (numberrecs == 0)
		{
			for (k=0; k<NCRO; k++)	EE[k]=0;
		}
		else
		{
			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (m=0; m<NLOCI; m++)  pointrec[k][m] = 0;
			}
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
			}	

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (m=0; m<NLOCI; m++)
			      		{
						if (pointrec[k][m] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[m];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[m];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}
		}
	}
	/* ************ Selected locus with free recombination *********** */;

	EE[0] = ( EE[0] & (~RM[0]) );
	if (uniform() < 0.5)		EE[0] = (EE[0] | RM[0]);

	/* *************************** Mask FF ************************** */;

	for (k=0; k<NCRO; k++)	FF[k] = ~EE[k];
}

/* ************************************************************************** */

printing_for_check()
{
	for (i=1; i<=pop; i++)
	for (l=1; l<=nind; l++)
	{
		if (l<=5)
		{
			fprintf(fsal,"\ngm[%d][%d][0]=%d\n", i, l, gm[i][l][0]);

			fprintf (fsal,"i=%d l=%d   ", i, l);

			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			{
				if ((gm[i][l][k] & RM[m]) == RM[m])	fprintf (fsal,"1 ");
				else						fprintf (fsal,"0 ");
			}
		}
	}
	return(0);
}

/* ************************************************************************** */

frequencies()
{
	for (k=0; k<NCRO; k++)
	for (m=0; m<NLOCI; m++)
	{
		freqRB1[k][m] = 0.0;
		freqRB2[k][m] = 0.0;
		freqSU1[k][m] = 0.0;
		freqSU2[k][m] = 0.0;
	}

	i=1; // RB1

		for (l=1; l<=NSAM; l++)
		{
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			{
				if ( (k==0) && (m==0) )	/****/;
				else if ((gm[i][l][k] & RM[m]) == RM[m])	freqRB1[k][m] += 1.0;
			}
		}
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	freqRB1[k][m] = freqRB1[k][m] / (double)NSAM;

	i=pop-1; // RB2

		for (l=1; l<=NSAM; l++)
		{
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			{
				if ( (k==0) && (m==0) )	/****/;
				else if ((gm[i][l][k] & RM[m]) == RM[m])	freqRB2[k][m] += 1.0;
			}
		}
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	freqRB2[k][m] = freqRB2[k][m] / (double)NSAM;


	i=2; // SU1

		for (l=1; l<=NSAM; l++)
		{
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			{
				if ( (k==0) && (m==0) )	/****/;
				else if ((gm[i][l][k] & RM[m]) == RM[m])	freqSU1[k][m] += 1.0;
			}
		}
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	freqSU1[k][m] = freqSU1[k][m] /(double)NSAM;

	i=pop; // SU2

		for (l=1; l<=NSAM; l++)
		{
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			{
				if ( (k==0) && (m==0) )	/****/;
				else if ((gm[i][l][k] & RM[m]) == RM[m])	freqSU2[k][m] += 1.0;
			}
		}
		for (k=0; k<NCRO; k++)
		for (m=0; m<NLOCI; m++)	freqSU2[k][m] = freqSU2[k][m] / (double)NSAM;
}


/* ************************************************************************** */

writing_phylfile()
{
	char outname[50]="phylfile";
	char id[10];

	sprintf (id, "%d", gen);
	strcat(outname,id);
	strcat(outname,".txt");
	printf("Crea fichero %s\n",outname);
	fphylH = fopen( outname, "w" );
	memset(outname,'\0',strlen(outname)); // limpia el array outname, se podría limpiar solo hasta el número a cambiar..
	strcat(outname,"phylfreqfile"); // reinicializa el array con el nombre base, no haría falta si se hubiera limpiado solo la parte que cambia

//	fprintf (fphylH,"gen=%d\n", gen);

//	fprintf (fphylH, "%d	%d\n", 4*NSAM, NCRO*NLOCI-2);

	for (l=1; l<=NSAM; l++)
	{
		fprintf (fphylH, "%d	%d\n", 4, NCRO*NLOCI-2);

		for (i=1; i<=pop; i++)
		{
			if ((i == 1)||(i == 2)||(i == pop-1)||(i == pop))
			{
				/***************** Population sampled *****************/

				if (i==1) 		fprintf (fphylH, "RB1       ");
				else if (i==2) 	fprintf (fphylH, "SU1       ");
				else if (i==pop-1) 	fprintf (fphylH, "RB2       ");
				else if (i==pop) 	fprintf (fphylH, "SU2       ");

				/************** Neutral linked haplotype **************/

				for (k=0; k<NCRO; k++)
				for (m=0; m<NLOCI; m++)
				{
					if ( (k==0) && (m==0) )	/****/;
					else if ((gm[i][l][k] & RM[m]) == RM[m])	fprintf (fphylH,"T");
					else						fprintf (fphylH,"A");
				}
				fprintf (fphylH,"\n");
			}
		}
	}
	return(0);
}

/* ************************************************************************** */

writing_phylfreqfile()
{
	char outname[50]="phylfreqfile";
	char id[10];

	sprintf (id, "%d", gen);
	strcat(outname,id);
	strcat(outname,".txt");
	printf("Crea fichero %s\n",outname);
	fphylF = fopen( outname, "w" );
	memset(outname,'\0',strlen(outname)); // limpia el array outname, se podría limpiar solo hasta el número a cambiar..
	strcat(outname,"phylfreqfile"); // reinicializa el array con el nombre base, no haría falta si se hubiera limpiado solo la parte que cambia

//	fprintf (fphylF,"gen=%d\n", gen);

	fprintf (fphylF, "%d	%d\n", 4, NCRO*NLOCI-2);

	for (i=1; i<=(NCRO*NLOCI-2); i++)	fprintf (fphylF, "2 ");
	fprintf (fphylF, "\n");

	for (i=1; i<=pop; i++)
	{
		if ((i == 1)||(i == 2)||(i == pop-1)||(i == pop))
		{

			if (i==1) 		fprintf (fphylF, "RB1       ");
			else if (i==2) 	fprintf (fphylF, "SU1       ");
			else if (i==pop-1) 	fprintf (fphylF, "RB2       ");
			else if (i==pop) 	fprintf (fphylF, "SU2       ");

			if (i==1)
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			if ( (k==0) && (m==0) )	/****/;
			else fprintf (fphylF,"%f ", freqRB1[k][m]);
						
			if (i==2)
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			if ( (k==0) && (m==0) )	/****/;
			else fprintf (fphylF,"%f ", freqSU1[k][m]);
		
			if (i==pop-1)
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			if ( (k==0) && (m==0) )	/****/;
			else fprintf (fphylF,"%f ", freqRB2[k][m]);
		
			if (i==pop)
			for (k=0; k<NCRO; k++)
			for (m=0; m<NLOCI; m++)
			if ( (k==0) && (m==0) )	/****/;
			else fprintf (fphylF,"%f ", freqSU2[k][m]);

			fprintf (fphylF,"\n");
		}
	}

	return(0);
}

/* ************************************************************************** */

output_file_gen0()
{
	// Print genotypes at generation 0 of pop/2 subpopulations of RB, and pop/2 of SU

	for (i=1; i<=pop; i++)
	{
		if ((i == 1)||(i == 2)||(i == pop-1)||(i == pop))
		{
			for (l=1; l<=NSAM; l++)
			{
				/************** Population sampled **************/

				if (i==1) 		fprintf (fdat0, "P1.%d", l);
				else if (i==2) 	fprintf (fdat0, "P2.%d", l);
				else if (i==pop-1) 	fprintf (fdat0, "P3.%d", l);
				else if (i==pop) 	fprintf (fdat0, "P4.%d", l);

				/************** Genotype RB (0) or SU (1) **************/

				if ((gm[i][l][0] & RM[0]) == RM[0])	fprintf (fdat0, "SU");
				else						fprintf (fdat0, "RB");

				/************** Original haplotype RB (0) or SU (1) **************/

				if ((gm[i][l][0] & RM[1]) == RM[1])	fprintf (fdat0, "hSU     ");
				else						fprintf (fdat0, "hRB     ");

				for (k=0; k<NCRO; k++)
				for (m=0; m<NLOCI; m++)
				{
//					if ( (k==0) && (m==0) )	/*******/;
//					else
//					{
						if ((gm[i][l][k] & RM[m]) == RM[m])	fprintf (fdat0,"1 ");
						else						fprintf (fdat0,"0 ");
//					}
				}
				fprintf (fdat0,"\n");
			}
		}
	}

	return(0);
}

/* ************************************************************************** */

finding_haplotypes()
{
	// Keep haplotypes of extreme populations only

	for (i=1; i<=pop; i++)
	{
		if ((i == 1)||(i == 2)||(i == pop-1)||(i == pop))
		{
			for (l=1; l<=NSAM; l++)
			{
				n ++;
				for (k=0; k<NCRO; k++)
				{
					haplotype[n][k] = gm[i][l][k];
				}
			}
		}
	}

	// Remove the first locus
	for (i=1; i<=n; i++)
	{
		haplotype[i][0] = ( haplotype[i][0] & (~RM[0]) );
	}

	// Find different haplotypes
	for (i=1; i<=n; i++)
	{
		if (typehaplotype[i] == 0)
		{
			type ++;
			typehaplotype[i] = type;
		}
		for (j=i+1; j<=n; j++)
		{
			if (typehaplotype[j] == 0)
			{
				for (k=0; k<NCRO; k++)	if (haplotype[j][k] != haplotype[i][k])	goto label;
				typehaplotype[j] = type;
			}
			label: /***/;		
		}
	}
	nhap = type;

	return(0);
}

/* ************************************************************************** */

output_file()
{
	// Print genotypes at last generation of extreme subpopulations

	for (i=1; i<=pop; i++)
	{
		if ((i == 1)||(i == 2)||(i == pop-1)||(i == pop))
		{
			for (l=1; l<=NSAM; l++)
			{
				h ++;

				/************** Population sampled **************/

				if (i==1) 		fprintf (fdat, "P1.%d", l);
				else if (i==2) 	fprintf (fdat, "P2.%d", l);
				else if (i==pop-1) 	fprintf (fdat, "P3.%d", l);
				else if (i==pop) 	fprintf (fdat, "P4.%d", l);

				/************** Genotype RB (0) or SU (1) **************/

				if ((gm[i][l][0] & RM[0]) == RM[0])	fprintf (fdat, "SU");
				else						fprintf (fdat, "RB");

				/************** Original haplotype RB (0) or SU (1) **************/

				if ((gm[i][l][0] & RM[1]) == RM[1])	fprintf (fdat, "hSU     ");
				else						fprintf (fdat, "hRB     ");

				for (k=0; k<NCRO; k++)
				for (m=0; m<NLOCI; m++)
				{
//					if ( (k==0) && (m==0) )  /*****/;
//					else
//					{
						if ((gm[i][l][k] & RM[m]) == RM[m])	fprintf (fdat,"1  ");
						else						fprintf (fdat,"0  ");
//					}
				}
				fprintf (fdat,"    %d\n", typehaplotype[h]);
			}
		}
	}

	return(0);
}

/* ************************************************************************** */

