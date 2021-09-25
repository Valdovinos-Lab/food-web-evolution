

//*********************************************************************************************
//****************************************OVERVIEW*********************************************
/* 	
This program calculates evolutionary food webs with invasions. 
This is the program version that we (JR Morris, KT Allhoff, F Valdovinos) used for our manuscript 
"Strange invaders increase disturbance and promote generalists in an evolving food web" submitted to Scientific Reports in 2021. 

The output is written into several files: 				
- before simulation start: "parameters.txt" summarises all necessary parameter values to reproduce a given simulatino run
- data for producing time series: "EvoInv_XXX.txt" (output every tmod time step, see output_timeseries)
- data to produce heatmaps (averaged throughout the simulation, output at the end of the simulation, see output_heatmap)
- data for turnover and niche overlap analysis (output at every tmod_turnover time step, see output_turnover)

*************************************************************************
Compile this program in your terminal using the following command:
	gcc MODEL_CLEAN.c -std=c99 -lm -lgsl -lgslcblas -ffast-math -o MODEL_CLEAN.out

Execute it locally via: 
	./MODEL_CLEAN.out [arguments]
	
	For example: ./MODEL_CLEAN.out res_CLEAN 12345 0.05 1 0.5 0.2 0.4 0.4 0.3 1
	with [res_CLEAN = name of subfolder for storing the results] and [12345 = seed for random numbers].
	Standard parameter set is: c0f=0.05, inlfow=1, loss=0.5, p=0.2, inv_z=0.4, s0=0.4, lowbound=b=0.3, sim_id = 1

Debug the code in your terminal using valgrind: 
	1) Rerun the compiler with additional flag -g
	2) Execute locally via: valgrind --leak-check=full -v ./MODEL_CLEAN.out [arguments]
*/


/*  ************#####JRM EDITS#####************** ORIGINAL WITH SPECIES HISTORY DATA 
 This version outputs traits and lifespan data for ALL species generated during simulation
 NOTE: this version edited by Michael Egan
 *****************************************************************************************************
 */



//****************************************************************************************************
//****************************************LOAD GSL PACKAGES*******************************************

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>			    // random number generator
#include <gsl/gsl_randist.h>
#include <gsl/gsl_odeiv2.h>			  // solver for the differential equations
#include <gsl/gsl_blas.h>			    // linear algebra with vectors and matrices
#include <gsl/gsl_linalg.h>			  // linear algebra for calculations of trophic positions

//*****************************************************************************************************
//************************************ASSIGN GLOBAL VARIABLES*******************************************

// ******************** global variables (some of them will be overwritten via input-variables!)
// ***** for the mutations
double s0=0.30;					// starting value for width of feeding range s
double lowbound=0.1;		             		 // lower boundary for the feeding ranges
double mut_z=0.1;				          // inheritance factor for the mutations
double inv_z;				              // inheritance factor for invaders
double p=0;				  	     // frequency of invasions: (number of invasions)/(number of mutations)
double delta_Tm=100;	     	            	 	 // A new species is added after every delta_Tm time steps
double extinct = 0.00000001;				// extinction threshold = initial biomass of new morphs     

// ***** interaction parameters:
double a0=1.0;					    	// attack rate factor 
double efficiency=0.85;		        	// efficiency factor 
double h0=0.4;					    	// handling time factor
double d0=0.3;				        	// mortality rate factor
double c0f = 0.05;			        	// competition parameter (for food)
double c0i = 0.00;				        // competition parameter (intraspecific offset)

// ***** resource:
double inflow=10; 				        // Input of inorganic nutrient
double loss=0.1;				  	// loss rate of inorganic nutrient
                                   
// ***** for the EvoVsInv project:
double* heatmap_S;   			        	// number of morphs (data for heatmap)
double* heatmap_B;				        // total biomass (data for heatmap) 		
int S_old_tmod_turnover=0;                   		// data for turnover analysis: number of species tmod_turnover time steps ago
     
// ***** program parameter:
int Smax = 100;			            // maximum species number (termination condition 1)
double controll1 = 1e-10;			    // relative accuracy of solver
double controll2 = 1e-12;			    // absolute accuracy of solver
int seed = 1234567;				    // seed for random numbers
char path[80] = "";                               // location of folder for output files
int sim_id=5;
double tend =  25000000;			    // maximum number of time steps (termination condition 2)
double tmod = 	  50000;			    // basic output for time series after everey tmod time step
double tmod_turnover=     10000; 			// output for turnover after every tmod_turnover time step 


// ******************** global pointers to important vectors and matrices: 
double* m_vec;   			// body masses of all species 
double* m_vec_old; 		// memory for turnover analysis
double* s_vec;					// feeding ranges of all species 
double* c_vec;					// feeding center of all species 
double* TL; 				  	// trophic positions of all species 
double** a_mat;					// attack rate matrix
double** c_mat;					// competition matrix


// #####****************************************************************#####
// #####********** species lifespan record (CODE EDITED BY MICHAEL EGAN)#####
enum SpeciesType {
	Ancestor = 0,
	Mutant,
	Invader
};

const char *speciesTypeNames[] = {"Ancestor", "Mutant", "Invader"};

typedef struct {
	int id;
	enum SpeciesType speciesType;
	double bodyMass;
	double feedingCenter;
	double feedingRange;
	double introduction;
	double extinction;
} Lifespan;

int *speciesId;
Lifespan *speciesLifespanRecord;

void output_speciesLifespanRecord(FILE *res_lifespan, Lifespan *speciesLifeSpanRecord, size_t nLifespanRecords);
// #####****************************************************************#####

//*****************************************************************************************************
//****************************************ASSIGN FUNCTIONS*********************************************

// ******************** functions
// ***** calculates one simulation:
double *web_calc(gsl_rng *r,FILE *res_res,FILE *res_pop,FILE *res_tp,FILE *res_m,FILE *res_c,FILE *res_s,FILE *res_heatmap, FILE *res_turnover, FILE *res_lifespan);

// ***** Introducing new species:
void ChooseParent(gsl_rng *r, int S, int* m, double y[]);
void ChooseParent2(gsl_rng *r, int S, int* m, double y[]);
void CreateNewSpecies(gsl_rng *r, int m, double z, double* mut_m, double* mut_c, double* mut_s);

// ***** evolve until next mutation (or immigration) event: 
int evolve_until_mut(double* t, double Tm, int S, FILE *res_res, FILE *res_pop,FILE *res_tp, FILE *res_m, FILE *res_c,FILE *res_s,double y[]);

// ***** calculate competition between two species:
double competition(double ci, double cj, double si, double sj);		
	
// ***** the function that contains the population dynamics:
int dynamics(double t, const double y[], double ydot[], void *params);	
	
// ***** calculates biomass flow from prey j to predator i:
double gij(const double y[], int S, int predator, int prey);
	
// ***** calculates the trophic positions of all species: 
double *tropos(const double y[], int S, FILE *res_tp);	

// ***** fills the attack rate and competition matrices: 
void calc_mat(int S);					
			
// ***** removes extinct species from the system: 
int Extinct(int S, double y[], int **removedSpeciesId, size_t *nRemovedSpecies);

// ***** calculates the functional diversity: 
// (as in Allhoff & Drossel (2016): "Biodiversity and ecosystem functioning in evolving food webs" (Phil Trans Roy Soc B))
double functional_diversity(int S, const double y[]);

// ***** output functions: (not all of them are always needed, so make sure to comment them in/out as you need them)
int output_timseries(FILE *res_res, FILE *res_pop, FILE *res_tp, FILE *res_m, FILE *res_c, FILE *res_s, const double y[], int S, double t);
int output_heatmap(FILE *res_heatmap);  			
int output_turnover(double t, int S, FILE *res_turnover, const double y[]);	


//*********************************************************************************************
//*********************************************************************************************
//***************************************MAIN FUNCTION*****************************************
//*********************************************************************************************
//*********************************************************************************************

int main(int argc, char *argv[])
{
	strcat(path, argv[1]);	              // first get the path to the output folder for this simulation (from global command when initiating simulation)

                                             // then get all other model parameters that are not yet defined above (this overwrites default values!):

        sscanf(argv[2], "%d", &seed);        // get seed for random numbers
        sscanf(argv[3], "%lf", &c0f);        // get interference competition parameter value
        sscanf(argv[4], "%lf", &inflow);     // get base resource inflow rate parameter
        sscanf(argv[5], "%lf", &loss);       // get base resource outflow rate parameter
        sscanf(argv[6], "%lf", &p);          // get probability of invasion/mutation parameter
        sscanf(argv[7], "%lf", &inv_z);      // get invader strangeness parameter
        sscanf(argv[8], "%lf", &s0);         // get starting feeding range (ANCESTOR) parameter
        sscanf(argv[9], "%lf", &lowbound);   // get feeding range lowerbound parameter
        sscanf(argv[10], "%d", &sim_id);     // get simulation id parameter


// ******************** create output files for simulation data:

	FILE *res_pop;								
	char format_pop[] = "/EvoInv_pop.txt";
	char filename_pop[80] = "";
	strcat(filename_pop, path);
	strcat(filename_pop, format_pop);
	res_pop = fopen(filename_pop, "w");

	FILE *res_res;								
	char format_res[] = "/EvoInv_res.txt";
	char filename_res[80] = "";
	strcat(filename_res, path);
	strcat(filename_res, format_res);
	res_res = fopen(filename_res, "w");

	FILE *res_tp;								
	char format_tp[] = "/EvoInv_tp.txt";
	char filename_tp[80] = "";
	strcat(filename_tp, path);
	strcat(filename_tp, format_tp);
	res_tp = fopen(filename_tp, "w");

	FILE *res_m;								
	char format_m[] = "/EvoInv_m.txt";
	char filename_m[80] = "";
	strcat(filename_m, path);
	strcat(filename_m, format_m);
	res_m = fopen(filename_m, "w");

	FILE *res_c;								
	char format_c[] = "/EvoInv_c.txt";
	char filename_c[80] = "";
	strcat(filename_c, path);
	strcat(filename_c, format_c);
	res_c = fopen(filename_c, "w");
	
	FILE *res_s;								
	char format_s[] = "/EvoInv_s.txt";
	char filename_s[80] = "";
	strcat(filename_s, path);
	strcat(filename_s, format_s);
	res_s = fopen(filename_s, "w");

	FILE *res_heatmap;								
	char format_heatmap[] = "%s/heatmap_%d.txt";
	char filename_heatmap[80] = "";
	sprintf(filename_heatmap, format_heatmap, path, sim_id);
	res_heatmap = fopen(filename_heatmap, "w");
	
	FILE *res_turnover;								
	char format_turnover[] = "%s/turnover.txt";
	char filename_turnover[80] = "";
	sprintf(filename_turnover, format_turnover, path);
	res_turnover = fopen(filename_turnover, "w");

	FILE *res_lifespan;
	char format_lifespan[] = "%s/lifespans.csv";
	char filename_lifespan[80] = "";
	sprintf(filename_lifespan, format_lifespan, path);
	res_lifespan = fopen(filename_lifespan, "w");

// ******************** Initialize random number generator
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T=gsl_rng_default;   			          // default random number generator (so called mt19937)
	gsl_rng_default_seed=seed;		      // a predefined seed-value makes simulations reproducable!
	r=gsl_rng_alloc(T);


// ********************  before doing anything else: print out all parameters for this simulation!
// ***** (So, if the simulation crashes, at least you know how to repeat it...)

	FILE *res_par;								
	char format_par[] = "/parameters.txt";
	char filename_par[80] = "";
	strcat(filename_par, path);
	strcat(filename_par, format_par);
	res_par = fopen(filename_par, "w");

	fprintf(res_par, "Mutation parameter: \ndelta_Tm: %lf\nextinct: %1.8f\n",delta_Tm,extinct);
	fprintf(res_par, "s0= %lf\nlb: %1.2f\n", s0, lowbound);
	fprintf(res_par, "mut_z: %lf\ninv_z: %lf\np: %1.2f\n\n",mut_z, inv_z, p);	
	fprintf(res_par, "interaction factors:\na0= %lf\nefficiency= %lf\n",a0,efficiency);
	fprintf(res_par, "h0= %lf\nd0= %lf\nc0f= %lf\nc0i= %lf\n\n",h0, d0, c0f,c0i);
	fprintf(res_par, "Resource inflow rate= %lf\nresource loss rate= %lf\n\n",inflow,loss);
	fprintf(res_par, "Program parameter:\nSmax =%d\n",Smax);
	fprintf(res_par, "controll1 =%1.12f\ncontroll2 =%1.12f\nseed:%d\n\n", controll1, controll2, seed);
	fprintf(res_par, "Important points in time:\ntend= %12.0f\ntmod=%lf\ntmod_turnover=%lf\n",tend, tmod, tmod_turnover);
	fclose(res_par);

// ******************** calculate simulation
	web_calc(r, res_res, res_pop, res_tp, res_m, res_c, res_s, res_heatmap, res_turnover, res_lifespan);

	
// ******************** close all files and free memory (THIS IS IMPORTANT TO AVOID SEGMENTATION FAULTS!!)

	fclose(res_pop);
	fclose(res_res);
	fclose(res_tp);
	fclose(res_m);
	fclose(res_c);
	fclose(res_s);
	fclose(res_heatmap);
	fclose(res_turnover);
	fclose(res_lifespan);
	
	gsl_rng_free(r);
	
	return(0);
}


//**************************************EVOLVE NETWORK******************************************
//**********************************************************************************************

double *web_calc(gsl_rng *r,FILE *res_res,FILE *res_pop,FILE *res_tp,FILE *res_m,FILE *res_c,FILE *res_s,FILE *res_heatmap, FILE *res_turnover, FILE *res_lifespan)
{
// ******************** Preparations...
	int S=1;							// current species number
	int i,j;						
 	double temp;				

	double *y=(double *) calloc((Smax+1),sizeof(double));	// memory for population sizes (format: resource,S1,S2,...)

	m_vec = (double *) calloc(Smax+1, sizeof(double));   	// body masses 
	m_vec_old = (double *) calloc(Smax+1, sizeof(double));   	// memory for turnover analysis 
	s_vec = (double *) calloc(Smax+1, sizeof(double));		// feeding ranges
	c_vec = (double *) calloc(Smax+1, sizeof(double));		// feeding center
	TL=(double *) calloc((1+Smax), sizeof(double));		// memory for trophic positions 

	a_mat = (double **) calloc(Smax+1, sizeof(double*));
	for(int i = 0; i <= Smax; i++)
	  a_mat[i] = (double *) calloc(Smax+1, sizeof(double));	// attack rate matrix

	c_mat = (double **) calloc(Smax+1, sizeof(double*));
	for(int i = 0; i <= Smax; i++)
	  c_mat[i] = (double *) calloc((Smax+1), sizeof(double)); 	// competition matrix

  	heatmap_S = (double *) calloc(1000, sizeof(double));   	// number of morphs (data for heatmap) 	(collect 100 data points during this simulation)
	heatmap_B = (double *) calloc(1000, sizeof(double));		// total biomass (data for heatmap)		(collect 100 data points during this simulation)

	// species lifespan record
	speciesId = (int *) calloc(Smax+1, sizeof(int));

	size_t nLifespanRecords = tend / delta_Tm;
	speciesLifespanRecord = (Lifespan *) calloc(nLifespanRecords, sizeof(Lifespan));

	// ********** Initial situation: One species + one resource 
	m_vec[0] = 1.0;				// resource's body mass
	speciesId[1] = 0; 				// ancestor species id number
	m_vec[1] = 100;				// ancestor species' body mass
	m_vec_old[0] = 1.0;				// memory for turnover calculation
	m_vec_old[1] = 100;				// memory for turnover calculation
	s_vec[1] = s0;					// ancestor species' feeding range
	c_vec[1] = 1.0;				// ancestor species' feeding centre
	y[0] = inflow/loss;				// resource biomass density = carrying capacity
	y[1] = extinct; 				// ancestor biomass density = extinction threshold

	calc_mat(S);					// attack rates + competition 

	// species lifespan record
	// ancestor
	speciesLifespanRecord[0].id = 0;
	speciesLifespanRecord[0].speciesType = Ancestor;
	speciesLifespanRecord[0].bodyMass = 100;
	speciesLifespanRecord[0].feedingCenter = 1.0;
	speciesLifespanRecord[0].feedingRange = s0;
	speciesLifespanRecord[0].introduction = 0;
	speciesLifespanRecord[0].extinction = INFINITY;

	// ********** for the mutation events
	double mut_m;					// mutant's body mass
	double mut_c;					// mutant's feeding center
	double mut_s;					// mutant's feeding range
	double mut_y;					// mutant's initial biomass density
	double z;							
	int m;						// index of the parent species 
	double Tm; 					// next mutation event 

	// species lifespan record
	int aSpeciesId = 0;

	// ********** Population dynamcis of ancestor species + resource until first mutation event:
	double t=0.0;					// Start time =0
	Tm = delta_Tm + t;				// End time =Tm
	S=evolve_until_mut(&t,Tm,S,res_res,res_pop,res_tp,res_m,res_c,res_s,y);


// ******************** Until either tend or Smax is reached: Add new species, calculate pop. dynamcs, add new species, ... 
	while (t<tend && S<Smax && S>0)
	{
		// ***** Before we go on: Save the data for the heatmap!
		int temp2=t;
		if ( fmod(temp2, tend/1000)==0.0  )

		{
			heatmap_S[(int) (temp2/(tend/1000)) -1] = S;				// data for number of morphs
				
			for(i=1;i<=S;i++)
				heatmap_B[(int) (temp2/(tend/1000)) -1] += y[i];	
		}
		
		
		// ***** ...and analyse the turnover!
		if( fmod(t,tmod_turnover)==0.0  )					
		{
			output_turnover(t, S, res_turnover,y);		// analyse turnover after every tmod_turnover time units
			S_old_tmod_turnover=S; 							      // remember number of species for next time
			
			for (int l=0; l<= Smax; l++)			        // remember species body masses for next time
				m_vec_old[l] = m_vec[l];
		}		
		
		enum SpeciesType aSpeciesType;
		
		// ***** Now start species introduction
		temp = gsl_rng_uniform_pos(r); 			
		if (temp>=p)			// ***** A mutation event takes place!
		{
			aSpeciesType = Mutant;

			z=mut_z;
			int success=0;
			while (success ==0)
			{
				ChooseParent(r, S, &m, y);		// select parent population proportional to individual density
				mut_y = extinct;			// initial biomass of the mutant

				if (y[m]- mut_y > 0) 			// mutants biomass is taken from the parent species
				{			 
					y[m]=y[m]-mut_y;
					success=1;
				}

				else if (y[m]- mut_y <= 0 )		// parent population is too small, mutation can not take place
				{					// (paranoia - this situaton should actually never occur on isolated patches!)
					mut_y=0;
					success=0;
				}
			}

			CreateNewSpecies(r, m, z, &mut_m, &mut_c, &mut_s);
		}

		else if (temp<p)		// ***** An invasion event takes place!
		{
			aSpeciesType = Invader;

			z = inv_z;
			ChooseParent2(r, S, &m, y);			// select parent population at random
			mut_y = extinct;				// initial biomass of the mutant is just added to the system

			CreateNewSpecies(r, m, z, &mut_m, &mut_c, &mut_s);
		}

		// add new species to lifespan record
		aSpeciesId++;

		int speciesLifespanIndex = aSpeciesId;
		speciesLifespanRecord[speciesLifespanIndex].id = aSpeciesId;
		speciesLifespanRecord[speciesLifespanIndex].speciesType = aSpeciesType;
		speciesLifespanRecord[speciesLifespanIndex].bodyMass = mut_m;
		speciesLifespanRecord[speciesLifespanIndex].feedingCenter = mut_c;
		speciesLifespanRecord[speciesLifespanIndex].feedingRange = mut_s;
		speciesLifespanRecord[speciesLifespanIndex].introduction = Tm;
		speciesLifespanRecord[speciesLifespanIndex].extinction = INFINITY;

		// ********** include the new species into the system (last position in all vectors)
		S=S+1;							// increase number of species	
		speciesId[S] = aSpeciesId;
		m_vec[S] = mut_m;					// body mass
		c_vec[S] = mut_c;					// feeding center
		s_vec[S] = mut_s;					// feeding range
		y[S]=mut_y;				 		// initial population density

		// ******************** calculate the population dynamics until the next mutation event takes place:
		int *removedSpeciesId;
		size_t nRemovedSpecies;

		S=Extinct(S,y, &removedSpeciesId, &nRemovedSpecies);						// exclude parent species if necessary
		assert(removedSpeciesId != NULL);

		for (int iRemovedSpecies = 0; iRemovedSpecies < nRemovedSpecies; iRemovedSpecies++)
		{
			int aRemovedSpeciesId = removedSpeciesId[iRemovedSpecies];

			int aRemovedSpeciesIndex = aRemovedSpeciesId;
			assert(speciesLifespanRecord[aRemovedSpeciesIndex].id == aRemovedSpeciesId);
			assert(speciesLifespanRecord[aRemovedSpeciesIndex].extinction == INFINITY);

			speciesLifespanRecord[aRemovedSpeciesIndex].extinction = Tm;
		}

		free(removedSpeciesId);

		calc_mat(S);						// up-date interaction matrices
		
		Tm = delta_Tm + t;					// Time of next species addition 	
		S=evolve_until_mut(&t,Tm,S,res_res,res_pop,res_tp,res_m,res_c, res_s,y);

	} // end of simulation (either the maximum time or the maximum species number is reached)

// ******************** final output
	output_timseries(res_res,res_pop,res_tp,res_m,res_c, res_s,y,S,t);
	output_turnover(t,S, res_turnover,y);

		// ***** Last time: Save the data for heatmap!
		if ( t==tend)
		{
			heatmap_S[999] = S;				
							
			for(i=1;i<=S;i++)
			 	heatmap_B[999]=heatmap_B[999] + y[i];	

		}  // Now we have two vectors: one for the number of morphs and another for the total biomass. Both contain 1000 data points.

//	output_heatmap(res_heatmap); 		// Calculate mean and CV for these two vectors and create one single line of output!

	output_speciesLifespanRecord(res_lifespan, speciesLifespanRecord, nLifespanRecords);

// ******************** free memory (THIS IS IMPORTANT TO AVOID SEGMENTATION FAULTS!!)
	free(y);
	free(TL);
	free(c_vec);
	free(s_vec);
	free(m_vec);
	free(m_vec_old);
	
	free(heatmap_S);
	free(heatmap_B);

	for(int i = 0; i <= Smax; i++)
	  free(a_mat[i]);
	free(a_mat);
	for(int i = 0; i <= Smax; i++)
	  free(c_mat[i]);
	free(c_mat);

	// species lifespan
	free(speciesId);
	free(speciesLifespanRecord);

	return 0;
}


//************************************************************************************************
// WHO is the next parent species, and WHERE does the mutation take place?
//************************************************************************************************
void ChooseParent(gsl_rng *r, int S, int* m, double y[])
{
	// ******************** the mutation probability of a population is proportional to its individual density!
	double AnzahlIndividuen=0;						
	double temp;
	double *Individuen=(double *) calloc((Smax+1),sizeof(double));	


	Individuen[0]= 0;
 	for (int i = 1; i <= S; i++)
	{
		Individuen[i]= y[i] / m_vec[i];	// individual density per population
		AnzahlIndividuen= AnzahlIndividuen + Individuen[i];	// number of all individuals in the system
	}


	temp = gsl_rng_uniform_pos(r) * AnzahlIndividuen;		// choose a mutating individual and...
	double Count_Individuen=0;
	int i=0;

	while (temp>Count_Individuen)					// ...find out to which population it belongs!
	{
		i+=1;
		Count_Individuen = Count_Individuen + Individuen[i];
	}

	*m=i;  // return parent species

	free(Individuen);
  	return;
}

//************************************************************************************************
// WHO is the next "parent" species for invaders? 
//************************************************************************************************
void ChooseParent2(gsl_rng *r, int S, int* m, double y[])
{
	// ******************** the invasion probability is just random (instead of being proportional to the individual density)

	double temp;
	int i=0;

	temp = gsl_rng_uniform_pos(r) * S;		// choose a random number between 0 and S (=number of species)

	for (i=0; temp>i; i++)
		;

	*m=i;						// return parent species

  	return;
}

//************************************************************************************************
// Create the new morph as a modification of its parent
//************************************************************************************************
void CreateNewSpecies(gsl_rng *r, int m, double z, double* mut_m, double* mut_c, double* mut_s)
{
	double temp;

	// ********** body masses
	temp = gsl_ran_gaussian(r, z);
	*mut_m = m_vec[m] * pow(10,temp);

	while(*mut_m < 1)								// (all species should be bigger than the external resource) 
	{
		temp = gsl_ran_gaussian(r, z);
		*mut_m = m_vec[m] * pow(10,temp);		
	}

	// ********** feeding ranges
	temp = gsl_ran_gaussian(r, z);
	*mut_s = s_vec[m] + temp;

	while(*mut_s < lowbound)							// (no extreme specialists alowed) 
	{
		temp = gsl_ran_gaussian(r, z);
		*mut_s = s_vec[m] + temp;
	}		

	// ********** feeding centers
	temp = gsl_ran_gaussian(r, z);			
	*mut_c = c_vec[m] * pow(10,temp);

	return;
}

//************************************************************************************************
// Calculate the population dynamics until the next species emerges
//************************************************************************************************
int evolve_until_mut(double* t, double Tm, int S, FILE *res_res,FILE *res_pop,FILE *res_tp, FILE *res_m, FILE *res_c, FILE *res_s,double y[])
{
	double t2=*t;
	double h=1e-2;

	// The solver expects a vector that contains all parameters needed to solve the differential equations.
	// Here, only the species number is needed, because all other parameters are global variables.
	double *params=(double *) calloc(1,sizeof(double));
	params[0]=S;					

	const gsl_odeiv2_step_type *Solv=gsl_odeiv2_step_rkf45;			// use the Runge-Kutta-method to solve...
	gsl_odeiv2_step *st=gsl_odeiv2_step_alloc(Solv,(1+S));			// ...a set of 1+S differential equations 
	gsl_odeiv2_evolve *e=gsl_odeiv2_evolve_alloc((1+S));
	gsl_odeiv2_control *con=gsl_odeiv2_control_y_new(controll1,controll2); 	// absolute error, relative error
	gsl_odeiv2_system sys={dynamics,NULL,(1+S),params};				// sys contains the dynamics of the system!

	// ******************** calculate population dynamics until the next species emerges:
	while(t2<Tm && t2<=tend && S<Smax && S>0)
	{
		// ***** output of current situation
		if( fmod(t2,tmod)==0.0  )							
			output_timseries(res_res,res_pop,res_tp,res_m,res_c,res_s,y,S,t2);

		// ***** solve differential equaton 
		int status=gsl_odeiv2_evolve_apply(e,con,st,&sys,&t2,Tm,&h,y); 	

		if(status!= GSL_SUCCESS)
			break;
	}

	assert(t2==Tm);	// The solver uses adaptive step sizes...
				// ...The end time, t2, is thus slightly above Tm and has to be set back to integer values.
	*t=t2;			// return current time
	
	int *removedSpeciesId;
	size_t nRemovedSpecies;

	S=Extinct(S,y, &removedSpeciesId, &nRemovedSpecies);		// In case some species went extinct
	assert(removedSpeciesId != NULL);

	for (int iRemovedSpecies = 0; iRemovedSpecies < nRemovedSpecies; iRemovedSpecies++)
	{
		int aRemovedSpeciesId = removedSpeciesId[iRemovedSpecies];

		int aRemovedSpeciesIndex = aRemovedSpeciesId;
		assert(speciesLifespanRecord[aRemovedSpeciesIndex].id == aRemovedSpeciesId);
		assert(speciesLifespanRecord[aRemovedSpeciesIndex].extinction == INFINITY);

		speciesLifespanRecord[aRemovedSpeciesIndex].extinction = Tm;
	}

	free(removedSpeciesId);

	// ******************** free memory (THIS IS IMPORTANT TO AVOID SEGMENTATION FAULTS!!)
	gsl_odeiv2_step_free(st);
	gsl_odeiv2_control_free(con);
	gsl_odeiv2_evolve_free(e);
	free(params);

	return S;
}



//************************************************************************************************
// Remove extinct species
//************************************************************************************************
int Extinct(int S, double y[], int **removedSpeciesId, size_t *nRemovedSpecies)
{
	int i,j,k;
	double temp=0;

	int iRemovedSpecies = 0;
	int *removedSpecies = calloc(S, sizeof(int));
	assert(removedSpecies != NULL);

	for (i=1; i<S+1; i++)
	{
		temp=0;

		if (y[i] < extinct)			// Is species i still alive?
			y[i] = 0;
		if (y[i] >= extinct)		
			temp = 1;			// If yes, then the species i is saved

		if (temp==0)				// If not, then remove this species from the system!
		{
			removedSpecies[iRemovedSpecies] = speciesId[i];
			iRemovedSpecies++;

			for (j=i; j<S; j++)
				{
					speciesId[j] = speciesId[j+1];
					m_vec[j] = m_vec[j+1];
					s_vec[j] = s_vec[j+1];
					c_vec[j] = c_vec[j+1];
				}

			speciesId[S] = 0;
			m_vec[S] = 0.0;			
			s_vec[S] = 0.0;
			c_vec[S] = 0.0;

			for (j=i; j<(S+1); j++)
				y[j]=y[j+1];
			
			y[S]=0;

			S=S-1;					// reduce species number
			i=i-1;
		}
	}

	*removedSpeciesId = removedSpecies;
	*nRemovedSpecies = iRemovedSpecies;

	return(S);
}


//************************************************************************************************
// Competition between species i and j = overlap of Gaussian feeding kernels 
//************************************************************************************************

double competition(double ci, double cj, double si, double sj)
{
  double temp1 = log10(ci)-log10(cj);
  double temp2 = si*si + sj*sj;
  return exp(-0.5*(temp1*temp1)/temp2) / sqrt(2.0*3.14159265358979323846*temp2);
}

//************************************************************************************************
// Calculate the attack rate and competition matrices
//************************************************************************************************
void calc_mat(int S)
{
	  const double pi = 3.14159265358979323846;
	  
	  // ********** Attack rate (Gaussian feeding kernels)
	  for(int j = 0; j <= S; j++)
	  	a_mat[0][j] = 0.0;
	  for(int i = 1; i <= S; i++)
	  {
		double temp1 = c_vec[i];
		for(int j = 0; j <= S; j++)
		{
			double temp2 = m_vec[j];
			double temp3 = log10(temp1/temp2);
			temp3 = (temp3*temp3) / (2.0*s_vec[i]*s_vec[i]);
			temp3 = a0/sqrt(2.0*pi)*exp(-temp3);
			a_mat[i][j] = temp3 *pow(m_vec[i],0.75)/s_vec[i];
		}
	  }

//   printf("Aij =");
//   for(int i = 0; i <= S; i++)
//   {
//     printf("\n");
//     for(int j = 0; j <= S; j++)
//       printf("%1.8f \t", a_mat[i][j]);
//   }
//   printf("\n");

	  // ********** Competition (overlap of Gaussian feeding kernels)
	  for(int j = 0; j <= S; j++)
	  	c_mat[0][j] = 0.0;

	  for(int i = 1; i <= S; i++)
	  {
		double normi = competition(c_vec[i], c_vec[i], s_vec[i], s_vec[i]);
		c_mat[i][0] = 0.0;
		for(int j = 1; j <= S; j++)
			c_mat[i][j] = c0f/normi*competition(c_vec[i], c_vec[j], s_vec[i], s_vec[j]);

		c_mat[i][i] = c0f + c0i;			// special case: intraspecific offset
	  }

//   printf("Cij =");
//   for(int i = 0; i <= S; i++)
//   {
//     printf("\n");
//     for(int j = 0; j <= S; j++)
//       printf("%1.8f \t", c_mat[i][j]);
//   }
//   printf("\n");

  	return;
}


//************************************************************************************************
// Population dynamics
//************************************************************************************************
int dynamics(double t, const double y[], double ydot[], void *params)
{
	double* pparams = (double *)params;
	int S = (int)*pparams;
	double emig_p0, emig_p1;
	
	double d0_k, inflow_k, loss_k;
 	
	// ***** denominator for functional response
	double denom[S+1];
	for(int i = 0; i <= S; i++)
	{
		denom[i] = 1.0;
		
		// 	for Holling II functional response:
		double hi = h0 * pow(m_vec[i],-0.75); 
		for(int j = 0; j <= S; j++)
			denom[i] += (hi * a_mat[i][j]) * y[j] ; 
		denom[i] *= m_vec[i];
	}

	// ***** resource dynamics 
   
    	ydot[0] = inflow - y[0] * loss;     // Chemostat growth

	for(int j = 1; j <= S; j++)
		ydot[0] -= a_mat[j][0]*y[j]*y[0]/denom[j];

    // ***** dynamics of the consumer species
   	for(int i = 1; i <= S; i++)
   	{
	   	ydot[i] = (efficiency*a_mat[i][0]/denom[i]) * y[0];  		
	   	for(int j = 1; j <= S; j++)
	    		ydot[i] += ( efficiency*a_mat[i][j]/denom[i] - a_mat[j][i]/denom[j] - c_mat[i][j] ) * y[j];

		 	ydot[i] -= d0 * pow(m_vec[i],-0.25);
			ydot[i] *= y[i];
	 	}

	return GSL_SUCCESS;
}


//*******************************************************************************************
// Calculate the consumption rate Gij 
// (is needed for the calculation of the trophic positions)
//*******************************************************************************************
double gij(const double y[], int S, int predator, int prey)
{
	double gij, denom;
	int j;

		denom = 1.0;
		double hi = h0 * pow(m_vec[predator],-0.75);

		for(j=0; j <= S; j++)
			denom += (hi* a_mat[predator][j]) *  y[j]; 

		denom = denom * m_vec[predator];

	gij = a_mat[predator][prey];
	gij = gij * y[prey];
	gij = gij / denom;

	return (gij);
}


//*******************************************************************************************
// Calculate the trophic positions (via the weighted trophic positions of the prey species)
//*******************************************************************************************
double *tropos(const double y[], int S, FILE *res_tp)
{
// ******************** We have to solve a system of the form (Matrix)Ta*(Vector)xa=(Vector)b
	int i,j;
	double norm;
	int sig=1;

	gsl_vector *b=gsl_vector_calloc(1+S);			// Inhomogeneity
	gsl_vector *xa=gsl_vector_calloc(1+S);		// The variables that we look for
	gsl_permutation *p=gsl_permutation_calloc(1+S);	// Permutation vector
	gsl_matrix *Ta=gsl_matrix_calloc(1+S,1+S);		// Matrix of the system

// ******************** Fill matrix Ta and Vector b:
	for(i=0;i<=S;i++)
	{
		gsl_matrix_set(Ta,i,0,gij(y,S,i,0)*efficiency);
		for (j=1; j<=S; j++)
			gsl_matrix_set(Ta, i, j, gij(y,S,i,j)*efficiency);
	}

	for(i=0;i<=S;i++)
	{
		gsl_vector_view tempp=gsl_matrix_row(Ta,i);
		norm=gsl_blas_dasum(&tempp.vector);

		if(norm!=0)
			for(j=0;j<=S;j++)
				gsl_matrix_set(Ta,i,j,gsl_matrix_get(Ta,i,j)/norm);

		gsl_matrix_set(Ta,i,i,gsl_matrix_get(Ta,i,i)-1);
	}

	for(i=1;i<=S;i++)
		gsl_vector_set(b,i,-1);

// ******************** Now solve the system!
	gsl_linalg_LU_decomp(Ta,p,&sig);
	gsl_linalg_LU_solve(Ta,p,b,xa);

	for(i=1;i<=S;i++)
		if (y[i]<extinct)
			gsl_vector_set(xa,i,0);		// species below the extinct threshold do not count for the evaluation!

// ******************** Output
	for(i=0;i<=S;i++)
	{
		TL[i]=gsl_vector_get(xa,i);
		 	fprintf(res_tp,"%1.6f \t", TL[i]);
	}
	fprintf(res_tp,"\n");

// ******************** free memory (THIS IS IMPORTANT TO AVOID SEGMENTATION FAULTS!!)
	gsl_matrix_free(Ta);
	gsl_vector_free(xa);
	gsl_permutation_free(p);
	gsl_vector_free(b);
	return 0;
}


//***************************************************************************************
//  Calculate the functional diversity (as in Allhoff&Drossel 2016 PhilTransB)
//***************************************************************************************
double functional_diversity(int S, const double y[])
{
	int i,j;
	double max;				// attack rate maximum 
	double C;				// feeding center of one species
	double M_prey;				// body mass of prey species
	double ais;				// attack rate of predator i on prey s

	double L=10;				// Length of logarithmic body mass axis
	double N=1000;				// number of sampling points
	double l=L/N;				// Length of intervals
	double s;				// current sampling point

	double I=0;				// Value of the integral over the envelope of all feeding kernels

	// ***** We calculate the intergral I via approximation with rectangles:

	for(j=0;j<N;j++)
	{
		max=0;
		s=(j)*l-(2);
		M_prey=s;

		for(i=0;i<(S+1);i++)
		{
			ais=0;
			if (y[i]>extinct) 
			{
				C=log10(c_vec[i]);
				ais = (C-M_prey)*(C-M_prey) / (2*s_vec[i]*s_vec[i]);
				ais= a0 / (s_vec[i]*sqrt(2*3.1415926535))*exp(-ais);
			}

			if(ais>max)
				max=ais;
		}

		I+=max*l;
	}

	return(I);				// I is our measure of functional diversity!
}


//*******************************************************************************************
// output of network properties and other main results 
// after every tmod time step
//*******************************************************************************************
int output_timseries(FILE *res_res,FILE *res_pop,FILE *res_tp, FILE *res_m,FILE *res_c, FILE *res_s,const double y[], int S, double t)
{
	int i,j;
	double *species=(double *) calloc(6,sizeof(double));		// number of species and biomass densities... 
	double *biomass=(double *) calloc(6,sizeof(double));		// ... per trophic level


	for (i=0; i<6; i++)					
	{
	species[i]=0;
	biomass[i]=0;
	}

	// ***** How many species and how much biomass in total? 
	for(i=1;i<=S;i++)
		if(y[i]>=extinct) 			
		{
			species[0]=species[0]+1;		
			biomass[0]=biomass[0]+y[i];	
		}

	// ***** How many species and how much biomass per trophic level

	fprintf(res_tp,"%1.0f \t %1.0f\t",species[0],t);
	tropos(y, S, res_tp);					// calculation + output of trophic positions


	for(i=1;i<=S;i++)
		if(y[i]>=extinct)				// only populations above the extinction threshold count for the evaluation!
		{
			for (j=1; j<=4; j++)
				if(TL[i]>=j-0.5 && TL[i]<j+0.5)
				{
					species[j]++;		// Species and biomass on level j 
					biomass[j]=biomass[j]+y[i];
				}
			if (TL[i]>=4.5)
			{
				species[5]++;			// Species and biomass on levels higher than 4.5 above the resource
				biomass[5]=biomass[5]+y[i];
			}
		}

	
	// ***** the file res_res contains lots of data of the network: 
	fprintf(res_res,"%1.0f\t", t);				// 1) time

	fprintf(res_res,"%1.3f \t", functional_diversity(S, y));	// 2) functional diversity
	fprintf(res_res,"%d \t", S);					// 3) global number of species 
	for(i=0;i<6;i++)
		fprintf(res_res,"%1.0f \t",species[i]);		// 4-9) number of species per trophic level

	fprintf(res_res, "%1.3f \t", y[0]);				// 10) Ressource biomass
	for(i=0;i<6;i++)
		fprintf(res_res,"%1.3f \t",biomass[i]);		// 11-16) local biomass (total and per trophic level)


	// ***** looking at the evolution of feeding ranges over time: 
	double mean_s=0, std_s=0;
	for (j=1; j<=S; j++)
		if(y[j]>=extinct)
			mean_s=mean_s + s_vec[j];   
	
	mean_s=mean_s/species[0];
	fprintf(res_res, "%1.3f \t", mean_s);				// 17) mean feeding range of all species>extinct 
	
	for (j=1; j<=S; j++)
		if(y[j]>=extinct)
			std_s = std_s + ((mean_s-s_vec[j])*(mean_s-s_vec[j]));
	
	std_s=std_s/species[0];
	std_s = sqrt(std_s);

	fprintf(res_res, "%1.3f \t", std_s);				// 18) standard deviation of the feeding ranges 
	fprintf(res_res,"\n");
	

	// ***** other data
	fprintf(res_m,"%1.0f \t %1.0f\t",species[0],t);		// body masses (of species>extinct) 
	for(i=0;i<=S;i++)
		if (y[i]>=extinct)
		  	fprintf(res_m,"%1.6f \t",log10(m_vec[i]));
		else
			fprintf(res_m,"0.000000 \t");
	fprintf(res_m,"\n");

	fprintf(res_c,"%1.0f \t %1.0f\t",species[0],t);		// feeding centres (of species>extinct) 
	for(i=0;i<=S;i++)
		if (y[i]>=extinct)
		  	fprintf(res_c,"%1.6f \t",log10(c_vec[i]));
		else
			fprintf(res_c,"0.000000 \t");
	fprintf(res_c,"\n");
	
	fprintf(res_s,"%1.0f \t %1.0f\t",species[0],t);		// feeding ranges (of species>extinct) 
	for(i=0;i<=S;i++)
		if (y[i]>=extinct)
		  	fprintf(res_s,"%1.6f \t",(s_vec[i]));
		else
			fprintf(res_s,"0.000000 \t");
	fprintf(res_s,"\n");

	fprintf(res_pop,"%1.0f \t %1.0f\t",species[0], t);		// population sizes (of species>extinct) 
	for(i=0;i<=S;i++)
		if (y[i]>=extinct)
	  		fprintf(res_pop,"%1.8f \t",y[i]);
		else
			fprintf(res_pop,"0.00000000 \t");
	fprintf(res_pop,"\n");		


// ******************** free memory (THIS IS IMPORTANT TO AVOID SEGMENTATION FAULTS!!)
	free(species);
	free(biomass);
	return 0;
}



//*******************************************************************************************
// output for the heatmaps
// (average over 1000 snap shots during the simulation runtime)
//*******************************************************************************************
int output_heatmap(FILE *res_heatmap)
{
	int i; 
	double S_mean=0, B_mean=0;
	double S_std=0, B_std=0;
	
	for (i=0; i<1000; i++)					// calculate mean of number of morphs and total biomass
	{
		S_mean = S_mean + heatmap_S[i];
		B_mean = B_mean + heatmap_B[i];
	}
	S_mean=S_mean/1000;
	B_mean=B_mean/1000;
	
	for (i=0; i<1000; i++)					// calculate standard deviation
	{
		S_std = S_std + (S_mean - heatmap_S[i]) * (S_mean - heatmap_S[i]);
		B_std = B_std + (B_mean - heatmap_B[i]) * (B_mean - heatmap_B[i]);
	}
	S_std=S_std/1000;
	S_std = sqrt(S_std);
	B_std=B_std/1000;
	B_std = sqrt(B_std);
	
	
	fprintf(res_heatmap, "%g\t %g\t %g\t %g\t %g\t %g\t %g\t \n", p, inv_z, lowbound, S_mean, B_mean, S_std, B_std);
	
	return 0;
}


//*******************************************************************************************
// output for turnover analysis
// (at every tmod_turnover time step or after every tmod_turnover/delta_Tm species introduction)
//*******************************************************************************************
int output_turnover(double t, int S, FILE *res_turnover,const double y[] )
{
	// ********** Calculate trophic similarity index of the two networks at time t and time t-tmod_turnover:
	int denominator= 0;
	int numerator=0; 
	
	for (int i=1; i<=S_old_tmod_turnover; i++) 		// for each of the species present in the old web...
	{
		for (int j=1; j<=i; j++)			// check whether it is still present in the new web!
			if (m_vec_old[i] == m_vec[j])
					numerator+=1;
	}
	
	denominator= S_old_tmod_turnover + S - numerator;
	fprintf(res_turnover, "%g\t %d\t %d \t %d \t", t, numerator, denominator, S);
	
	// ********** Calculate niche overlap index of the network at time t

	double *input=(double *) calloc((S+1),sizeof(double));	// vector input[i]: total feeding input of species i (per unit biomass)
	double *loss=(double *) calloc((S+1),sizeof(double));	// vector loss[i]: total feeding loss of species i (per unit biomass)

	for(int i=1;i<=S;i++)
	{
		input[i] = input[i] + efficiency*gij(y, S, i,0)* y[i];
		for(int j=1;j<=S;j++)			
		{
			input[i]=input[i] + efficiency*gij(y, S, i,j) *y[i];
			loss[i] =loss[i]  + efficiency*gij(y, S, j,i) *y[j];
		}
	}

	double** p_mat = (double **) calloc(S+1, sizeof(double*));
	for(int i=0; i <= S; i++)
	  	p_mat[i] = (double *) calloc(S+1, sizeof(double));	// matrix p[i][j]: proportion of energy that i gains via predation from j 
	double** q_mat = (double **) calloc(S+1, sizeof(double*));
	for(int i=0; i <= S; i++)
	  	q_mat[i] = (double *) calloc(S+1, sizeof(double));	// matrix q[i][j]: proportion of energy that i looses due to predation from j 

	for(int i=1;i<=S;i++)
	{
		p_mat[i][0]= (efficiency*gij(y, S, i,0)*y[i]) / input[i];
		for(int j=1;j<=S;j++)			
		{
			p_mat[i][j] = (efficiency*gij(y, S, i,j)*y[i]) / input[i]; 
			q_mat[i][j] = (efficiency*gij(y, S, j,i)*y[j]) / loss[i];
		}
	}


	double** overlap_mat = (double **) calloc(S+1, sizeof(double*));
	for(int i=0; i <= S; i++)
	  	overlap_mat[i] = (double *) calloc(S+1, sizeof(double));	// matrix overlap[ij]: niche overlap of species i and j 

	double denominator1= 0;
	double denominator2= 0;
	double numerator1=0; 
	double numerator2=0; 

	for(int i=1;i<=S;i++)
	{
		for(int j=1;j<=S;j++)
		{
			numerator1= p_mat[i][0] * p_mat[j][0];
			numerator2= 0;
			denominator1= p_mat[i][0] * p_mat[i][0]; 
			denominator2= p_mat[j][0] * p_mat[j][0]; 

			for(int l=1;l<=S;l++)
			{
				// See eq. 3 in Yodzis 1999 "In search of operational trophospecies in a tropical aquatic food web" for explanation!
				numerator1 = numerator1 + (p_mat[i][l] * p_mat[j][l]) ; 
				numerator2 = numerator2 + (q_mat[i][l] * q_mat[j][l]) ; 
				overlap_mat[i][j] = (numerator1 + numerator2);

				denominator1 = denominator1 + p_mat[i][l] * p_mat[i][l] + q_mat[i][l] * q_mat[i][l];
				denominator2 = denominator2 + p_mat[j][l] * p_mat[j][l] + q_mat[j][l] * q_mat[j][l];
				overlap_mat[i][j] = overlap_mat[i][j] / (sqrt(denominator1*denominator2));
			}
		}
	}

	double nicheoverlap =0;
	for(int i=1;i<=S;i++)
		for(int j=1;j<=S;j++)
			nicheoverlap = nicheoverlap + overlap_mat[i][j];

	nicheoverlap = nicheoverlap / (S*S); 
	fprintf(res_turnover, "%g \n", nicheoverlap);

	free(input);
	free(loss);

	for(int i = 0; i <= S; i++)
	  free(p_mat[i]);
	free(p_mat);
	for(int i = 0; i <= S; i++)
	  free(q_mat[i]);
	free(q_mat);
	for(int i = 0; i <= S; i++)
	  free(overlap_mat[i]);
	free(overlap_mat);
	
	return 0;
}

void output_speciesLifespanRecord(FILE *res_lifespan, Lifespan *speciesLifespanRecord, size_t nLifespanRecords)
{
	fprintf(res_lifespan, \
		"id,speciesType,bodyMass,feedingCenter,feedingRange,introduction,extinction\n");

	for (int i = 0; i < nLifespanRecords; i++)
	{
		Lifespan lifespanRecord = speciesLifespanRecord[i];

		fprintf(res_lifespan, \
			"%d,%s,%.6f,%.6f,%.6f,%.0f,%.0f\n", \
			lifespanRecord.id+1, \
			speciesTypeNames[lifespanRecord.speciesType], \
			lifespanRecord.bodyMass, \
			lifespanRecord.feedingCenter, \
			lifespanRecord.feedingRange, \
			lifespanRecord.introduction, \
			lifespanRecord.extinction);
	}
}
