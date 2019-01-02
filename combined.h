#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#include <time.h>

#define	D_max 200  // Max number of dimensions of the search space
#define	S_max 200 // Max swarm size
#define R_max 200 // Max number of runs
#define PI 3.14
int D;
int LocalLimit; // Spider Monkey
int GlobalLimit; // Spider Monkey
int limit;// ABC
double Population[S_max][D_max];
double Population_new[S_max*2][D_max+3];// this has d pe //crowding distt and d+1 pe unke f1 value and d+2 pe f2 valued
double L_Population_new[S_max*2][D_max+3];// this has d pe //crowding distt and d+1 pe unke f1 value and d+2 pe f2 valued
double fun_val[S_max ][2];  
double fun_val_new[S_max *2 ][4]; 
double L_fun_val_new[S_max *2 ][4];  
//double fitness[S_max ][2]; 
double prob[S_max ]; 
double new_position[D_max]; 
double* ObjValSol; 
double FitnessSol; 
int neighbour, param2change; 
double GlobalMin[2]; 
double GlobalLeaderPosition[D_max]; 
double LocalMin[S_max /2][2]; //
double LocalLeaderPosition[S_max /2][D_max]; 
int LocalLimitCount[S_max /2];
double GlobalMins[R_max]; 
int GlobalLimitCount;
int gpoint[S_max ][2];
double r,r1,r2; 
FILE * f_run,*f_run1;
int Pr;
double part;
double acc_err;
double lb[D_max],ub[D_max];
double feval;
double mean_feval,total_feval,mean_error;
int lo,group_size,g,iter;
double obj_val;
double cr;
double Foods[S_max ][D_max]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
double f[S_max ];  /*f is a vector holding objective function values associated with food sources */
double fitness_abc[S_max ]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
double trial[S_max ]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
double prob_abc[S_max ]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
double solution [D_max]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
double GlobalParams[D_max]; /*Parameters of the optimum solution*/
double GlobalMins_abc[R_max]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
//double pi;
//double E;
long double E;			// exp(1). Useful for some test functions
struct landscape funct; // When (x,f(x)) is read from a file
//int iter;
struct archive memPos;
double nEval;			// Number of fitness evaluations
long double pi;			// Useful for some test functions
int run;

double fv[2]={0};

double  * fun(double sol[]);

void initilize_params(int Pr)
{
 int d;
 if (Pr==0)
        {
           D=30;
           obj_val=0;
           acc_err=1.0e-5;
           for ( d = 0; d < D; d++ )
           {
               lb[d] = 1; ub[d] =50;
           }

        }
        
        
}


double * fun(double x[])
   {
     feval++;
    
      

 switch (Pr)
    {
      case 0:
//void UF1( double *f,int k)
	{       int nx=30;
		unsigned int j, i,count1, count2;
		double sum1, sum2, yj;
                //double fv[2];
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
                
                printf("\n%s\n","in func uf1");


		for(j = 2; j <= nx; j++)
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0)
			{
				sum2 += yj;
				count2++;
			}
			else
			{
				sum1 += yj;
				count1++;
			}
		}
		fv[0] = x[0]+ 2.0 * sum1 / (double)count1;
		fv[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
                //fprintf(f_run1,"%d \t %lf \t %lf \t  \t  \n",z,f[0],f[1]);
                //printf("\n\n%lf %lf\n",fv[0],fv[1]);//trying to print function values for obj function problem in generating random number....also ask which form he want output.
	}
return fv;
}
}

