
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "combined.h"
#include "main.h"
#define Pop_size 20
#define Max_iterations 5
#define Total_Run 5
#define max_part 5  /*this represents the maximum group size*/
#define PI 3.14


int try,c_popnew,L_c_popnew;
double* ftemp;
double* tempf;
//ftemp= malloc(2*sizeof(double));
//tempf= malloc(2*sizeof(double));
double fun_val[200][2];
double fun_val_new[400 ][4];
double fun_val_copy[400][2];
double L_fun_val_new[400 ][4];
double L_fun_val_copy[400][2];  //first for sorting index on basis of obj1 values ,2nd is obj1 value,3rd is obj2 value and 4th is sorting index on basis of obj2
int tempf1,tempf2,tempf1index,tempf2index;
double sumofcd; 



//int check_dominance(int *a, int *b);


/*double CalculateFitness(double fun)
 {
	 double result=0;
	 if(fun>=0)
	 {result=1/(fun+1);
	 }
	 else
	 {
		 result=1+fabs(fun);
	 }
	 return result;
 }*/





void create_g()
{
   printf("\n%s\n", "i am going to create grp");
int g=0;
  for(lo=0;lo<Pop_size;)
    {

          group_size=(lo+(int)(Pop_size/part)-1);
           gpoint[g][0]=lo;
           gpoint[g][1]=group_size;
           if((Pop_size-group_size)<(Pop_size/part))
                {
   	         gpoint[g][1]=(Pop_size-1);
		break;
                }
           
     //printf("\n groups are...\n %d %d \n",lo,group_size);
     printf("\n%s\n", "i am creating grps.. ");
                g = g+1;
		lo = group_size+1;
                //printf("\n lo value is...%d",lo);
    }  
   //printf("\n");
	try=g;
      printf("this is trial%d \n",g);
	//try = g;
	//printf("thisoooooo is trial%d \n",try);
	//printf("\n %s ","okslll"); 
	//printf("\n %s "," okspp"); 
 }
void initialize()
{
   int i,j,k;
   	for(i=0;i<Pop_size;i++)
	{
        // printf("\nindividual is %d \n",i);
       for (j=0;j<D;j++)
		{
                 Population[i][j]=(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(ub[j]-lb[j])+lb[j];
		 new_position[j]=Population[i][j];
                //printf("poplulation is %lf",Population[i][j]);
		}// for j
		   Population[i][D]=0;//for global non domin
           Population[i][D+1]=0;//for local non domin
         ftemp=fun(new_position);
         for(k=0;k<2;k++)
          {
	   fun_val[i][k]=ftemp[k];
           printf("fun value after initializing %lf",fun_val[i][k]);
	  //fitness[i][k]=CalculateFitness(fun_val[i][k]);
          }
      }// for i
    // Initialize Global Leader Learning
    GlobalMin[0]=fun_val[0][0];
    GlobalMin[1]=fun_val[0][1];
    for(i=0;i<D;i++)
        GlobalLeaderPosition[i]=Population[0][i];
    GlobalLimitCount=0;
    // Initialize Local Leader Learning
    for(k=0;k<try;k++)
	{
        LocalMin[k][0]=fun_val[gpoint[k][0]][0];
        LocalMin[k][1]=fun_val[gpoint[k][0]][1]; 
        LocalLimitCount[k]=0;   
         for(i=0;i<D;i++)
           LocalLeaderPosition[k][i]=Population[gpoint[k][0]][i];            
        }
   //printf("%s","ok"); 
}

void GlobalLearning()
{
   int i,j,rnum;
    int k,m,p,flag=0;
    double G_trial[2]={0};
    int numdoms=0;
    double valcheck[2]={0};
    double valcheck1[2]={0};
     G_trial[0]=GlobalMin[0];
     G_trial[1]=GlobalMin[1];
     
   /*for(i=0;i<Pop_size;i++)//here find nondominant sol
	{
	if (fun_val[i]<GlobalMin)
		{
        GlobalMin=fun_val[i];
        for(j=0;j<D;j++)
           GlobalLeaderPosition[j]=Population[i][j];
        }
        //printf("\n %d ",i); 
	}// for i*/

      
    
    for(k=0;k<Pop_size;k++)
    {   flag=0;
        for(m=0;m<Pop_size;m++)
        {
            if(m!=k)
            {

              valcheck[0]=fun_val[k][0];
              valcheck[1]=fun_val[k][1];
              valcheck1[0]=fun_val[m][0];
              valcheck1[1]=fun_val[m][1];
             p=check_dominance(&valcheck[0],&valcheck1[0]);
             if(p==-1)
             {
                 flag=1;
                 break;
             }
            }
        }
        if(flag==0)
        Population[k][D]=1;
        else
        Population[k][D]=0;
        

        
         
    }
    /*for(i=0;i<5;i++)
    if(obj_array[i][2]==1)
    printf("\nthis is non dominant %d",i);*/
    for(i=0;i<Pop_size;i++)
    if(Population[i][D]==1)
    {
        numdoms=numdoms+1;
    printf("\nthis is non dominant %d and func values are %lf \t %lf",i,fun_val[i][0],fun_val[i][1]);
    
    }
    printf("\nnumber of non dominant in global learning  %d",numdoms);
    if(numdoms >= 1)
    { 
        do
        {    srand(time(NULL));
             rnum = (rand() % (Pop_size - 0 + 1)) + 0;
             //printf("\nrandoms in process\n %d",rnum);
        }while(Population[rnum][D]!=1);
        printf("\nrandom dominant is %d",rnum);
        
    
    GlobalMin[0]=fun_val[rnum][0];
    GlobalMin[0]=fun_val[rnum][0];
        for(j=0;j<D;j++)
           GlobalLeaderPosition[j]=Population[rnum][j];
    
    }


	if(fabs(G_trial[0]-GlobalMin[0])<acc_err && fabs(G_trial[1]-GlobalMin[1])<acc_err)
	       GlobalLimitCount=GlobalLimitCount+1;
       else
          GlobalLimitCount=0;
	
 }

void LocalLearning()
{  
   printf("\n%s\n", "i am going to local learnig");
   int i,j,k,m,p,flag=0;
   int numdoms=0,rnum;
   
    double OldMin[Pop_size/2][2];
    for(k=0;k<try;k++)
	{
        OldMin[k][0]=LocalMin[k][0];
        OldMin[k][1]=LocalMin[k][1];
	   //printf("\n %d", k);          
    }// end for k
    
  for(k=0;k<try;k++)
   {   
      numdoms=0;
      printf("\n%s\n", "i am going to try grp local learnig");
	for(i=gpoint[k][0];i<=gpoint[k][1];i++)
	{
          //printf("\n%d  is the i value",i);
          flag=0;
          for(m=gpoint[k][0];m<=gpoint[k][1];m++)
	  {
            //printf("\n%d  is the m value",m);
            if(m!=i)
            {
              p=check_dominance(&fun_val[i][0],&fun_val[m][0]);
              //printf("\n%d  is the p value",p);
              if(p==-1)
               {
                 flag=1;
                 break;
                printf("%s","we are breaking ");
               }
            }
	  }
          if(flag==0)
          {
          Population[i][D+1]=1;
          printf("\n%d  is the dominant in this grp ",i);
          }
          else
          Population[i][D+1]=0;
          
          
    }//for i =gpoint..
         
        



       for(i=gpoint[k][0];i<=gpoint[k][1];i++)
        if(Population[i][D+1]==1)
       {
        numdoms=numdoms+1;
       printf("\nthis is non dominant %d",i);
       }
       printf("\n%s\n","hello");
       //printf("\nnumber of non dominant in this group %d",numdoms);
       if(numdoms >=1)
         { 
        do
        {    srand(time(NULL));
             rnum = (rand() % (gpoint[k][1] -gpoint[k][0] + 1)) + gpoint[k][0];
             //printf("\nrandoms in process\n %d",rnum);
        }while(Population[rnum][D+1]!=1);
        printf("\nrandom dominant is %d",rnum);
        
       
       LocalMin[k][0]=fun_val[rnum][0];
       LocalMin[k][1]=fun_val[rnum][1];
        for(j=0;j<D;j++)
            LocalLeaderPosition[k][j]=Population[rnum][j];
       }//numdoms>=1


   }// for k
   
    for(k=0;k<try;k++)
	{
        if (fabs(OldMin[k][0]-LocalMin[k][0])<acc_err && fabs(OldMin[k][1]-LocalMin[k][1])<acc_err)
             LocalLimitCount[k]=LocalLimitCount[k]+1;    
        else
              LocalLimitCount[k]=0;        
       }// for k
	//printf("\n %s "," oksss"); 	
 }// end of fun


void LocalLeaderPhase(int k) //for now i am changing the posn to new position without comparing fitness here.
{
  int i,j;
  lo=gpoint[k][0];
  group_size=gpoint[k][1];
  for (i=lo;i<=group_size;i++)
  {
        int PopRand;
        do{
           PopRand=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(group_size-lo)+lo);
        }while(PopRand==i);
        
        for(j=0;j<D;j++){
        if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )>=cr)
		{
            new_position[j]=Population[i][j]+(LocalLeaderPosition[k][j]-Population[i][j])*((double)rand() / ((double)(RAND_MAX)+(double)(1)))+
                            (Population[PopRand][j]-Population[i][j])*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5)*2;
		     
	         }
	  
        else
        {
            new_position[j]=Population[i][j];    
        }
        if (new_position[j]<lb[j])
           new_position[j]=lb[j];
        if (new_position[j]>ub[j])
           new_position[j]=ub[j];
	//printf("\n %d ",i); 
  }// end of for j loop
        ObjValSol=fun(new_position);  //here i need to change....
        int p=check_dominance(&ObjValSol[0],&fun_val[i][0]);//here error can come...due to pointer
        //FitnessSol=CalculateFitness(ObjValSol);
          if (p==1)
           {
              for(j=0;j<D;j++)
               Population[i][j]=new_position[j];
              fun_val[i][0]=ObjValSol[0];
              fun_val[i][1]=ObjValSol[1];
         //fitness[i]=FitnessSol;
          }//for if p==1

          if (p==0)
		{
			for(j=0;j<D;j++)
                            L_Population_new[L_c_popnew][j]=new_position[j];
			L_fun_val_new[L_c_popnew][0]=L_c_popnew;
			L_fun_val_new[L_c_popnew][3]=L_c_popnew;
			L_fun_val_new[L_c_popnew][1]=ObjValSol[0];
			L_fun_val_new[L_c_popnew][2]=ObjValSol[1];
			L_c_popnew++;
		}//for p==0
		
         
	 }// for i
}// end of fun
void GlobalLeaderPhase(int k)
{

   int i,j,l;
  
  lo=gpoint[k][0];
  group_size=gpoint[k][1];
  i=lo;
  l=lo; /* l is the count*/
 while(l<group_size)
  {
         //printf("inside gl phase 
     
       if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )<prob[i]) 
        {  
              l++;      
         int PopRand;
        do{
       
        PopRand=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(group_size-lo)+lo);
        }while(PopRand==i);
        param2change=(int)((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*D);
        for(j=0;j<D;j++)
            new_position[j]=Population[i][j];
        new_position[param2change]=Population[i][param2change]+(GlobalLeaderPosition[param2change]-Population[i][param2change])*
(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )+(Population[PopRand][param2change]-Population[i][param2change])*((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )-0.5)*2;
        if (new_position[param2change]<lb[param2change])
           new_position[param2change]=lb[param2change];
        if (new_position[param2change]>ub[param2change])
           new_position[param2change]=ub[param2change];
        ObjValSol=fun(new_position);
         int p=check_dominance(&ObjValSol[0],&fun_val[i][0]);//here error can come...due to pointer
        //FitnessSol=CalculateFitness(ObjValSol);
        if (p==1)
        {  for(j=0;j<D;j++)
                   Population[i][j]=new_position[j];
              fun_val[i][0]=ObjValSol[0];
              fun_val[i][1]=ObjValSol[1];
              //fitness[i]=FitnessSol;
        }//for if p==1
		if (p==0)
		{
			for(j=0;j<D;j++)
                   Population_new[c_popnew][j]=new_position[j];
			fun_val_new[c_popnew][0]=c_popnew;
			fun_val_new[c_popnew][3]=c_popnew;
			fun_val_new[c_popnew][1]=ObjValSol[0];
			fun_val_new[c_popnew][2]=ObjValSol[1];
			c_popnew++;
		}//for p==0
		
		
		
    }// if prob[i]
     i++;
     if (i==(group_size))
        i=lo;
 }//while l
}//fun

void GlobalLeaderDecision()
{
if(GlobalLimitCount>GlobalLimit)
{
        
             
             GlobalLimitCount=0;
            
         if(part<max_part) 
        {     
              //printf("Global ... Partition\n");
             part=part+1;
              create_g(); 
              LocalLearning();
        }
        else
           {     //printf("Global ... Partition\n");
             part=1;
              create_g(); 
              LocalLearning();
        }
                        
                              
}
}
void LocalLeaderDecision()
{
     int i,j,k;
for(k=0;k<try;k++)
   {
          if(LocalLimitCount[k]>LocalLimit)
          {
            // printf("Local Limit Cross\n");
             for(i=gpoint[k][0];i<=gpoint[k][1];i++)
              {
                   for(j=0;j<D;j++){
                        if((   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )>=cr){
                        Population[i][j]=Population[i][j]+(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*(ub[j]-lb[j]);
                    }// end of if
                    else
                    {
                        Population[i][j]=Population[i][j]+(GlobalLeaderPosition[j]-Population[i][j])*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )+(Population[i][j]-LocalLeaderPosition[k][j])*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );    
                    }// end of else
                     if (Population[i][j]<lb[j])
                       Population[i][j]=lb[j];
                    if (Population[i][j]>ub[j])
                       Population[i][j]=ub[j];
                    }// end of for j
                    
	                tempf=fun( Population[i]);
                        fun_val[i][0]=tempf[0];
                        fun_val[i][1]=tempf[1];
                        //
	                //fitness[i]=CalculateFitness(fun_val[i]);
	           
	           
             }//For i
             
              LocalLimitCount[k]=0;
              //printf("\n %s", "hello");
         }//if
         
          
   }//k
  

}// End of Function


void CalculateProbabilities()
{
     int i;
     double maxfit;

    for (i=0;i<Pop_size;i++)
        {
         prob[i]=Population[i][D+3]/sumofcd;
         //prob[i]=1;
        }

}

int check_dominance(double *a, double *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;


    for (i=0; i<2; i++)
                {
                    if (a[i] < b[i])
                       flag1 = 1;
                    else if (a[i] > b[i])
                        flag2 = 1;
                        
                    
                }
                if (flag1==1 && flag2==0)
                {
                    return (1);
                }
                else
                {
                    if (flag1==0 && flag2==1)
                    {
                        return (-1);
                    }
                    else
                    {
                        return (0);
                    }
                }

}



int main()
{
int run,j;
double mean;
mean=0;
srand(time(NULL));
 E = exp( 1 );
        pi = acos( -1 );
//f_run = fopen( "SpiderMonkey_General1_SBCRor.ods", "w" );
//f_run1 = fopen( "SpiderMonkey_General1_SBCR_100or.ods", "w" );    
//fprintf(f_run," Prob \t Dim \t   Mean Fun Val   \t SD \t   Error \t Succ Mean Fun Eval \t  RUN \t Succ Rate \t Acceptable Error\n  ");
//fprintf(f_run1," Prob  \t No. of RUN  \t Fun Val \t   Error \n  ");
int d;


for ( Pr =0;Pr<1;Pr++)
{
    
       initilize_params(Pr);
double mean_error=0,error=0,total_feval=0;
double mean=0.0,var=0.0,sd=0.0;
mean_feval=0;
int succ_rate=0;
LocalLimit=D*Pop_size;
GlobalLimit=2;//global limit ko kamm kar rhi hu....so that groups will be formed( pehle = to Pop_size tha)
for(run=0;run<Total_Run;run++)
{

printf("\n\nThis is run number %d\n\n",run);


initialize();
GlobalLearning();
LocalLearning();
feval=0;
part=1;
create_g();

iter=0;
error=0;
cr=0.1;

for (iter=0;iter<Max_iterations;iter++)
 {                
    c_popnew=0;
    L_c_popnew=0;
    //copying initial population to new one for sorting 
	for(int k=0;k<Pop_size;k++)
	{
		for(int j=0;j<D;j++)
                        {
                         Population_new[c_popnew][j]=Population[k][j];
			L_Population_new[L_c_popnew][j]=Population[k][j];
                        }

			
		fun_val_new[c_popnew][0]=c_popnew;
		fun_val_new[c_popnew][3]=c_popnew;
		fun_val_new[c_popnew][1]=fun_val[k][0];
		fun_val_new[c_popnew][2]=fun_val[k][1];
		c_popnew++;


		L_fun_val_new[L_c_popnew][0]=L_c_popnew;
		L_fun_val_new[L_c_popnew][3]=L_c_popnew;
		L_fun_val_new[L_c_popnew][1]=fun_val[k][0];
		L_fun_val_new[L_c_popnew][2]=fun_val[k][1];
		L_c_popnew++;
             
		
	}


    
    printf("\n\nThis is iteration  number %d\n\n",iter);
  
    for(int k=0;k<try;k++)
    {
          LocalLeaderPhase(k);
		
    }
	
/////////////////////////////////////////////////////////



      
	//after this i have to choose Population based on crowding distance from new population whose size is c_popnew formed..sorting is in inc order
	
	//now we are sorting the fun val new array on basis of fun values to calculate crowding dist.
	
	//taking fun_val_new(0) sorted index for obj1
	
       printf("L_c_popnew after Lll is %d",L_c_popnew);
       
	for (int i = 0; i <= L_c_popnew; i++) //here i am storing the fun values of pop as a copy to update after acc to crowding distance...
        {
		L_Population_new[i][D+1]=L_fun_val_new[i][1];
		L_Population_new[i][D+2]=L_fun_val_new[i][2];
	}
		
	
	printf("\n%s\n","now sorting function value ");
    for (int i = 0; i <=L_c_popnew; i++) 
    {		
       //printf("\n%s\n","now sorting function value ");
       for (int j = 0; j <= L_c_popnew-i-1; j++)
	   {
           if (L_fun_val_new[j][1] > L_fun_val_new[j+1][1])
		   {
              
			  tempf1 = L_fun_val_new[j][1];
              L_fun_val_new[j][1]=L_fun_val_new[j+1][1];
              L_fun_val_new[j+1][1]= tempf1;
			  
			  tempf1index = L_fun_val_new[j][0];
              L_fun_val_new[j][0]=L_fun_val_new[j+1][0];
              L_fun_val_new[j+1][0]= tempf1index;
			  
		   }
		   
		   if (L_fun_val_new[j][2] > L_fun_val_new[j+1][2])

		   {
              
			  tempf2 = L_fun_val_new[j][2];
              L_fun_val_new[j][2]=L_fun_val_new[j+1][2];
              L_fun_val_new[j+1][2]= tempf2;
			  
			  tempf2index = L_fun_val_new[j][3];
              L_fun_val_new[j][3]=L_fun_val_new[j+1][3];
              L_fun_val_new[j+1][3]= tempf2index;
			  
		   }
	   }
		   
	}//for i =0 to popnew
	   //after this the fun_val_new will have indices sorted on basis of obj values..
	   // calculating crowding distances...
	
    for(int k=0;k<=L_c_popnew;k++)//crowding dist on basis of f1...
	{
		int index=0;
		for(int i=0;i<=L_c_popnew;i++)
		{
		if(L_fun_val_new[i][0]==k)
                index=i;
                }	
		
		if(index!=0 && index!=L_c_popnew)
			L_Population_new[k][D]=(L_fun_val_new[index+1][1]-L_fun_val_new[index-1][1])/1000;
			
               else 
		   L_Population_new[k][D]=999999;
        
        //printf("\n%s\n","now sorting cd ");			
		
	}//for k=0 to popsize
	
	
	
	for(int k=0;k<=L_c_popnew;k++)//crowding dist on basis of f2...
	{
		int index=0;
		for(int i=0;i<=L_c_popnew;i++)
		{
		if(L_fun_val_new[i][3]==k)
                 index=i;
                }	
		
		if(index!=0 && index!=L_c_popnew)
			L_Population_new[k][D]=L_Population_new[k][D]+(L_fun_val_new[index+1][2]-L_fun_val_new[index-1][2])/1000;
			
                else 
			L_Population_new[k][D]=999999;
        			
		
       printf("\n cd is......%lf for ....%d \n",L_Population_new[k][D],k);

	}//for k=0 to popsize
	
	//now i have to update the old pop and old fun_val on basis of least crowding distance...
	

        sumofcd=0;
	for(int i=0;i<Pop_size;i++)//in this step i can also copy cd to original for calculating probabity..
	{
		double mincd=9999999;
		int indexmin;
		for(int j=0;j<=L_c_popnew;j++)
		{
			if(L_Population_new[j][D]<mincd)
			{
				mincd=L_Population_new[j][D];
				indexmin=j;
			}
		}
		printf(" %d  minm cd is %lf ",i,mincd); 
		L_Population_new[indexmin][D]=9999999;
		
	    for(int k=0;k<D;k++)
			Population[i][k]=L_Population_new[indexmin][k];
		fun_val[i][0]=L_Population_new[indexmin][D+1];
		fun_val[i][1]=L_Population_new[indexmin][D+2];
                Population[i][D+3]=mincd;//storing the cd of the pop in d+3 row for calculating probability..
                sumofcd=sumofcd+mincd;
		
	}
		
	//after this the population and fun_val will be sorted acc to crowding dist 		
		
	






//////////////////////////////////////////////// ////////////////////////////////////         	
	
    CalculateProbabilities();//from here i have to change.....
/////////////////////////////////////////////////////////////////
	
		
     for(int k=0;k<try;k++)
    {
          GlobalLeaderPhase(k);
    }
	
	//after this i have to choose Population based on crowding distance from new population whose size is c_popnew formed..sorting is in inc order
	
	//now we are sorting the fun val new array on basis of fun values to calculate crowding dist.
	
	//taking fun_val_new(0) sorted index for obj1
	
       printf("c_popnew after gll is %d",c_popnew);
       
	for (int i = 0; i <= c_popnew; i++) //here i am storing the fun values of pop as a copy to update after acc to crowding distance...
    {
		Population_new[i][D+1]=fun_val_new[i][1];
		Population_new[i][D+2]=fun_val_new[i][2];
	}
		
	
	printf("\n%s\n","now sorting function value ");
	for (int i = 0; i <= c_popnew; i++) 
    {		
       //printf("\n%s\n","now sorting function value ");
       for (int j = 0; j <= c_popnew-i-1; j++)
	   {
           if (fun_val_new[j][1] > fun_val_new[j+1][1])
		   {
              
			  tempf1 = fun_val_new[j][1];
              fun_val_new[j][1]=fun_val_new[j+1][1];
              fun_val_new[j+1][1]= tempf1;
			  
			  tempf1index = fun_val_new[j][0];
              fun_val_new[j][0]=fun_val_new[j+1][0];
              fun_val_new[j+1][0]= tempf1index;
			  
		   }
		   
		   if (fun_val_new[j][2] > fun_val_new[j+1][2])
		   {
              
			  tempf2 = fun_val_new[j][2];
              fun_val_new[j][2]=fun_val_new[j+1][2];
              fun_val_new[j+1][2]= tempf2;
			  
			  tempf2index = fun_val_new[j][3];
              fun_val_new[j][3]=fun_val_new[j+1][3];
              fun_val_new[j+1][3]= tempf2index;
			  
		   }
	   }
		   
	}//for i =0 to popnew
	   //after this the fun_val_new will have indices sorted on basis of obj values..
	   // calculating crowding distances...
	
    for(int k=0;k<=c_popnew;k++)//crowding dist on basis of f1...
	{
		int index=0;
		for(int i=0;i<=c_popnew;i++)
		{
		if(fun_val_new[i][0]==k)
        index=i;
        }	
		
		if(index!=0 && index!=c_popnew)
			Population_new[k][D]=(fun_val_new[index+1][1]-fun_val_new[index-1][1])/1000;
			
        else 
			Population_new[k][D]=999999;
        
        //printf("\n%s\n","now sorting cd ");			
		
	}//for k=0 to popsize
	
	
	
	for(int k=0;k<=c_popnew;k++)//crowding dist on basis of f1...
	{
		int index=0;
		for(int i=0;i<=c_popnew;i++)
		{
		if(fun_val_new[i][3]==k)
                index=i;
                }	
		
		if(index!=0 && index!=c_popnew)
			Population_new[k][D]=Population_new[k][D]+(fun_val_new[index+1][2]-fun_val_new[index-1][2])/1000;
			
        else 
			Population_new[k][D]=999999;
        			
		
       printf("\n cd is......%lf for ....%d \n",Population_new[k][D],k);

	}//for k=0 to popsize
	
	//now i have to update the old pop and old fun_val on basis of least crowding distance...
	
	for(int i=0;i<Pop_size;i++)//in this step i can also copy cd to original for calculating probabity..
	{
		double mincd=9999999;
		int indexmin;
		for(int j=0;j<=c_popnew;j++)
		{
			if(Population_new[j][D]<mincd)
			{
				mincd=Population_new[j][D];
				indexmin=j;
			}
		}
		printf(" %d  minm cd is %lf ",i,mincd); 
		Population_new[indexmin][D]=9999999;
		
	    for(int k=0;k<D;k++)
			Population[i][k]=Population_new[indexmin][k];
		fun_val[i][0]=Population_new[indexmin][D+1];
		fun_val[i][1]=Population_new[indexmin][D+2];
		
	}
		
	//after this the population and fun_val will be sorted acc to crowding dist 		
		
	
			
		
	
	
	
	
            
    GlobalLearning();
    LocalLearning();
    LocalLeaderDecision();
    GlobalLeaderDecision ();

    

   	


}//iter
}//run


for(int i=0;i<Pop_size;i++)//for printing last pareto front
{

if(Population[i][D]==1)
printf("\nthe pareto front is  %lf and  %lf \n",fun_val[i][0],fun_val[i][1]);


}



}//pr
}//main



                 
/* 
   
    /*if(fabs(GlobalMin-obj_val)<=acc_err)
    {
        succ_rate+=1;
        mean_feval=mean_feval+feval;
	 
        break;
        
    }
   
    cr=cr+(0.4/Max_iterations);*/
   
   
//}//for iter
/*error=fabs(GlobalMin-obj_val);
mean_error=mean_error+error;
fprintf(f_run1,"%d \t %d  \t %e \t %e \n",Pr,run+1,GlobalMin,error);
printf("Pr %d Dim %d  Run %d  F_val %e Error %e \n",Pr, D, run+1,GlobalMin,error);
GlobalMins[run]=GlobalMin;
mean=mean+GlobalMin;
total_feval=total_feval+feval;*/

//}// end of run

/*mean=mean/Total_Run;
mean_error=mean_error/Total_Run;
if(succ_rate>0)
 mean_feval=mean_feval/(double)(succ_rate);
total_feval=total_feval/Total_Run;
 for(int k=0;k<Total_Run;k++)
    var=var+pow(GlobalMins[k]-mean,2);
    var=var/Total_Run;
    sd=sqrt(var);

printf("Means of %d runs: %e \n",Total_Run,mean);
fprintf(f_run," %d \t %d \t  %e \t %e \t %e \t  %f \t %d \t %d  \t %e \n",Pr,D, mean,sd, mean_error,total_feval,Total_Run, succ_rate,acc_err );*/



/*}// end of pr
fclose(f_run);
fclose(f_run1);

/return 0;

}


//}
/*FILE * f_synth;
FILE * f_expl;
FILE * fLandscape;
FILE * f_init;
FILE * f_init_save;

FILE * f_coeff;

FILE *f_swarm;*/ // For information about the variable swarm size*/




