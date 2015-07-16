#include "headers.h"
#include "init_2.h"
#include "declarations.h"

void initialize(vector<vertex>& node,vector<fval>& fvar)
{  
	int ind; 

	for(int j=0;j<=ny;j++)
	{
		for(int i=0;i<=nx;i++) 
		{   
			ind = i + j*str_x;   			
    						    
			fvar[ind].F = 0.0; 
			fvar[ind].dp = 0.0;   			 			
			fvar[ind].res = 0.0; 
			fvar[ind].err = 0.0; 
			fvar[ind].rhs = 0.0; 
			fvar[ind].div = 0.0; 
			fvar[ind].f_rhs = 0.0;
			fvar[ind].cor_rhs = 0.0; 
			
			for(int var=0;var<NUM_VAR;var++)
			{
			  fvar[ind].u[var] = 0.0; 
			  fvar[ind].ux[var] = 0.0; 
			  fvar[ind].uy[var] = 0.0; 
			  fvar[ind].uxx[var] = 0.0; 
			  fvar[ind].uyy[var] = 0.0; 
			  fvar[ind].u0[var] = 0.0; 
			}		
					 		
		} 
	}				
}

/**********************************************************************/
//Initializing the values of the RHS, i.e, F and others for various levels of multigrid 

void mg_initialize(vector<mg_grid>& level, vector<vertex>& node)
{
	int ind,lev;
	int sp; 

	for(lev=0;lev<mg_levels;lev++)
	{	
		sp = pow(2,lev); 

		for(int j=0;j<=ny;j=j+sp)
		{
	  		for(int i=0;i<=nx;i=i+sp)
	  		{
			  	ind = i+j*str_x; 	  		  			
	  						     			
       				level[lev].phi_s[ind] = 0.0; 	       				       			
	       			level[lev].F[ind] = 0.0; 
	       			level[lev].rhs[ind] =0.0;
	  			level[lev].res[ind] = 0.0;
	  			level[lev].cor_rhs[ind] = 0.0;
	  			level[lev].coeff[ind]=0.0; 
	  		}	
		}
	} 
}

/**********************************************************************/
//Initializing the isotropic random turbulence flow field

void turbulent_initialize(vector<fval>& fvar, vector<vertex>& node)
{
	int ind;
	double x_co,y_co; 
		
	for(int j=0; j<=ny;j++)
	{
		for(int i=0;i<=nx;i++)
		{
			ind = i + j*str_x;
			
			x_co = i*dx;
			y_co = j*dy;

			fvar[ind].u[0] = velocity_component(x_co,y_co,0);   //0--->x-component of velocity 			
			fvar[ind].u[1] = velocity_component(x_co,y_co,1);   //1--->y-component of velocity			
			fvar[ind].u[2] =  0.0;
		}
	}		
}

/**********************************************************************/
//Function to compute the velocity component 

double velocity_component(double x_co,double y_co,int comp)
{
	double a,b; //Real and imaginary parts of the streamfunction DFT
	double sum_y = 0.0, sum_x=0.0,k,E,tan_t,res_A,res_B;  //Running variables for summation in inverse DFT
				
	for(int k2=ky_min;k2<=ky_max;k2++)
	{
		for(int k1=kx_min;k1<=kx_max;k1++)
		{			
			k = sqrt( pow(k1,2) + pow(k2,2) );
			
			E = Q*pow(k,8)*exp(-4.*(k/kp)*(k/kp));
			
			tan_t = tan(random_num()); 	
			
			a = (E/M_PI*pow(k,3))/(1.+tan_t*tan_t); 
			b = a*tan_t; 
			
			res_A = (sqrt(pow(a,2) + pow(b,2)))/(tot_p); 
			res_B = random_num();
			
			if(comp==0)
			{
				sum_x = sum_x - k1*res_A*cos(res_B + k1*x_co + k2*y_co); 
				sum_y = sum_y - k1*res_A*sin(res_B + k1*x_co + k2*y_co); 
			}
			
			if(comp==1)
			{
				sum_x = sum_x + k2*res_A*cos(res_B + k1*x_co + k2*y_co); 
				sum_y = sum_y + k2*res_A*sin(res_B + k1*x_co + k2*y_co);
			}						
		}	
	}
	
	return sqrt(pow(sum_x,2) + pow(sum_y,2));
}

/**********************************************************************/
//Function that returns a random number between 0 and twopi

double random_num()
{		
	double MIN_RAND = 0.0, MAX_RAND = twopi/8.;
	double range = MAX_RAND - MIN_RAND;
	
	return range*((((double) rand()) / (double) RAND_MAX)) + MIN_RAND ;
}

/**********************************************************************/
