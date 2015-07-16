#include "headers.h"
#include "init_2.h"
#include "declarations.h"

double div_calc(vector<vertex>& node, vector<fval>& fvar)
{
	double max_div=0.0; 
	int ind; 
	
	per_tdma1x(fvar,0);
	per_tdma1y(fvar,1); 
	
	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			ind = i + j*str_x; 		
			
			fvar[ind].div = fvar[ind].ux[0] + fvar[ind].uy[1]; 
			
		 	if( abs(fvar[ind].div)>max_div ) max_div = abs(fvar[ind].div); 		
		}
	}
	
	return max_div; 
}

/*********************************************************************************************/
/****The L2 norm of the residuals from time step to time step for all the basic variables******/
/****L_inf norm for U-velocity and V-velocity****/ 
/****The file pointer is passed to the function from the solver function*/

void global_residuals(vector<vertex>& node, vector<fval>& fvar,ofstream& fp_res,int ite)
{
  int ind; 
  double time = ite*dt; //Computing the time at which these are computed
  double u_res=0.0,v_res=0.0, p_res=0.0; 
  double u_exact=0.0, v_exact=0.0, max_u_err=0.0, max_v_err=0.0; 
  double u_error =0.0, v_error=0.0; 
  
  for(int j=0;j<ny;j++)
  {
    for(int i=0;i<nx;i++)
    {
	ind = i + j*str_x;   
	
	u_res = u_res + (fvar[ind].u[0] - fvar[ind].u0[0])*(fvar[ind].u[0] - fvar[ind].u0[0]); 
	v_res = v_res + (fvar[ind].u[1] - fvar[ind].u0[1])*(fvar[ind].u[1] - fvar[ind].u0[1]); 
	p_res = p_res + (fvar[ind].u[2] - fvar[ind].u0[2])*(fvar[ind].u[2] - fvar[ind].u0[2]);		
	
	/*u_exact = -(cos(k*node[ind].x[0])*sin(k*node[ind].x[1])*exp(-2.*k*k*time/Re)); 
	v_exact =  ( sin(k*node[ind].x[0])*cos(k*node[ind].x[1])*exp(-2.*k*k*time/Re) );
	
	u_error = abs( u_exact - fvar[ind].u[0] ); 
	v_error = abs( v_exact - fvar[ind].u[1] );   					
	
	if(u_error>=max_u_err)
	{
		max_u_err = u_error; 
	}
	
	if(v_error>=max_v_err)
	{
		max_v_err = v_error; 
	}*/	
    }  
  }
  
  fp_res<<time<<"		"<<u_res<<"	"<<v_res<<"	"<<p_res<<"	"<<max_u_err<<"		"<<max_v_err<<"\n"; 
}
/*********************************************************************************************/
