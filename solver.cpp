#include "headers.h"
#include "init_2.h"
#include "declarations.h"

///////////////////Solver function////////////////////

void solver(vector<vertex>& node, vector<fval>& fvar, vector<mg_grid>& level, pbcs& pbc, int restart_count)
{		
	double max_div; 
	ofstream fp_res; //Filestream to write the residuals as per time for the 3 flow variables
	
	fp_res.open("global_residuals.dat"); 
	
	cout<<"The time step is "<<dt<<"\n"; 
	cout<<"The no.of timesteps is "<<ite<<"\n"; 
	cout<<"The spacing is dx = "<<dx<< " dy= "<<dy<<"\n"; 

	for(int t=restart_count+1;t<=ite;t++)
	{
		explicit_solver(node,fvar,level,pbc);
		per_bcs(node,fvar); 		
		max_div = div_calc(node,fvar);
		  
		cout<<"----------- "<<"ite "<<t<<" "<<max_div<<"------------"<<"\n"; 	
		global_residuals(node,fvar,fp_res,t);
		
		if(t%file_freq==0) write_to_file(node,fvar,t); 		
		if(t%r_file_freq==0) write_restart(fvar,t);
	}
	
	fp_res.close(); 
}

///////////////////////////////////////////////////////////////////////////////////

void explicit_solver(vector<vertex>& node, vector<fval>& fvar, vector<mg_grid>& level, pbcs& pbc)
{ 
	int ind; 
	
	for(int j=0;j<ny;j++)
	{		
		for(int i=0;i<nx;i++)
		{
			ind = i + j*str_x; 
			
			fvar[ind].u0[0] = fvar[ind].u[0];                                    //U1  stage-1 mod-RK scheme 
			fvar[ind].u0[1] = fvar[ind].u[1]; 		 		     //V1 
			fvar[ind].u0[2] = fvar[ind].u[2]; 				     //Pressure for the time residual
		}
	}
	
	mg_clear_levels(level); 
	poisson_source(fvar,level,4.0); 
	mg_poisson_solver(level,fvar,pbc);
	per_tdma1x(fvar,2); 
	per_tdma1y(fvar,2); 

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/4.)*( fvar[ind].u[3] - fvar[ind].ux[2]);     //U2 stage-2 mod-RK scheme 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/4.)*( fvar[ind].u[4] - fvar[ind].uy[2]);     //V2 
		}
	}												

	mg_clear_levels(level);
	poisson_source(fvar,level,3.0);
	mg_poisson_solver(level,fvar,pbc);
	per_tdma1x(fvar,2); 
	per_tdma1y(fvar,2); 

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/3.)*( fvar[ind].u[3] - fvar[ind].ux[2]);    //U3 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/3.)*( fvar[ind].u[4] - fvar[ind].uy[2]);	   //V3 
		}
	}

	mg_clear_levels(level);
	poisson_source(fvar,level,2.0);
	mg_poisson_solver(level,fvar,pbc);	
	per_tdma1x(fvar,2); 
	per_tdma1y(fvar,2); 

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++) 
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + (dt/2.)*( fvar[ind].u[3] - fvar[ind].ux[2]);    //U4 
			fvar[ind].u[1] = fvar[ind].u0[1] + (dt/2.)*( fvar[ind].u[4] - fvar[ind].uy[2]);	   //V4 
		}
	}

	mg_clear_levels(level);
	poisson_source(fvar,level,1.0);
	mg_poisson_solver(level,fvar,pbc); 
	per_tdma1x(fvar,2); 
	per_tdma1y(fvar,2); 

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			ind = i + j*str_x; 
	
			fvar[ind].u[0] = fvar[ind].u0[0] + dt*( fvar[ind].u[3] - fvar[ind].ux[2]);    	//U(n+1)
			fvar[ind].u[1] = fvar[ind].u0[1] + dt*( fvar[ind].u[4] - fvar[ind].uy[2]);    	//V(n+1) 
		}
	}	
}
