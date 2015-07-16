#include "headers.h"
#include "init_2.h"
#include "declarations.h"

///////////////////Function overloaded of Writing to a file////////////////////
//Writing data at frequency specified by file_freq in the "init" file

void write_to_file(vector<vertex>& node, vector<fval>& fvar, int t)
{
	ofstream out_put; 
	
	int p;
	double flo_time=t*dt; 

	string file("data_2D_turb");
	
	stringstream tag; 

	tag<<t;

	file = file + "_" + tag.str() + ".vtk";   	

	out_put.open(file.c_str());

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
	out_put<<"SCALARS"<<" divergence"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" div_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].div<<"\n"; 
		}
	}

	/*out_put<<"SCALARS"<<" error_pr"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" er_pr_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2] +  ( 0.25*( cos(2.*k*node[p].x[0]) + cos(2.*k*node[p].x[1]) )*exp(-4.*k*k*flo_time/Re) )<<"\n";	
		}
	}
	
	out_put<<"SCALARS"<<" error_u"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" er_u_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0] + (cos(k*node[p].x[0])*sin(k*node[p].x[1])*exp(-2.*k*k*flo_time/Re)) <<"\n"; 
		}
	}

	out_put<<"SCALARS"<<" error_v"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" er_v_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[1] - ( sin(k*node[p].x[0])*cos(k*node[p].x[1])*exp(-2.*k*k*flo_time/Re) ) <<"\n"; 
		}
	}*/	
	
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].res<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" F"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" F_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].F<<"\n"; 
		}
	}

	out_put<<"SCALARS"<<" Vorticity"<<" double"<<"\n";
        out_put<<"LOOKUP_TABLE"<<" Vorticity_lp"<<"\n";

        for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        out_put<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }


	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	out_put.close(); 
}

///////////////////Writing to a file////////////////////

void write_to_file(vector<vertex>& node, vector<fval>& fvar)
{	

	ofstream out_put;  
	int p; 

	out_put.open("2D_turb.vtk");

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			//out_put<<i<<"	"<<j<<"	   "<<fvar[p].u[2]<<"\n"; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
		
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 

			out_put<<fvar[p].res<<"\n"; 
		}
	}
		
	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	out_put.close();
}
