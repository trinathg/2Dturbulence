#include "headers.h"
#include "init_2.h"
#include "declarations.h"

/**************************************************************************************************************************************/
/*Specifying and initializing the coefficients globally*/
/*Heavily modified the coefficients for doubly periodic condition on Feb 2, 2014*/
/******Heavily modified to accomodate the periodic BCs******/ 

vector<double> hx(mg_levels), hy(mg_levels), hX2(mg_levels), hY2(mg_levels), hxx(mg_levels), hyy(mg_levels); //The grid sizes at various levels of the multigrid cycle
vector<double> XA(mg_levels),YA(mg_levels),XH(mg_levels),YH(mg_levels),YB(mg_levels);

/*Case-1*/
vector<double> dc_b1_11(mg_levels), dc_b1_12(mg_levels), dc_b1_13(mg_levels), dc_nb1_21(mg_levels), dc_nb1_22(mg_levels), dc_nb1_23(mg_levels), dc_nb1_24(mg_levels), dc_i1_32(mg_levels), dc_i1_33(mg_levels), dc_i1_34(mg_levels), dc_i1_35(mg_levels), dc_i1_36(mg_levels); 
vector<double> dc_b2_11(mg_levels), dc_b2_12(mg_levels), dc_b2_13(mg_levels), dc_nb2_21(mg_levels), dc_nb2_22(mg_levels), dc_nb2_23(mg_levels), dc_nb2_24(mg_levels), dc_i2_32(mg_levels), dc_i2_33(mg_levels), dc_i2_34(mg_levels), dc_i2_35(mg_levels), dc_i2_36(mg_levels); 
vector<double> dc_b3_11(mg_levels), dc_b3_12(mg_levels), dc_b3_13(mg_levels), dc_nb3_21(mg_levels), dc_nb3_22(mg_levels), dc_nb3_23(mg_levels), dc_nb3_24(mg_levels), dc_i3_32(mg_levels), dc_i3_33(mg_levels), dc_i3_34(mg_levels), dc_i3_35(mg_levels), dc_i3_36(mg_levels); 

vector<double> dc_ny1_11(mg_levels), dc_ny1_12(mg_levels), dc_ny1_13(mg_levels), dc_ny1_21(mg_levels), dc_ny1_22(mg_levels), dc_ny1_23(mg_levels), dc_ny1_24(mg_levels), dc_ny1_32(mg_levels), dc_ny1_33(mg_levels), dc_ny1_34(mg_levels), dc_ny1_35(mg_levels), dc_ny1_36(mg_levels);
vector<double> dc_ny_11(mg_levels), dc_ny_12(mg_levels), dc_ny_13(mg_levels), dc_ny_21(mg_levels), dc_ny_22(mg_levels), dc_ny_23(mg_levels), dc_ny_24(mg_levels), dc_ny_32(mg_levels), dc_ny_33(mg_levels), dc_ny_34(mg_levels), dc_ny_35(mg_levels), dc_ny_36(mg_levels);
 
/*Case-2*/
vector<double> dc_b4_11(mg_levels), dc_b4_12(mg_levels), dc_b4_13(mg_levels), dc_nb4_21(mg_levels), dc_nb4_22(mg_levels), dc_nb4_23(mg_levels), dc_nb4_24(mg_levels), dc_i4_32(mg_levels), dc_i4_33(mg_levels), dc_i4_34(mg_levels), dc_i4_35(mg_levels), dc_i4_36(mg_levels); 
vector<double> dc_b5_11(mg_levels), dc_b5_12(mg_levels), dc_b5_13(mg_levels), dc_nb5_21(mg_levels), dc_nb5_22(mg_levels), dc_nb5_23(mg_levels), dc_nb5_24(mg_levels), dc_i5_32(mg_levels), dc_i5_33(mg_levels), dc_i5_34(mg_levels), dc_i5_35(mg_levels), dc_i5_36(mg_levels); 
vector<double> dc_b6_11(mg_levels), dc_b6_12(mg_levels), dc_b6_13(mg_levels), dc_nb6_21(mg_levels), dc_nb6_22(mg_levels), dc_nb6_23(mg_levels), dc_nb6_24(mg_levels), dc_i6_32(mg_levels), dc_i6_33(mg_levels), dc_i6_34(mg_levels), dc_i6_35(mg_levels), dc_i6_36(mg_levels); 
vector<double> dc_b7_11(mg_levels), dc_b7_12(mg_levels), dc_b7_13(mg_levels), dc_nb7_21(mg_levels), dc_nb7_22(mg_levels), dc_nb7_23(mg_levels), dc_nb7_24(mg_levels), dc_i7_32(mg_levels), dc_i7_33(mg_levels), dc_i7_34(mg_levels), dc_i7_35(mg_levels), dc_i7_36(mg_levels); 
vector<double> dc_pb7_11(mg_levels), dc_pb7_12(mg_levels), dc_pb7_13(mg_levels), dc_pb7_21(mg_levels), dc_pb7_22(mg_levels), dc_pb7_23(mg_levels), dc_pb7_24(mg_levels), dc_pb7_32(mg_levels), dc_pb7_33(mg_levels), dc_pb7_34(mg_levels), dc_pb7_35(mg_levels), dc_pb7_36(mg_levels);

/*Case-3*/
vector<double> dc_b8_11(mg_levels), dc_b8_12(mg_levels), dc_b8_13(mg_levels), dc_nb8_21(mg_levels), dc_nb8_22(mg_levels), dc_nb8_23(mg_levels), dc_nb8_24(mg_levels), dc_i8_32(mg_levels), dc_i8_33(mg_levels), dc_i8_34(mg_levels), dc_i8_35(mg_levels), dc_i8_36(mg_levels); 
vector<double> dc_b9_11(mg_levels), dc_b9_12(mg_levels), dc_b9_13(mg_levels), dc_nb9_21(mg_levels), dc_nb9_22(mg_levels), dc_nb9_23(mg_levels), dc_nb9_24(mg_levels), dc_i9_32(mg_levels), dc_i9_33(mg_levels), dc_i9_34(mg_levels), dc_i9_35(mg_levels), dc_i9_36(mg_levels); 
vector<double> dc_b10_11(mg_levels), dc_b10_12(mg_levels), dc_b10_13(mg_levels), dc_nb10_21(mg_levels), dc_nb10_22(mg_levels), dc_nb10_23(mg_levels), dc_nb10_24(mg_levels), dc_i10_32(mg_levels), dc_i10_33(mg_levels), dc_i10_34(mg_levels), dc_i10_35(mg_levels), dc_i10_36(mg_levels); 
vector<double> dc_b11_11(mg_levels), dc_b11_12(mg_levels), dc_b11_13(mg_levels), dc_nb11_21(mg_levels), dc_nb11_22(mg_levels), dc_nb11_23(mg_levels), dc_nb11_24(mg_levels), dc_i11_32(mg_levels), dc_i11_33(mg_levels), dc_i11_34(mg_levels), dc_i11_35(mg_levels), dc_i11_36(mg_levels); 
vector<double> dc_b12_11(mg_levels), dc_b12_12(mg_levels), dc_b12_13(mg_levels), dc_nb12_21(mg_levels), dc_nb12_22(mg_levels), dc_nb12_23(mg_levels), dc_nb12_24(mg_levels), dc_i12_32(mg_levels), dc_i12_33(mg_levels), dc_i12_34(mg_levels), dc_i12_35(mg_levels), dc_i12_36(mg_levels); 

/*Coefficients for the RHS side*/

vector<double> dc_b13_11(mg_levels), dc_b13_12(mg_levels), dc_b13_13(mg_levels), dc_nb13_21(mg_levels), dc_nb13_22(mg_levels), dc_nb13_23(mg_levels), dc_nb13_24(mg_levels), dc_nb13_25(mg_levels); 
vector<double> dc_b14_11(mg_levels), dc_b14_12(mg_levels), dc_b14_13(mg_levels), dc_nb14_21(mg_levels), dc_nb14_22(mg_levels), dc_nb14_23(mg_levels), dc_nb14_24(mg_levels), dc_nb14_25(mg_levels); 

/****************************************************************************************************************************************/
/*******************************************Multi grid coefficients**********************************************************************/

void mg_coeff()
{
	/*ofstream coeff_c, coeff_d; 
	
	coeff_c.open("coeff_c.dat");
	coeff_d.open("coeff_d.dat");*/
	
	double sum_coeffs=0.0; 
	
	for(int lev=0;lev<mg_levels;lev++)
	{	
		hx[lev] = dx*pow(2,lev); 
		hy[lev] = dy*pow(2,lev);
		
		hX2[lev] = hx[lev]*hx[lev];
		hY2[lev] = hy[lev]*hy[lev]; 
		
		hxx[lev] = 1./( hx[lev]*hx[lev] ); 
		hyy[lev] = 1./( hy[lev]*hy[lev] ); 
						
		/********Neumann Coefficients after BCs have been accommodated*********/
		/********Heavily modified coefficients for doubly periodic domain******/ 
				
		XA[lev] = alpha*hxx[lev];	//Immediate coefficients 
		YA[lev] = a*hyy[lev];
		
		XH[lev] = hxx[lev];		//The central coefficients in the difference formulation 
		YH[lev] = -2.*(a+b)*hyy[lev];
		
		YB[lev] = b*hyy[lev]; 		//Out on the far edge
	
	
		/************************************Case-1 (j=1)********************************/ 
		//j 
		
		dc_b1_11[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 		
		dc_b1_12[lev] = a*XH[lev] + YH[lev]*alpha; 
		dc_b1_13[lev] = b*XH[lev]; //Two more entries are similar to the previous ones  		
		
		dc_nb1_21[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_nb1_22[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_nb1_23[lev] = a*XH[lev] + alpha*YH[lev]; 
		dc_nb1_24[lev] = b*XH[lev]; //One more entry same as this owing to periodicity				

		dc_i1_32[lev] = XH[lev]*b;    
		dc_i1_33[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_i1_34[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_i1_35[lev] = XH[lev]*a + YH[lev]*alpha;
		dc_i1_36[lev] = XH[lev]*b; 		
				
		/*coeff_d<<dc_b1_11[lev]<<"\n"; 
		coeff_d<<dc_b1_12[lev]<<"\n"; 
		coeff_d<<dc_b1_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb1_21[lev]<<"\n"; 
		coeff_d<<dc_nb1_22[lev]<<"\n"; 
		coeff_d<<dc_nb1_23[lev]<<"\n"; 
		coeff_d<<dc_nb1_24[lev]<<"\n"; 
		
		coeff_d<<dc_i1_32[lev]<<"\n"; 
		coeff_d<<dc_i1_33[lev]<<"\n"; 
		coeff_d<<dc_i1_34[lev]<<"\n"; 
		coeff_d<<dc_i1_35[lev]<<"\n";
		coeff_d<<dc_i1_36[lev]<<"\n";*/
		
		//j+1
		
		dc_b2_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_b2_12[lev] = a*XA[lev] + YA[lev]*alpha; 
		dc_b2_13[lev] = b*XA[lev]; //Two more entries are similar to the previous ones  		

		dc_nb2_21[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_nb2_22[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_nb2_23[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_nb2_24[lev] = b*XA[lev]; 	//One more entry due to periodicity 			 
		
		dc_i2_32[lev] = XA[lev]*b; 	
		dc_i2_33[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i2_34[lev] = -2.*(a+b)*XA[lev] + YA[lev];
		dc_i2_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i2_36[lev] = XA[lev]*b;
						
		/*coeff_d<<dc_b2_11[lev]<<"\n"; 
		coeff_d<<dc_b2_12[lev]<<"\n"; 
		coeff_d<<dc_b2_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb2_21[lev]<<"\n"; 
		coeff_d<<dc_nb2_22[lev]<<"\n"; 
		coeff_d<<dc_nb2_23[lev]<<"\n"; 
		coeff_d<<dc_nb2_24[lev]<<"\n"; 
		
		coeff_d<<dc_i2_32[lev]<<"\n"; 
		coeff_d<<dc_i2_33[lev]<<"\n"; 
		coeff_d<<dc_i2_34[lev]<<"\n"; 
		coeff_d<<dc_i2_35[lev]<<"\n";
		coeff_d<<dc_i2_36[lev]<<"\n";*/

		//j+2 
		
		dc_b3_11[lev] = YB[lev]; 
		dc_b3_12[lev] = YB[lev]*alpha; 
		dc_b3_13[lev] = 0.0; //Two more entries similar to the previous ones and one among them is zero.
				     //One non-zero periodic entry  		

		dc_nb3_21[lev] = YB[lev]*alpha; 
		dc_nb3_22[lev] = YB[lev]; 	
		dc_nb3_23[lev] = YB[lev]*alpha; 
		dc_nb3_24[lev] = 0.0; 			

		dc_i3_32[lev] = 0.0;
		dc_i3_33[lev] = YB[lev]*alpha;
		dc_i3_34[lev] = YB[lev];
		dc_i3_35[lev] = YB[lev]*alpha;
		dc_i3_36[lev] = 0.0;	
						
		//j = ny-2 
		
		dc_ny1_11[lev] = YB[lev]; 
		dc_ny1_12[lev] = alpha*YB[lev]; 
		dc_ny1_13[lev] = alpha*YB[lev];  //This is a periodic entry at the start/end of the row
		
		dc_ny1_21[lev] = YB[lev]*alpha; 
		dc_ny1_22[lev] = YB[lev]; 
		dc_ny1_23[lev] = YB[lev]*alpha; 
		dc_ny1_24[lev] = 0.0; 		//No periodic entries on this line
		
		dc_ny1_32[lev] = 0.0; 
		dc_ny1_33[lev] = YB[lev]*alpha; 
		dc_ny1_34[lev] = YB[lev]; 
		dc_ny1_35[lev] = YB[lev]*alpha; 
		dc_ny1_36[lev] = 0.0; 		//No periodic entries on this line
		
		//j=ny-1 
		
		dc_ny_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_ny_12[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_ny_13[lev] = b*XA[lev]; //Two more periodic entries at the end of the row. dc_ny_13 and dc_ny_12 strictly in that order
		
		dc_ny_21[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_ny_22[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_ny_23[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_ny_24[lev] = b*XA[lev]; //One periodic entry at the end of the row same as dc_ny_24 
		
		dc_ny_32[lev] = XA[lev]*b; 
		dc_ny_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_ny_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_ny_35[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_ny_36[lev] = XA[lev]*b; 
									
		/*coeff_d<<dc_b3_11[lev]<<"\n"; 
		coeff_d<<dc_b3_12[lev]<<"\n"; 
		coeff_d<<dc_b3_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb3_21[lev]<<"\n"; 
		coeff_d<<dc_nb3_22[lev]<<"\n"; 
		coeff_d<<dc_nb3_23[lev]<<"\n"; 
		coeff_d<<dc_nb3_24[lev]<<"\n"; 
		
		coeff_d<<dc_i3_32[lev]<<"\n"; 
		coeff_d<<dc_i3_33[lev]<<"\n"; 
		coeff_d<<dc_i3_34[lev]<<"\n"; 
		coeff_d<<dc_i3_35[lev]<<"\n";
		coeff_d<<dc_i3_36[lev]<<"\n";*/
		
		/*if(lev==0)
		{
			cout<<" ----------------Case-1 checks---------------"<<"\n";
			
			cout<<" Row-1: "<<dc_b1_11[lev] + 2.*( dc_b1_12[lev] + dc_b1_13[lev] ) + 2.*( dc_b2_11[lev] + 2.*(dc_b2_12[lev] + dc_b2_13[lev]) ) + 2.*( dc_b3_11[lev] + 2.*(dc_b3_12[lev] + dc_b3_13[lev])) <<"\n"; 		
			
			cout<<" Row-2: "<<dc_nb1_21[lev] + dc_nb1_22[lev] + dc_nb1_23[lev] + 2.*dc_nb1_24[lev] + dc_nb2_21[lev] + dc_nb2_22[lev] + dc_nb2_23[lev] + 2.*dc_nb2_24[lev] + dc_nb3_21[lev] + dc_nb3_22[lev] + dc_nb3_23[lev] + 2.*dc_nb3_24[lev] + dc_ny1_21[lev] + dc_ny1_22[lev] + dc_ny1_23[lev] + 2.*dc_ny1_24[lev] + dc_ny_21[lev] + dc_ny_22[lev] + dc_ny_23[lev] + 2.*dc_ny_24[lev]<<"\n"; 
			
			cout<<" Row-3: "<<dc_i1_32[lev] + dc_i1_33[lev] + dc_i1_34[lev] + dc_i1_35[lev] + dc_i1_36[lev] + dc_i2_32[lev] + dc_i2_33[lev] + dc_i2_34[lev] + dc_i2_35[lev] + dc_i2_36[lev] + dc_i3_32[lev] + dc_i3_33[lev] + dc_i3_34[lev] + dc_i3_35[lev] + dc_i3_36[lev] + dc_ny1_32[lev] + dc_ny1_33[lev] + dc_ny1_34[lev] + dc_ny1_35[lev] + dc_ny1_36[lev] + dc_ny_32[lev] + dc_ny_33[lev] + dc_ny_34[lev] + dc_ny_35[lev] + dc_ny_36[lev]<<"\n"; 
		}*/
	
						
		/********************************Case-2 (j=2 in code here)********************************/ 
		
		//j-1 

		dc_b4_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_b4_12[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_b4_13[lev] = b*XA[lev]; //Two more entries owing to periodicity

		dc_nb4_21[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_nb4_22[lev] = -2.*(a+b)*XA[lev] + YA[lev];		
		dc_nb4_23[lev] = XA[lev]*a + YA[lev]*alpha; 		
		dc_nb4_24[lev] = XA[lev]*b; //One more entry owing to periodicity

		dc_i4_32[lev] = XA[lev]*b;    
		dc_i4_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i4_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i4_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i4_36[lev] = XA[lev]*b;
		
		/*coeff_d<<dc_b4_11[lev]<<"\n"; 
		coeff_d<<dc_b4_12[lev]<<"\n"; 
		coeff_d<<dc_b4_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb4_21[lev]<<"\n"; 
		coeff_d<<dc_nb4_22[lev]<<"\n"; 
		coeff_d<<dc_nb4_23[lev]<<"\n"; 
		coeff_d<<dc_nb4_24[lev]<<"\n"; 
		
		coeff_d<<dc_i4_32[lev]<<"\n"; 
		coeff_d<<dc_i4_33[lev]<<"\n"; 
		coeff_d<<dc_i4_34[lev]<<"\n"; 
		coeff_d<<dc_i4_35[lev]<<"\n";
		coeff_d<<dc_i4_36[lev]<<"\n";
		
		if(lev==0)  
		{
			sum_coeffs = sum_coeffs + dc_b4_11[lev] + 2.*dc_b4_12[lev] + 2.*dc_b4_13[lev] + dc_nb4_21[lev] + dc_nb4_22[lev] + dc_nb4_23[lev] + 2.*dc_nb4_24[lev] + dc_i4_32[lev] + dc_i4_33[lev] + dc_i4_34[lev] + dc_i4_35[lev] + dc_i4_36[lev];			
		}*/

		//j 

		dc_b5_11[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_b5_12[lev] = a*XH[lev] + YH[lev]*alpha; 
		dc_b5_13[lev] = b*XH[lev]; 

		dc_nb5_21[lev] = a*XH[lev] + alpha*YH[lev]; 
		dc_nb5_22[lev] = -2.*(a+b)*XH[lev] + YH[lev];		
		dc_nb5_23[lev] = a*XH[lev] + alpha*YH[lev]; 		
		dc_nb5_24[lev] = XH[lev]*b;  //One more entry which is same as this one owing to periodicity 

		dc_i5_32[lev] = XH[lev]*b;    
		dc_i5_33[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_i5_34[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_i5_35[lev] = XH[lev]*a + YH[lev]*alpha;
		dc_i5_36[lev] = XH[lev]*b;							
		
		/*coeff_d<<dc_b5_11[lev]<<"\n"; 
		coeff_d<<dc_b5_12[lev]<<"\n"; 
		coeff_d<<dc_b5_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb5_21[lev]<<"\n"; 
		coeff_d<<dc_nb5_22[lev]<<"\n"; 
		coeff_d<<dc_nb5_23[lev]<<"\n"; 
		coeff_d<<dc_nb5_24[lev]<<"\n"; 
		
		coeff_d<<dc_i5_32[lev]<<"\n"; 
		coeff_d<<dc_i5_33[lev]<<"\n"; 
		coeff_d<<dc_i5_34[lev]<<"\n"; 
		coeff_d<<dc_i5_35[lev]<<"\n";
		coeff_d<<dc_i5_36[lev]<<"\n";
		
		if(lev==0)
		{
		  sum_coeffs = sum_coeffs + dc_b5_11[lev] + 2.*dc_b5_12[lev] + 2.*dc_b5_13[lev] + dc_nb5_21[lev] + dc_nb5_22[lev] + dc_nb5_23[lev] + 2.*dc_nb5_24[lev] + dc_i5_32[lev] + dc_i5_33[lev] + dc_i5_34[lev] + dc_i5_35[lev] + dc_i5_36[lev];		  
		}*/

		//j+1

		dc_b6_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_b6_12[lev] = a*XA[lev] + YA[lev]*alpha; 
		dc_b6_13[lev] = b*XA[lev]; //Two more periodic entries at the end of this row

		dc_nb6_21[lev] = a*XA[lev] + YA[lev]*alpha; 
		dc_nb6_22[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_nb6_23[lev] = a*XA[lev] + YA[lev]*alpha; 
		dc_nb6_24[lev] = XA[lev]*b; //One periodic entry at the end of this row

		dc_i6_32[lev] = XA[lev]*b;    
		dc_i6_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i6_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i6_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i6_36[lev] = XA[lev]*b;
				
		/*coeff_d<<dc_b6_11[lev]<<"\n"; 
		coeff_d<<dc_b6_12[lev]<<"\n"; 
		coeff_d<<dc_b6_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb6_21[lev]<<"\n"; 
		coeff_d<<dc_nb6_22[lev]<<"\n"; 
		coeff_d<<dc_nb6_23[lev]<<"\n"; 
		coeff_d<<dc_nb6_24[lev]<<"\n"; 
		
		coeff_d<<dc_i6_32[lev]<<"\n"; 
		coeff_d<<dc_i6_33[lev]<<"\n"; 
		coeff_d<<dc_i6_34[lev]<<"\n"; 
		coeff_d<<dc_i6_35[lev]<<"\n";
		coeff_d<<dc_i6_36[lev]<<"\n";
		
		if(lev==0)
		{
			sum_coeffs = sum_coeffs + dc_b6_11[lev] + 2.*dc_b6_12[lev] + 2.*dc_b6_13[lev] + dc_nb6_21[lev] + dc_nb6_22[lev] + dc_nb6_23[lev] + 2.*dc_nb6_24[lev] + dc_i6_32[lev] + dc_i6_33[lev] + dc_i6_34[lev] + dc_i6_35[lev] + dc_i6_36[lev];			
		}*/

		//j+2 

		dc_b7_11[lev] =  YB[lev] ; 
		dc_b7_12[lev] =  YB[lev]*alpha; // One periodic entry at the end of the row with entry same as this

		dc_nb7_21[lev] = YB[lev]*alpha; 
		dc_nb7_22[lev] = YB[lev]; 
		dc_nb7_23[lev] = YB[lev]*alpha; 
		dc_nb7_24[lev] = 0.0; 

		dc_i7_32[lev] = 0.0;    
		dc_i7_33[lev] = YB[lev]*alpha; 
		dc_i7_34[lev] = YB[lev]; 
		dc_i7_35[lev] = YB[lev]*alpha; 
		dc_i7_36[lev] = 0.0; 
		
		/*coeff_d<<dc_b7_11[lev]<<"\n"; 
		coeff_d<<dc_b7_12[lev]<<"\n"; 
		coeff_d<<dc_b7_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb7_21[lev]<<"\n"; 
		coeff_d<<dc_nb7_22[lev]<<"\n"; 
		coeff_d<<dc_nb7_23[lev]<<"\n"; 
		coeff_d<<dc_nb7_24[lev]<<"\n"; 
		
		coeff_d<<dc_i7_32[lev]<<"\n"; 
		coeff_d<<dc_i7_33[lev]<<"\n"; 
		coeff_d<<dc_i7_34[lev]<<"\n"; 
		coeff_d<<dc_i7_35[lev]<<"\n";
		coeff_d<<dc_i7_36[lev]<<"\n";*/
		
		
		//Periodic entry 
		
		dc_pb7_11[lev] = YB[lev];  
		dc_pb7_12[lev] = YB[lev]*alpha; 
		dc_pb7_13[lev] = YB[lev]*alpha; 
						
		dc_pb7_21[lev] = YB[lev]*alpha; 
		dc_pb7_22[lev] = YB[lev]; 
		dc_pb7_23[lev] = YB[lev]*alpha; 
		dc_pb7_24[lev] = 0.0; 
		
		dc_pb7_32[lev] = 0.0; 
		dc_pb7_33[lev] = YB[lev]*alpha;
		dc_pb7_34[lev] = YB[lev]; 
		dc_pb7_35[lev] = YB[lev]*alpha; 
		dc_pb7_36[lev] = 0.0; 
				
		
											
   	      /******************************Case-3 (j=3 and the interior rows)***************************/
   	      //These coefficients were not modified for periodic case earlier 
   	      //Modifying them now on Feb 15, 2014 
   	         	 

	        //j-2 
	        
		dc_b8_11[lev] =  YB[lev]; 
		dc_b8_12[lev] =  YB[lev]*alpha; 
		dc_b8_13[lev] =  YB[lev]*alpha;  //This is the periodic entry at the end of the row

		dc_nb8_21[lev] = YB[lev]*alpha ; 
		dc_nb8_22[lev] = YB[lev]; 
		dc_nb8_23[lev] = YB[lev]*alpha; 
		dc_nb8_24[lev] = 0.0; 

		dc_i8_32[lev] = 0.0;    
		dc_i8_33[lev] = YB[lev]*alpha; 
		dc_i8_34[lev] = YB[lev]; 
		dc_i8_35[lev] = YB[lev]*alpha;
		dc_i8_36[lev] = 0.0;
		
		/*coeff_d<<dc_b8_11[lev]<<"\n"; 
		coeff_d<<dc_b8_12[lev]<<"\n"; 
		coeff_d<<dc_b8_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb8_21[lev]<<"\n"; 
		coeff_d<<dc_nb8_22[lev]<<"\n"; 
		coeff_d<<dc_nb8_23[lev]<<"\n"; 
		coeff_d<<dc_nb8_24[lev]<<"\n"; 
		
		coeff_d<<dc_i8_32[lev]<<"\n"; 
		coeff_d<<dc_i8_33[lev]<<"\n"; 
		coeff_d<<dc_i8_34[lev]<<"\n"; 
		coeff_d<<dc_i8_35[lev]<<"\n";
		coeff_d<<dc_i8_36[lev]<<"\n";*/
	
		//j-1

		dc_b9_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_b9_12[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_b9_13[lev] = b*XA[lev]; //Two more periodic entries  in this row

		dc_nb9_21[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_nb9_22[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_nb9_23[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_nb9_24[lev] = XA[lev]*b; //One more periodic entry in this row

		dc_i9_32[lev] = XA[lev]*b;    
		dc_i9_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i9_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i9_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i9_36[lev] = XA[lev]*b;
		
		/*coeff_d<<dc_b9_11[lev]<<"\n"; 
		coeff_d<<dc_b9_12[lev]<<"\n"; 
		coeff_d<<dc_b9_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb9_21[lev]<<"\n"; 
		coeff_d<<dc_nb9_22[lev]<<"\n"; 
		coeff_d<<dc_nb9_23[lev]<<"\n"; 
		coeff_d<<dc_nb9_24[lev]<<"\n"; 
		
		coeff_d<<dc_i9_32[lev]<<"\n"; 
		coeff_d<<dc_i9_33[lev]<<"\n"; 
		coeff_d<<dc_i9_34[lev]<<"\n"; 
		coeff_d<<dc_i9_35[lev]<<"\n";
		coeff_d<<dc_i9_36[lev]<<"\n";*/
	
		//j 

		dc_b10_11[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_b10_12[lev] = a*XH[lev] + alpha*YH[lev]; 
		dc_b10_13[lev] = XH[lev]*b; //Two more periodic entries in this row 

		dc_nb10_21[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_nb10_22[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_nb10_23[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_nb10_24[lev] = XH[lev]*b; //One more periodic entry in this row

		dc_i10_32[lev] = XH[lev]*b;    
		dc_i10_33[lev] = XH[lev]*a + YH[lev]*alpha; 
		dc_i10_34[lev] = -2.*(a+b)*XH[lev] + YH[lev]; 
		dc_i10_35[lev] = XH[lev]*a + YH[lev]*alpha;
		dc_i10_36[lev] = XH[lev]*b;
		
		/*coeff_d<<dc_b10_11[lev]<<"\n"; 
		coeff_d<<dc_b10_12[lev]<<"\n"; 
		coeff_d<<dc_b10_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb10_21[lev]<<"\n"; 
		coeff_d<<dc_nb10_22[lev]<<"\n"; 
		coeff_d<<dc_nb10_23[lev]<<"\n"; 
		coeff_d<<dc_nb10_24[lev]<<"\n"; 
		
		coeff_d<<dc_i10_32[lev]<<"\n"; 
		coeff_d<<dc_i10_33[lev]<<"\n"; 
		coeff_d<<dc_i10_34[lev]<<"\n"; 
		coeff_d<<dc_i10_35[lev]<<"\n";
		coeff_d<<dc_i10_36[lev]<<"\n";*/
		
		//j+1 

		dc_b11_11[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_b11_12[lev] = a*XA[lev] + alpha*YA[lev]; 
		dc_b11_13[lev] = b*XA[lev];  //Two more periodic entries in this row

		dc_nb11_21[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_nb11_22[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_nb11_23[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_nb11_24[lev] = XA[lev]*b; //One more periodic entry in this row

	        dc_i11_32[lev] = XA[lev]*b;    
		dc_i11_33[lev] = XA[lev]*a + YA[lev]*alpha; 
		dc_i11_34[lev] = -2.*(a+b)*XA[lev] + YA[lev]; 
		dc_i11_35[lev] = XA[lev]*a + YA[lev]*alpha;
		dc_i11_36[lev] = XA[lev]*b;
		
		/*coeff_d<<dc_b11_11[lev]<<"\n"; 
		coeff_d<<dc_b11_12[lev]<<"\n"; 
		coeff_d<<dc_b11_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb11_21[lev]<<"\n"; 
		coeff_d<<dc_nb11_22[lev]<<"\n"; 
		coeff_d<<dc_nb11_23[lev]<<"\n"; 
		coeff_d<<dc_nb11_24[lev]<<"\n"; 
		
		coeff_d<<dc_i11_32[lev]<<"\n"; 
		coeff_d<<dc_i11_33[lev]<<"\n"; 
		coeff_d<<dc_i11_34[lev]<<"\n"; 
		coeff_d<<dc_i11_35[lev]<<"\n";
		coeff_d<<dc_i11_36[lev]<<"\n";*/
		
		//j+2 

		dc_b12_11[lev] =  YB[lev];
		dc_b12_12[lev] =  YB[lev]*alpha; 
		dc_b12_13[lev] =  YB[lev]*alpha; //This is the periodic entry

		dc_nb12_21[lev] = YB[lev]*alpha ; 
		dc_nb12_22[lev] = YB[lev]; 
		dc_nb12_23[lev] = YB[lev]*alpha; 
		dc_nb12_24[lev] = 0.0; 

		dc_i12_32[lev] = 0.0;    
		dc_i12_33[lev] = YB[lev]*alpha; 
		dc_i12_34[lev] = YB[lev]; 
		dc_i12_35[lev] = YB[lev]*alpha;
		dc_i12_36[lev] = 0.0;		
		
		/*coeff_d<<dc_b12_11[lev]<<"\n"; 
		coeff_d<<dc_b12_12[lev]<<"\n"; 
		coeff_d<<dc_b12_13[lev]<<"\n"; 
		
		coeff_d<<dc_nb12_21[lev]<<"\n"; 
		coeff_d<<dc_nb12_22[lev]<<"\n"; 
		coeff_d<<dc_nb12_23[lev]<<"\n"; 
		coeff_d<<dc_nb12_24[lev]<<"\n"; 
		
		coeff_d<<dc_i12_32[lev]<<"\n"; 
		coeff_d<<dc_i12_33[lev]<<"\n"; 
		coeff_d<<dc_i12_34[lev]<<"\n"; 
		coeff_d<<dc_i12_35[lev]<<"\n";
		coeff_d<<dc_i12_36[lev]<<"\n";*/
			
		//coeff_c.close(); 
		//coeff_d.close();
		
		}
}

/**************************************************************************************************/
/******************************Multigrid interpolating function************************************/
/**************************************************************************************************/

//The levels go like this  level[0]:the finest , level[mg_levels-1]: the coarsest 
//This is always acted upon from coarse to fine grids. So when the integer is passed it is an indication to pass the data onto the finer grid 
void mg_interpolate(vector<mg_grid>& level,int lev)
{
	//cout<<"In the interpolate function\n";  
	int ind; 

	int sp = pow(2,lev);     // Spacing for the level "lev"
	int sp1 = pow(2,lev-1);  
	
	int en_inx = nx-sp, en_iny = ny-sp; 

        /*Actual loop*/
        
	for(int j=0;j<=en_iny;j=j+sp1)
	{
		for(int i=0;i<=en_inx;i=i+sp1)
		{
			ind = i + j*str_x; 

			if( (i%sp!=0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.25*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] + level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] ) ;  //Checked for PBC

			if( (i%sp!=0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.5*( level[lev].phi_s[ind+sp1] + level[lev].phi_s[ind-sp1] ); //Checked for PBC

			if( (i%sp==0)&&(j%sp!=0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + 0.5*( level[lev].phi_s[ind+sp1*str_x] + level[lev].phi_s[ind-sp1*str_x] ); //Checked for PBC

			if( (i%sp==0)&&(j%sp==0) ) level[lev-1].phi_s[ind] = level[lev-1].phi_s[ind] + level[lev].phi_s[ind];  //Checked for PBC
		}
	}	
}

/**************************************************************************************************/
/********************************Multigrid Restricting function************************************/
//Restricts field values from fine grid to coarse grid 
void mg_restrict(vector<mg_grid>& level, int lev)
{
	//cout<<"Transfering data from grid level "<<lev<<" to "<<lev+1<<"\n"; 
	
	int sp = pow(2,lev); 
	int sp1 = pow(2,lev+1); //Restricting the field data on a finer grid to a coarser grid (higher level no.) 
	
	int st_inx = 0, st_iny = 0; 
	int en_inx = nx-sp, en_iny = ny-sp; 		
	int ind; 
	
	//Full weighted operator for the interior cells of the coarse grid 
	//Modified for the boundary periodicity 
	
	for(int j=0;j<=en_iny;j=j+sp)
	{
		for(int i=0;i<=en_inx;i=i+sp)
		{
			ind = i+j*str_x; 
			
			if( (i%sp1==0) && (j%sp1==0) )
			{				
				if( i!=0 && j!=0 && i!=en_inx && j!=en_iny )
				{					
					level[lev+1].rhs[ind] = ( level[lev].res[ind-sp-sp*str_x] + level[lev].res[ind-sp+sp*str_x] + level[lev].res[ind+sp-sp*str_x] + level[lev].res[ind+sp+sp*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[ind-sp] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator 					
				}
				
				if( i==0 && j!=0 && i!=en_inx && j!=en_iny)
				{
					level[lev+1].rhs[ind] = ( level[lev].res[en_inx+j*str_x-sp*str_x] + level[lev].res[en_inx+j*str_x+sp*str_x] + level[lev].res[ind+sp-sp*str_x] + level[lev].res[ind+sp+sp*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[en_inx+j*str_x] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator			
				}
				
				if( i==en_inx && i!=0 && j!=0 && j!=en_iny)
				{
					level[lev+1].rhs[ind] = ( level[lev].res[ind-sp-sp*str_x] + level[lev].res[ind-sp+sp*str_x] + level[lev].res[st_inx+j*str_x-sp*str_x] + level[lev].res[st_inx+j*str_x+sp*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[ind-sp] + level[lev].res[st_inx+j*str_x] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator			
				}
				
				if( j==0 && j!=en_iny && i!=0 && i!=en_inx)
				{
					level[lev+1].rhs[ind] = ( level[lev].res[i-sp+en_iny*str_x] + level[lev].res[ind-sp+sp*str_x] + level[lev].res[i+sp+en_iny*str_x] + level[lev].res[ind+sp+sp*str_x] + 2.*( level[lev].res[i+en_iny*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[ind-sp] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator				
				}
				
				if( j==en_iny && j!=0 && i!=0 && i!=en_inx)
				{
					level[lev+1].rhs[ind] = ( level[lev].res[ind-sp-sp*str_x] + level[lev].res[i-sp] + level[lev].res[ind+sp-sp*str_x] + level[lev].res[i+sp] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[i] + level[lev].res[ind-sp] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator				
				}
				
				if( (i==0) && (j==0) )
				{
					level[lev+1].rhs[ind] = ( level[lev].res[en_inx+en_iny*str_x] + level[lev].res[en_inx+j*str_x+sp*str_x] + level[lev].res[st_inx+sp+en_iny*str_x] + level[lev].res[ind+sp+sp*str_x] + 2.*( level[lev].res[st_inx+en_iny*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[en_inx+st_iny*str_x] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator
				}
				
				if( (i==0) && (j==en_iny) )
				{
					level[lev+1].rhs[ind] = ( level[lev].res[en_inx+j*str_x-sp*str_x] + level[lev].res[en_inx+st_iny*str_x] + level[lev].res[ind+sp-sp*str_x] + level[lev].res[st_inx+sp+st_iny*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[st_inx+st_iny*str_x] + level[lev].res[en_inx+j*str_x] + level[lev].res[ind+sp] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator	
				}
				
				if( (i==en_inx) && (j==0) )
				{
					level[lev+1].rhs[ind] = ( level[lev].res[en_inx-sp+en_iny*str_x] + level[lev].res[ind-sp+sp*str_x] + level[lev].res[st_inx+en_iny*str_x] + level[lev].res[st_inx+sp*str_x] + 2.*( level[lev].res[en_inx+en_iny*str_x] + level[lev].res[ind+sp*str_x] + level[lev].res[ind-sp] + level[lev].res[st_inx+st_iny*str_x] ) + 4.*( level[lev].res[ind] ) )*(1./16.) ; //Full weighted operator			
				}
				
				if( (i==en_inx) && (j==en_iny) )
				{
					level[lev+1].rhs[ind] = ( level[lev].res[ind-sp-sp*str_x] + level[lev].res[en_inx-sp+st_iny*str_x] + level[lev].res[st_inx+(en_iny-sp)*str_x] + level[lev].res[st_inx+st_iny*str_x] + 2.*( level[lev].res[ind-sp*str_x] + level[lev].res[en_inx+st_iny*str_x] + level[lev].res[ind-sp] + level[lev].res[st_inx + en_iny*str_x] ) + 4.*(level[lev].res[ind]) )*(1./16.) ; //Full weighted operator
				}				
			}				
		}	
	}
							
}

/********************************************************************/
/*Updating the pressures after relaxations are done in a multi-grid fashion*/
/********************************************************************/
void mg_final(vector<mg_grid>& level, vector<fval>& fvar, int lev)
{	
	for(int i=0;i<tot_p;i++)
	{
		fvar[i].u[2] = level[lev].phi_s[i]; 				
	}
}

/********************************************************************/
/*Reading correction coefficients for various grids 
  from the null space vector of MATLAB*/ 
/********************************************************************/
void mg_read_coeff(vector<mg_grid>& level)
{
  int ind,st_inx,st_iny,step_x,step_y;   
  int en_inx, en_iny, sp; 
  ifstream coeff_input[mg_levels];
  //ofstream coeff_check;  
  //coeff_check.open("coeff_check.dat"); 
  
  string base_file("coeff_from_lab_"),file; 
  stringstream tag;             
      
  for(int lev=0;lev<mg_levels;lev++)
  {
     	//cout<<"Reading coefficients at level "<<lev<<"\n";
     
  	st_inx = 0; 
        st_iny = 0;
        
        sp = pow(2,lev); 
        
        en_inx = nx-sp; 
        en_iny = ny-sp; 
        
        step_x = pow(2,lev); 
  	step_y = pow(2,lev);
  	
  	tag<<lev+1;  	
  	file = base_file + tag.str() + ".dat";   
  		
	coeff_input[lev].open(file.c_str()); 
	
	for(int j=st_inx;j<=en_iny;j=j+step_x)
  	{
	    for(int i=st_iny;i<=en_inx;i=i+step_y)
    	    {
		ind = i + j*str_x;       
		coeff_input[lev]>>level[lev].coeff[ind];      		      		   				      		
    	    }  
  	}
  	
  	tag.str("");
  	coeff_input[lev].close();   
  }
    	//coeff_check.close(); 
}

/********************************************************************************************/
/*************************Computing residual for just the interior grid**********************/

void mg_residual_neu(vector<mg_grid>& level, int lev)
{
   int sp = pow(2,lev), ind; 
   
   int st_iny = 0, en_iny = ny-sp; 
   int st_inx = 0, en_inx = nx-sp; 
   
   double LHS, RHS;    	
   double RHS1, RHS2, RHS3, RHS4, RHS5, RHS6; 
      
    for(int j=st_iny;j<=en_iny;j=j+sp)
    {
      for(int i=st_inx;i<=en_inx;i=i+sp)
      {	        
	ind = i + j*str_x; 
	
       /*******************Case-1 j=1***********************/ 
	if(j==st_iny)
	{
	  if(i==st_inx) //Checked 
	  {	        	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b1_13[lev]*level[lev].phi_s[en_inx-sp + j*str_x] + dc_b1_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[en_inx+j*str_x+sp*str_x];	     	    
	   
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[en_inx + (j+2*sp)*str_x]; 	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i+sp+(en_iny-sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[en_inx+(en_iny-sp)*str_x]; 
	    
	    RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[en_inx-sp+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[en_inx+en_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 	    	   	    
	    
	    level[lev].res[ind] =  level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc -  dc_b1_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind];    
	    
	    //Checking the indexing error now 
	    
	    /*
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 	    	    
	    cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+j*str_x+sp*str_x<<"\n"<<en_inx+j*str_x+sp*str_x<<"\n";	    	    
	    cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+j*str_x+2*sp*str_x<<"\n";	    
	    cout<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"<<en_inx+(en_iny-sp)*str_x<<"\n";	    
	    cout<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<i+2*sp+en_iny*str_x<<"\n"<<en_inx-sp+en_iny*str_x<<"\n"<<en_inx+en_iny*str_x<<"\n"; 	    
	    */
	    
	    //cout<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b2_11[lev]<<" \n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n";      s	    
	  }
	  else if(i==st_inx+sp) //checked 
	  {    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb1_24[lev]*level[lev].phi_s[en_inx+j*str_x];
	       	    
   	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
   	    
   	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x]; 
   	       	      	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x-sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x+sp]; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(en_iny)*str_x-sp] + dc_ny_22[lev]*level[lev].phi_s[i+(en_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(en_iny)*str_x+sp] + dc_ny_24[lev]*level[lev].phi_s[i+(en_iny)*str_x+2*sp] + dc_ny_24[lev]*level[lev].phi_s[en_inx+(en_iny)*str_x]; 
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
	    	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind] ; 	     
	    
	    /*cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 	    
	   
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+j*str_x+sp*str_x<<"\n"; 
	   
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	   
	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"; 
	   
	    cout<<i-sp+(en_iny)*str_x<<"\n"<<i+(en_iny)*str_x<<"\n"<<i+sp+(en_iny)*str_x<<"\n"<<i+2*sp+(en_iny)*str_x<<"\n"<<en_inx+(en_iny)*str_x<<"\n";*/
	    
	    //cout<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n";	 
	  }
	 else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS2 = dc_i2_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS3 = dc_i3_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x]; 
	        
	    RHS4 = dc_ny1_33[lev]*level[lev].phi_s[i-sp+(en_iny-sp)*str_x] + dc_ny1_34[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_35[lev]*level[lev].phi_s[i+sp+(en_iny-sp)*str_x];
	    
	    RHS5 = dc_ny_32[lev]*level[lev].phi_s[i-2*sp+en_iny*str_x] + dc_ny_33[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_ny_34[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_35[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_ny_36[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 	    	    	   
	    
	    //cout<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i3_32[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i3_36[lev]<<"\n"<<dc_ny1_32[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny1_36[lev]<<"\n"<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n";         
	    
	    /*cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"; 
	    cout<<i-2*sp+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<i+2*sp+en_iny*str_x<<"\n"; 
	    cout<<"------------------------------------------------------------------------------------\n"; */ 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] -RHS ) - level[lev].point_correc  - level[lev].phi_s[ind]*dc_i1_34[lev]*level[lev].coeff[ind] ; 	    
	  }
	  else if(i==en_inx-sp) //Checked 
	  {	    	    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb1_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x];
	    	    	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x+sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x-sp]; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(en_iny)*str_x+sp] + dc_ny_22[lev]*level[lev].phi_s[i+(en_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(en_iny)*str_x-sp] + dc_ny_24[lev]*level[lev].phi_s[i+(en_iny)*str_x-2*sp] + dc_ny_24[lev]*level[lev].phi_s[st_inx+(en_iny)*str_x];
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	   	
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind];
	    	    
 	    /*cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
 	    cout<<st_inx+j*str_x+sp*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
 	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";  
 	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n";
 	    cout<<st_inx+en_iny*str_x<<"\n"<<i+(en_iny)*str_x-2*sp<<"\n"<<i+(en_iny)*str_x-sp<<"\n"<<i+(en_iny)*str_x<<"\n"<<i+(en_iny)*str_x+sp<<"\n"; */
	    	    
	    //cout<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb2_24[lev]<< "\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb3_24[lev]<< "\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_ny1_24[lev]<< "\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny_24[lev]<< "\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev];    	        	    
	  }
	  else //Checked 
	  {	    	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b1_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x]; 	   	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i-sp+(en_iny-sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[st_inx+(en_iny-sp)*str_x]; 
	    
	    RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i-2*sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[st_inx+sp+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[st_inx+en_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind]; 	    
	    
	   /* cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"; 
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    cout<<st_inx+(en_iny-sp)*str_x<<"\n"<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"; 
	    cout<<st_inx+en_iny*str_x<<"\n"<<st_inx+sp+en_iny*str_x<<"\n"<<i-2*sp+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n";*/	     
	    
	    //cout<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b2_12[lev]<<" \n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_11[lev]<<"\n";  		    	    	    	     	    
	  }	  	
	}	
	
	/*****************Case-2 (j=2)**********************/ 
	
	if(j==st_iny+sp)
	{
	  if(i==st_inx) //Checked 
	  {	    	    	    	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];  
	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b5_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	     
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[en_inx+(j+2*sp)*str_x]; 
	    	    	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_13[lev]*level[lev].phi_s[en_inx + en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	      	  	    	     
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";
	    cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+(j+2*sp)*str_x<<"\n";
	    cout<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<en_inx+en_iny*str_x<<"\n";  */ 
	    	    
	    //cout<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b5_11[lev]<<" \n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_13[lev]<<"\n"<<dc_b7_13[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_13[lev]<<"\n";  	    
	  }
	  else if(i==st_inx+sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb5_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_24[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x] + dc_pb7_24[lev]*level[lev].phi_s[en_inx+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	       	
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind]; 
	    
          /*cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";*/	            
	    
	    //cout<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n "<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n";	    	        
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
	  {	    	    	    
	    RHS1 = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS2 = dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS3 = dc_i6_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS4 = dc_i7_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS5 = dc_pb7_33[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_34[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_35[lev]*level[lev].phi_s[i+sp+en_iny*str_x]; 
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    	  	    	    	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_i5_34[lev]*level[lev].coeff[ind]; 
	    
	   /* cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	    cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n";
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n";  
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";    
	    
	    cout<<"------------------------------------------------------------------------------\n";*/ 	    
	    
	   //cout<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n "<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i7_32[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i7_36[lev]<<"\n"<<dc_pb7_32[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_pb7_36[lev]<<"\n";    	    
	  }
	  else if(i==en_inx-sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb5_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x]; 
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i-sp+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    	    	  	        
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";*/      	        
	    	    	    	    
	   //cout<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n "<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n "<<dc_nb7_21[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n";	       
	  }
	  else //checked 
	  {	    	     	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b5_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x];	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x]; 
	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[st_inx+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    	
	    
	    level[lev].res[ind] =  level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind]; 	    
	    
          /*cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n";
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n";
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    cout<<st_inx+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n";*/  	    
	    
	    //cout<<"The index in super matrix is "<<ind<<"\n";
	    	    
	    //cout<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n "<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n";  	    	    
	  }  
	}
	
	/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
	if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	{
	  if(i==st_inx) //Checked 
	  {	    	     	    
	    RHS1 = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];  
	    
	    RHS2 = dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b9_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];
	    
	    RHS3 = dc_b10_12[lev]*level[lev].phi_s[ind+sp] + dc_b10_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b10_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b10_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS4 = dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b11_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 	    	    	    
	    RHS5 = dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[en_inx+(j+2*sp)*str_x];	   
	     
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    		           
	    
	   //cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<" \n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n "<<dc_b11_12[lev]<<"\n"q<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n";	     
	   
	   /*cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n"; 
	   cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n";
	   cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	   cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	   cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+(j+2*sp)*str_x<<"\n";   
	   
	   cout<<"------------------------------------------------------------------------------\n";*/	      
	   	   
	   level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind]; 	      	       	   	      	    
	   //cout<<"--------------------------------------------------------"<<"\n";	    
	  }
	  else if(i==st_inx+sp)  //Checked 
	  {	    	    	    
	    RHS1 = dc_nb8_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS2 = dc_nb9_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];
	    
	    RHS3 = dc_nb10_21[lev]*level[lev].phi_s[ind-sp] + dc_nb10_23[lev]*level[lev].phi_s[ind+sp] + dc_nb10_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb10_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS4 = dc_nb11_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS5 = dc_nb12_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb10_22[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    
	    cout<<"-------------------------------------------------\n"; */
	    
	    //cout<<dc_nb8_21[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb9_21[lev]<<" \n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n "<<dc_nb12_21[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_24[lev]<<"\n"<<dc_nb12_24[lev]<<"\n";	       	    
	   //cout<<"--------------------------------------------------------"<<"\n";	    
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)  //Checked 
	  {    	    
	    RHS1 = dc_i8_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i8_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i8_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS2 = dc_i9_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i9_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i9_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i9_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i9_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS3 = dc_i10_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i10_33[lev]*level[lev].phi_s[ind-sp] +  dc_i10_35[lev]*level[lev].phi_s[ind+sp] + dc_i10_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS4 = dc_i11_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i11_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i11_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i11_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i11_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS5 = dc_i12_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i12_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i12_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_i10_34[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	    cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    
	    cout<<"-------------------------------------------------\n";*/
	    
	    //cout<<dc_i8_33[lev]<<"\n"<<dc_i8_34[lev]<<"\n"<<dc_i8_35[lev]<<"\n"<<dc_i9_32[lev]<<"\n"<<dc_i9_33[lev]<<"\n"<<dc_i9_34[lev]<<"\n"<<dc_i9_35[lev]<<"\n"<<dc_i9_36[lev]<<"\n"<<dc_i10_32[lev]<<"\n"<<dc_i10_33[lev]<<"\n"<<dc_i10_34[lev]<<"\n"<<dc_i10_35[lev]<<"\n"<<dc_i10_36[lev]<<"\n"<<dc_i11_32[lev]<<"\n"<<dc_i11_33[lev]<<"\n"<<dc_i11_34[lev]<<"\n"<<dc_i11_35[lev]<<"\n"<<dc_i11_36[lev]<<"\n"<<dc_i12_33[lev]<<"\n"<<dc_i12_34[lev]<<"\n"<<dc_i12_35[lev]<<"\n";
	    
	    //cout<<"-------------------------------------------------\n"; 	      	    	   	    	     
	  }
	  else if(i==en_inx-sp) //Checked
	  {	    	     	    
	    RHS1 = dc_nb8_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] ; 
	    
	    RHS2 = dc_nb9_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb10_21[lev]*level[lev].phi_s[ind+sp] + dc_nb10_23[lev]*level[lev].phi_s[ind-sp] + dc_nb10_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb10_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS4 = dc_nb11_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS5 = dc_nb12_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 	  
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind]  - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb10_22[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";   
	    
	    cout<<"--------------------------------------------------------"<<"\n";*/
	    
	    //cout<<dc_nb8_24[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_21[lev]<<"\n"<<dc_nb9_24[lev]<<" \n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_21[lev]<<"\n "<<dc_nb12_24[lev]<<"\n"<<dc_nb12_24[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_21[lev]<<"\n";	       	    
	    
	    //cout<<"--------------------------------------------------------"<<"\n";	    	    	    	    
	  }
	  else //Checked 
	  {    	    	    
	    RHS1 = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x]; 
	    
	    RHS2 = dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b9_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_b10_12[lev]*level[lev].phi_s[ind-sp] + dc_b10_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b10_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b10_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS4 = dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b11_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 	    
	    RHS5 = dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind]; 
	    
	    /*cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+(j)*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"; 
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    
	    cout<<"--------------------------------------------------------"<<"\n";*/
	    
	   //cout<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_11[lev]<<"\n"<<dc_b9_12[lev]<<" \n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n "<<dc_b11_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n";	       
	    
	   //cout<<"--------------------------------------------------------"<<"\n";	      	    
	  }	  	  
	}
	
	/*****************Case-4(j=ny-2)**********************/ 
	
	if(j==en_iny-sp)
	{
	  if(i==st_inx)  //checked 
	  {	     	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 

	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b5_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 

	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];
	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[en_inx+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	       	 
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind];   	    	    
	    
	   //cout<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b5_11[lev]<<" \n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"; 	     
	    
	   /*cout<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<en_inx+st_iny*str_x<<"\n"; 
	   cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n";
	   cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	   cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	   cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	   
	   cout<<"-------------------------------------------------------------------------------\n";*/	   	    	   	            
	  }
	  else if(i==st_inx+sp) //Checked
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb5_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i+sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	       
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind]; 
	     	    
	    //cout<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n "<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n";	    
	    	    	        
	    /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n";
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";*/	    
	  }
	  else if(i>st_inx+sp && i<en_inx-sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS2 = dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS3 = dc_i6_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS4 = dc_i7_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x];
	    
	    RHS5 = dc_pb7_33[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_34[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_35[lev]*level[lev].phi_s[i+sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_i5_34[lev]*level[lev].coeff[ind]; 
		    
	    //cout<<dc_pb7_32[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_pb7_36[lev]<<"\n"<<dc_i7_32[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i7_36[lev]<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n "<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n";	
	   	   
	   /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	   cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	   cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	   cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	   cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n";  
	   	   
	   cout<<"--------------------------------------------------------------------\n";*/    	  	    	    	    	   
	  }
	  else if(i==en_inx-sp) //checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb5_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x];
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i-sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    	    	
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind]; 
	    	    	              	    	    
	   //cout<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n "<<dc_nb7_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n";	    
	   
	   /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	   cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	   cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	   cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"; 
	   cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";   
	   
	   cout<<"--------------------------------------------------------------------\n";*/   
	  }
	  else //Checked 
	  {	    	    	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b5_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x];  
	    	
	    RHS5 =  dc_pb7_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[st_inx+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc - level[lev].phi_s[ind]*(dc_b5_11[lev]*level[lev].coeff[ind]); 	    
	    
	    //cout<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n "<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"; 
	    
	    /*cout<<st_inx+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"; 
	    cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n";
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n";*/  
	    	    
	  }  
	}
		
	/*****************Case-5(j=ny-1)**********************/ 
		
	if(j==en_iny) 
	{
	  if(i==st_inx) 
	  {	  	    	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b1_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x]; 	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i+sp+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[en_inx+(st_iny+sp)*str_x]; 
	    
	    RHS5 =  dc_ny_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i+2*sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[en_inx-sp+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[en_inx+st_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*(dc_b1_11[lev]*level[lev].coeff[ind]);    
	   
	    /*cout<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<i+2*sp+st_iny*str_x<<"\n"<<en_inx-sp+st_iny*str_x<<"\n"<<en_inx+st_iny*str_x<<"\n";
	    cout<<i+(st_iny+sp)*str_x<<"\n"<<i+sp+(st_iny+sp)*str_x<<"\n"<<en_inx+(st_iny+sp)*str_x<<"\n";
	    cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n";
	    cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n";*/ 	    	     
	  	    
	    //cout<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b2_11[lev]<<" \n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n";	    		       	    
	  }
	  else if(i==st_inx+sp) //Checked 
	  {	    	   	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb1_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];	       	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x-sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x+sp] ; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(st_iny)*str_x-sp] + dc_ny_22[lev]*level[lev].phi_s[i+(st_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(st_iny)*str_x+sp] ;    	    
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	     
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind]; 	 	  
	    
	    /*cout<<i+(st_iny)*str_x-sp<<"\n"<<i+(st_iny)*str_x<<"\n"<<i+(st_iny)*str_x+sp<<"\n"<<i+(st_iny)*str_x+2*sp<<"\n"<<en_inx+(st_iny)*str_x<<"\n"; 
	    cout<<i+(st_iny+sp)*str_x-sp<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x+sp<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n";*/ 	    
	  	    
	    //cout<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"; 	    	     	    	    
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
	  {	    	    	    
	     RHS1 = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp]; 
	     
	     RHS2 = dc_i2_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	     
	     RHS3 =  dc_i3_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x];   
	    
	     RHS4 = dc_ny1_33[lev]*level[lev].phi_s[i-sp+(st_iny+sp)*str_x] + dc_ny1_34[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_35[lev]*level[lev].phi_s[i+sp+(st_iny+sp)*str_x];
	     
	     RHS5 = dc_ny_32[lev]*level[lev].phi_s[i-2*sp+st_iny*str_x] + dc_ny_33[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_ny_34[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_35[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_ny_36[lev]*level[lev].phi_s[i+2*sp+st_iny*str_x];
	     
	     RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	         
	     
	     level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*(dc_i1_34[lev]*level[lev].coeff[ind]); 
	     
	     //cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_32[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny1_36[lev]<<"\n"<<dc_i3_32[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i3_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"; 
	     
	     
	    /* cout<<i-2*sp+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<i+2*sp+st_iny*str_x   <<"\n"; 
	     cout<<i-sp+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+sp+(st_iny+sp)*str_x<<"\n"; 
	     cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	     cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	     cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n";    	     
	     
	     cout<<"--------------------------------------------------------------\n";*/
	  }
	  else if(i==en_inx-sp)  //Checked
	  {    	    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb1_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] ;	    
   	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x+sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x-sp] ; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+st_iny*str_x+sp] + dc_ny_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+st_iny*str_x-sp] + dc_ny_24[lev]*level[lev].phi_s[i+st_iny*str_x-2*sp] + dc_ny_24[lev]*level[lev].phi_s[st_inx+st_iny*str_x];
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	    	    	    	    
	    
	    level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*(dc_nb1_22[lev]*level[lev].coeff[ind]) ;  	    
	    	        
	    //cout<<dc_ny_24[lev]<< "\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny1_24[lev]<< "\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<"\n"<<dc_nb3_24[lev]<< "\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb2_24[lev]<< "\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n";	    
	    
	   /*cout<<st_inx+st_iny*str_x<<"\n"<<i+st_iny*str_x-2*sp<<"\n"<<i+st_iny*str_x-sp<<"\n"<<i+st_iny*str_x<<"\n"<<i+st_iny*str_x+sp<<"\n"; 
	    cout<<i+(st_iny+sp)*str_x-sp<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x+sp<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";*/ 	    
	  } 
	  else //Checked 
	  {	  	     	     
	     RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b1_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	     
	     RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	     
	     RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x]  + dc_b3_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x];    	    
	     
	     RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i-sp+(st_iny+sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[st_inx+(st_iny+sp)*str_x]; 
	    
	     RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i-2*sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[st_inx+sp+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[st_inx+st_iny*str_x];	    

	     RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    	     
	     
	     level[lev].res[ind] = level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc - level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind];
	     
	     //level[lev].phi_s[ind] = sin(2.*M_PI*i*dx)*sin(2.*M_PI*j*dy); 	 	    
	     	     	    
	     //cout<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b2_12[lev]<<" \n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"; 	
	     
	     /*cout<<st_inx+st_iny*str_x<<"\n"<<st_inx+sp+st_iny*str_x<<"\n"<<i-2*sp+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"; 
	     cout<<st_inx+(st_iny+sp)*str_x<<"\n"<<i-sp+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x<<"\n"; 
	     cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"; 
	     cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	     cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n";*/                                           
	  }	  	
	}
		
      }    
    }

}

/****************************************************************************/
/************Computing the total RHS sum for Neumann correction**************/

void mg_compute_rhs_tot(vector<mg_grid>& level, int lev)
{
     int ind; 

     int sp = pow(2,lev);      
     
     int st_inx = 0; 
     int en_inx = nx-sp;      
     
     int st_iny = 0; 
     int en_iny = ny-sp;   
     
     /*ofstream mg_rhs_tot; 
     mg_rhs_tot.open("mg_rhs_tot.dat");
     
     cout<<"In rhs to for level "<<lev<<"\n"; */     			
     
     for(int j=st_iny;j<=en_iny;j=j+sp)
     {
    	for(int i=st_inx;i<=en_inx;i=i+sp)
    	{
	      ind = i + j*str_x; 
	      
	      if(j==st_iny)
	      {
		if(i==st_inx)
		{	  
		  level[lev].rhs[ind] =  ( level[lev].F[ind] + alpha*level[lev].F[ind+sp] + alpha*level[lev].F[en_inx+j*str_x] ) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] + alpha*level[lev].F[en_inx+j*str_x+sp*str_x] ) + alpha*( level[lev].F[i+en_iny*str_x] + alpha*level[lev].F[i+sp+en_iny*str_x] + alpha*level[lev].F[en_inx+en_iny*str_x] ); 		  
		}
		else if(i>=st_inx+sp && i<=en_inx-sp)
		{	  
		  level[lev].rhs[ind] = ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) + alpha*( alpha*level[lev].F[i-sp+en_iny*str_x] + level[lev].F[i+en_iny*str_x] + alpha*level[lev].F[i+sp+en_iny*str_x] ); 
		}
		else 
		{
		  level[lev].rhs[ind] = ( level[lev].F[ind] + alpha*level[lev].F[ind-sp] + alpha*level[lev].F[st_inx+j*str_x] ) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind-sp+sp*str_x] + alpha*level[lev].F[st_inx+j*str_x+sp*str_x] ) + alpha*( level[lev].F[i+en_iny*str_x] + alpha*level[lev].F[i-sp+en_iny*str_x] + alpha*level[lev].F[st_inx+en_iny*str_x] );			  
		}	      
	      }      
	      else if(j==st_iny+sp)
	      {
		if(i==st_inx)
		{
		  level[lev].rhs[ind] = alpha*( level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] + alpha*level[lev].F[en_inx+j*str_x-sp*str_x] ) + ( level[lev].F[ind] + alpha*level[lev].F[ind+sp] + alpha*level[lev].F[en_inx+j*str_x] ) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] + alpha*level[lev].F[en_inx+j*str_x+sp*str_x] ); 			  
		}
		else if(i>=st_inx+sp && i<=en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) ; 			  
		}		
		else 
		{
		  level[lev].rhs[ind] = alpha*( alpha*level[lev].F[st_inx+j*str_x-sp*str_x] + alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] ) + ( alpha*level[lev].F[st_inx+j*str_x] + alpha*level[lev].F[ind-sp]  + level[lev].F[ind] ) + alpha*( alpha*level[lev].F[st_inx+j*str_x+sp*str_x] + alpha*level[lev].F[ind-sp+sp*str_x]  + level[lev].F[ind+sp*str_x] );			  
		}      
	      }
	      else if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	      {
		if(i==st_inx)
		{
		  level[lev].rhs[ind] = alpha*(level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] + alpha*level[lev].F[en_inx+j*str_x-sp*str_x]) + (level[lev].F[ind] + alpha*level[lev].F[ind+sp] + alpha*level[lev].F[en_inx+j*str_x]) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] + alpha*level[lev].F[en_inx+j*str_x+sp*str_x] );		   			  
		}
		else if(i>=st_inx+sp && i<=en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] );		   			  
		}
		else 
		{
		  level[lev].rhs[ind] = alpha*( alpha*level[lev].F[st_inx+j*str_x-sp*str_x] + alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] ) + ( alpha*level[lev].F[st_inx+j*str_x] + alpha*level[lev].F[ind-sp]  + level[lev].F[ind] ) + alpha*( alpha*level[lev].F[st_inx+sp*str_x+j*str_x] + alpha*level[lev].F[ind-sp+sp*str_x]  + level[lev].F[ind+sp*str_x] ); 		  
		}	      
	      }
	      else if(j==en_iny-sp)
	      {
		if(i==st_inx)
		{	    
		  level[lev].rhs[ind] = alpha*( level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] + alpha*level[lev].F[en_inx+j*str_x-sp*str_x] ) + ( level[lev].F[ind] + alpha*level[lev].F[ind+sp] + alpha*level[lev].F[en_inx+j*str_x] ) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] + alpha*level[lev].F[en_inx+j*str_x+sp*str_x] );		  
		}
		else if(i>=st_inx+sp && i<=en_inx-sp)
		{
		  level[lev].rhs[ind] = alpha*( alpha*level[lev].F[ind-sp-sp*str_x]  + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha*( alpha*level[lev].F[ind-sp+sp*str_x] + level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind+sp+sp*str_x] ) ;			  
		}
		else 
		{		  
  		  level[lev].rhs[ind] = alpha*( level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind-sp-sp*str_x] + alpha*level[lev].F[st_inx+j*str_x-sp*str_x] ) + ( level[lev].F[ind] + alpha*level[lev].F[ind-sp] + alpha*level[lev].F[st_inx+j*str_x] ) + alpha*( level[lev].F[ind+sp*str_x] + alpha*level[lev].F[ind-sp+sp*str_x] + alpha*level[lev].F[st_inx+j*str_x+sp*str_x] ) ;		    		  
		}     
	      }
	      else 
	      {		
		if(i==st_inx)
		{	  
		  level[lev].rhs[ind] =  ( level[lev].F[ind] + alpha*level[lev].F[ind+sp] + alpha*level[lev].F[en_inx+j*str_x] ) + alpha*( level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] + alpha*level[lev].F[en_inx+j*str_x-sp*str_x] ) + alpha*( level[lev].F[i+st_iny*str_x] + alpha*level[lev].F[i+sp+st_iny*str_x] + alpha*level[lev].F[en_inx+st_iny*str_x] ); 		  
		}
		else if(i>=st_inx+sp && i<=en_inx-sp)
		{	  
		  level[lev].rhs[ind] = ( alpha*level[lev].F[ind-sp] + level[lev].F[ind] + alpha*level[lev].F[ind+sp] ) + alpha*( alpha*level[lev].F[ind-sp-sp*str_x] + level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind+sp-sp*str_x] ) + alpha*( alpha*level[lev].F[i-sp+st_iny*str_x] + level[lev].F[i+st_iny*str_x] + alpha*level[lev].F[i+sp+st_iny*str_x] ); 			  
		}
		else 
		{		  			  
		  level[lev].rhs[ind] =  ( level[lev].F[ind] + alpha*level[lev].F[ind-sp] + alpha*level[lev].F[st_inx+j*str_x] ) + alpha*( level[lev].F[ind-sp*str_x] + alpha*level[lev].F[ind-sp-sp*str_x] + alpha*level[lev].F[st_inx+j*str_x-sp*str_x] ) +  alpha*( level[lev].F[i+st_iny*str_x] + alpha*level[lev].F[i-sp+st_iny*str_x] + alpha*level[lev].F[st_inx+st_iny*str_x] );		  		  
		}		
	      }            
	      
	      //cout<<"i= "<<i<<"j= "<<j<<"	"<<level[lev].rhs[ind]<<"\n";
	    		    	
	    	//mg_rhs_tot<<i<<"	"<<j<<"		"<<level[lev].rhs[ind]<<"\n";  
	    	
	    	//cout<<ind<<"	"<<level[lev].rhs[ind]<<"\n";
	    } 
	  }
	  
	  //mg_rhs_tot.close(); 	        
}
/*************************************************************************************************/
//Gauss Seidel SOR smoothing at level "lev" and with no. of iterations specified as "ite_num"
void mg_gauss_seidel6(vector<mg_grid>& level, int lev, int ite_num, int up, int v_level, pbcs& pbc)
{
  int ind;       
  int nx_level,ny_level; 
  
  double RHS, point_correc,res_norm,err_norm=1.0; 
  double rhs_sum=0.0,tol=1.0e-2; 
  double omega=1.0; //Relaxation factor 
  double RHS1, RHS2, RHS3, RHS4, RHS5, RHS6; 
  
  int i,j;
  int sp = pow(2,lev);      
  
  int st_inx = 0;      //Changed for the periodic case
  int en_inx = nx-sp;  //Changed for the periodic case    //The right most line is left out 
  
  int st_iny = 0;      //Changed for the periodic case 
  int en_iny = ny-sp;  //Changed for the periodic case    //The top most line is left out
  
  nx_level = (nx/pow(2,lev))+1; //No. of points in X-direction at each level 
  ny_level = (ny/pow(2,lev))+1; //No. of points in Y-direction at each level
    
 // cout<<"nx_level= "<<nx_level<<"\nny_level= "<<ny_level<<"\n"; 
  
  mg_coeff();
  mg_read_coeff(level); 
  
  if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  {       
        mg_compute_rhs_tot(level,lev); 
        rhs_sum = mg_eval_rhs_sum(level,lev); 
  }  
     
  if( (lev!=v_level) && (up==0) )
  {     	
      for(int j=sp;j<ny;j=j+sp)
      {
      	for(int i=sp;i<nx;i=i+sp)
      	{
      		ind = i + j*str_x;
      		
      		level[lev].phi_s[ind] = 0.0;       		      	 	      	
      	}      
      }
   
      rhs_sum = mg_eval_rhs_sum(level,lev);               
  }
    
  level[lev].point_correc = rhs_sum/((nx_level-1)*(ny_level-1));       
 
  /*      
  cout<<"*********************************************"<<"\n"; 
  cout<<"The RHS sum is "<<rhs_sum<<"\n"; 
  //cout<<"The point correc is "<<level[lev].point_correc<<"\n";  
  cout<<"*********************************************"<<"\n";
  */
     
  for(int ite=0;ite<ite_num;ite++)
  { 
    #pragma omp parallel for default(shared) private(i,j,ind,RHS) schedule(static)
    for(j=st_iny;j<=en_iny;j=j+sp)
    {
      for(i=st_inx;i<=en_inx;i=i+sp)
      {	        	
	ind = i + j*str_x; 
	        	
	//cout<<i<<"	"<<j<<"		"<<level[lev].F[ind]<<"\n";
	
	/*******************Case-1 j=1***********************/ 
	if(j==st_iny)
	{
	  if(i==st_inx) //Checked 
	  {	        	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b1_13[lev]*level[lev].phi_s[en_inx-sp + j*str_x] + dc_b1_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[en_inx+j*str_x+sp*str_x];	     	    
	   
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[en_inx + (j+2*sp)*str_x]; 	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i+sp+(en_iny-sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[en_inx+(en_iny-sp)*str_x]; 
	    
	    RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[en_inx-sp+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[en_inx+en_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
	    	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/( dc_b1_11[lev]*level[lev].coeff[ind] ); 	
	    
	    level[lev].res[ind] = dc_b1_11[lev]*level[lev].coeff[ind]*level[lev].phi_s[ind] - (1.-omega)*level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc );    
	    
	    //Checking the indexing error now 
	    
	    /*
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 	    	    
	    cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+j*str_x+sp*str_x<<"\n"<<en_inx+j*str_x+sp*str_x<<"\n";	    	    
	    cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+j*str_x+2*sp*str_x<<"\n";	    
	    cout<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"<<en_inx+(en_iny-sp)*str_x<<"\n";	    
	    cout<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<i+2*sp+en_iny*str_x<<"\n"<<en_inx-sp+en_iny*str_x<<"\n"<<en_inx+en_iny*str_x<<"\n"; 	    
	    */
	    
	    //cout<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b2_11[lev]<<" \n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n";      s	    
	  }
	  else if(i==st_inx+sp) //checked 
	  {    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb1_24[lev]*level[lev].phi_s[en_inx+j*str_x];
	       	    
   	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
   	    
   	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x]; 
   	       	      	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x-sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x+sp]; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(en_iny)*str_x-sp] + dc_ny_22[lev]*level[lev].phi_s[i+(en_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(en_iny)*str_x+sp] + dc_ny_24[lev]*level[lev].phi_s[i+(en_iny)*str_x+2*sp] + dc_ny_24[lev]*level[lev].phi_s[en_inx+(en_iny)*str_x]; 
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb1_22[lev]*level[lev].coeff[ind] ); 	   
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind] - (1.-omega)*level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind] -  omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 	     
	    
	    /*cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 	    
	   
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+j*str_x+sp*str_x<<"\n"; 
	   
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	   
	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"; 
	   
	    cout<<i-sp+(en_iny)*str_x<<"\n"<<i+(en_iny)*str_x<<"\n"<<i+sp+(en_iny)*str_x<<"\n"<<i+2*sp+(en_iny)*str_x<<"\n"<<en_inx+(en_iny)*str_x<<"\n";*/
	    
	    //cout<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n";	 
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS2 = dc_i2_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS3 = dc_i3_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x]; 
	        
	    RHS4 = dc_ny1_33[lev]*level[lev].phi_s[i-sp+(en_iny-sp)*str_x] + dc_ny1_34[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_35[lev]*level[lev].phi_s[i+sp+(en_iny-sp)*str_x];
	    
	    RHS5 = dc_ny_32[lev]*level[lev].phi_s[i-2*sp+en_iny*str_x] + dc_ny_33[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_ny_34[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_35[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_ny_36[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
	    	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] +  omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] -RHS ) - level[lev].point_correc )/( dc_i1_34[lev]*level[lev].coeff[ind] );
	    
	    //cout<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i3_32[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i3_36[lev]<<"\n"<<dc_ny1_32[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny1_36[lev]<<"\n"<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n";         
	    
	    /*cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n"; 
	    cout<<i-2*sp+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<i+2*sp+en_iny*str_x<<"\n"; 
	    cout<<"------------------------------------------------------------------------------------\n"; */ 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_i1_34[lev]*level[lev].coeff[ind] - (1.-omega)*level[lev].phi_s[ind]*dc_i1_34[lev]*level[lev].coeff[ind] -  omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] -RHS ) - level[lev].point_correc ); 	    
	  }
	  else if(i==en_inx-sp) //Checked 
	  {	    	    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb1_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x];
	    	    	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x+sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x-sp]; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(en_iny)*str_x+sp] + dc_ny_22[lev]*level[lev].phi_s[i+(en_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(en_iny)*str_x-sp] + dc_ny_24[lev]*level[lev].phi_s[i+(en_iny)*str_x-2*sp] + dc_ny_24[lev]*level[lev].phi_s[st_inx+(en_iny)*str_x];
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/( dc_nb1_22[lev]*level[lev].coeff[ind] );	 	
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_nb1_22[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc );
	    	    
 	    /*cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
 	    cout<<st_inx+j*str_x+sp*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
 	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";  
 	    cout<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"<<i+sp+(en_iny-sp)*str_x<<"\n";
 	    cout<<st_inx+en_iny*str_x<<"\n"<<i+(en_iny)*str_x-2*sp<<"\n"<<i+(en_iny)*str_x-sp<<"\n"<<i+(en_iny)*str_x<<"\n"<<i+(en_iny)*str_x+sp<<"\n"; */
	    	    
	    //cout<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb2_24[lev]<< "\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb3_24[lev]<< "\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_ny1_24[lev]<< "\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny_24[lev]<< "\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev];    	        	    
	  }
	  else //Checked 
	  {	    	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b1_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x]; 	   	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(en_iny-sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i-sp+(en_iny-sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[st_inx+(en_iny-sp)*str_x]; 
	    
	    RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i-2*sp+en_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[st_inx+sp+en_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[st_inx+en_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind]  = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b1_11[lev]*level[lev].coeff[ind]); 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_b1_11[lev]*level[lev].coeff[ind] - omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 
	    
	    
	   /* cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"; 
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    cout<<st_inx+(en_iny-sp)*str_x<<"\n"<<i-sp+(en_iny-sp)*str_x<<"\n"<<i+(en_iny-sp)*str_x<<"\n"; 
	    cout<<st_inx+en_iny*str_x<<"\n"<<st_inx+sp+en_iny*str_x<<"\n"<<i-2*sp+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n";*/	     
	    
	    //cout<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b2_12[lev]<<" \n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_11[lev]<<"\n";  		    	    	    	     	    
	  }	  	
	}	
	
	/*****************Case-2 (j=2)**********************/ 
	
	if(j==st_iny+sp)
	{
	  if(i==st_inx) //Checked 
	  {	    	    	    	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];  
	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b5_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	     
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[en_inx+(j+2*sp)*str_x]; 
	    	    	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_13[lev]*level[lev].phi_s[en_inx + en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_b5_11[lev]*level[lev].coeff[ind] );   	  	    	     
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
	    
	    /*cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";
	    cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+(j+2*sp)*str_x<<"\n";
	    cout<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n"<<en_inx+en_iny*str_x<<"\n";  */ 
	    	    
	    //cout<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b5_11[lev]<<" \n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_13[lev]<<"\n"<<dc_b7_13[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_13[lev]<<"\n";  	    
	  }
	  else if(i==st_inx+sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb5_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_24[lev]*level[lev].phi_s[i+2*sp+en_iny*str_x] + dc_pb7_24[lev]*level[lev].phi_s[en_inx+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind]  =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_nb5_22[lev]*level[lev].coeff[ind]);   	
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind] - omega*(level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc); 
	    
          /*cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";*/	            
	    
	    //cout<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n "<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n";	    	        
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
	  {	    	    	    
	    RHS1 = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS2 = dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS3 = dc_i6_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS4 = dc_i7_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS5 = dc_pb7_33[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_34[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_35[lev]*level[lev].phi_s[i+sp+en_iny*str_x]; 
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].phi_s[ind] =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc )/( dc_i5_34[lev]*level[lev].coeff[ind] ); 	  	    	    
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_i5_34[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_i5_34[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] -RHS) - level[lev].point_correc ); 
	    
	   /* cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	    cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n";
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n";  
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";    
	    
	    cout<<"------------------------------------------------------------------------------\n";*/ 	    
	    
	   //cout<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n "<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i7_32[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i7_36[lev]<<"\n"<<dc_pb7_32[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_pb7_36[lev]<<"\n";    	    
	  }
	  else if(i==en_inx-sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb5_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x]; 
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i+sp+en_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i-sp+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].phi_s[ind]  =  (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb5_22[lev]*level[lev].coeff[ind] );	  	        
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_nb5_22[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
	    
	    /*cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";
	    cout<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n"<<i+sp+en_iny*str_x<<"\n";*/      	        
	    	    	    	    
	   //cout<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n "<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n "<<dc_nb7_21[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n";	       
	  }
	  else //checked 
	  {	    	     	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b5_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x];	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x]; 
	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i-sp+en_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[st_inx+en_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind]  = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b5_11[lev]*level[lev].coeff[ind]);	
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
	    
	    
          /*cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n";
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n";
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    cout<<st_inx+en_iny*str_x<<"\n"<<i-sp+en_iny*str_x<<"\n"<<i+en_iny*str_x<<"\n";*/  	    
	    
	    //cout<<"The index in super matrix is "<<ind<<"\n";
	    	    
	    //cout<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n "<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n";  	    	    
	  }  
	}
	
	/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
	if(j>=st_iny+2*sp && j<=en_iny-2*sp)
	{
	  if(i==st_inx) //Checked 
	  {	    	     	    
	    RHS1 = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];  
	    
	    RHS2 = dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b9_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];
	    
	    RHS3 = dc_b10_12[lev]*level[lev].phi_s[ind+sp] + dc_b10_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b10_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b10_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS4 = dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b11_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 	    	    	    
	    RHS5 = dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[en_inx+(j+2*sp)*str_x];	   
	     
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_b10_11[lev]*level[lev].coeff[ind]);   		           
	    
	   //cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<" \n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n "<<dc_b11_12[lev]<<"\n"q<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n";	     
	   
	   /*cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n"; 
	   cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n";
	   cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	   cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	   cout<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"<<en_inx+(j+2*sp)*str_x<<"\n";   
	   
	   cout<<"------------------------------------------------------------------------------\n";*/	      
	   	   
	   level[lev].res[ind] = level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind] - (1.-omega)*level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind] - + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 	      	       	   	      	    
	   //cout<<"--------------------------------------------------------"<<"\n";	    
	  }
	  else if(i==st_inx+sp)  //Checked 
	  {	    	    	    
	    RHS1 = dc_nb8_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS2 = dc_nb9_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x];
	    
	    RHS3 = dc_nb10_21[lev]*level[lev].phi_s[ind-sp] + dc_nb10_23[lev]*level[lev].phi_s[ind+sp] + dc_nb10_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb10_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS4 = dc_nb11_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS5 = dc_nb12_21[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb10_22[lev]*level[lev].coeff[ind] ); 	    	    
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb10_22[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_nb10_22[lev]*level[lev].coeff[ind] - omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    
	    cout<<"-------------------------------------------------\n"; */
	    
	    //cout<<dc_nb8_21[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb9_21[lev]<<" \n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n "<<dc_nb12_21[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_24[lev]<<"\n"<<dc_nb12_24[lev]<<"\n";	       	    
	   //cout<<"--------------------------------------------------------"<<"\n";	    
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)  //Checked 
	  {    	    
	    RHS1 = dc_i8_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i8_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i8_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS2 = dc_i9_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i9_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i9_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i9_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i9_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS3 = dc_i10_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i10_33[lev]*level[lev].phi_s[ind-sp] +  dc_i10_35[lev]*level[lev].phi_s[ind+sp] + dc_i10_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS4 = dc_i11_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i11_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i11_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i11_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i11_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS5 = dc_i12_33[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_i12_34[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_i12_35[lev]*level[lev].phi_s[ind+sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_i10_34[lev]*level[lev].coeff[ind] );  
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_i10_34[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_i10_34[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	    cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	    cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"; 
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n"; 
	    
	    cout<<"-------------------------------------------------\n";*/
	    
	    //cout<<dc_i8_33[lev]<<"\n"<<dc_i8_34[lev]<<"\n"<<dc_i8_35[lev]<<"\n"<<dc_i9_32[lev]<<"\n"<<dc_i9_33[lev]<<"\n"<<dc_i9_34[lev]<<"\n"<<dc_i9_35[lev]<<"\n"<<dc_i9_36[lev]<<"\n"<<dc_i10_32[lev]<<"\n"<<dc_i10_33[lev]<<"\n"<<dc_i10_34[lev]<<"\n"<<dc_i10_35[lev]<<"\n"<<dc_i10_36[lev]<<"\n"<<dc_i11_32[lev]<<"\n"<<dc_i11_33[lev]<<"\n"<<dc_i11_34[lev]<<"\n"<<dc_i11_35[lev]<<"\n"<<dc_i11_36[lev]<<"\n"<<dc_i12_33[lev]<<"\n"<<dc_i12_34[lev]<<"\n"<<dc_i12_35[lev]<<"\n";
	    
	    //cout<<"-------------------------------------------------\n"; 	      	    	   	    	     
	  }
	  else if(i==en_inx-sp) //Checked
	  {	    	     	    
	    RHS1 = dc_nb8_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb8_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb8_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] ; 
	    
	    RHS2 = dc_nb9_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb9_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb9_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb9_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb10_21[lev]*level[lev].phi_s[ind+sp] + dc_nb10_23[lev]*level[lev].phi_s[ind-sp] + dc_nb10_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb10_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS4 = dc_nb11_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb11_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb11_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb11_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS5 = dc_nb12_21[lev]*level[lev].phi_s[ind+sp+2*sp*str_x] + dc_nb12_22[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_nb12_23[lev]*level[lev].phi_s[ind-sp+2*sp*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind]  - RHS) - level[lev].point_correc )/(dc_nb10_22[lev]*level[lev].coeff[ind]); 	  
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_nb10_22[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*(dc_nb10_22[lev]*level[lev].coeff[ind]) - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind]  - RHS) - level[lev].point_correc ); 
	    
	    /*cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";
	    cout<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"<<ind+sp+2*sp*str_x<<"\n";   
	    
	    cout<<"--------------------------------------------------------"<<"\n";*/
	    
	    //cout<<dc_nb8_24[lev]<<"\n"<<dc_nb8_24[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_21[lev]<<"\n"<<dc_nb9_24[lev]<<" \n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_21[lev]<<"\n "<<dc_nb12_24[lev]<<"\n"<<dc_nb12_24[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_21[lev]<<"\n";	       	    
	    
	    //cout<<"--------------------------------------------------------"<<"\n";	    	    	    	    
	  }
	  else //Checked 
	  {    	    	    
	    RHS1 = dc_b8_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b8_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x]; 
	    
	    RHS2 = dc_b9_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b9_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b9_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b9_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_b10_12[lev]*level[lev].phi_s[ind-sp] + dc_b10_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b10_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b10_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS4 = dc_b11_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b11_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b11_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b11_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 	    
	    RHS5 = dc_b12_11[lev]*level[lev].phi_s[ind+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[ind-sp+2*sp*str_x] + dc_b12_12[lev]*level[lev].phi_s[st_inx+(j+2*sp)*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b10_11[lev]*level[lev].coeff[ind]);	 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind] -  (1.-omega)*level[lev].phi_s[ind]*dc_b10_11[lev]*level[lev].coeff[ind] - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
	    
	    /*cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+(j)*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"; 
	    cout<<st_inx+(j+2*sp)*str_x<<"\n"<<ind-sp+2*sp*str_x<<"\n"<<ind+2*sp*str_x<<"\n"; 
	    
	    cout<<"--------------------------------------------------------"<<"\n";*/
	    
	   //cout<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_11[lev]<<"\n"<<dc_b9_12[lev]<<" \n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n "<<dc_b11_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n";	       
	    
	   //cout<<"--------------------------------------------------------"<<"\n";	      	    
	  }	  	  
	}
	
	/*****************Case-4(j=ny-2)**********************/ 
	
	if(j==en_iny-sp)
	{
	  if(i==st_inx)  //checked 
	  {	     	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[en_inx-sp+(j+sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 

	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind+sp] + dc_b5_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b5_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 

	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];
	    
	    RHS5 = dc_pb7_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[en_inx+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_b5_11[lev]*level[lev].coeff[ind]);   	 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] - (1.-omega)*level[lev].phi_s[ind]*dc_b5_11[lev]*level[lev].coeff[ind] - omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc);   	    	    
	    
	   //cout<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b5_11[lev]<<" \n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"; 	     
	    
	   /*cout<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<en_inx+st_iny*str_x<<"\n"; 
	   cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n";
	   cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	   cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n"; 
	   cout<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx-sp+(j+sp)*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n"; 
	   
	   cout<<"-------------------------------------------------------------------------------\n";*/	   	    	   	            
	  }
	  else if(i==st_inx+sp) //Checked
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind+2*sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[en_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind-sp] + dc_nb5_23[lev]*level[lev].phi_s[ind+sp] + dc_nb5_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb5_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x]; 
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i+sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/( dc_nb5_22[lev]*level[lev].coeff[ind] );    
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*( dc_nb5_22[lev]*level[lev].coeff[ind] ) -  (1.-omega)*level[lev].phi_s[ind]*( dc_nb5_22[lev]*level[lev].coeff[ind] ) - omega*(level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 
	     	    
	    //cout<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n "<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n";	    
	    	    	        
	    /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n";
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n"; 
	    cout<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n"<<en_inx+(j+sp)*str_x<<"\n";*/	    
	  }
	  else if(i>st_inx+sp && i<en_inx-sp) //Checked 
	  {	    	    	    
	    RHS1 = dc_i4_32[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_i4_33[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_i4_34[lev]*level[lev].phi_s[ind+sp*str_x] + dc_i4_35[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_i4_36[lev]*level[lev].phi_s[ind+2*sp+sp*str_x]; 
	    
	    RHS2 = dc_i5_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i5_33[lev]*level[lev].phi_s[ind-sp] + dc_i5_35[lev]*level[lev].phi_s[ind+sp] + dc_i5_36[lev]*level[lev].phi_s[ind+2*sp]; 
	    
	    RHS3 = dc_i6_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i6_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i6_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i6_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i6_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	    
	    RHS4 = dc_i7_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i7_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i7_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x];
	    
	    RHS5 = dc_pb7_33[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_34[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_35[lev]*level[lev].phi_s[i+sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_i5_34[lev]*level[lev].coeff[ind]); 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*(dc_i5_34[lev]*level[lev].coeff[ind]) -  (1.-omega)*level[lev].phi_s[ind]*(dc_i5_34[lev]*level[lev].coeff[ind]) - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 
		    
	    //cout<<dc_pb7_32[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_pb7_36[lev]<<"\n"<<dc_i7_32[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i7_36[lev]<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n "<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n";	
	   	   
	   /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	   cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	   cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	   cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"; 
	   cout<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n"<<ind+2*sp+sp*str_x<<"\n";  
	   	   
	   cout<<"--------------------------------------------------------------------\n";*/    	  	    	    	    	   
	  }
	  else if(i==en_inx-sp) //checked 
	  {	    	    	    
	    RHS1 = dc_nb4_21[lev]*level[lev].phi_s[ind+sp+sp*str_x] + dc_nb4_22[lev]*level[lev].phi_s[ind+sp*str_x] + dc_nb4_23[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_nb4_24[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_nb5_21[lev]*level[lev].phi_s[ind+sp] + dc_nb5_23[lev]*level[lev].phi_s[ind-sp] + dc_nb5_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb5_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_nb6_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb6_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb6_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb6_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_nb7_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb7_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb7_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x];
	    
	    RHS5 = dc_pb7_21[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_pb7_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_23[lev]*level[lev].phi_s[i-sp+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb5_22[lev]*level[lev].coeff[ind] );	
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*( dc_nb5_22[lev]*level[lev].coeff[ind] ) - (1.-omega)*level[lev].phi_s[ind]*( dc_nb5_22[lev]*level[lev].coeff[ind] ) - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 	    	              	    	    
	   //cout<<dc_pb7_24[lev]<<"\n"<<dc_pb7_24[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_24[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n "<<dc_nb7_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n";	    
	   
	   /*cout<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"; 
	   cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	   cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	   cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"; 
	   cout<<st_inx+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n"<<ind+sp+sp*str_x<<"\n";   
	   
	   cout<<"--------------------------------------------------------------------\n";*/   
	  }
	  else //Checked 
	  {	    	    	    
	    RHS1 = dc_b4_11[lev]*level[lev].phi_s[ind+sp*str_x] + dc_b4_12[lev]*level[lev].phi_s[ind-sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[ind-2*sp+sp*str_x] + dc_b4_13[lev]*level[lev].phi_s[st_inx+sp+(j+sp)*str_x] + dc_b4_12[lev]*level[lev].phi_s[st_inx+(j+sp)*str_x]; 
	    
	    RHS2 = dc_b5_12[lev]*level[lev].phi_s[ind-sp] + dc_b5_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b5_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b5_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS3 = dc_b6_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b6_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b6_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b6_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS4 = dc_b7_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_b7_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x];  
	    	
	    RHS5 =  dc_pb7_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_pb7_12[lev]*level[lev].phi_s[st_inx+st_iny*str_x];
	    
	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc )/(dc_b5_11[lev]*level[lev].coeff[ind]);    	    	    	    	    	        
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*(dc_b5_11[lev]*level[lev].coeff[ind]) - (1.-omega)*level[lev].phi_s[ind]*(dc_b5_11[lev]*level[lev].coeff[ind]) -  omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 	    
	    
	    //cout<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n "<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"; 
	    
	    /*cout<<st_inx+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"; 
	    cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n";
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"; 
	    cout<<st_inx+(j+sp)*str_x<<"\n"<<st_inx+sp+(j+sp)*str_x<<"\n"<<ind-2*sp+sp*str_x<<"\n"<<ind-sp+sp*str_x<<"\n"<<ind+sp*str_x<<"\n";*/  
	    	    
	  }  
	}
		
	/***********************************Case-5(j=en_iny)******************************************/ 
		
	if(j==en_iny)
	{
	  if(i==st_inx) 
	  {
	  	    	    	    
	    RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind+sp] + dc_b1_13[lev]*level[lev].phi_s[ind+2*sp] + dc_b1_13[lev]*level[lev].phi_s[en_inx-sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[en_inx-sp+(j-sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x]; 	        
	    
	    RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i+sp+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[en_inx+(st_iny+sp)*str_x]; 
	    
	    RHS5 =  dc_ny_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i+2*sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[en_inx-sp+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[en_inx+st_iny*str_x];	    

	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_b1_11[lev]*level[lev].coeff[ind]);
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*(dc_b1_11[lev]*level[lev].coeff[ind]) -  (1.-omega)*level[lev].phi_s[ind]*(dc_b1_11[lev]*level[lev].coeff[ind]) - omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc );    
	   
	    /*cout<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<i+2*sp+st_iny*str_x<<"\n"<<en_inx-sp+st_iny*str_x<<"\n"<<en_inx+st_iny*str_x<<"\n";
	    cout<<i+(st_iny+sp)*str_x<<"\n"<<i+sp+(st_iny+sp)*str_x<<"\n"<<en_inx+(st_iny+sp)*str_x<<"\n";
	    cout<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"<<en_inx+(j-2*sp)*str_x<<"\n";
	    cout<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx-sp+(j-sp)*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx-sp+j*str_x<<"\n"<<en_inx+j*str_x<<"\n";*/ 	    	     
	  	    
	    //cout<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b2_11[lev]<<" \n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n";	    		       	    
	  }
	  else if(i==st_inx+sp) //Checked 
	  {	    	   	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind-sp] + dc_nb1_23[lev]*level[lev].phi_s[ind+sp] + dc_nb1_24[lev]*level[lev].phi_s[ind+2*sp] + dc_nb1_24[lev]*level[lev].phi_s[en_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind+2*sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[en_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[ind+2*sp-2*sp*str_x] + dc_nb3_24[lev]*level[lev].phi_s[en_inx+(j-2*sp)*str_x];	       	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x-sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x+sp] ; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+(st_iny)*str_x-sp] + dc_ny_22[lev]*level[lev].phi_s[i+(st_iny)*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+(st_iny)*str_x+sp] ;    	    
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS) - level[lev].point_correc )/( dc_nb1_22[lev]*level[lev].coeff[ind] ); 
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*( dc_nb1_22[lev]*level[lev].coeff[ind] ) - (1.-omega)*level[lev].phi_s[ind]*( dc_nb1_22[lev]*level[lev].coeff[ind] ) -  omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS) - level[lev].point_correc ); 	 	  
	    
	    /*cout<<i+(st_iny)*str_x-sp<<"\n"<<i+(st_iny)*str_x<<"\n"<<i+(st_iny)*str_x+sp<<"\n"<<i+(st_iny)*str_x+2*sp<<"\n"<<en_inx+(st_iny)*str_x<<"\n"; 
	    cout<<i+(st_iny+sp)*str_x-sp<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x+sp<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"<<en_inx+(j-sp)*str_x<<"\n"; 
	    cout<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n"<<en_inx+j*str_x<<"\n";*/ 	    
	  	    
	    //cout<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"; 	    	     	    	    
	  }
	  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
	  {	    	    	    
	     RHS1 = dc_i1_32[lev]*level[lev].phi_s[ind-2*sp] + dc_i1_33[lev]*level[lev].phi_s[ind-sp] + dc_i1_35[lev]*level[lev].phi_s[ind+sp] + dc_i1_36[lev]*level[lev].phi_s[ind+2*sp]; 
	     
	     RHS2 = dc_i2_32[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_i2_33[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_i2_34[lev]*level[lev].phi_s[ind-sp*str_x] + dc_i2_35[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_i2_36[lev]*level[lev].phi_s[ind+2*sp-sp*str_x]; 
	     
	     RHS3 =  dc_i3_33[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] + dc_i3_34[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_i3_35[lev]*level[lev].phi_s[ind+sp-2*sp*str_x];   
	    
	     RHS4 = dc_ny1_33[lev]*level[lev].phi_s[i-sp+(st_iny+sp)*str_x] + dc_ny1_34[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_35[lev]*level[lev].phi_s[i+sp+(st_iny+sp)*str_x];
	     
	     RHS5 = dc_ny_32[lev]*level[lev].phi_s[i-2*sp+st_iny*str_x] + dc_ny_33[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_ny_34[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_35[lev]*level[lev].phi_s[i+sp+st_iny*str_x] + dc_ny_36[lev]*level[lev].phi_s[i+2*sp+st_iny*str_x];
	     
	     RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	     level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc)/(dc_i1_34[lev]*level[lev].coeff[ind]);	         
	     
	     level[lev].res[ind] = level[lev].phi_s[ind]*(dc_i1_34[lev]*level[lev].coeff[ind]) -  (1.-omega)*level[lev].phi_s[ind]*(dc_i1_34[lev]*level[lev].coeff[ind]) - omega*( level[lev].coeff[ind]*(level[lev].rhs[ind] - RHS) - level[lev].point_correc); 
	     
	     //cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_32[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny1_36[lev]<<"\n"<<dc_i3_32[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i3_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"; 
	     
	     
	    /* cout<<i-2*sp+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"<<i+sp+st_iny*str_x<<"\n"<<i+2*sp+st_iny*str_x   <<"\n"; 
	     cout<<i-sp+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+sp+(st_iny+sp)*str_x<<"\n"; 
	     cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	     cout<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"<<ind+2*sp-sp*str_x<<"\n"; 
	     cout<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n"<<ind+2*sp<<"\n";    	     
	     
	     cout<<"--------------------------------------------------------------\n";*/
	  }
	  else if(i==en_inx-sp)  //Checked
	  {    	    	    	    
	    RHS1 = dc_nb1_21[lev]*level[lev].phi_s[ind+sp] + dc_nb1_23[lev]*level[lev].phi_s[ind-sp] + dc_nb1_24[lev]*level[lev].phi_s[ind-2*sp] + dc_nb1_24[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	    
	    RHS2 = dc_nb2_21[lev]*level[lev].phi_s[ind+sp-sp*str_x] + dc_nb2_22[lev]*level[lev].phi_s[ind-sp*str_x] + dc_nb2_23[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_nb2_24[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	    
	    RHS3 = dc_nb3_21[lev]*level[lev].phi_s[ind+sp-2*sp*str_x] + dc_nb3_22[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_nb3_23[lev]*level[lev].phi_s[ind-sp-2*sp*str_x] ;	    
   	    
   	    RHS4 = dc_ny1_21[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x+sp] + dc_ny1_22[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_23[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x-sp] ; 
   	    
   	    RHS5 = dc_ny_21[lev]*level[lev].phi_s[i+st_iny*str_x+sp] + dc_ny_22[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_23[lev]*level[lev].phi_s[i+st_iny*str_x-sp] + dc_ny_24[lev]*level[lev].phi_s[i+st_iny*str_x-2*sp] + dc_ny_24[lev]*level[lev].phi_s[st_inx+st_iny*str_x];
   	    
   	    RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    	    
	    level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/(dc_nb1_22[lev]*level[lev].coeff[ind]);	    	    	    
	    
	    level[lev].res[ind] = level[lev].phi_s[ind]*(dc_nb1_22[lev]*level[lev].coeff[ind]) - (1.-omega)*level[lev].phi_s[ind]*(dc_nb1_22[lev]*level[lev].coeff[ind]) - omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc );  	    
	    	        
	    //cout<<dc_ny_24[lev]<< "\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny1_24[lev]<< "\n"<<dc_ny1_24[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<"\n"<<dc_nb3_24[lev]<< "\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb2_24[lev]<< "\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n";	    
	    
	   /*cout<<st_inx+st_iny*str_x<<"\n"<<i+st_iny*str_x-2*sp<<"\n"<<i+st_iny*str_x-sp<<"\n"<<i+st_iny*str_x<<"\n"<<i+st_iny*str_x+sp<<"\n"; 
	    cout<<i+(st_iny+sp)*str_x-sp<<"\n"<<i+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x+sp<<"\n"; 
	    cout<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"<<ind+sp-2*sp*str_x<<"\n"; 
	    cout<<st_inx+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"<<ind+sp-sp*str_x<<"\n"; 
	    cout<<st_inx+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n"<<ind+sp<<"\n";*/ 	    
	  }
	  else //Checked 
	  {	  	     	     
	     RHS1 = dc_b1_12[lev]*level[lev].phi_s[ind-sp] + dc_b1_13[lev]*level[lev].phi_s[ind-2*sp] + dc_b1_13[lev]*level[lev].phi_s[st_inx+sp+j*str_x] + dc_b1_12[lev]*level[lev].phi_s[st_inx+j*str_x]; 
	     
	     RHS2 = dc_b2_11[lev]*level[lev].phi_s[ind-sp*str_x] + dc_b2_12[lev]*level[lev].phi_s[ind-sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[ind-2*sp-sp*str_x] + dc_b2_13[lev]*level[lev].phi_s[st_inx+sp+(j-sp)*str_x] + dc_b2_12[lev]*level[lev].phi_s[st_inx+(j-sp)*str_x]; 
	     
	     RHS3 = dc_b3_11[lev]*level[lev].phi_s[ind-2*sp*str_x] + dc_b3_12[lev]*level[lev].phi_s[ind-sp-2*sp*str_x]  + dc_b3_12[lev]*level[lev].phi_s[st_inx+(j-2*sp)*str_x];    	    
	     
	     RHS4 = dc_ny1_11[lev]*level[lev].phi_s[i+(st_iny+sp)*str_x] + dc_ny1_12[lev]*level[lev].phi_s[i-sp+(st_iny+sp)*str_x] + dc_ny1_13[lev]*level[lev].phi_s[st_inx+(st_iny+sp)*str_x]; 
	    
	     RHS5 = dc_ny_11[lev]*level[lev].phi_s[i+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[i-sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[i-2*sp+st_iny*str_x] + dc_ny_13[lev]*level[lev].phi_s[st_inx+sp+st_iny*str_x] + dc_ny_12[lev]*level[lev].phi_s[st_inx+st_iny*str_x];	    

	     RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
	    
	     level[lev].phi_s[ind] = (1.-omega)*level[lev].phi_s[ind] + omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc )/( dc_b1_11[lev]*level[lev].coeff[ind] );
	     
	     level[lev].res[ind] = level[lev].phi_s[ind]*( dc_b1_11[lev]*level[lev].coeff[ind] ) -  (1.-omega)*level[lev].phi_s[ind]*( dc_b1_11[lev]*level[lev].coeff[ind] ) - omega*( level[lev].coeff[ind]*( level[lev].rhs[ind] - RHS ) - level[lev].point_correc );
	     
	     //level[lev].phi_s[ind] = sin(2.*M_PI*i*dx)*sin(2.*M_PI*j*dy); 	 	    
	     	     	    
	     //cout<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_13[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b2_12[lev]<<" \n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"; 	
	     
	     /*cout<<st_inx+st_iny*str_x<<"\n"<<st_inx+sp+st_iny*str_x<<"\n"<<i-2*sp+st_iny*str_x<<"\n"<<i-sp+st_iny*str_x<<"\n"<<i+st_iny*str_x<<"\n"; 
	     cout<<st_inx+(st_iny+sp)*str_x<<"\n"<<i-sp+(st_iny+sp)*str_x<<"\n"<<i+(st_iny+sp)*str_x<<"\n"; 
	     cout<<st_inx+(j-2*sp)*str_x<<"\n"<<ind-sp-2*sp*str_x<<"\n"<<ind-2*sp*str_x<<"\n"; 
	     cout<<st_inx+(j-sp)*str_x<<"\n"<<st_inx+sp+(j-sp)*str_x<<"\n"<<ind-2*sp-sp*str_x<<"\n"<<ind-sp-sp*str_x<<"\n"<<ind-sp*str_x<<"\n"; 
	     cout<<st_inx+j*str_x<<"\n"<<st_inx+sp+j*str_x<<"\n"<<ind-2*sp<<"\n"<<ind-sp<<"\n"<<ind<<"\n";*/                                           
	  }	  	
	}
	
	//level[lev].res[ind] = level[lev].point_correc;
	
	// cout<<ind<<"	"<<level[lev].phi_s[ind]<<"\n";
	// cout<<ind<<"	"<<level[lev].coeff[ind]<<"\n"; 			
      }    
        //cout<<"------------------------------------------\n"; 
      
    }                    
    
    //	cout<<ite<<"\n";        
  }    
  
  //impose_zero_mean(level,lev);
  
  if(lev==0)
  {  	
  	mg_bcs_neu(level,lev,pbc); 
  }
   
}

/******************************************************************************/
/*******************Computing residual norm at level lev***********************/
double mg_res_norm(vector<mg_grid>& level, int lev)
{
	int sp = pow(2,lev); 
	int ind; 
	double abs_res = 0.0; 
	
	for(int j=0;j<ny;j=j+sp)
	{
		for(int i=0;i<nx;i=i+sp)
		{
			ind = i + j*str_x; 
			
			abs_res = abs_res + level[lev].res[ind]*level[lev].res[ind]; 		
		}	
	}
			
	return sqrt(abs_res);  
}

/*******************************************************************************************/
//Multigrid Poisson Solver 
void mg_poisson_solver(vector<mg_grid>& level,vector<fval>& fvar, pbcs& pbc)
{		
	//Assigning the field variable pressure to the multi-level variable at level-0 	
	double res_norm=0.0; 
	int num_v_c=1; 

	//ofstream res_out;	
	//res_out.open("mg_v_residual.dat");
	
	mg_coeff();  	      // Coefficients for smoothing operator at all levels
	mg_read_coeff(level); // Reading the null space vectors for all grid sizes
			       	
	/******************************************/		       	
		
        /*for(int v_c=0;v_c<num_v_c;v_c++)  //The basic V-cycling 
	{	
		v_cycle(level,0,pbc); 		
		res_norm = mg_res_norm(level,0); 	
		res_out<<v_c<<"		"<<log10(res_norm)<<"\n"; 			
	}*/
			
	/******************************************/
		
	fmg(level,pbc); 		
	mg_final(level,fvar,0);	
 	//res_out.close(); 
}
/*******************************************************************************************/
//Full Multi Grid (FMG) Technique 
void fmg(vector<mg_grid>& level, pbcs& pbc)
{
	/***********************/
      	//Exact computation at the coarsest level	
	
	int coar_lev = mg_levels-1; //The level numbering starts from "0"
	
	int nx_levelp = ((nx)/pow(2,coar_lev)); //Points in X-direction at the corasest level after dicarding a periodic line of values
        int ny_levelp = ((ny)/pow(2,coar_lev)); //Points in Y-direction at the corasest level after dicarding a periodic line of values
	
      //mg_evaluate_rhs(level,coar_lev,pbc);    
        mg_compute_rhs_tot(level,coar_lev);
         
        double rhs_sum = mg_eval_rhs_sum(level,coar_lev);        
        //cout<<"The rhs_sum at the start of fmg is "<<rhs_sum<<"\n";		
	level[coar_lev].point_correc = rhs_sum/( (nx_levelp)*(ny_levelp) ); 	
	arma_direct_solve(level,coar_lev);  
	mg_interpolate(level,coar_lev); 
	
	/***********************/
	int v_start = coar_lev-1;  
	
	for(int i=v_start;i>=0;i--) 
	{		
		v_cycle(level,i,pbc);
		if(i>0) mg_interpolate(level,i);	//"If" condition to make sure that the values are not interpolated from the finest grid
	}	
}	

/*******************************************************************************************/
//V-cycle required for Full Multi Grid (FMG)

void v_cycle(vector<mg_grid>& level, int v_level, pbcs& pbc)
{
	int up=0,neu=600; 
	double rhs_sum, nx_level, ny_level; 
	
	//cout<<"No. of mg levels is "<<mg_levels<<"\n"; 
	
	for(int i=v_level;i<mg_levels;i++)	//Descending the stair case //Decremented mg_levels by 1 for testing purposes 
	{
		if(i<mg_levels-1)
		{		
			//mg_gauss_seidel6(level,i,neu,up,v_level,pbc); 
			mg_conjugate_gradient(level,i,neu,up,v_level);
			mg_residual_neu(level,i);		
			mg_restrict(level,i); 				
		}
		else
		{			      
		      	nx_level = ((nx)/pow(2,i)); //No. of points in X-direction at each level after discarding a line of periodic values
			ny_level = ((ny)/pow(2,i)); //No. of points in Y-direction at each level after discarding a line of periodic values
				
		      	rhs_sum = mg_eval_rhs_sum(level,i);		      	
		      	level[i].point_correc = rhs_sum/( (nx_level)*(ny_level) );		      	 
		      	arma_direct_solve(level,i); 	
		      	//impose_zero_mean(level,i);		      			 					
		}
	}	
			
	up = 1;  //Variable indicating if the cycle is going up or down. 
			
	for(int i=mg_levels-1; i>v_level;i--)  //Ascending the stair case  //Decrementing the upper bound of i by 1 to make it mg_levels-2
	{
		mg_interpolate(level,i);		
		mg_conjugate_gradient(level,i-1,2*neu,up,v_level);						
		//mg_gauss_seidel6(level,i-1,2*neu,up,v_level,pbc); 							
	}
	
	mg_residual_neu(level,v_level);
	//cout<<mg_res_norm(level,v_level)<<"\n";  	
}

/******************************************************************************************************/
void impose_zero_mean(vector<mg_grid>& level, int lev)
{
	int ind; 
	double sum=0.0, mean=0.0;
	int sp = pow(2,lev); 
	
	int spacin_x = nx/sp, spacin_y = ny/sp;
	int tot_inp = (spacin_x)*(spacin_y); //Total no. of interior points required for imposing zero mean. Grid points on the far right and top row are removed.   	
	
	for(int j=0;j<=ny-sp;j=j+sp)
	{
	   for(int i=0;i<=nx-sp;i=i+sp)
	   {
	   	ind = i + j*str_x; 
	   	sum = sum + level[lev].phi_s[ind]; 	   
	   }
	}	
		
	mean = sum/tot_inp; 	
				
	for(int j=0;j<=ny-sp;j=j+sp)
	{
	   for(int i=0;i<=nx-sp;i=i+sp)
	   {	
	   	ind = i + j*str_x; 
	   	level[lev].phi_s[ind] = level[lev].phi_s[ind] - mean; 
	   }
	}	
}

/******************************************************************************************************/
/****Evaluating the RHS sum to further get the point correction for applying consistency condition*****/

double mg_eval_rhs_sum(vector<mg_grid>& level, int lev)
{
	int ind;
	double rhs_sum=0.0; 	
	int sp = pow(2,lev); 
	
	int st_inx = 0; 
	int en_inx = nx-sp; 
	int st_iny = 0; 
	int en_iny = ny-sp; 
	
	for(int j=st_iny;j<=en_iny;j=j+sp)
	{
		for(int i=st_inx;i<=en_inx;i=i+sp)
		{
			ind = i + j*str_x; 
			
			rhs_sum = rhs_sum + level[lev].coeff[ind] * level[lev].rhs[ind]; 							
		}	
	}
	
	return rhs_sum; 
}

/********************************************************************/
//Modifying the direct solver coefficients for the matrix inversion in case of a doubly periodic domain
//All the indices are incremented by 1 as the periodic case has all the points for computation

void arma_direct_solve(vector<mg_grid>& level, int lev)
{
	int ind, ind_m, i_m,j_m; 
				
	int sp = pow(2,lev);
	int st_inx=0, st_iny=0;  	        //True indices of the periodic computational domain
	int en_inx = nx_per, en_iny = ny_per;   //True indices of the periodic computational domain
		
	int nx_sol = (nx_per/pow(2,lev)) + 1, ny_sol = (ny_per/pow(2,lev)) + 1; //No of points of the coarsest grid
	int tot_p_sol = (nx_sol)*(ny_sol); //No.of points that are solved for directly. The top and the farthest right line of points are removed due to periodicity
		 
	int str_m = nx_sol; //Stride for the coarsest grid walk (also reduces by 1 for the periodic case)
		
	/*cout<<"nx_sol = "<<nx_sol<<"\n";
	cout<<"ny_sol = "<<ny_sol<<"\n";
	
	cout<<"hxx[lev] at level"<<lev<<"is"<<hxx[lev]<<"\n"; 
	cout<<"hyy[lev] at level"<<lev<<"is"<<hyy[lev]<<"\n";*/
	
	mat A = zeros<mat>(tot_p_sol,tot_p_sol);           //The elliptical operator matrix in full
	mat B = zeros<mat>(tot_p_sol,1);	
		
	mat trim_A = zeros<mat>(tot_p_sol-1,tot_p_sol-1);  //Trimmed matrix which is non-singular 
	mat trim_B = zeros<mat>(tot_p_sol-1,1);	
		
	/*cout<<"In the direct solver\n";*/
		
	for(int j=0;j<=ny_per;j=j+sp)
	{
		for(int i=0;i<=nx_per;i=i+sp)
		{
			ind = i + j*str_x; 	
			
			i_m = (i/sp); 
			j_m = (j/sp); 		

			ind_m = i_m + j_m*str_m; 
			
			//cout<<"i_m= "<<i_m<<"j_m= "<<j_m<<"\n"; 
			//cout<<"Act-ind-m:   "<<ind_m<<"\n"; 
			
			/*The various cases would start now*/
			
			/*******************Case-1*********************************/
			
			/*******************Case-1 j=st_iny***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx)
			  { 			   			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m) = dc_b1_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_b1_12[lev]; 			    			    			
			       
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b2_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m+str_m) = dc_b2_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_b2_12[lev];
			      
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b3_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m+2*str_m) = dc_b3_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_b3_12[lev];  
			      	    
			    A(ind_m,i_m+(ny_sol-2)*str_m) = dc_ny1_11[lev]; 
			    A(ind_m,i_m+1+(ny_sol-2)*str_m) = dc_ny1_12[lev]; 
			    A(ind_m,i_m+2+(ny_sol-2)*str_m) = dc_ny1_13[lev];
			    A(ind_m,nx_sol-2+(ny_sol-2)*str_m) = dc_ny1_13[lev]; 
			    A(ind_m,nx_sol-1+(ny_sol-2)*str_m) = dc_ny1_12[lev];
			    
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_ny_11[lev]; 
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_ny_12[lev]; 
			    A(ind_m,i_m+2+(ny_sol-1)*str_m) = dc_ny_13[lev];
			    A(ind_m,nx_sol-2+(ny_sol-1)*str_m) = dc_ny_13[lev]; 
			    A(ind_m,nx_sol-1+(ny_sol-1)*str_m) = dc_ny_12[lev];
			  }
			  else if(i==st_inx+sp)
			  {     			    			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_nb1_24[lev];
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb2_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_nb2_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb3_24[lev]; 			    
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_nb3_24[lev];	    			      			    
				
			    A(ind_m,i_m-1+(ny_sol-2)*str_m) = dc_ny1_21[lev]; 
			    A(ind_m,i_m+(ny_sol-2)*str_m) = dc_ny1_22[lev]; 
			    A(ind_m,i_m+1+(ny_sol-2)*str_m) = dc_ny1_23[lev]; 
			    A(ind_m,i_m+2+(ny_sol-2)*str_m) = dc_ny1_24[lev]; 			    
			    A(ind_m,nx_sol-1+(ny_sol-2)*str_m) = dc_ny1_24[lev];			   			    
			    
    			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_ny_21[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_ny_22[lev]; 
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_ny_23[lev]; 
			    A(ind_m,i_m+2+(ny_sol-1)*str_m) = dc_ny_24[lev]; 			    
			    A(ind_m,nx_sol-1+(ny_sol-1)*str_m) = dc_ny_24[lev];
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-3*sp)
			  {			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1+str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i3_36[lev]; 
			    			    
			    A(ind_m,i_m-2+(ny_sol-2)*str_m) = dc_ny1_32[lev]; 
			    A(ind_m,i_m-1+(ny_sol-2)*str_m) = dc_ny1_33[lev]; 
			    A(ind_m,i_m+(ny_sol-2)*str_m) = dc_ny1_34[lev];
			    A(ind_m,i_m+1+(ny_sol-2)*str_m) = dc_ny1_35[lev]; 
			    A(ind_m,i_m+2+(ny_sol-2)*str_m) = dc_ny1_36[lev];
			    
    			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_ny_32[lev]; 
			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_ny_33[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_ny_34[lev];
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_ny_35[lev]; 
			    A(ind_m,i_m+2+(ny_sol-1)*str_m) = dc_ny_36[lev];			    			   		    
			  }
			  else if(i==en_inx-2*sp)
			  { 
			    A(ind_m,j_m*str_m) = dc_nb1_24[lev];   			     			    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,j_m*str_m+str_m) = dc_nb2_24[lev];
			    A(ind_m,ind_m-2+str_m) = dc_nb2_24[lev]; 			    
			    A(ind_m,ind_m-1+str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb2_21[lev]; 			    
			    
			    A(ind_m,j_m*str_m+2*str_m) = dc_nb3_24[lev];
			    A(ind_m,ind_m-2+2*str_m) = dc_nb3_24[lev]; 			    			  
			    A(ind_m,ind_m-1+2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb3_21[lev];
			    			    		    
			    A(ind_m,(ny_sol-2)*str_m) = dc_ny1_24[lev];
			    A(ind_m,i_m-2+(ny_sol-2)*str_m) = dc_ny1_24[lev]; 			    			  
			    A(ind_m,i_m-1+(ny_sol-2)*str_m) = dc_ny1_23[lev]; 
			    A(ind_m,i_m+(ny_sol-2)*str_m) = dc_ny1_22[lev]; 
			    A(ind_m,i_m+1+(ny_sol-2)*str_m) = dc_ny1_21[lev];
			    
   			    A(ind_m,(ny_sol-1)*str_m) = dc_ny_24[lev];
			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_ny_24[lev]; 			    			  
			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_ny_23[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_ny_22[lev]; 
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_ny_21[lev]; 		 					    
			  }
			  else 
			  {		    				    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    A(ind_m,j_m*str_m) = dc_b1_12[lev];
			    A(ind_m,1+j_m*str_m) = dc_b1_13[lev]; 			     
			    
			    A(ind_m,ind_m+str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b2_13[lev]; 			    
			    A(ind_m,j_m*str_m+str_m) = dc_b2_12[lev];
			    A(ind_m,1+j_m*str_m+str_m) = dc_b2_13[lev];
			    
			    A(ind_m,ind_m+2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b3_13[lev];			    
			    A(ind_m,j_m*str_m+2*str_m) = dc_b3_12[lev];			    
			    A(ind_m,1+j_m*str_m+2*str_m) = dc_b3_13[lev];
			    			    			    
			    A(ind_m,i_m+(ny_sol-2)*str_m) = dc_ny1_11[lev]; 
			    A(ind_m,i_m-1+(ny_sol-2)*str_m) = dc_ny1_12[lev]; 
			    A(ind_m,i_m-2+(ny_sol-2)*str_m) = dc_ny1_13[lev];			    
			    A(ind_m,(ny_sol-2)*str_m) = dc_ny1_12[lev];			    
			    A(ind_m,1+(ny_sol-2)*str_m) = dc_ny1_13[lev];
			    
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_ny_11[lev]; 
			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_ny_12[lev]; 
			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_ny_13[lev];			    
			    A(ind_m,(ny_sol-1)*str_m) = dc_ny_12[lev];			    
			    A(ind_m,1+(ny_sol-1)*str_m) = dc_ny_13[lev];			    
			  }	  	
			}
			
			/*******************Case-2 j=st_iny+sp***********************/ 
			
			if(j==st_iny+sp)
			{
			  if(i==st_inx)
			  {			    			    			    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b4_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m-str_m) = dc_b4_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_b4_12[lev];
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m) = dc_b5_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_b5_12[lev];			    
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b6_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m+str_m) = dc_b6_13[lev]; 			    
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_b6_12[lev];
			     
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b7_13[lev];			    
			    A(ind_m,nx_sol-2+j_m*str_m+2*str_m) = dc_b7_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_b7_12[lev];
			    
			    A(ind_m,i_m + (ny_sol-1)*str_m) = dc_pb7_11[lev]; 
			    A(ind_m,i_m+1 + (ny_sol-1)*str_m) = dc_pb7_12[lev]; 
			    A(ind_m,i_m+2 + (ny_sol-1)*str_m) = dc_pb7_13[lev];			    
			    A(ind_m,nx_sol-2 + (ny_sol-1)*str_m) = dc_pb7_13[lev]; 
			    A(ind_m,nx_sol-1 + (ny_sol-1)*str_m) = dc_pb7_12[lev];			     			
			  }
			  else if(i==st_inx+sp) 
			  {	    			    		    
			    A(ind_m,ind_m-1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb4_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_nb4_24[lev];			    
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_nb5_24[lev];
			    	    	    
			    A(ind_m,ind_m-1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb6_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb7_24[lev];			    
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_nb7_24[lev];
			    
  		            A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_pb7_21[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_pb7_22[lev]; 
			    A(ind_m,i_m+1+(ny_sol-1)*str_m)= dc_pb7_23[lev]; 
			    A(ind_m,i_m+2+(ny_sol-1)*str_m) = dc_pb7_24[lev];			    
			    A(ind_m,nx_sol-1+(ny_sol-1)*str_m) = dc_pb7_24[lev];
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-3*sp)
			  {			    
			    A(ind_m,ind_m-2-str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]; 
			    
			    A(ind_m,ind_m-2+str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m+str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2+2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m+2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_i7_36[lev]; 
			    
    			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_pb7_32[lev]; 
			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_pb7_33[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m)   = dc_pb7_34[lev]; 
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_pb7_35[lev]; 
			    A(ind_m,i_m+2+(ny_sol-1)*str_m) = dc_pb7_36[lev];
			  }
			  else if(i==en_inx-2*sp)
			  {	    			       
			    A(ind_m,ind_m+1-str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb4_24[lev]; 
			    A(ind_m,j_m*str_m-str_m) = dc_nb4_24[lev];
			    
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    A(ind_m,j_m*str_m) = dc_nb5_24[lev];
			    	    	    
			    A(ind_m,ind_m+1+str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb6_24[lev];
			    A(ind_m,j_m*str_m+str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m+1+2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb7_24[lev];		      
			    A(ind_m,j_m*str_m+2*str_m) = dc_nb7_24[lev];
			     			    
			    A(ind_m,i_m+1+(ny_sol-1)*str_m) = dc_pb7_21[lev]; 
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_pb7_22[lev]; 
			    A(ind_m,i_m-1+(ny_sol-1)*str_m)= dc_pb7_23[lev]; 
			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_pb7_24[lev];		      
			    A(ind_m,(ny_sol-1)*str_m) = dc_pb7_24[lev];
			  }
			  else 
			  {	    
			    A(ind_m,ind_m-str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b4_13[lev]; 
			    A(ind_m,1+j_m*str_m-str_m) = dc_b4_13[lev];
			    A(ind_m,j_m*str_m-str_m) = dc_b4_12[lev];
			    			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    A(ind_m,1+j_m*str_m) = dc_b5_13[lev];
			    A(ind_m,j_m*str_m) = dc_b5_12[lev];
			    
			    A(ind_m,ind_m+str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b6_13[lev]; 
			    A(ind_m,1+j_m*str_m+str_m) = dc_b6_13[lev];
			    A(ind_m,j_m*str_m+str_m) = dc_b6_12[lev];
			    			    
			    A(ind_m,ind_m+2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b7_13[lev];
			    A(ind_m,1+j_m*str_m+2*str_m) = dc_b7_13[lev];
			    A(ind_m,j_m*str_m+2*str_m) = dc_b7_12[lev];
			    
			    A(ind_m,i_m+(ny_sol-1)*str_m) = dc_pb7_11[lev]; 
			    A(ind_m,i_m-1+(ny_sol-1)*str_m) = dc_pb7_12[lev]; 
			    A(ind_m,i_m-2+(ny_sol-1)*str_m) = dc_pb7_13[lev];
			    A(ind_m,1+(ny_sol-1)*str_m) = dc_pb7_13[lev];
			    A(ind_m,(ny_sol-1)*str_m) = dc_pb7_12[lev];
			  }  
			}
			
			/*******************Case-3 j>=st_iny+2*sp && j<=en_iny-2*sp*************/
			
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx)
			  {
			    			    		    			    	    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b8_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m-2*str_m) = dc_b8_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_b8_12[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b9_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m-str_m) = dc_b9_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_b9_12[lev]; 
			     				    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b10_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m) = dc_b10_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_b10_12[lev]; 
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b11_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m+str_m) = dc_b11_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_b11_12[lev];
			    			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_b12_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m+2*str_m) = dc_b12_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_b12_12[lev];
			  }
			  else if(i==st_inx+sp) 
			  {		    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb8_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_nb8_24[lev];
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb9_24[lev];	
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_nb9_24[lev];
			    
			    A(ind_m,ind_m-1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb10_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_nb10_24[lev];
			    
			    A(ind_m,ind_m-1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb11_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m+2+2*str_m) = dc_nb12_24[lev];   		    
			    A(ind_m,nx_sol-1+j_m*str_m+2*str_m) = dc_nb12_24[lev];
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-3*sp)
			  {    			   			    
			    A(ind_m,ind_m-1-2*str_m) = dc_i8_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i8_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i8_35[lev];
			    
			    A(ind_m,ind_m-2-str_m) = dc_i9_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i9_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i9_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i9_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i9_36[lev];
			    
			    A(ind_m,ind_m-2) = dc_i10_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i10_33[lev]; 
			    A(ind_m,ind_m) = dc_i10_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i10_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i10_36[lev];
			    
			    A(ind_m,ind_m-2+str_m) = dc_i11_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i11_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i11_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i11_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i11_36[lev];
			    
			    A(ind_m,ind_m-1+2*str_m) = dc_i12_33[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_i12_34[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_i12_35[lev];			    
			  }
			  else if(i==en_inx-2*sp)
			  {	     	  			    			    
			    A(ind_m,j_m*str_m-2*str_m) = dc_nb8_24[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_nb8_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb8_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_nb8_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb8_24[lev]; 
			    
			    A(ind_m,j_m*str_m-str_m) = dc_nb9_24[lev];
			    A(ind_m,ind_m+1-str_m) = dc_nb9_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb9_22[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_nb9_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb9_24[lev];	
			    
			    A(ind_m,j_m*str_m) = dc_nb10_24[lev];
			    A(ind_m,ind_m+1) = dc_nb10_21[lev]; 
			    A(ind_m,ind_m) = dc_nb10_22[lev]; 
			    A(ind_m,ind_m-1) = dc_nb10_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb10_24[lev];
			    
			    A(ind_m,j_m*str_m+str_m) = dc_nb11_24[lev];
			    A(ind_m,ind_m+1+str_m) = dc_nb11_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb11_22[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_nb11_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb11_24[lev];
			    
			    A(ind_m,j_m*str_m+2*str_m) = dc_nb12_24[lev];
			    A(ind_m,ind_m+1+2*str_m) = dc_nb12_21[lev]; 
			    A(ind_m,ind_m+2*str_m) = dc_nb12_22[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_nb12_23[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_nb12_24[lev];			    
			  }
			  else 
			  {    			    		    
			    A(ind_m,ind_m-2*str_m) = dc_b8_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b8_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b8_13[lev];
			    A(ind_m,1+j_m*str_m-2*str_m) = dc_b8_13[lev];
			    A(ind_m,j_m*str_m-2*str_m) = dc_b8_12[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b9_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b9_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b9_13[lev];
			    A(ind_m,1+j_m*str_m-str_m) = dc_b9_13[lev];
			    A(ind_m,j_m*str_m-str_m) = dc_b9_12[lev];
			    
			    A(ind_m,ind_m) = dc_b10_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b10_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b10_13[lev];
			    A(ind_m,1+j_m*str_m-str_m) = dc_b10_13[lev];
	  	            A(ind_m,j_m*str_m-str_m) = dc_b10_12[lev];
			    
			    A(ind_m,ind_m+str_m) = dc_b11_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b11_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b11_13[lev]; 
			    A(ind_m,1+j_m*str_m+str_m) = dc_b11_13[lev];
			    A(ind_m,j_m*str_m+str_m) = dc_b11_12[lev];
			    
			    A(ind_m,ind_m+2*str_m) = dc_b12_11[lev]; 
			    A(ind_m,ind_m-1+2*str_m) = dc_b12_12[lev]; 
			    A(ind_m,ind_m-2+2*str_m) = dc_b12_13[lev];
			    A(ind_m,1+j_m*str_m+2*str_m) = dc_b12_13[lev];
			    A(ind_m,j_m*str_m+2*str_m) = dc_b12_12[lev];		      
			  }
			}
					
			/************************Case-4 j==en_iny-sp*************************************/
			if(j==en_iny-sp)
			{
			  if(i==st_inx)
			  {	    			    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_b4_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m+str_m) = dc_b4_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_b4_12[lev];
			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b5_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m) = dc_b5_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_b5_12[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b6_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m-str_m) = dc_b6_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_b6_12[lev];
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b7_13[lev];			      
			    A(ind_m,nx_sol-2+j_m*str_m-2*str_m) = dc_b7_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_b7_12[lev];	    			    	    			    			    
			    A(ind_m,i_m+st_iny*str_m) = dc_pb7_11[lev]; 
			    A(ind_m,i_m+1+st_iny*str_m) = dc_pb7_12[lev]; 
			    A(ind_m,i_m+2+st_iny*str_m) = dc_pb7_13[lev];			      
			    A(ind_m,nx_sol-2+st_iny*str_m) = dc_pb7_13[lev];
			    A(ind_m,nx_sol-1+st_iny*str_m) = dc_pb7_12[lev];			    			    
			  }
			  else if(i==st_inx+sp) 
			  {			       			    
			    A(ind_m,ind_m-1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m+1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_nb4_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m+str_m) = dc_nb4_24[lev];
			    
			    A(ind_m,ind_m-1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m+1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb5_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_nb5_24[lev];
			    	    	    
			    A(ind_m,ind_m-1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m+1-str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb6_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb7_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_nb7_24[lev];
			    
			    A(ind_m,i_m-1+st_iny*str_m) = dc_pb7_21[lev]; 
			    A(ind_m,i_m+st_iny*str_m) = dc_pb7_22[lev]; 
			    A(ind_m,i_m+1+st_iny*str_m)= dc_pb7_23[lev]; 
			    A(ind_m,i_m+2+st_iny*str_m) = dc_pb7_24[lev];
			    A(ind_m,nx_sol-1+st_iny*str_m) = dc_pb7_24[lev];    		    
			  }
			  else if(i>st_inx+sp && i<en_inx-2*sp)
			  {  	  			    
			    A(ind_m,ind_m-2+str_m) = dc_i4_32[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_i4_33[lev]; 
			    A(ind_m,ind_m+str_m) = dc_i4_34[lev]; 
			    A(ind_m,ind_m+1+str_m) = dc_i4_35[lev]; 
			    A(ind_m,ind_m+2+str_m) = dc_i4_36[lev];			    	  
			    
			    A(ind_m,ind_m-2) = dc_i5_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i5_33[lev]; 
			    A(ind_m,ind_m)   = dc_i5_34[lev]; 
			    A(ind_m,ind_m+1) = dc_i5_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i5_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i6_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i6_33[lev]; 
			    A(ind_m,ind_m-str_m)   = dc_i6_34[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_i6_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i6_36[lev]; 
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i7_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i7_33[lev]; 
			    A(ind_m,ind_m-2*str_m)   = dc_i7_34[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_i7_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i7_36[lev];
			    
   			    A(ind_m,i_m-2+st_iny*str_m) = dc_pb7_32[lev]; 
			    A(ind_m,i_m-1+st_iny*str_m) = dc_pb7_33[lev]; 
			    A(ind_m,i_m+st_iny*str_m)   = dc_pb7_34[lev]; 
			    A(ind_m,i_m+1+st_iny*str_m) = dc_pb7_35[lev]; 
			    A(ind_m,i_m+2+st_iny*str_m) = dc_pb7_36[lev];			    
			  }
			  else if(i==en_inx-2*sp)
			  {
			    A(ind_m,j_m*str_m+str_m) = dc_nb4_24[lev];     	  	    			    			    
			    A(ind_m,ind_m+1+str_m) = dc_nb4_21[lev]; 
			    A(ind_m,ind_m+str_m) = dc_nb4_22[lev]; 
			    A(ind_m,ind_m-1+str_m)= dc_nb4_23[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_nb4_24[lev]; 
			    
			    A(ind_m,j_m*str_m) = dc_nb5_24[lev];
			    A(ind_m,ind_m+1) = dc_nb5_21[lev]; 
			    A(ind_m,ind_m) = dc_nb5_22[lev]; 
			    A(ind_m,ind_m-1)= dc_nb5_23[lev]; 
			    A(ind_m,ind_m-2) = dc_nb5_24[lev]; 
			    
			    A(ind_m,j_m*str_m-str_m) = dc_nb6_24[lev];	    	    
			    A(ind_m,ind_m+1-str_m) = dc_nb6_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb6_22[lev]; 
			    A(ind_m,ind_m-1-str_m)= dc_nb6_23[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_nb6_24[lev];
			    
			    A(ind_m,j_m*str_m-2*str_m) = dc_nb7_24[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_nb7_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb7_22[lev]; 
			    A(ind_m,ind_m-1-2*str_m)= dc_nb7_23[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_nb7_24[lev];
			    
			    A(ind_m,st_iny*str_m) = dc_pb7_24[lev];
			    A(ind_m,i_m+1+st_iny*str_m) = dc_pb7_21[lev]; 
			    A(ind_m,i_m+st_iny*str_m) = dc_pb7_22[lev]; 
			    A(ind_m,i_m-1+st_iny*str_m) = dc_pb7_23[lev]; 
			    A(ind_m,i_m-2+st_iny*str_m) = dc_pb7_24[lev];			    
			  }
			  else 
			  {    			      	    			    
			    A(ind_m,ind_m+str_m) = dc_b4_11[lev]; 
			    A(ind_m,ind_m-1+str_m) = dc_b4_12[lev]; 
			    A(ind_m,ind_m-2+str_m) = dc_b4_13[lev];
			    A(ind_m,1+j_m*str_m+str_m) = dc_b4_13[lev];
			    A(ind_m,j_m*str_m+str_m) = dc_b4_12[lev];
			     			    
			    A(ind_m,ind_m) = dc_b5_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b5_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b5_13[lev]; 
			    A(ind_m,1+j_m*str_m) = dc_b5_13[lev];
			    A(ind_m,j_m*str_m) = dc_b5_12[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b6_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b6_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b6_13[lev]; 
			    A(ind_m,1+j_m*str_m-str_m) = dc_b6_13[lev];
			    A(ind_m,j_m*str_m-str_m) = dc_b6_12[lev];
			    
			    A(ind_m,ind_m-2*str_m) = dc_b7_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b7_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b7_13[lev];		    
			    A(ind_m,1+j_m*str_m-2*str_m) = dc_b7_13[lev];
			    A(ind_m,j_m*str_m-2*str_m) = dc_b7_12[lev];	
			    
   			    A(ind_m,i_m+st_iny*str_m) = dc_pb7_11[lev]; 
			    A(ind_m,i_m-1+st_iny*str_m) = dc_pb7_12[lev]; 
			    A(ind_m,i_m-2+st_iny*str_m) = dc_pb7_13[lev];		    
			    A(ind_m,1+st_iny*str_m) = dc_pb7_13[lev];
			    A(ind_m,st_iny*str_m) = dc_pb7_12[lev];				    
			  }  
			}
			
			/**************************Case-5(j=en_iny)*****************************/ 
	
			if(j==en_iny)
			{
			  if(i==st_inx)
			  {	    			        			    
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m+1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m+2) = dc_b1_13[lev]; 
			    A(ind_m,nx_sol-2+j_m*str_m) = dc_b1_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_b1_12[lev];
			    
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_b2_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m-str_m) = dc_b2_13[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_b2_12[lev];
			    
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_b3_13[lev];
			    A(ind_m,nx_sol-2+j_m*str_m-2*str_m) = dc_b3_13[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_b3_12[lev];
			    
			    A(ind_m,i_m+(st_iny+1)*str_m) = dc_ny1_11[lev]; 
			    A(ind_m,i_m+1+(st_iny+1)*str_m) = dc_ny1_12[lev]; 
			    A(ind_m,i_m+2+(st_iny+1)*str_m) = dc_ny1_13[lev];
			    A(ind_m,nx_sol-2+(st_iny+1)*str_m) = dc_ny1_13[lev];
			    A(ind_m,nx_sol-1+(st_iny+1)*str_m) = dc_ny1_12[lev];
			    
    			    A(ind_m,i_m+(st_iny)*str_m) = dc_ny_11[lev]; 
			    A(ind_m,i_m+1+(st_iny)*str_m) = dc_ny_12[lev]; 
			    A(ind_m,i_m+2+(st_iny)*str_m) = dc_ny_13[lev];
			    A(ind_m,nx_sol-2+(st_iny)*str_m) = dc_ny_13[lev];
			    A(ind_m,nx_sol-1+(st_iny)*str_m) = dc_ny_12[lev];			    			    
			  }
			  else if(i==st_inx+sp)
			  {	    			     	  			    
			    A(ind_m,ind_m-1) = dc_nb1_21[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m+2) = dc_nb1_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m) = dc_nb1_24[lev];
			    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_21[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_nb2_24[lev]; 
			    A(ind_m,nx_sol-1+j_m*str_m-str_m) = dc_nb2_24[lev]; 
			    
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_21[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_nb3_24[lev];
			    A(ind_m,nx_sol-1+j_m*str_m-2*str_m) = dc_nb3_24[lev];			    
			    			    			    
			    A(ind_m,i_m-1+(st_iny+1)*str_m) = dc_ny1_21[lev]; 
			    A(ind_m,i_m+(st_iny+1)*str_m) = dc_ny1_22[lev]; 
			    A(ind_m,i_m+1+(st_iny+1)*str_m) = dc_ny1_23[lev]; 
			    A(ind_m,i_m+2+(st_iny+1)*str_m) = dc_ny1_24[lev];
			    A(ind_m,nx_sol-1+(st_iny+1)*str_m) = dc_ny1_24[lev];
			    
			    A(ind_m,i_m-1+(st_iny)*str_m) = dc_ny_21[lev]; 
			    A(ind_m,i_m+(st_iny)*str_m) = dc_ny_22[lev]; 
			    A(ind_m,i_m+1+(st_iny)*str_m) = dc_ny_23[lev]; 
			    A(ind_m,i_m+2+(st_iny)*str_m) = dc_ny_24[lev];
			    A(ind_m,nx_sol-1+(st_iny)*str_m) = dc_ny_24[lev];
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-3*sp)
			  {	    			    			    			    
			    A(ind_m,ind_m-2) = dc_i1_32[lev]; 
			    A(ind_m,ind_m-1) = dc_i1_33[lev]; 
			    A(ind_m,ind_m) = dc_i1_34[lev];
			    A(ind_m,ind_m+1) = dc_i1_35[lev]; 
			    A(ind_m,ind_m+2) = dc_i1_36[lev]; 
			    
			    A(ind_m,ind_m-2-str_m) = dc_i2_32[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_i2_33[lev]; 
			    A(ind_m,ind_m-str_m) = dc_i2_34[lev];
			    A(ind_m,ind_m+1-str_m) = dc_i2_35[lev]; 
			    A(ind_m,ind_m+2-str_m) = dc_i2_36[lev];
			    
			    A(ind_m,ind_m-2-2*str_m) = dc_i3_32[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_i3_33[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_i3_34[lev];
			    A(ind_m,ind_m+1-2*str_m) = dc_i3_35[lev]; 
			    A(ind_m,ind_m+2-2*str_m) = dc_i3_36[lev];			    			    
			    
			    A(ind_m,i_m-2+(st_iny+1)*str_m) = dc_ny1_32[lev]; 
			    A(ind_m,i_m-1+(st_iny+1)*str_m) = dc_ny1_33[lev]; 
			    A(ind_m,i_m+(st_iny+1)*str_m) = dc_ny1_34[lev];
			    A(ind_m,i_m+1+(st_iny+1)*str_m) = dc_ny1_35[lev]; 
			    A(ind_m,i_m+2+(st_iny+1)*str_m) = dc_ny1_36[lev];
			    
			    A(ind_m,i_m-2+(st_iny)*str_m) = dc_ny_32[lev]; 
			    A(ind_m,i_m-1+(st_iny)*str_m) = dc_ny_33[lev]; 
			    A(ind_m,i_m+(st_iny)*str_m) = dc_ny_34[lev];
			    A(ind_m,i_m+1+(st_iny)*str_m) = dc_ny_35[lev]; 
			    A(ind_m,i_m+2+(st_iny)*str_m) = dc_ny_36[lev];			    
			  }
			  else if(i==en_inx-2*sp)
			  { 
			    A(ind_m,j_m*str_m) = dc_nb1_24[lev];			     			    		    
			    A(ind_m,ind_m-2) = dc_nb1_24[lev];
    			    A(ind_m,ind_m-1) = dc_nb1_23[lev]; 
			    A(ind_m,ind_m) = dc_nb1_22[lev]; 
			    A(ind_m,ind_m+1) = dc_nb1_21[lev]; 
			    
			    A(ind_m,j_m*str_m-str_m) = dc_nb2_24[lev];
			    A(ind_m,ind_m-2-str_m) = dc_nb2_24[lev]; 			    
			    A(ind_m,ind_m-1-str_m) = dc_nb2_23[lev]; 
			    A(ind_m,ind_m-str_m) = dc_nb2_22[lev]; 
			    A(ind_m,ind_m+1-str_m) = dc_nb2_21[lev]; 			    
			    
			    A(ind_m,j_m*str_m-2*str_m) = dc_nb3_24[lev];
			    A(ind_m,ind_m-2-2*str_m) = dc_nb3_24[lev]; 			    			  
			    A(ind_m,ind_m-1-2*str_m) = dc_nb3_23[lev]; 
			    A(ind_m,ind_m-2*str_m) = dc_nb3_22[lev]; 
			    A(ind_m,ind_m+1-2*str_m) = dc_nb3_21[lev];			    			     
			    			    
			    A(ind_m,(st_iny+1)*str_m) = dc_ny1_24[lev];
			    A(ind_m,i_m-2+(st_iny+1)*str_m) = dc_ny1_24[lev]; 			    			  
			    A(ind_m,i_m-1+(st_iny+1)*str_m) = dc_ny1_23[lev]; 
			    A(ind_m,i_m+(st_iny+1)*str_m) = dc_ny1_22[lev]; 
			    A(ind_m,i_m+1+(st_iny+1)*str_m) = dc_ny1_21[lev];
			    
			    A(ind_m,(st_iny)*str_m) = dc_ny_24[lev];
			    A(ind_m,i_m-2+(st_iny)*str_m) = dc_ny_24[lev]; 			    			  
			    A(ind_m,i_m-1+(st_iny)*str_m) = dc_ny_23[lev]; 
			    A(ind_m,i_m+(st_iny)*str_m) = dc_ny_22[lev]; 
			    A(ind_m,i_m+1+(st_iny)*str_m) = dc_ny_21[lev];
			  }
			  else 
			  {   				
			    A(ind_m,j_m*str_m) = dc_b1_12[lev]; 
			    A(ind_m,1+j_m*str_m) = dc_b1_13[lev];      	    	    		    		    			     
			    A(ind_m,ind_m) = dc_b1_11[lev]; 
			    A(ind_m,ind_m-1) = dc_b1_12[lev]; 
			    A(ind_m,ind_m-2) = dc_b1_13[lev]; 
			    
			    A(ind_m,j_m*str_m) = dc_b2_12[lev]; 
			    A(ind_m,1+j_m*str_m) = dc_b2_13[lev]; 
			    A(ind_m,ind_m-str_m) = dc_b2_11[lev]; 
			    A(ind_m,ind_m-1-str_m) = dc_b2_12[lev]; 
			    A(ind_m,ind_m-2-str_m) = dc_b2_13[lev]; 
			    
			    A(ind_m,j_m*str_m) = dc_b3_12[lev];
			    A(ind_m,1+j_m*str_m) = dc_b3_13[lev];
			    A(ind_m,ind_m-2*str_m) = dc_b3_11[lev]; 
			    A(ind_m,ind_m-1-2*str_m) = dc_b3_12[lev]; 
			    A(ind_m,ind_m-2-2*str_m) = dc_b3_13[lev];
			    		    
			    A(ind_m,(st_iny+1)*str_m) = dc_ny1_12[lev];
			    A(ind_m,1+(st_iny+1)*str_m) = dc_ny1_13[lev];
			    A(ind_m,i_m+(st_iny+1)*str_m) = dc_ny1_11[lev]; 
			    A(ind_m,i_m-1+(st_iny+1)*str_m) = dc_ny1_12[lev]; 
			    A(ind_m,i_m-2+(st_iny+1)*str_m) = dc_ny1_13[lev];
			    
   			    A(ind_m,(st_iny)*str_m) = dc_ny_12[lev];
			    A(ind_m,1+(st_iny)*str_m) = dc_ny_13[lev];
			    A(ind_m,i_m+(st_iny)*str_m) = dc_ny_11[lev]; 
			    A(ind_m,i_m-1+(st_iny)*str_m) = dc_ny_12[lev]; 
			    A(ind_m,i_m-2+(st_iny)*str_m) = dc_ny_13[lev]; 			     
			  }	  	
			}
			
			/**********************************************/
			//Populating the full rhs matrix 
			
			B[ind_m] = ( level[lev].rhs[ind] - ( level[lev].point_correc/level[lev].coeff[ind] ) );
		}	
	}
		
	//A.save("Full_matrix",raw_ascii);
	//trim_A.save("trim_A.dat", raw_ascii);
	
	/**********************************************/
	//Reducing the matrix to a non-singular matrix 
	
        trim_A = A.submat(0,0,tot_p_sol-1,tot_p_sol-1); 
	trim_B = B.submat(0,0,tot_p_sol-1,0);
	
	/*cout<<"The determinant of the trimmed matrix is   "<<det(trim_A)<<"\n";*/
	
	mat X = solve(trim_A,trim_B);
	
	/*trim_B.save("trim_B.dat", raw_ascii);
	X.save("Solution_X.dat", raw_ascii); */
	
	for(int j=0;j<=ny_per;j=j+sp)
	{
		for(int i=0;i<=nx_per;i=i+sp)
		{
			i_m = (i/sp); 
			j_m = (j/sp); 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 
						
			//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
					
			if(ind_m<=tot_p_sol-2)
			{
				level[lev].phi_s[ind] = X(ind_m);	
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 			
			}
		}
	}
}

/*********************************************************************************************/
//Conjugate gradient solver for Periodic Neumann-Poisson problem
void mg_conjugate_gradient(vector<mg_grid>& level, int lev, int sm_ite, int up, int v_level)
{
	int sp = pow(2,lev); 
	
	int nx_sol = (nx_per/pow(2,lev)), ny_sol = (ny_per/pow(2,lev)); //No of spacings on the grid to be solved for
	int tot_p_sol = (nx_sol+1)*(ny_sol+1); //No.of points that are solved for directly. 
					       //The top most and the farthest right lines of points are removed due to periodicity
	int ind; 				   
	int str_m = nx_sol+1; 
	int i_m, j_m, ind_m; 	
	
	int st_inx = 0, en_inx = nx_sol; 
	int st_iny = 0, en_iny = ny_sol; 					   

	vec ri(tot_p_sol-1), rip(tot_p_sol-1), pi(tot_p_sol-1), pip(tot_p_sol-1); 
	
	vec ritil(tot_p_sol-1), riptil(tot_p_sol-1), ub(tot_p_sol), ux(tot_p_sol), b(tot_p_sol-1); 
	
	vec x=zeros<vec>(tot_p_sol-1), Api(tot_p_sol-1);
	
	vec error(tot_p_sol-1);  	

	double alphai,betai; 
	double tol=1.0e-10, rhs_sum=0; 	
	int count=0; 
	double res_norm=100.0,max_norm=100.0; 
	
	//ofstream res_out; 
	//res_out.open("nonsingular_cg_residual.dat");	
	
	/********Point correction for solvability**********/ 
	if(lev==v_level) //As both part of the ascending and descending parts of the cycle 
  	{     
      		mg_compute_rhs_tot(level,lev); 
      		rhs_sum = mg_eval_rhs_sum(level,lev);             
  	}  
  
	if( (lev!=v_level) && (up==0))
  	{     	
      		for(int j=st_iny;j<=en_iny;j=j+sp)
      		{
		      	for(int i=st_inx;i<=en_inx;i=i+sp)
      			{
		      		ind = i + j*str_x;
      		
      				level[lev].phi_s[ind] = 0.0;       		      	 	      	
      			}      
      		}
   	
      		rhs_sum = mg_eval_rhs_sum(level,lev);               
  	}
  	
   
  	level[lev].point_correc = rhs_sum/(tot_p_sol); 		
  	
  	/********************************************************************/ 
  	
		
	for(int j=0;j<=ny_sol;j=j+sp)
	{
		for(int i=0;i<=nx_sol;i=i+sp)
		{
			ind = i + j*str_x; 							
			
			i_m = (i/sp);
			j_m = (j/sp); 
			
			ind_m = i_m + j_m*str_m; 
			
			ub(ind_m) = level[lev].rhs[ind] - (level[lev].point_correc/level[lev].coeff[ind]); ;	//Untrimmed RHS	(Includes the last point whose value is to be trimmed) 
			ux(ind_m) = level[lev].phi_s[ind]; 							
		}
	}
		
	/****Initialization of 'b' done*******/ 
	
	b = ub.submat(0,0,tot_p_sol-2,0);   //Trimming the column matrix so that the pinned point is eliminated
	x = ux.submat(0,0,tot_p_sol-2,0); 
	
	ri = b - mulvec(x,lev); 
	pi = ri; 
			
	while(count<sm_ite)
	//while(res_norm>tol)
	{		
		Api = mulvec(pi,lev); 		
		alphai = dot(ri,ri)/dot(pi,Api);
								
		x = x + alphai*pi; 
		
		rip = ri - alphai*Api; 		
		
		betai  = dot(rip,rip)/dot(ri,ri);		
		pip = rip + betai*pi; 	
		
		pi = pip;
		ri = rip; 	
						
	        res_norm = norm(ri,2); 		
		max_norm = norm(ri,"inf"); 			

		//cout<<count<<"	"<<res_norm<<"  "<<norm(ri,"-inf")<<"  "<<norm(ri,"inf")<<"\n";
		//res_out<<count<<"  "<<res_norm<<"  "<<norm(ri,"-inf")<<"  "<<norm(ri,"inf")<<"\n";

		count++; 
	}
		
	//cout<<"No. of iterations is "<<count<<"\n"; 

	for(int j=0;j<=ny_per;j=j+sp)
	{
		for(int i=0;i<=nx_per;i=i+sp)
		{
			i_m = (i/sp); 
			j_m = (j/sp); 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 					
					
			if(ind_m<=tot_p_sol-2)
			{
				level[lev].phi_s[ind] = x(ind_m);
				level[lev].res[ind] = ri(ind_m);			
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 
				level[lev].res[ind] = 0.0; 						
			}			
		}
	}

	//res_out.close();	 
}

/*************************************************************/ 
//Functions linking up Momentum and Pressure Poisson's equation 

void poisson_source(vector<fval>& fvar, vector<mg_grid>& level, double Q)
{
	int ind;

        per_tdma1x(fvar,0); //du/dx
        per_tdma1y(fvar,0); //du/dy 

        per_tdma1x(fvar,1); //dv/dx 
        per_tdma1y(fvar,1); //dv/dy 

        per_tdma2x(fvar,0);//d2u/dx2
        per_tdma2y(fvar,0);//d2u/dy2 

        per_tdma2x(fvar,1);//d2v/dx2
        per_tdma2y(fvar,1);//d2v/dy2 

        for(int j=0;j<ny;j++)
        {
                for(int i=0;i<nx;i++)
                {
                        ind = i + j*str_x;

                        fvar[ind].u[3] = -( fvar[ind].u[0]*fvar[ind].ux[0] + fvar[ind].u[1]*fvar[ind].uy[0] ) +  (1./Re)*( fvar[ind].uxx[0] + fvar[ind].uyy[0] );

                        fvar[ind].u[4] = -( fvar[ind].u[0]*fvar[ind].ux[1] + fvar[ind].u[1]*fvar[ind].uy[1] ) +  (1./Re)*( fvar[ind].uxx[1] + fvar[ind].uyy[1] );
                }
        }

        per_tdma1x(fvar,3);
        per_tdma1y(fvar,4);

        for(int j=0;j<ny;j++)
        {
                for(int i=0;i<nx;i++)
                {
                        ind = i + j*str_x;

                        fvar[ind].F = fvar[ind].ux[3] + fvar[ind].uy[4] + (Q/dt)*( fvar[ind].ux[0] + fvar[ind].uy[1] );
                        level[0].F[ind] = fvar[ind].F;                         
                }
        }
}

/**************************************************************/
//Function to compute Pressure BCs based on momentum equation 
void compute_pbc(vector<fval>& fvar, pbcs& pbc)
{
          int ind, ind1, ind2; 
  
	  /*Modifying evaluation of RHS in line with the B.Cs for the driven cavity problem*/
	  
	  tdma2x(fvar,0);  //Evaluate d2u/dx2 at all points in the domain
	  tdma2y(fvar,1);  //Evlauate d2v/dy2 at all points in the domain
	   
	  for(int i=0;i<=nx;i++) //Getting the pressure gradient in the top and bottom walls of the domain
	  {    
	     ind1 = i; 					//Bottom wall
	     ind2 = i + ny*str_x;      			//Top wall
	     
	     pbc.p1y[i] = (1./Re)*fvar[ind1].uyy[1];      //Bottom wall
	     pbc.p2y[i] = (1./Re)*fvar[ind2].uyy[1];      //Top wall
	  }
	  
	  for(int j=0;j<=ny;j++) //Getting the pressure gradient at the left and right faces of the domain 
	  {
	     ind1 = j*str_x; 				//Left wall
	     ind2 = nx + j*str_x;			//Right wall
	     
	     pbc.p1x[j] = (1./Re)*fvar[ind1].uxx[0];     //Left wall
	     pbc.p2x[j] = (1./Re)*fvar[ind2].uxx[0];     //Right wall 
	  } 
}

/*****************************************************************/

void mg_bcs_neu(vector<mg_grid>& level,int lev, pbcs& pbc)
{
	int ind;
	int sp = pow(2,lev); 
	  
	 //Top boundary
 	 for(int i=0;i<=nx;i=i+sp)
	 {
		ind = i + ny*str_x; 								
	
		level[lev].phi_s[ind] = level[lev].phi_s[i]; 	//Copying the periodic boundary condition
	 } 

	  //Right boundary 
	 for(int j=0;j<=ny;j=j+sp)
	 {
		ind = nx + j*str_x; 						
		level[lev].phi_s[ind] = level[lev].phi_s[j*str_x];  		//Copying the periodic boundary condition
	 }  
}

/*****************************************************************/

void mg_clear_levels(vector<mg_grid>& level)
{
	int ind; 
	
	for(int lev=1;lev<mg_levels;lev++)
	{
		for(int j=0;j<=ny;j++)
		{
			for(int i=0;i<=nx;i++)
			{
				ind = i + j*str_x; 
				level[lev].phi_s[ind] = 0.0; 
				level[lev].F[ind] = 0.0; 				
			}		
		}	
	}	
}

/*****************************************************************/
