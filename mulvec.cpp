#include "headers.h"
#include "init_2.h"
#include "declarations.h"

vec mulvec(vec X, int lev)
{ 	
	int sp=pow(2,lev), ind; 
	
	int nx_sol = (nx_per/pow(2,lev)), ny_sol = (ny_per/pow(2,lev)); //No of spacings on the grid to be solved for
	int tot_p_sol = (nx_sol+1)*(ny_sol+1); //No.of points that are solved for directly. The top and the farthest right lines of points are  	
					   //removed due to periodicity	
					   
	vec ans(tot_p_sol), subans(tot_p_sol-1); 					   
			
	int i_m,j_m,str_m=nx_sol+1,ind_m; 
	int m_inx, m_iny; 
	
	int st_inx = 0;       //Changed for the periodic case
        int en_inx = nx_per;  //Changed for the periodic case    
  
  	int st_iny = 0;       //Changed for the periodic case 
  	int en_iny = ny_per;  //Changed for the periodic case   
  	  				
	m_inx = nx_sol; 
	m_iny = ny_sol; 
  	
	double RHS1, RHS2, RHS3, RHS4, RHS5, RHS6; 
	
	for(int j=0;j<=ny_per;j=j+sp)
	{
		for(int i=0;i<=nx_per;i=i+sp)
		{
			ind = i + j*str_x; 						
			
			i_m = (i/sp); 
			j_m = (j/sp); 
			
			ind_m = i_m + j_m*str_m; 
			
			/*******************Case-1 j=1***********************/ 
			if(j==st_iny)
			{
			  if(i==st_inx) //Checked 
			  {	        	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+1] + dc_b1_13[lev]*X[ind_m+2] + dc_b1_13[lev]*X[m_inx-1 + j_m*str_m] + dc_b1_12[lev]*X[m_inx+j_m*str_m]; 
			    	    
			    RHS2 = dc_b2_11[lev]*X[ind_m+str_m] + dc_b2_12[lev]*X[ind_m+1+str_m] + dc_b2_13[lev]*X[ind_m+2+str_m] + dc_b2_13[lev]*X[m_inx-1+j_m*str_m+str_m] + dc_b2_12[lev]*X[m_inx+j_m*str_m+str_m];	     	    
			   
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*str_m] + dc_b3_12[lev]*X[ind_m+1+2*str_m] + dc_b3_12[lev]*X[m_inx+(j_m+2)*str_m]; 	        
			    
			    RHS4 = dc_ny1_11[lev]*X[i_m+(m_iny-1)*str_m] + dc_ny1_12[lev]*X[i_m+1+(m_iny-1)*str_m] + dc_ny1_13[lev]*X[m_inx+(m_iny-1)*str_m]; 
			    
			    RHS5 = dc_ny_11[lev]*X[i_m + m_iny*str_m] + dc_ny_12[lev]*X[i_m+1+m_iny*str_m] + dc_ny_13[lev]*X[i_m+2+m_iny*str_m] + dc_ny_13[lev]*X[m_inx-1+m_iny*str_m] ;	    

			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
			    
			    //cout<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"; 			    
			    			    	    			    
			  }
			  else if(i==st_inx+sp) //checked 
			  {    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+1] + dc_nb1_24[lev]*X[ind_m+2] + dc_nb1_24[lev]*X[m_inx+j_m*str_m];
			       	    
		   	    RHS2 = dc_nb2_21[lev]*X[ind_m-1+str_m] + dc_nb2_22[lev]*X[ind_m+str_m] + dc_nb2_23[lev]*X[ind_m+1+str_m] + dc_nb2_24[lev]*X[ind_m+2+str_m] + dc_nb2_24[lev]*X[m_inx+(j_m+1)*str_m]; 
		   	    
		   	    RHS3 = dc_nb3_21[lev]*X[ind_m-1+2*str_m] + dc_nb3_22[lev]*X[ind_m+2*str_m] + dc_nb3_23[lev]*X[ind_m+1+2*str_m]; 
		   	       	      	    
		   	    RHS4 = dc_ny1_21[lev]*X[i_m+(m_iny-1)*str_m-1] + dc_ny1_22[lev]*X[i_m+(m_iny-1)*str_m] + dc_ny1_23[lev]*X[i_m+(m_iny-1)*str_m+1]; 
		   	    
		   	    RHS5 = dc_ny_21[lev]*X[i_m+(m_iny)*str_m-1] + dc_ny_22[lev]*X[i_m+(m_iny)*str_m] + dc_ny_23[lev]*X[i_m+(m_iny)*str_m+1] + dc_ny_24[lev]*X[i_m+(m_iny)*str_m + 2] ; //+ dc_ny_24[lev]*X[en_inx-sp+(en_iny)*str_m ]; 
		   	    
		   	    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 			  
		   	    
		   	   // cout<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"; 			    
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1] + dc_i1_36[lev]*X[ind_m+2]; 
			    
			    RHS2 = dc_i2_32[lev]*X[ind_m-2+str_m] + dc_i2_33[lev]*X[ind_m-1+str_m] + dc_i2_34[lev]*X[ind_m+str_m] + dc_i2_35[lev]*X[ind_m+1+str_m] + dc_i2_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS3 = dc_i3_33[lev]*X[ind_m-1+2*str_m] + dc_i3_34[lev]*X[ind_m+2*str_m] + dc_i3_35[lev]*X[ind_m+1+2*str_m]; 
			
			    RHS4 = dc_ny1_33[lev]*X[i_m-1+(m_iny-1)*str_m] + dc_ny1_34[lev]*X[i_m+(m_iny-1)*str_m] + dc_ny1_35[lev]*X[i_m+1+(m_iny-1)*str_m];
			    
			    if(i+2*sp==en_inx)
			    {
			    	RHS5 = dc_ny_32[lev]*X[i_m-2+m_iny*str_m] + dc_ny_33[lev]*X[i_m-1+m_iny*str_m] + dc_ny_34[lev]*X[i_m+m_iny*str_m] + dc_ny_35[lev]*X[i_m+1+m_iny*str_m];
			    	
			    	/*cout<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n";
			    	
			    	cout<<"--------------------------------------";*/
			    }
			    else
			    {
			    	RHS5 = dc_ny_32[lev]*X[i_m-2+m_iny*str_m] + dc_ny_33[lev]*X[i_m-1+m_iny*str_m] + dc_ny_34[lev]*X[i_m+m_iny*str_m] + dc_ny_35[lev]*X[i_m+1+m_iny*str_m] + dc_ny_36[lev]*X[i_m+2+m_iny*str_m];
			    	
			    	/*cout<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"; 	
			    	
			    	cout<<"--------------------------------------";*/		    				    	
			    }
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 			    	     	    
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m+1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-1] + dc_nb1_24[lev]*X[ind_m-2] + dc_nb1_24[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+1+str_m] + dc_nb2_22[lev]*X[ind_m+str_m] + dc_nb2_23[lev]*X[ind_m-1+str_m] + dc_nb2_24[lev]*X[ind_m-2+str_m] + dc_nb2_24[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+1+2*str_m] + dc_nb3_22[lev]*X[ind_m+2*str_m] + dc_nb3_23[lev]*X[ind_m-1+2*str_m];
			    	    	    
		   	    RHS4 = dc_ny1_21[lev]*X[i_m+(m_iny-1)*str_m+1] + dc_ny1_22[lev]*X[i_m+(m_iny-1)*str_m] + dc_ny1_23[lev]*X[i_m+(m_iny-1)*str_m-1]; 
		   			   	    
		   	    RHS5 = dc_ny_22[lev]*X[i_m+(m_iny)*str_m] + dc_ny_23[lev]*X[i_m+(m_iny)*str_m-1] + dc_ny_24[lev]*X[i_m+(m_iny)*str_m-2] + dc_ny_24[lev]*X[st_inx+(m_iny)*str_m];
		   	    
		   	    //cout<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"; 
		   			   	    
		   	    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			   
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m-1] + dc_b1_13[lev]*X[ind_m-2] + dc_b1_13[lev]*X[st_inx+1+j_m*str_m] + dc_b1_12[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m+str_m] + dc_b2_12[lev]*X[ind_m-1+str_m] + dc_b2_13[lev]*X[ind_m-2+str_m] + dc_b2_13[lev]*X[st_inx+1+(j_m+1)*str_m] + dc_b2_12[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m+2*str_m] + dc_b3_12[lev]*X[ind_m-1+2*str_m] + dc_b3_12[lev]*X[st_inx+(j_m+2)*str_m]; 	   	        
			    
			    RHS4 = dc_ny1_11[lev]*X[i_m+(m_iny-1)*str_m] + dc_ny1_12[lev]*X[i_m-1+(m_iny-1)*str_m] + dc_ny1_13[lev]*X[st_inx+(m_iny-1)*str_m]; 
			    
			    RHS5 = dc_ny_12[lev]*X[i_m-1+m_iny*str_m] + dc_ny_13[lev]*X[i_m-2+m_iny*str_m] + dc_ny_13[lev]*X[st_inx+1+m_iny*str_m] + dc_ny_12[lev]*X[st_inx+m_iny*str_m];	    

			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			    //cout<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_ny1_13[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n";
			    		    			      		     
			  }	  	
			}	
			
			/*****************Case-2 (j=2)**********************/ 
	
			if(j==st_iny+sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	    	    	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-str_m] + dc_b4_12[lev]*X[ind_m+1-str_m] + dc_b4_13[lev]*X[ind_m+2-str_m] + dc_b4_13[lev]*X[m_inx-1+(j_m-1)*str_m] + dc_b4_12[lev]*X[m_inx+(j_m-1)*str_m];  
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+1] + dc_b5_13[lev]*X[ind_m+2] + dc_b5_13[lev]*X[m_inx-1+j_m*str_m] + dc_b5_12[lev]*X[m_inx+j_m*str_m]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+str_m] + dc_b6_12[lev]*X[ind_m+1+str_m] + dc_b6_13[lev]*X[ind_m+2+str_m] + dc_b6_13[lev]*X[m_inx-1+(j_m+1)*str_m] + dc_b6_12[lev]*X[m_inx+(j_m+1)*str_m]; 
			     
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*str_m] + dc_b7_12[lev]*X[ind_m+1+2*str_m] + dc_b7_12[lev]*X[m_inx+(j_m+2)*str_m]; 
			    	    	    
			    RHS5 = dc_pb7_11[lev]*X[i_m+m_iny*str_m] + dc_pb7_12[lev]*X[i_m+1+m_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	 
			    
			    //cout<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";			       			    			  
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-1-str_m] + dc_nb4_22[lev]*X[ind_m-str_m] + dc_nb4_23[lev]*X[ind_m+1-str_m] + dc_nb4_24[lev]*X[ind_m+2-str_m] + dc_nb4_24[lev]*X[m_inx+(j_m-1)*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+1] + dc_nb5_24[lev]*X[ind_m+2] + dc_nb5_24[lev]*X[m_inx+j_m*str_m]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-1+str_m] + dc_nb6_22[lev]*X[ind_m+str_m] + dc_nb6_23[lev]*X[ind_m+1+str_m] + dc_nb6_24[lev]*X[ind_m+2+str_m] + dc_nb6_24[lev]*X[m_inx+(j_m+1)*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-1+2*str_m] + dc_nb7_22[lev]*X[ind_m+2*str_m] + dc_nb7_23[lev]*X[ind_m+1+2*str_m];
			    
			    RHS5 = dc_pb7_21[lev]*X[i_m-1+m_iny*str_m] + dc_pb7_22[lev]*X[i_m+m_iny*str_m] + dc_pb7_23[lev]*X[i_m+1+m_iny*str_m] + dc_pb7_24[lev]*X[i_m+2+m_iny*str_m] ;
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
			    
			    //cout<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_24[lev]<<"\n";			    			    			    	   			  	 
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked
			  {	    	    	    
			    RHS1 = dc_i4_32[lev]*X[ind_m-2-str_m] + dc_i4_33[lev]*X[ind_m-1-str_m] + dc_i4_34[lev]*X[ind_m-str_m] + dc_i4_35[lev]*X[ind_m+1-str_m] + dc_i4_36[lev]*X[ind_m+2-str_m]; 
			    
			    RHS2 = dc_i5_32[lev]*X[ind_m-2] + dc_i5_33[lev]*X[ind_m-1] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+1] + dc_i5_36[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_i6_32[lev]*X[ind_m-2+str_m] + dc_i6_33[lev]*X[ind_m-1+str_m] + dc_i6_34[lev]*X[ind_m+str_m] + dc_i6_35[lev]*X[ind_m+1+str_m] + dc_i6_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS4 = dc_i7_33[lev]*X[ind_m-1+2*str_m] + dc_i7_34[lev]*X[ind_m+2*str_m] + dc_i7_35[lev]*X[ind_m+1+2*str_m];
			    
					    
		   	    RHS5 = dc_pb7_33[lev]*X[i_m-1+m_iny*str_m] + dc_pb7_34[lev]*X[i_m+m_iny*str_m] + dc_pb7_35[lev]*X[i_m+1+m_iny*str_m]; 
				    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			    /*cout<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"; 	
			    	
			    cout<<"--------------------------------------";*/ 
			    		    	        	    
			  }
			  else if(i==en_inx-sp) //Checked 
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m+1-str_m] + dc_nb4_22[lev]*X[ind_m-str_m] + dc_nb4_23[lev]*X[ind_m-1-str_m] + dc_nb4_24[lev]*X[ind_m-2-str_m] + dc_nb4_24[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m+1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-1] + dc_nb5_24[lev]*X[ind_m-2] + dc_nb5_24[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m+1+str_m] + dc_nb6_22[lev]*X[ind_m+str_m] + dc_nb6_23[lev]*X[ind_m-1+str_m] + dc_nb6_24[lev]*X[ind_m-2+str_m] + dc_nb6_24[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m+1+2*str_m] + dc_nb7_22[lev]*X[ind_m+2*str_m] + dc_nb7_23[lev]*X[ind_m-1+2*str_m]; 
			    
			    RHS5 = dc_pb7_22[lev]*X[i_m+m_iny*str_m] + dc_pb7_23[lev]*X[i_m-1+m_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			    //cout<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"; 		    	    			    	       
			  }
			  else //checked 
			  {	    	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m-str_m] + dc_b4_12[lev]*X[ind_m-1-str_m] + dc_b4_13[lev]*X[ind_m-2-str_m] + dc_b4_13[lev]*X[st_inx+1+(j_m-1)*str_m] + dc_b4_12[lev]*X[st_inx+(j_m-1)*str_m]; 	    
			    
			    RHS2 =  dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-1] + dc_b5_13[lev]*X[ind_m-2] + dc_b5_13[lev]*X[st_inx+1+j_m*str_m] + dc_b5_12[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m+str_m] + dc_b6_12[lev]*X[ind_m-1+str_m] + dc_b6_13[lev]*X[ind_m-2+str_m] + dc_b6_13[lev]*X[st_inx+1+(j_m+1)*str_m] + dc_b6_12[lev]*X[st_inx+(j_m+1)*str_m];	    
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m+2*str_m] + dc_b7_12[lev]*X[ind_m-1+2*str_m] + dc_b7_12[lev]*X[st_inx+(j_m+2)*str_m]; 
			    
			    RHS5 = dc_pb7_12[lev]*X[i_m-1+m_iny*str_m] + dc_pb7_12[lev]*X[st_inx+m_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			   // cout<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n";		    
			  }  
			}
	
			/*****************Case-3 (j>=3 and j<=ny-3)**********************/ 
	
			if(j>=st_iny+2*sp && j<=en_iny-2*sp)
			{
			  if(i==st_inx) //Checked 
			  {	    	     	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*str_m] + dc_b8_12[lev]*X[ind_m+1-2*str_m] + dc_b8_12[lev]*X[m_inx+(j_m-2)*str_m];  
			    			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-str_m] + dc_b9_12[lev]*X[ind_m+1-str_m] + dc_b9_13[lev]*X[ind_m+2-str_m] + dc_b9_13[lev]*X[m_inx-1+(j_m-1)*str_m] + dc_b9_12[lev]*X[m_inx+(j_m-1)*str_m];
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m+1] + dc_b10_13[lev]*X[ind_m+2] + dc_b10_13[lev]*X[m_inx-1+j_m*str_m] + dc_b10_12[lev]*X[m_inx+j_m*str_m]; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+str_m] + dc_b11_12[lev]*X[ind_m+1+str_m] + dc_b11_13[lev]*X[ind_m+2+str_m] + dc_b11_13[lev]*X[m_inx-1+(j_m+1)*str_m] + dc_b11_12[lev]*X[m_inx+(j_m+1)*str_m]; 	    	    
			    
			    if(j+2*sp==en_iny)
			    {			    
			    	RHS5 = dc_b12_11[lev]*X[ind_m+2*str_m] + dc_b12_12[lev]*X[ind_m+1+2*str_m];		
			    	
			    	//cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"; 			    		   		    	
			    }
			    else
			    {	    
			    	RHS5 = dc_b12_11[lev]*X[ind_m+2*str_m] + dc_b12_12[lev]*X[ind_m+1+2*str_m] + dc_b12_12[lev]*X[m_inx+(j_m+2)*str_m];	   	
			    	//cout<<dc_b8_11[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n";
			    	
			    	//cout<<"--------------------------------------------";
			    }
			     
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	    			    			   	    
			  }
			  else if(i==st_inx+sp)  //Checked 
			  {	    	    	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m-1-2*str_m] + dc_nb8_22[lev]*X[ind_m-2*str_m] + dc_nb8_23[lev]*X[ind_m+1-2*str_m]; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m-1-str_m] + dc_nb9_22[lev]*X[ind_m-str_m] + dc_nb9_23[lev]*X[ind_m+1-str_m] + dc_nb9_24[lev]*X[ind_m+2-str_m] + dc_nb9_24[lev]*X[m_inx+(j_m-1)*str_m];
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m-1] + dc_nb10_22[lev]*X[ind_m] +dc_nb10_23[lev]*X[ind_m+1] + dc_nb10_24[lev]*X[ind_m+2] + dc_nb10_24[lev]*X[m_inx+j_m*str_m]; 
			    			    			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m-1+str_m] + dc_nb11_22[lev]*X[ind_m+str_m] + dc_nb11_23[lev]*X[ind_m+1+str_m] + dc_nb11_24[lev]*X[ind_m+2+str_m] + dc_nb11_24[lev]*X[m_inx+(j_m+1)*str_m]; 
			    
			    RHS5 = dc_nb12_21[lev]*X[ind_m-1+2*str_m] + dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m+1+2*str_m];
			    			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			    
			    
			   /* cout<<dc_nb8_21[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_23[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_23[lev]<<"\n";
			    	
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp)  //Checked 
			  {    	    
			    RHS1 = dc_i8_33[lev]*X[ind_m-1-2*str_m] + dc_i8_34[lev]*X[ind_m-2*str_m] + dc_i8_35[lev]*X[ind_m+1-2*str_m]; 
			    
			    RHS2 = dc_i9_32[lev]*X[ind_m-2-str_m] + dc_i9_33[lev]*X[ind_m-1-str_m] + dc_i9_34[lev]*X[ind_m-str_m] + dc_i9_35[lev]*X[ind_m+1-str_m] + dc_i9_36[lev]*X[ind_m+2-str_m]; 
			    
			    RHS3 = dc_i10_32[lev]*X[ind_m-2] + dc_i10_33[lev]*X[ind_m-1] + dc_i10_34[lev]*X[ind_m] + dc_i10_35[lev]*X[ind_m+1] + dc_i10_36[lev]*X[ind_m+2]; 
			    
			    RHS4 = dc_i11_32[lev]*X[ind_m-2+str_m] + dc_i11_33[lev]*X[ind_m-1+str_m] + dc_i11_34[lev]*X[ind_m+str_m] + dc_i11_35[lev]*X[ind_m+1+str_m] + dc_i11_36[lev]*X[ind_m+2+str_m]; 
			    
			    RHS5 = dc_i12_33[lev]*X[ind_m-1+2*str_m] + dc_i12_34[lev]*X[ind_m+2*str_m] + dc_i12_35[lev]*X[ind_m+1+2*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			    /*cout<<dc_i8_33[lev]<<"\n"<<dc_i8_34[lev]<<"\n"<<dc_i8_35[lev]<<"\n"<<dc_i9_32[lev]<<"\n"<<dc_i9_33[lev]<<"\n"<<dc_i9_34[lev]<<"\n"<<dc_i9_35[lev]<<"\n"<<dc_i9_36[lev]<<"\n"<<dc_i10_32[lev]<<"\n"<<dc_i10_33[lev]<<"\n"<<dc_i10_34[lev]<<"\n"<<dc_i10_35[lev]<<"\n"<<dc_i10_36[lev]<<"\n"<<dc_i11_32[lev]<<"\n"<<dc_i11_33[lev]<<"\n"<<dc_i11_34[lev]<<"\n"<<dc_i11_35[lev]<<"\n"<<dc_i11_36[lev]<<"\n"<<dc_i12_33[lev]<<"\n"<<dc_i12_34[lev]<<"\n"<<dc_i12_35[lev]<<"\n"; 			    		       
			    
			    cout<<"--------------------------------------------";*/
			  }
			  else if(i==en_inx-sp) //Checked
			  {	    	     	    
			    RHS1 = dc_nb8_21[lev]*X[ind_m+1-2*str_m] + dc_nb8_22[lev]*X[ind_m-2*str_m] + dc_nb8_23[lev]*X[ind_m-1-2*str_m] ; 
			    
			    RHS2 = dc_nb9_21[lev]*X[ind_m+1-str_m] + dc_nb9_22[lev]*X[ind_m-str_m] + dc_nb9_23[lev]*X[ind_m-1-str_m] + dc_nb9_24[lev]*X[ind_m-2-str_m] + dc_nb9_24[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS3 = dc_nb10_21[lev]*X[ind_m+1] + dc_nb10_22[lev]*X[ind_m] + dc_nb10_23[lev]*X[ind_m-1] + dc_nb10_24[lev]*X[ind_m-2] + dc_nb10_24[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS4 = dc_nb11_21[lev]*X[ind_m+1+str_m] + dc_nb11_22[lev]*X[ind_m+str_m] + dc_nb11_23[lev]*X[ind_m-1+str_m] + dc_nb11_24[lev]*X[ind_m-2+str_m] + dc_nb11_24[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    if(j+2*sp==en_iny)
			    {			    
			    	RHS5 = dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m-1+2*str_m];			    	
			    }
			    else
			    {
				RHS5 = dc_nb12_21[lev]*X[ind_m+1+2*str_m] + dc_nb12_22[lev]*X[ind_m+2*str_m] + dc_nb12_23[lev]*X[ind_m-1+2*str_m];	    		    
				
				//cout<<dc_nb8_23[lev]<<"\n"<<dc_nb8_22[lev]<<"\n"<<dc_nb8_21[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_24[lev]<<"\n"<<dc_nb9_23[lev]<<"\n"<<dc_nb9_22[lev]<<"\n"<<dc_nb9_21[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_24[lev]<<"\n"<<dc_nb10_23[lev]<<"\n"<<dc_nb10_22[lev]<<"\n"<<dc_nb10_21[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_24[lev]<<"\n"<<dc_nb11_23[lev]<<"\n"<<dc_nb11_22[lev]<<"\n"<<dc_nb11_21[lev]<<"\n"<<dc_nb12_23[lev]<<"\n"<<dc_nb12_22[lev]<<"\n"<<dc_nb12_21[lev]<<"\n"; 	
			    	
			    	//cout<<"----------------------------------------------\n";		    	
			    }
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 			    			    
			  }
			  else //Checked 
			  {    	    	    
			    RHS1 = dc_b8_11[lev]*X[ind_m-2*str_m] + dc_b8_12[lev]*X[ind_m-1-2*str_m] + dc_b8_12[lev]*X[st_inx+(j_m-2)*str_m]; 
			    
			    RHS2 = dc_b9_11[lev]*X[ind_m-str_m] + dc_b9_12[lev]*X[ind_m-1-str_m] + dc_b9_13[lev]*X[ind_m-2-str_m] + dc_b9_13[lev]*X[st_inx+1+(j_m-1)*str_m] + dc_b9_12[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS3 = dc_b10_11[lev]*X[ind_m] + dc_b10_12[lev]*X[ind_m-1] + dc_b10_13[lev]*X[ind_m-2] + dc_b10_13[lev]*X[st_inx+1+j_m*str_m] + dc_b10_12[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS4 = dc_b11_11[lev]*X[ind_m+str_m] + dc_b11_12[lev]*X[ind_m-1+str_m] + dc_b11_13[lev]*X[ind_m-2+str_m] + dc_b11_13[lev]*X[st_inx+1+(j_m+1)*str_m] + dc_b11_12[lev]*X[st_inx+(j_m+1)*str_m]; 	    
			    
			    if(j+2*sp==en_iny)
			    {
			    	RHS5 = dc_b12_12[lev]*X[ind_m-1+2*str_m] + dc_b12_12[lev]*X[st_inx+(j_m+2)*str_m]; 	    	
			    }
			    else
			    {
			    	RHS5 = dc_b12_11[lev]*X[ind_m+2*str_m] + dc_b12_12[lev]*X[ind_m-1+2*str_m] + dc_b12_12[lev]*X[st_inx+(j_m+2)*str_m];
			    	//cout<<dc_b8_12[lev]<<"\n"<<dc_b8_12[lev]<<"\n"<<dc_b8_11[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_13[lev]<<"\n"<<dc_b9_12[lev]<<"\n"<<dc_b9_11[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_13[lev]<<"\n"<<dc_b10_12[lev]<<"\n"<<dc_b10_11[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_13[lev]<<"\n"<<dc_b11_12[lev]<<"\n"<<dc_b11_11[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_12[lev]<<"\n"<<dc_b12_11[lev]<<"\n";
			    	
			    	//cout<<"---------------------------------------\n";
			    }
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			    			  		  
			  }	  	  
			  
			}
	
			/*****************Case-4(j=ny-2)**********************/ 
	
			if(j==en_iny-sp)
			{
			  if(i==st_inx)  //checked 
			  {	     	    
			    RHS1 = dc_b4_11[lev]*X[ind_m+str_m] + dc_b4_12[lev]*X[ind_m+1+str_m] + dc_b4_13[lev]*X[ind_m+2+str_m] + dc_b4_13[lev]*X[m_inx-1+(j_m+1)*str_m] ; 

			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m+1] + dc_b5_13[lev]*X[ind_m+2] + dc_b5_13[lev]*X[m_inx-1+j_m*str_m] + dc_b5_12[lev]*X[m_inx+j_m*str_m]; 

			    RHS3 = dc_b6_11[lev]*X[ind_m-str_m] + dc_b6_12[lev]*X[ind_m+1-str_m] + dc_b6_13[lev]*X[ind_m+2-str_m] + dc_b6_13[lev]*X[m_inx-1+(j_m-1)*str_m] + dc_b6_12[lev]*X[m_inx+(j_m-1)*str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*str_m] + dc_b7_12[lev]*X[ind_m+1-2*str_m] + dc_b7_12[lev]*X[m_inx+(j_m-2)*str_m];
			    
			    RHS5 = dc_pb7_11[lev]*X[i_m+st_iny*str_m] + dc_pb7_12[lev]*X[i_m+1+st_iny*str_m] + dc_pb7_12[lev]*X[m_inx+st_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			     	            
			    			   			    
			    //cout<<dc_pb7_11[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b4_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"; 
			    
			  }
			  else if(i==st_inx+sp) //Checked
			  {	    	    	    
			    RHS1 = dc_nb4_21[lev]*X[ind_m-1+str_m] + dc_nb4_22[lev]*X[ind_m+str_m] + dc_nb4_23[lev]*X[ind_m+1+str_m] + dc_nb4_24[lev]*X[ind_m+2+str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m-1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m+1] + dc_nb5_24[lev]*X[ind_m+2] + dc_nb5_24[lev]*X[m_inx+j_m*str_m]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m-1-str_m] + dc_nb6_22[lev]*X[ind_m-str_m] + dc_nb6_23[lev]*X[ind_m+1-str_m] + dc_nb6_24[lev]*X[ind_m+2-str_m] + dc_nb6_24[lev]*X[m_inx+(j_m-1)*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m-1-2*str_m] + dc_nb7_22[lev]*X[ind_m-2*str_m] + dc_nb7_23[lev]*X[ind_m+1-2*str_m]; 
			    
			    RHS5 = dc_pb7_21[lev]*X[i_m-1+st_iny*str_m] + dc_pb7_22[lev]*X[i_m+st_iny*str_m] + dc_pb7_23[lev]*X[i_m+1+st_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5; 
			    
			    //cout<<dc_pb7_21[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_23[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb4_21[lev]<<"\n"<<dc_nb4_22[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"; 			    			     	    	    
			  }
			  else if(i>st_inx+sp && i<en_inx-sp) //Checked 
			  {	    
			  	if(i+2*sp==en_inx)
			  	{
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2+str_m] + dc_i4_33[lev]*X[ind_m-1+str_m] + dc_i4_34[lev]*X[ind_m+str_m] + dc_i4_35[lev]*X[ind_m+1+str_m]; 
			  		
			  		//cout<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n";		  					   		
			  	}
			  	else
			  	{
			  		RHS1 = dc_i4_32[lev]*X[ind_m-2+str_m] + dc_i4_33[lev]*X[ind_m-1+str_m] + dc_i4_34[lev]*X[ind_m+str_m] + dc_i4_35[lev]*X[ind_m+1+str_m] + dc_i4_36[lev]*X[ind_m+2+str_m]; 	
			  		
			  		/*cout<<dc_pb7_33[lev]<<"\n"<<dc_pb7_34[lev]<<"\n"<<dc_pb7_35[lev]<<"\n"<<dc_i7_33[lev]<<"\n"<<dc_i7_34[lev]<<"\n"<<dc_i7_35[lev]<<"\n"<<dc_i6_32[lev]<<"\n"<<dc_i6_33[lev]<<"\n"<<dc_i6_34[lev]<<"\n"<<dc_i6_35[lev]<<"\n"<<dc_i6_36[lev]<<"\n"<<dc_i5_32[lev]<<"\n"<<dc_i5_33[lev]<<"\n"<<dc_i5_34[lev]<<"\n"<<dc_i5_35[lev]<<"\n"<<dc_i5_36[lev]<<"\n"<<dc_i4_32[lev]<<"\n"<<dc_i4_33[lev]<<"\n"<<dc_i4_34[lev]<<"\n"<<dc_i4_35[lev]<<"\n"<<dc_i4_36[lev]<<"\n";		  					  	 	
			  		
			  		cout<<"ind_m= "<<ind_m<<"\n"; 
			  		
			  		cout<<"----------------------------------------------------\n"; */
			  	}			  	    	    			    
			    
			    RHS2 = dc_i5_32[lev]*X[ind_m-2] + dc_i5_33[lev]*X[ind_m-1] + dc_i5_34[lev]*X[ind_m] + dc_i5_35[lev]*X[ind_m+1] + dc_i5_36[lev]*X[ind_m+2]; 
			    
			    RHS3 = dc_i6_32[lev]*X[ind_m-2-str_m] + dc_i6_33[lev]*X[ind_m-1-str_m] + dc_i6_34[lev]*X[ind_m-str_m] + dc_i6_35[lev]*X[ind_m+1-str_m] + dc_i6_36[lev]*X[ind_m+2-str_m]; 
			    
			    RHS4 = dc_i7_33[lev]*X[ind_m-1-2*str_m] + dc_i7_34[lev]*X[ind_m-2*str_m] + dc_i7_35[lev]*X[ind_m+1-2*str_m];
			    
			    RHS5 = dc_pb7_33[lev]*X[i_m-1+st_iny*str_m] + dc_pb7_34[lev]*X[i_m+st_iny*str_m] + dc_pb7_35[lev]*X[i_m+1+st_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			     	    	 
			  }
			  else if(i==en_inx-sp) //checked 
			  {	    	    	    
			    RHS1 = dc_nb4_22[lev]*X[ind_m+str_m] + dc_nb4_23[lev]*X[ind_m-1+str_m] + dc_nb4_24[lev]*X[ind_m-2+str_m] + dc_nb4_24[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    RHS2 = dc_nb5_21[lev]*X[ind_m+1] + dc_nb5_22[lev]*X[ind_m] + dc_nb5_23[lev]*X[ind_m-1] + dc_nb5_24[lev]*X[ind_m-2] + dc_nb5_24[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS3 = dc_nb6_21[lev]*X[ind_m+1-str_m] + dc_nb6_22[lev]*X[ind_m-str_m] + dc_nb6_23[lev]*X[ind_m-1-str_m] + dc_nb6_24[lev]*X[ind_m-2-str_m] + dc_nb6_24[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS4 = dc_nb7_21[lev]*X[ind_m+1-2*str_m] + dc_nb7_22[lev]*X[ind_m-2*str_m] + dc_nb7_23[lev]*X[ind_m-1-2*str_m];
			    
			    RHS5 = dc_pb7_21[lev]*X[i_m+1+st_iny*str_m] + dc_pb7_22[lev]*X[i_m+st_iny*str_m] + dc_pb7_23[lev]*X[i_m-1+st_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			    //cout<<"ind_m= "<<ind_m<<"\n";
			    
			    //cout<<dc_pb7_23[lev]<<"\n"<<dc_pb7_22[lev]<<"\n"<<dc_pb7_21[lev]<<"\n"<<dc_nb7_23[lev]<<"\n"<<dc_nb7_22[lev]<<"\n"<<dc_nb7_21[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_24[lev]<<"\n"<<dc_nb6_23[lev]<<"\n"<<dc_nb6_22[lev]<<"\n"<<dc_nb6_21[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_24[lev]<<"\n"<<dc_nb5_23[lev]<<"\n"<<dc_nb5_22[lev]<<"\n"<<dc_nb5_21[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_24[lev]<<"\n"<<dc_nb4_23[lev]<<"\n"<<dc_nb4_22[lev]<<"\n";
			    		    	    
			  }
			  else //Checked 
			  {	    	    	    
			    RHS1 = dc_b4_12[lev]*X[ind_m-1+str_m] + dc_b4_13[lev]*X[ind_m-2+str_m] + dc_b4_13[lev]*X[st_inx+1+(j_m+1)*str_m] + dc_b4_12[lev]*X[st_inx+(j_m+1)*str_m]; 
			    
			    RHS2 = dc_b5_11[lev]*X[ind_m] + dc_b5_12[lev]*X[ind_m-1] + dc_b5_13[lev]*X[ind_m-2] + dc_b5_13[lev]*X[st_inx+1+j_m*str_m] + dc_b5_12[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS3 = dc_b6_11[lev]*X[ind_m-str_m] + dc_b6_12[lev]*X[ind_m-1-str_m] + dc_b6_13[lev]*X[ind_m-2-str_m] + dc_b6_13[lev]*X[st_inx+1+(j_m-1)*str_m] + dc_b6_12[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS4 = dc_b7_11[lev]*X[ind_m-2*str_m] + dc_b7_12[lev]*X[ind_m-1-2*str_m] + dc_b7_12[lev]*X[st_inx+(j_m-2)*str_m];  
			    	
			    RHS5 =  dc_pb7_11[lev]*X[i_m+st_iny*str_m] + dc_pb7_12[lev]*X[i_m-1+st_iny*str_m] + dc_pb7_12[lev]*X[st_inx+st_iny*str_m];
			    
			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			    
			   // cout<<dc_pb7_12[lev]<<"\n"<<dc_pb7_12[lev]<<"\n"<<dc_pb7_11[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_12[lev]<<"\n"<<dc_b7_11[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_13[lev]<<"\n"<<dc_b6_12[lev]<<"\n"<<dc_b6_11[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_13[lev]<<"\n"<<dc_b5_12[lev]<<"\n"<<dc_b5_11[lev]<<"\n"<<dc_b4_12[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_13[lev]<<"\n"<<dc_b4_12[lev]<<"\n";			    		    			    	    
			  }  
			}
		
			/****************************************Case-5(j=en_iny)******************************************/ 
		
			if(j==en_iny)
			{
			  if(i==st_inx) 
			  {			  	    	    	    
			    RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m+1] + dc_b1_13[lev]*X[ind_m+2] + dc_b1_13[lev]*X[m_inx-1+j_m*str_m] ; 
			    
			    RHS2 = dc_b2_11[lev]*X[ind_m-str_m] + dc_b2_12[lev]*X[ind_m+1-str_m] + dc_b2_13[lev]*X[ind_m+2-str_m] + dc_b2_13[lev]*X[m_inx-1+(j_m-1)*str_m] + dc_b2_12[lev]*X[m_inx+(j_m-1)*str_m]; 
			    
			    RHS3 = dc_b3_11[lev]*X[ind_m-2*str_m] + dc_b3_12[lev]*X[ind_m+1-2*str_m] + dc_b3_12[lev]*X[m_inx+(j_m-2)*str_m]; 	        
			    
			    RHS4 = dc_ny1_11[lev]*X[i_m+(st_iny+1)*str_m] + dc_ny1_12[lev]*X[i_m+1+(st_iny+1)*str_m] + dc_ny1_12[lev]*X[m_inx+(st_iny+1)*str_m]; 
			    
			    RHS5 =  dc_ny_11[lev]*X[i_m+st_iny*str_m] + dc_ny_12[lev]*X[i_m+1+st_iny*str_m] + dc_ny_13[lev]*X[i_m+2+st_iny*str_m] + dc_ny_13[lev]*X[m_inx-1+st_iny*str_m] + dc_ny_12[lev]*X[m_inx+st_iny*str_m];	    

			    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			    	    		       	    
			    
			    //cout<<dc_ny_11[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_13[lev]<<"\n"<<dc_ny_12[lev]<<"\n"<<dc_ny1_11[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_ny1_12[lev]<<"\n"<<dc_b3_11[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b3_12[lev]<<"\n"<<dc_b2_11[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_13[lev]<<"\n"<<dc_b2_12[lev]<<"\n"<<dc_b1_11[lev]<<"\n"<<dc_b1_12[lev]<<"\n"<<dc_b1_13[lev]<<"\n"<<dc_b1_13[lev]<<"\n"; 
			    
			  }
			  else if(i==st_inx+sp) //Checked 
			  {	    	   	    
			    RHS1 = dc_nb1_21[lev]*X[ind_m-1] + dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m+1] + dc_nb1_24[lev]*X[ind_m+2] ; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m-1-str_m] + dc_nb2_22[lev]*X[ind_m-str_m] + dc_nb2_23[lev]*X[ind_m+1-str_m] + dc_nb2_24[lev]*X[ind_m+2-str_m] + dc_nb2_24[lev]*X[m_inx+(j_m-1)*str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m-1-2*str_m] + dc_nb3_22[lev]*X[ind_m-2*str_m] + dc_nb3_23[lev]*X[ind_m+1-2*str_m] + dc_nb3_24[lev]*X[ind_m+2-2*str_m] + dc_nb3_24[lev]*X[m_inx+(j_m-2)*str_m];
			    	       	    
		   	    RHS4 = dc_ny1_21[lev]*X[i_m+(st_iny+1)*str_m-1] + dc_ny1_22[lev]*X[i_m+(st_iny+1)*str_m] + dc_ny1_23[lev]*X[i_m+(st_iny+1)*str_m+1] ; 
		   	    
		   	    RHS5 = dc_ny_21[lev]*X[i_m+(st_iny)*str_m-1] + dc_ny_22[lev]*X[i_m+(st_iny)*str_m] + dc_ny_23[lev]*X[i_m+(st_iny)*str_m+1] + dc_ny_24[lev]*X[i_m+2+(st_iny)*str_m] + dc_ny_24[lev]*X[m_inx+(st_iny)*str_m] ;    	    
		   	    
		   	    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			     	    	     	    	    
		   	    
		   	    /*cout<<"ind_m = "<<ind_m<<"\n"; 
		   	       	    
		   	    
		   	    cout<<dc_ny_21[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb3_24[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb1_21[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_24[lev]<<"\n";*/
		   	    		   	        
			  }
			  else if(i>=st_inx+2*sp && i<=en_inx-2*sp) //Checked 
			  {	    	    	    
			     if(i+2*sp==en_inx)
			     {
			     	RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1];     
			     	/*cout<<"ind_m= "<<ind_m<<"\n"; 
			     		
				cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n";*/
			     }
			     else
			     {
			     	RHS1 = dc_i1_32[lev]*X[ind_m-2] + dc_i1_33[lev]*X[ind_m-1] + dc_i1_34[lev]*X[ind_m] + dc_i1_35[lev]*X[ind_m+1] + dc_i1_36[lev]*X[ind_m+2];
			     	
			     	/*cout<<"ind_m= "<<ind_m<<"\n"; 
			     		
				cout<<dc_ny_32[lev]<<"\n"<<dc_ny_33[lev]<<"\n"<<dc_ny_34[lev]<<"\n"<<dc_ny_35[lev]<<"\n"<<dc_ny_36[lev]<<"\n"<<dc_ny1_33[lev]<<"\n"<<dc_ny1_34[lev]<<"\n"<<dc_ny1_35[lev]<<"\n"<<dc_i3_33[lev]<<"\n"<<dc_i3_34[lev]<<"\n"<<dc_i3_35[lev]<<"\n"<<dc_i2_32[lev]<<"\n"<<dc_i2_33[lev]<<"\n"<<dc_i2_34[lev]<<"\n"<<dc_i2_35[lev]<<"\n"<<dc_i2_36[lev]<<"\n"<<dc_i1_32[lev]<<"\n"<<dc_i1_33[lev]<<"\n"<<dc_i1_34[lev]<<"\n"<<dc_i1_35[lev]<<"\n"<<dc_i1_36[lev]<<"\n";
				
				cout<<"---------------------------------------------------\n";*/ 
			     }						  			     
			     
			     RHS2 = dc_i2_32[lev]*X[ind_m-2-str_m] + dc_i2_33[lev]*X[ind_m-1-str_m] + dc_i2_34[lev]*X[ind_m-str_m] + dc_i2_35[lev]*X[ind_m+1-str_m] + dc_i2_36[lev]*X[ind_m+2-str_m]; 
			     
			     RHS3 =  dc_i3_33[lev]*X[ind_m-1-2*str_m] + dc_i3_34[lev]*X[ind_m-2*str_m] + dc_i3_35[lev]*X[ind_m+1-2*str_m];   
			    
			     RHS4 = dc_ny1_33[lev]*X[i_m-1+(st_iny+1)*str_m] + dc_ny1_34[lev]*X[i_m+(st_iny+1)*str_m] + dc_ny1_35[lev]*X[i_m+1+(st_iny+1)*str_m];
			     
			     RHS5 = dc_ny_32[lev]*X[i_m-2+st_iny*str_m] + dc_ny_33[lev]*X[i_m-1+st_iny*str_m] + dc_ny_34[lev]*X[i_m+st_iny*str_m] + dc_ny_35[lev]*X[i_m+1+st_iny*str_m] + dc_ny_36[lev]*X[i_m+2+st_iny*str_m];
			     
			     ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;			    
			  }
			  else if(i==en_inx-sp)  //Checked
			  {    	    	    	    
			    RHS1 = dc_nb1_22[lev]*X[ind_m] + dc_nb1_23[lev]*X[ind_m-1] + dc_nb1_24[lev]*X[ind_m-2] + dc_nb1_24[lev]*X[st_inx+j_m*str_m]; 
			    
			    RHS2 = dc_nb2_21[lev]*X[ind_m+1-str_m] + dc_nb2_22[lev]*X[ind_m-str_m] + dc_nb2_23[lev]*X[ind_m-1-str_m] + dc_nb2_24[lev]*X[ind_m-2-str_m] + dc_nb2_24[lev]*X[st_inx+(j_m-1)*str_m]; 
			    
			    RHS3 = dc_nb3_21[lev]*X[ind_m+1-2*str_m] + dc_nb3_22[lev]*X[ind_m-2*str_m] + dc_nb3_23[lev]*X[ind_m-1-2*str_m] ;	    
		   	    
		   	    RHS4 = dc_ny1_21[lev]*X[i_m+(st_iny+1)*str_m+1] + dc_ny1_22[lev]*X[i_m+(st_iny+1)*str_m] + dc_ny1_23[lev]*X[i_m+(st_iny+1)*str_m-1] ; 
		   	    
		   	    RHS5 = dc_ny_21[lev]*X[i_m+st_iny*str_m+1] + dc_ny_22[lev]*X[i_m+st_iny*str_m] + dc_ny_23[lev]*X[i_m+st_iny*str_m-1] + dc_ny_24[lev]*X[i_m+st_iny*str_m-2] + dc_ny_24[lev]*X[st_inx+st_iny*str_m];
		   	    
		   	    ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;
		   	    
		   	    /*cout<<"ind_m= "<<ind_m<<"\n"; 			    	    
		   	    
		   	    cout<<dc_ny_24[lev]<<"\n"<<dc_ny_24[lev]<<"\n"<<dc_ny_23[lev]<<"\n"<<dc_ny_22[lev]<<"\n"<<dc_ny_21[lev]<<"\n"<<dc_ny1_23[lev]<<"\n"<<dc_ny1_22[lev]<<"\n"<<dc_ny1_21[lev]<<"\n"<<dc_nb3_23[lev]<<"\n"<<dc_nb3_22[lev]<<"\n"<<dc_nb3_21[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_24[lev]<<"\n"<<dc_nb2_23[lev]<<"\n"<<dc_nb2_22[lev]<<"\n"<<dc_nb2_21[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_24[lev]<<"\n"<<dc_nb1_23[lev]<<"\n"<<dc_nb1_22[lev]<<"\n"; */
			  }
			  else //Checked 
			  {	  	     	     
			     RHS1 = dc_b1_11[lev]*X[ind_m] + dc_b1_12[lev]*X[ind_m-1] + dc_b1_13[lev]*X[ind_m-2] + dc_b1_13[lev]*X[st_inx+1+j_m*str_m] + dc_b1_12[lev]*X[st_inx+j_m*str_m]; 
			     
			     RHS2 = dc_b2_11[lev]*X[ind_m-str_m] + dc_b2_12[lev]*X[ind_m-1-str_m] + dc_b2_13[lev]*X[ind_m-2-str_m] + dc_b2_13[lev]*X[st_inx+1+(j_m-1)*str_m] + dc_b2_12[lev]*X[st_inx+(j_m-1)*str_m]; 
			     
			     RHS3 = dc_b3_11[lev]*X[ind_m-2*str_m] + dc_b3_12[lev]*X[ind_m-1-2*str_m]  + dc_b3_12[lev]*X[st_inx+(j_m-2)*str_m];    	    
			     
			     RHS4 = dc_ny1_11[lev]*X[i_m+(st_iny+1)*str_m] + dc_ny1_12[lev]*X[i_m-1+(st_iny+1)*str_m] + dc_ny1_13[lev]*X[st_inx+(st_iny+1)*str_m]; 
			    
			     RHS5 = dc_ny_11[lev]*X[i_m+st_iny*str_m] + dc_ny_12[lev]*X[i_m-1+st_iny*str_m] + dc_ny_13[lev]*X[i_m-2+st_iny*str_m] + dc_ny_13[lev]*X[st_inx+1+st_iny*str_m] + dc_ny_12[lev]*X[st_inx+st_iny*str_m];	    

			     ans(ind_m) = RHS1 + RHS2 + RHS3 + RHS4 + RHS5;	
			     //This mulvec part belongs to the pinned point. Hence, it is trimmed out and not involved in computation at all. 
			     			     
			  }	  	
			}											
		
		}	
	}

	subans = ans.submat(0,0,tot_p_sol-2,0); 
	return subans; 
}
