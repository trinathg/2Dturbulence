extern void make_mesh(void); 
extern void initialize(vector<vertex>&,vector<fval>&);
extern void write_to_file(vector<vertex>&,vector<fval>&);  
extern void write_to_file(vector<vertex>&,vector<fval>&,int); //Function overloading to write files at certain frequency 

extern void tdma1x(vector<fval>&,int);
extern void tdma1y(vector<fval>&,int); 
 
extern void tdma2x(vector<fval>&,int); 
extern void tdma2y(vector<fval>&,int); 

extern void per_tdma1x(vector<fval>&, int);
extern void per_tdma1y(vector<fval>&, int);
extern void per_tdma2x(vector<fval>&, int);
extern void per_tdma2y(vector<fval>&, int); 

extern void filter_gaitonde_8(vector<vertex>&, vector<fval>&);
 
extern void solver(vector<vertex>&,vector<fval>&,vector<mg_grid>&,pbcs&, int);

extern void explicit_solver(vector<vertex>&,vector<fval>&, vector<mg_grid>&, pbcs&);
extern void bcs(vector<vertex>&,vector<fval>&); 
extern void per_bcs(vector<vertex>&,vector<fval>&); 
extern void compute_pbc(vector<fval>&, pbcs&);
extern double div_calc(vector<vertex>&,vector<fval>&);

extern double residual_calc(vector<fval>&); 
extern void global_residuals(vector<vertex>&, vector<fval>&, ofstream& ,int );
extern double error_calc(vector<fval>&);

extern double get_wall_time();

extern void bcs_neu(vector<fval>&); 
extern void poisson_Gauss_Seidel6(vector<fval>&);
extern void poisson_source(vector<fval>& ,vector<mg_grid>&, double);
extern void read_coeff(vector<double>&, vector<double>&, double &);
extern void evaluate_rhs6(vector<fval>&);
extern void evaluate_rhs_total(vector<fval>&, vector<double>&); 
extern void evlauate_F(vector<fval>&); 

extern void write_restart(vector<fval>&,int); 
extern void read_restart(vector<fval>&, int &); 
extern void mg_restart(vector<fval>&, vector<mg_grid>&); 

extern vec mulvec(vec,int);

/****************************************************************************/
/***************************Multigrid functions******************************/ 

extern void mg_coeff(); 
extern void mg_interpolate(vector<mg_grid>&, int);

extern void mg_restrict(vector<mg_grid>&, int);
extern void mg_final(vector<mg_grid>&,vector<fval>&,int); 

extern void mg_read_coeff(vector<mg_grid>&); 
extern void mg_evaluate_rhs(vector<mg_grid>&, int);

extern void mg_residual_neu(vector<mg_grid>& , int);

extern void mg_compute_rhs_tot(vector<mg_grid>&, int); 
extern double mg_eval_rhs_sum(vector<mg_grid>&, int);

extern void mg_gauss_seidel6(vector<mg_grid>&, int, int, int, int,pbcs&); //Overloaded function for FMG
extern double mg_res_norm(vector<mg_grid>&, int);

extern void mg_poisson_solver(vector<mg_grid>&,vector<fval>&, pbcs&);
extern void v_cycle(vector<mg_grid>&, int, pbcs&); //Overloaded function for FMG
extern void fmg(vector<mg_grid>&, pbcs&); //Function for FMG 
extern void mg_bcs_neu(vector<mg_grid>& ,int , pbcs&);
extern void mg_conjugate_gradient(vector<mg_grid>&, int, int, int, int);
extern vec mulvec(vec, int);

extern void mg_bcs_neu(vector<mg_grid>&, int);
extern void mg_initialize(vector<mg_grid>&, vector<vertex>&);

extern void impose_zero_mean(vector<mg_grid>&, int);
extern void arma_direct_solve(vector<mg_grid>&, int);
extern void mg_clear_levels(vector<mg_grid>&);

extern void turbulent_initialize(vector<fval>&, vector<vertex>&);
extern double velocity_component(double,double,int);
extern double random_num();

/****************************************************************************/
