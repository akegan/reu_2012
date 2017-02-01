
/***********************************************************************************
 * The purpose of this code is to solve a pair of non-linear first order ODEs      *
 * derived from a one-dimensional solution of MHD equations in an solar            *
 * prominence using the 4th order Runge-Kutta solver from the GNU GSL Library 1.12 *
 *																				   *
 * The equations in question are:												   *
 * xi=theta^5/2																	   *
 * zeta=1/10*kappa_o*beta^2*[(beta^2-H^2)/(1+G_o^2+H^2) *d(xi)/dH				   *
 *																				   *
 * and the Runge-Kutta algorithm will be solving:								   *
 * d(xi)/dH=10[(1+G_o^2+H^2)/kappa_o*beta^2*(beta-H^2)]*zeta					   *	
 * d(zeta)/dH=(beta-H^2)*xi^(2*(n-1)/5)-gamma*beta								   *
 ***********************************************************************************/



#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_rng.h>
#include<math.h>
#include<gsl/gsl_odeiv.h>
#include<gsl/gsl_matrix.h>
#include<string.h>

/*********************************************
 **  Define all of the initial Conditions   **
 **                                         **
 *********************************************/
#define gamma 0.66667//0.69325 ////0.476787//
/*Set the Constants*/
double G_o=0;
double kappa_o=1.0; //1.2
double beta=5.0; //based on value from paper 1
//float gam=1.5; //will be set later by user, gam because gamma is system variable?
size_t dimensions=3; //# of dimensions for RK code

/*Set the order of the temperature term*/
double n=1.0; //n=1 for now, then move it up to n=2

/*Set the step size*/
int total_steps=1000000; //this will have to vary, h is then calculated below

/*SET INITIAL CONDITIONS*/

#define init_H 0.0
#define fin_H sqrt(beta)

double init_xi=1.0;
double init_zeta=0;

double init_dzeta=1;

/*define the structure that will hold all of the particle information*/
struct variables {
	double H; //"spatial" coordinate/time coordinate
	double xi;//xi
	double zeta;//zeta
	double dxi_dH;// deriv of xi wrt H
	double dzeta_dH;//deriv of zeta wrt H
};

/*Set up file output*/
char output_file[30]="prominence_out_gamma_";//prefix for output file
char FILE_NAME[30];



/******************************************
 **Func: The RK aux. function that specifies all of the physics
 **      Takes in the time, the variables y and the parameters
 **      y[ ]={xi, zeta}
 **      derivs[ ]={dxi_dH,dzeta_dH}
 ******************************************/
int func (double t, const double y[], double derivs[], void *params)
{
	double* constants= *(double **)params;//{G_o^2, kappa_o, beta,n, gamma}
	//printf("Constants: G-0 %lf,kappa %lf, beta %lf, n %lf, gamma %lf\n",constants[0],constants[1], constants[2], constants[3], constants[4]);//debug
	 //printf("derivs[o] %lf, derivs[1] %lf, y[0] %lf, y[1] %lf\n",derivs[0], derivs[1], y[0], y[1]);
	//printf("gamma passed to the equation: %lf\n", constants[4]);
	//double test_var=constants[4]*constants[2];
	//printf("test_var= %lf, regular multiplication= %lf\n", test_var, 1.6667*constants[2]);
	derivs[0]=10*(1+constants[0]+t*t)*y[1]/(constants[1]*constants[2]*constants[2]*(constants[2]-t*t));
	//derivs[1]=(constants[2]-t*t)*pow(y[0],2*(constants[3]-1)/5)-gamma*constants[2];
	derivs[1]=(constants[2]-t*t)*pow(y[0],2*(n-1)/5)-gamma*constants[2];
	//printf("is the power working correctly? %lf\n", pow(y[0], .4));
	//	printf("time is: %lf\n", t);
	return GSL_SUCCESS;
}



/**************************************
 ****       Main Function          ****
 **************************************/
int main(void){
	double init_dxi=(1-gamma)*beta;
	printf("Do to the processer statements work?\n beta= %lf, fin_H= %lf \n", beta, fin_H);

	/*Set the step size, h based on beta and total_steps*/
	double h=fabs(init_H-fin_H)/ (float) total_steps;//hopefully this works!
	printf("calculted h= %lf\n", h);
	
	/*Collect number of particles, threshold, from user*/
	/*printf("what would you like gamma to be?: ");
	char line [20]; 
	fgets(line, sizeof(line),stdin);
	sscanf(line, "%lf",&gam);*/
	
	//printf("Made it past scanning for gamma\n");
	
	/*Time how long the whole program takes*/
	time_t start, end;
	double dif;
	time(&start);

	/*Set the initial coordinates*/
	struct variables v;//struct where all of the data can be stored
	v.H=init_H; //"spatial" coordinate/time coordinate
	v.xi=init_xi ;//xi
	v.zeta=init_zeta;//zeta
	v.dxi_dH=init_dxi;// deriv of xi wrt H
	v.dzeta_dH=init_dzeta;//deriv of zeta wrt H
	
	
	/*******************************************************************
	 * MAIN LOOP********************************************************
	 *          ********************************************************
	 *******************************************************************/
	
    
	/*Open the file to write initial coordinates to*/
    char gamma_string[9];
    sprintf(gamma_string, "%1.4lf.txt",gamma); //name file with value of gamma
    strcpy(FILE_NAME, output_file);
    strcat(FILE_NAME, gamma_string);
	char* mode="w"; //specify file mode, write/append
	
    FILE*current_file=fopen(FILE_NAME,mode);
    printf("filename: %s\n",FILE_NAME);
    if(current_file==NULL){
		printf("Error opening file %s\n",FILE_NAME);
		return(0); 
	}
	//  fclose(current_file);
    

	
	printf("problem before printing?\n");
	//fprintf(current_file, "Total_steps: %d\n",total_steps/10000);
	fprintf(current_file, "%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH, v.dzeta_dH);
	printf("why are you not printing??\n");
	//fclose(current_file);
    /*********************************
     *** Start RK Stepping         ***
     *********************************/
    
    /*Set the physical constants*/
    
	double G_o2=G_o*G_o;
	printf("G_o2= %lf\n",G_o2);
    double const_array[]={G_o2, kappa_o, beta,n};
    double* constants=const_array;
	printf("constants: %lf\n",constants[3]);
    
    /*Setup the RK solver*/
    void *jac=NULL;//not using jacobian
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rk4;//set the solver to RK-4
    gsl_odeiv_system sys = {func, jac, dimensions, &constants};  
    
    
    double t, step;//t is for seconds, step is for counting steps between seconds
    double t_h;//for the smaller steps
    double input[]={0.0,0.0};//declare the input array
    double y_err[dimensions-1];//output error
    double dydt_in[dimensions-1], dydt_out[dimensions-1];//input and output derivs
	
    /***********************
     **Main RK solver Loop**
     ***********************/
    
	
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, dimensions);//allocate the solver  
	//	GSL_ODEIV_FN_EVAL(&sys, t_h, input, dydt_in); 
	// 	printf("i=%d\n",i);
	t=init_H;//-2*fabs(h);
	//double t_h=0;
	int steps=0;
	int inside_step;
	
	while(steps<total_steps){
		
		
		for(inside_step=0; inside_step<10000;inside_step++){
		input[0]=v.xi;
		input[1]=v.zeta;
		
		
		dydt_in[0]=v.dxi_dH;
		dydt_in[1]=v.dzeta_dH;
		//printf("time here is: %lf\n", t);
		t_h=t;
		int status = gsl_odeiv_step_apply (s, t_h, h, input, y_err, dydt_in,dydt_out, &sys);
		if(status!=GSL_SUCCESS){
		printf("GSL stepping not a success\n");}
		v.xi=input[0];
		v.zeta=input[1];
		
		v.dxi_dH=dydt_out[0];
		v.dzeta_dH=dydt_out[1];
		t=t+h;//next time step
		//printf("t=%lf, h=%lf\n",t, h); //to debug
		v.H=t;
			steps=steps+1;	
			//printf("inside step: %d, steps: %d\n", inside_step, steps);
		}
		/*print to file*/  
		//printf("steps=%d\n", steps);
		fprintf(current_file, "%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	  
		
			printf("%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	 
		
		
		
		
		
		
		
	}
	;//??
	fclose(current_file);
    printf("%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	  
	
	
	
	
	time(&end);
	dif=difftime (end,start);
	//printf("Total Particles: %lu, Runtime: %.8lf\n",original, dif);
	
	
	return(0);
}
