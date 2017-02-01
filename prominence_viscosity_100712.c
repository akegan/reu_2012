
/***********************************************************************************
 * The purpose of this code is to solve a pair of non-linear first order ODEs      *
 * derived from a one-dimensional solution of MHD equations in an solar            *
 * prominence using the 4th order Runge-Kutta solver from the GNU GSL Library 1.12 *
 *	
 * This particular version of the codes deals with the allowance of viscosity in 
 * the simpler isothermal case. The dependent variable, W is the de-dimensionalized
 * velocity of the plasma, a y-dependent vector in the z-direction. W=Lamda_0/eta*v
 * 
 * The equations in question are:												   *
 * xi=dW/DH																		   *
 * zeta=W																		   *
 *																				   *
 * and the Runge-Kutta algorithm will be solving:								   *
 * xi=d(zeta)/dH					   *	
 * dxi/dH=1/(3*zeta^2*mu_0) [(beta-H^2)-6*mu_0*zeta*xi^2-2*zeta					   *
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

/*Set the Constants*/
#define mu_0 1.0

#define beta_orig 5.0; //based on value from paper 1
double beta=5.0;//this is just due to bizarre 
size_t dimensions=3; //# of dimensions for RK code

/*Set the step size*/
int total_steps=1000000; //this will have to vary, h is then calculated below

/*SET INITIAL CONDITIONS*/
#define init_H 0.0
#define fin_H sqrt(beta)
double init_xi=0.0;
double init_zeta=-1.56;

/*define the structure that will hold all of the particle information*/
struct variables {
	double H; //"spatial" coordinate/time coordinate
	double xi;//xi
	double zeta;//zeta
	double dxi_dH;// deriv of xi wrt H
	double dzeta_dH;//deriv of zeta wrt H
};

/*Set up file output*/
char output_file[30]="prominence_visc_mu_";//prefix for output file
char FILE_NAME[30];



/******************************************
 **Func: The RK aux. function that specifies all of the physics
 **      Takes in the time, the variables y and the parameters
 **      y[ ]={u, xi}
 **      derivs[ ]={dzeta_dH,dxi_dH}
 ******************************************/
int func (double t, const double y[], double derivs[], void*params)
{
	double* constants= *(double **)params;//points to void? we'll see how this works
	//printf("PRE: t=%lf, y[0]= %lf, y[1]=%lf, derivs[0]=%lf, derivs[1]=%lf\n",t, y[0], y[1], derivs[0], derivs[1]);
	derivs[0]=y[1];//d(zeta)/dH=xi
	derivs[1]=(0.5*(beta-t*t)+y[0]-mu_0*y[0]*y[1]*y[1])/(mu_0*y[0]*y[0]);
	//derivs[1]=(0.5*(beta-t*t)-mu_0*y[0]*y[1]*y[1]+y[0])/(y[0]*y[0]*mu_0);//dxi/dH=1/(3*zeta^2*mu_0) [(beta-H^2)-6*mu_0*zeta*xi^2-2*zeta
	//printf("xi= %lf, zeta/u= %lf, time= %lf\n", y[1],y[0],t);
	//printf("POS: t=%lf, y[0]= %lf, y[1]=%lf, derivs[0]=%lf, derivs[1]=%lf\n",t, y[0], y[1], derivs[0], derivs[1]);
	//printf("time is: %lf\n", t);
	return GSL_SUCCESS;
}



/**************************************
 ****       Main Function          ****
 **************************************/
int main(void){
	
	/*Set the step size, h based on beta and total_steps*/
	double h=fabs(init_H-fin_H)/ (float) total_steps;//hopefully this works!
	printf("calculted h= %lf\n", h);
	
	/*Time how long the whole program takes*/
	time_t start, end;
	double dif;
	time(&start);

	/*Set the initial coordinates*/
	struct variables v;//struct where all of the data can be stored
	v.H=init_H; //"spatial" coordinate/time coordinate
	v.xi=init_xi ;//xi
	v.zeta=init_zeta;//zeta
	v.dxi_dH=0;// deriv of xi wrt H--doesn't matter
	v.dzeta_dH=0;//deriv of zeta wrt H--doesn't matter
	
	
	/*******************************************************************
	 * MAIN LOOP********************************************************
	 *          ********************************************************
	 *******************************************************************/
	
    
	/*Open the file to write initial coordinates to*/
    char mu_string[9];
    sprintf(mu_string, "%1.4lf.txt",mu_0); //name file with value of mu
    strcpy(FILE_NAME, output_file);
    strcat(FILE_NAME, mu_string);
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
	printf("%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH, v.dzeta_dH);
	printf("why are you not printing??\n");
	//fclose(current_file);
    /*********************************
     *** Start RK Stepping         ***
     *********************************/
    
    /*Set the physical constants*/
    double const_array[]={1,1,1};//for this version, don't need these
	printf("was that the problem?\n");
    double* constants=const_array;
	
    
    /*Setup the RK solver*/
    void *jac=NULL;//not using jacobian
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rk4;//set the solver to RK-4
    gsl_odeiv_system sys = {func, jac, dimensions, &constants};  
    
    
    double t, step;//t is for seconds, step is for counting steps between seconds
    double t_h;//for the smaller steps
    double input[dimensions-1];//declare the input array
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
		input[0]=v.zeta;
		input[1]=v.xi;
		
		
		dydt_in[0]=v.dzeta_dH;
		dydt_in[1]=v.dxi_dH;
		//printf("time here is: %lf\n", t);
		t_h=t;
		int status = gsl_odeiv_step_apply (s, t_h, h, input, y_err, dydt_in,dydt_out, &sys);
		if(status!=GSL_SUCCESS){
		printf("GSL stepping not a success\n");}
		v.zeta=input[0];
		v.xi=input[1];
			//printf("yerr=%lf, %lf\n", y_err[0], y_err[1]);
		v.dzeta_dH=dydt_out[0];
		v.dxi_dH=dydt_out[1];
		t=t+h;//next time step
		//printf("t=%lf, h=%lf\n",t, h); //to debug
		v.H=t;
			steps=steps+1;	
			//printf("inside step: %d, steps: %d\n", inside_step, steps);
		}
		/*print to file*/  
		//printf("steps=%d\n", steps);
		fprintf(current_file, "%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	  
		printf("%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	 //print out if you need/want to
		
		
		
		
		
		
		
	}
	;//??
	fclose(current_file);
    printf("%lf, %lf, %lf, %lf, %lf\n", v.H, v.xi, v.zeta, v.dxi_dH,v.dzeta_dH);	  
	
	
	
	
	time(&end);
	dif=difftime (end,start);
	//printf("Total Particles: %lu, Runtime: %.8lf\n",original, dif);
	
	
	return(0);
}
