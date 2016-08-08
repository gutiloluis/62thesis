#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;
const gsl_rng_type *T;
gsl_rng *R;

double k_r = ;
double k_a = ;
double k_b = ;

double kon_ab = ;
double koff_ab = ;

double kb_dab = ;
double ku_dab = ;

double g_r = ;
double g_a = ;
double g_ab = ;
double g_b = ;

int main(){



	return 0;

}

void step(double *t, double *d, double *r, double *a, double *b, double *ab, double *dab){

	double K = k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r + g_a**a + g_b**b + g_ab**ab;

	double Dt = -log(gsl_rng_uniform(R))/K;
  double which = gsl_rng_uniform(R);

	*t+=Dt;

	if(which < k_r**d/K){
		*r+=1;
	} else if(which >= k_r**d/K && which < (k_r**d + k_a**r)/K){
		*a+=1;
	} else if(which >= (k_r**d + k_a**r)/K && which < (k_r**d + k_a**r + k_b**r)/K){
		*b+=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b)/K){
		*a-=1;
		*b-=1;
		*ab+=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab)/K){
		*a+=1;
		*b+=1;
		*ab-=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab)/K){
		*d-=1;
		*ab-=1;
		*dab+=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab)/K){
		*d+=1;
		*ab+=1;
		*dab-=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r)/K){
		*r-=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r + g_a**a)/K){
		*a-=1;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r + g_a**a)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab + g_r**r + g_a**a + g_b**b)/K){
		*b-=1;
	} else {
		*ab-=1;
	}
}
