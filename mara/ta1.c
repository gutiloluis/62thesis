#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;

const gsl_rng_type *T;
gsl_rng *R;

double k_r = 1.0;
double k_r_l = 0.1; /*Leakage fraction of the total.*/
/*double k_a = 0.2;*/ /*0.2 default*/
double k_b = 1.0; /*1,0 default*/

double kon_ab = 1.0;
double koff_ab = 1.0;

double kb_dab = 25.0;
double ku_dab = 1.0;

double g_r = 0.1;
double g_a = 0.02;
double g_ab = 0.02;
double g_b = 0.3;

double t_total_1;
double n_cells;

double a_threshold = 1.0;

void step(double *t, double *d, double *r, double *a, double *b, double *ab, double *dab, double k_a);
void cell(double t_total_1, double *a_f, double *b_f, double *ab_f, double *t2_f, double k_a);
void cells(double t_total_1);

int main(){

	T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,2);

	*&t_total_1 = 30.0/g_b;

	cells(t_total_1);
	/*cells(t_total_1);*/

	return 0;

}

void step(double *t, double *d, double *r, double *a, double *b, double *ab, double *dab, double k_a){

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
		k_r*=k_r_l;
	} else if(which >= (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab)/K && which < (k_r**d + k_a**r + k_b**r + kon_ab**a**b + koff_ab**ab + kb_dab**d**ab + ku_dab**dab)/K){
		*d+=1;
		*ab+=1;
		*dab-=1;
		k_r/=k_r_l;
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

void cell(double t_total_1, double *a_f, double *b_f, double *ab_f, double *t2_f, double k_a){

	double t=0;
	double d=1;
	double r=0;
	double a=0;
	double b=0;
	double ab=0;
	double dab=0;

	/*FILE *out = fopen("one_cell.dat", "w");*/


	while(t<t_total_1){
		step(&t, &d, &r, &a, &b, &ab, &dab, k_a);
		/*fprintf(out, "%f %f %f %f %f %f %f\n", t, d, r, a, b, ab, dab);*/
	}
	*a_f = a;
	*b_f = b;
	*ab_f = ab;

	double t2 = 0;
	while(a>a_threshold){
		step(&t2, &d, &r, &a, &b, &ab, &dab, 0.0);
		/*fprintf(out, "%f %f %f %f %f %f %f\n", t+t2, d, r, a, b, ab, dab);*/
	}
	*t2_f = t2;
	/*fclose(out);*/
}

void cells(double t_total_1){
	double a_f;
	double b_f;
	double ab_f;
	double t2_f;

	FILE *outs = fopen("many_cells.dat", "w");

	double k_a=0.2;

	while(k_a<=2.0){
		cell(t_total_1, &a_f, &b_f, &ab_f, &t2_f, k_a);
		fprintf(outs, "%f %f %f %f %f\n", t2_f, a_f, b_f, ab_f, k_a);
		k_a+=0.005;
	}
	fclose(outs);
}
