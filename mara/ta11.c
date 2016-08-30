#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;

const gsl_rng_type *T;
gsl_rng *R;


double V = 1.0;

double t_total_1;
double n_cells;



void step(double *t, double *d, double *r, double *a, double *b, double *ab, double *dab);
void cell(double t_total_1);
void cells(double n_cells, double t_total_1);

int main(){

	T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,0);

	*&t_total_1 = 100.0/g_b;

	cell(t_total_1);

	return 0;

}

/*Species
States of the promoter, the first number represents the number of A (MarA) molecules bound and the second the number of MarR_2 molecules bound
p00
p01
p02
p10
p11
p12
m: mRNA
au: marA unfolded
a: marA
ru: marR unfolded
r: MarR
r2: MarR_2
s: salicylate
r2s: MarR_2 bound to salicylate.
*/

/*Reactions
1, (3) a dissociation from p
2. (3) a association to p
3. (4) r2 dissociation from p
4. (4) r2 association to p
5. (6) transcription
6. m degradation
7. a translation
8. r translation
9. a folding
10. r folding
11. r dimerization
12. r2 dimer disruption
13. au degradation
14. a degradation
15. ru degradation
16. r degradation
17. r2 degradation
18. r2 ans s binding
19. r2s unbinding
*/


void step(double *t, double *p00, double *p01, double *p02, double *p10, double *p11, double *p12, double *m, double *au, double *a, double *ru, double *r, double *r2, double *s, double *r2s){

	double K = ( ka_u**p10 + ka_u**p11 + ka_u**p12 )  +  ( ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p) )  +  ( kr_u**p01 + 2.0*kr_u**p02 + kr_u**p11 + 2.0*kr_u**p12 )  +  ( 2.0*kr**p00**r2/v + kr**p01**r2/v + 2.0*kr**p10**r2/(v*al)  +  kr**p11**r2/(v*al_p) )  +  ( al00**p00 + al01**p01 + al10**p10 + al11**p11 + al12**p12 + al02**p02 )  +  ( l_m**m )  +  ()

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

void cell(double t_total){

	double t=0;
	double d=1;
	double r=0;
	double a=0;
	double b=0;
	double ab=0;
	double dab=0;

	FILE *out = fopen("one_one_cell.dat", "w");

	fprintf(out, "%f %f %f %f %f %f %f %f\n",t,d,r,a,b,ab,dab,a+ab);
	while(t<t_total){
		step(&t, &d, &r, &a, &b, &ab, &dab);
		fprintf(out, "%f %f %f %f %f %f %f %f\n",t,d,r,a,b,ab,dab,a+ab);
	}
	fclose(out);
}
