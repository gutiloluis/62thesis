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
16. (2) r degradation (with r2)
18. r2 ans s binding
19. r2s unbinding
*/


void step(double *t, double *p00, double *p01, double *p02, double *p10, double *p11, double *p12, double *m, double *au, double *a, double *ru, double *r, double *r2, double *s, double *r2s){

	double K = ( ka_u**p10 + ka_u**p11 + ka_u**p12 )  +  ( ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p) )  +  ( kr_u**p01 + 2.0*kr_u**p02 + kr_u**p11 + 2.0*kr_u**p12 )  +  ( 2.0*kr**p00**r2/v + kr**p01**r2/v + 2.0*kr**p10**r2/(v*al)  +  kr**p11**r2/(v*al_p) )  +  ( al00**p00 + al01**p01 + al10**p10 + al11**p11 + al12**p12 + al02**p02 )  +  ( l_m**m )  +  ( b_a**m )  +  ( b_r**m )  +  ( kf_a**au )  +  ( kf_r**ru )  +  ( kd_r**r*(*r-1)/(v*2.0) )  +  ( kd_ru**r2 )  +  ( l_au**au )  +  ( l_a**a )  +  ( l_ru**ru )  +  ( l_r**r + l_r2**r2 )  +  ( ks**r2**s/v )  +  ( ks_u**r2s );

	double Dt = -log(gsl_rng_uniform(R))/K;
  double w = gsl_rng_uniform(R);

	*t+=Dt;

	/* a dissociation from p10 */
	if( w < ka_u**p10/K ){
		*p10 -= 1;
		*p00 += 1;
		*a += 1;
	}
	/* a dissociation from p11 */
	else if( w >= ka_u**p10/K && w < ( ka_u**p10 + ka_u**p11 )/K ){
		*p11 -= 1;
		*p01 += 1;
		*a += 1;
	}
	/* a dissociation from p12 */
	else if( w >= ( ka_u**p10 + ka_u**p11 )/K && w < ( ka_u**p10 + ka_u**p11 + ka_u**p12 )/K ){
		*p12 -= 1;
		*p02 += 1;
		*a += 1;
	}
	/*a association to p00 */
	else if( w >= ( ka_u**p10 + ka_u**p11 + ka_u**p12 )/K && w < ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v) )/K){
		*p00 -= 1;
		*p10 += 1;
		*a -= 1;
	}
	/*a association to p01 */
	else if( w >= ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v) )/K && w < ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v + ka**p01**a/(v*b)) )/K ){
		*p01 -= 1;
		*p11 += 1;
		*a -= 1;
	}
	/*a association to p02 */
	else if( w >= ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v + ka**p01**a/(v*b)) )/K && w < ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p)) )/K ){
		*p02 -= 1;
		*p12 += 1;
		*a -= 1;
	}
	/*r2 dissociation from p01 */
	else if( w >= ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p)) )/K && w < ( (ka_u**p10 + ka_u**p11 + ka_u**p12) + (ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p)) + (kr_u**p01) )/K ){
		*p01 -= 1;
		*p00 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p02 */
	else if(  ){
		*p02 -= 1;
		*p01 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p11 */
	else if(  ){
		*p11 -= 1;
		*p10 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p12 */
	else if(  ){
		*p12 -= 1;
		*p11 += 1;
		*r2 += 1;
	}
	/*r2 association to p00 */
	else if(  ){
		*p00 -= 1;
		*p01 += 1;
		*r2 -= 1;
	}
	/*r2 association to p01 */
	else if(  ){
		*p01 -= 1;
		*p02 += 1;
		*r2 -= 1;
	}
	/*r2 association to p10 */
	else if(  ){
		*p10 -= 1;
		*p11 += 1;
		*r2 -= 1;
	}
	/*r2 association to p11 */
	else if(  ){
		*p11 -= 1;
		*p12 += 1;
		*r2 -= 1;
	}
	/*transcription with p00 */
	else if(  ){

	}
	/*transcription with p00 */
	else if(  ){

	}
	/*transcription with p00 */
	else if(  ){

	}
	/*transcription with p00 */
	else if(  ){

	}
	/*transcription with p00 */
	else if(  ){

	}
	/*transcription with p00 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}
	/*r2 dissociation from p02 */
	else if(  ){

	}



	if(which < k_r**d/K){
		*r+=1;
	} else if(which >= k_r**d/K && which < (k_r**d + k_a**r)/K){
		*a+=1;
	} else if(which >= (k_r**d + k_a**r)/K && which < (k_r**d + k_a**r + k_b**r)/K){
		*b+=1;
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
