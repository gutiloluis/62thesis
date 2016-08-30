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

	/* a dissociation from p10 */
	K1 = ( ka_u**p10 );
	/* a dissociation from p11 */
	K2 = K1 + ka_u**p11;
	/* a dissociation from p12 */
	K3 = K2 + ka_u**p12;
	/*a association to p00 */
	K4 = K3 + ka**p00**a/v;
	/*a association to p01 */
	K5 = K4 + ka**p01**a/(v*b);
	/*a association to p02 */
	K6 = K5 + ka**p02**a/(v*b_p);
	/*r2 dissociation from p01 */
	K7 = K6 + kr_u**p01;
	/*r2 dissociation from p02 */
	K8 = K7 + 2.0*kr_u**p02;
	/*r2 dissociation from p11 */
	K9 = K8 + kr_u**p11;
	/*r2 dissociation from p12 */
	K10 = K9 + 2.0*kr_u**p12;
	/*r2 association to p00 */
	K11 = K10 + 2.0*kr**p00**r2/v;
	/*r2 association to p01 */
	K12 = K11 + kr**p01**r2/v;
	/*r2 association to p10 */
	K13 = K12 + 2.0*kr**p10**r2/(v*al);
	/*r2 association to p11 */
	K14 = K13 + kr**p11**r2/(v*al_p);
	/*transcription with p00 */
	K15 = K14 + al00**p00;
	/*transcription with p01*/
	K16 = K15 + al01**p01;
	/*transcription with p10 */
	K17 = K16 + al10**p10;
	/*transcription with p11 */
	K18 = K17 + al11**p11;
	/*transcription with p12 */
	K19 = K18 + al12**p12;
	/*transcription with p02 */
	K20 = K19 + al02**p02;
	/*m degradation */
	K21 = K20 + l_m**m;
	/*a translation */
	K22 = K21 + b_a**m;
	/*r translation */
	K23 = K22 + b_r**m;
	/*a folding */
	K24 = K23 + kf_a**au;
	/*r folding */
	K25 = K24 + kf_r**ru;
	/*r dimerization */
	K26 = K25 + kd_r**r*(*r-1)/(v*2.0);
	/*r2 dimer disruption */
	K27 = K26 + kd_ru**r2;
	/*au degradation */
	K28 = K27 + l_au**au;
	/*a degradation */
	K29 = K28 + l_a**a;
	/*ru degradation */
	K30 = K29 + l_ru**ru;
	/*r degradation */
	K31 = K30 + l_r**r;
	/*r2 degradation */
	K32 = K31 + l_r2**r2;
	/*r2 and s association */
	K33 = K32 + ks**r2**s/v;
	/*r2s dissociation*/
	K34 = K33 + ks_u**r2s;

	double Dt = -log(gsl_rng_uniform(R))/K;
  double w = gsl_rng_uniform(R);

	*t+=Dt;

	/* a dissociation from p10 */
	if(  ){
		*p10 -= 1;
		*p00 += 1;
		*a += 1;
	}
	/* a dissociation from p11 */
	else if(  ){
		*p11 -= 1;
		*p01 += 1;
		*a += 1;
	}
	/* a dissociation from p12 */
	else if(  ){
		*p12 -= 1;
		*p02 += 1;
		*a += 1;
	}
	/*a association to p00 */
	else if(  ){
		*p00 -= 1;
		*p10 += 1;
		*a -= 1;
	}
	/*a association to p01 */
	else if(  ){
		*p01 -= 1;
		*p11 += 1;
		*a -= 1;
	}
	/*a association to p02 */
	else if(  ){
		*p02 -= 1;
		*p12 += 1;
		*a -= 1;
	}
	/*r2 dissociation from p01 */
	else if(  ){
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
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p01 */
	else if(  ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p10 */
	else if(  ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p11 */
	else if(  ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p12 */
	else if(  ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p02 */
	else if(  ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*m degradation */
	else if(  ){
		*m -= 1;
	}
	/*a translation */
	else if(  ){
		*au += 1;
	}
	/*r translation */
	else if(  ){
		*ru += 1;
	}
	/*a folding */
	else if(  ){
		*au -= 1;
		*a += 1;
	}
	/*r folding */
	else if(  ){
		*ru -= 1;
		*r += 1;
	}
	/*r dimerization */
	else if(  ){
		*r2 += 1;
		*r -= 2;
	}
	/*r2 dimer disruption */
	else if(  ){
		*r2 -= 1;
		*r += 2;
	}
	/*au degradation */
	else if(  ){
		*au -= 1;
	}
	/*a degradation */
	else if(  ){
		*a -= 1;
	}
	/*ru degradation */
	else if(  ){
		*ru -= 1;
	}
	/*r degradation */
	else if(  ){
 		*r -= 1;
	}
	/*r2 degradation */
	else if(  ){
		*r2 -= 1;
	}
	/*r2 and s association */
	else if(  ){
		*s -= 1;
		*r2 -= 1;
		*r2s += 1;
	}
	/*r2s dissociation*/
	else {
		*s += 1;
		*r2 += 1;
		*r2s += 1;
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
