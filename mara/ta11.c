#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;

const gsl_rng_type *T;
gsl_rng *R;

double ka_u = 1.8;
double ka;
double kr_u = 1.8;
double kr;

double c_inh1 = 800;
double c_inh2 = 10;
double c_act = 80;

double al00 = 0.4;
double al01;
double al10;
double al11;
double al12;
double al02;

double l_m;

double b_a = 0.34*20.0;
double b_r = 0.044*20.0;

double kf_a = 5.0;
double kf_r = 5.0;

double kd_r = 0.01;
double kd_ru;

double l_r;
double l_r2;
double l_au;
double l_a;
double l_ru;

double ks = 20.0;
double ks_u = 0.5;

double al = 1000.0;
double al_p = 1.5;
double b = 1.5;
double b_p = 1.5;

double v = 1.0;

double t_total;

double K, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34;

void step(double *t, double *p00, double *p01, double *p02, double *p10, double *p11, double *p12, double *m, double *au, double *a, double *ru, double *r, double *r2, double *s, double *r2s);
void cell(double t_total);

int main(){

	T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,0);

	*&ka = ka_u/1500.0;
	*&kr = kr_u/150.0;

	*&al01 = al00/c_inh1;
	*&al10 = al00*c_act;
	*&al11 = al00*c_act/c_inh1;
	*&al12 = al00*c_act/(c_inh1*c_inh2);
	*&al02 = al00/(c_inh1*c_inh2);

	*&l_m = log(2.0)/24.0;
	*&l_r = log(2.0)/24.0;
	*&l_a = log(2);

	*&kd_ru = kd_r / 50.0;

	*&l_r2 = l_r;
	*&l_au = l_r;
	*&l_ru = l_r;

	*&t_total = 1000.0/(l_a + l_r);

	cell(t_total);

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

	K = ( ka_u**p10 + ka_u**p11 + ka_u**p12 )  +  ( ka**p00**a/v + ka**p01**a/(v*b) + ka**p02**a/(v*b_p) )  +  ( kr_u**p01 + 2.0*kr_u**p02 + kr_u**p11 + 2.0*kr_u**p12 )  +  ( 2.0*kr**p00**r2/v + kr**p01**r2/v + 2.0*kr**p10**r2/(v*al)  +  kr**p11**r2/(v*al_p) )  +  ( al00**p00 + al01**p01 + al10**p10 + al11**p11 + al12**p12 + al02**p02 )  +  ( l_m**m )  +  ( b_a**m )  +  ( b_r**m )  +  ( kf_a**au )  +  ( kf_r**ru )  +  ( kd_r**r*(*r-1)/(v*2.0) )  +  ( kd_ru**r2 )  +  ( l_au**au )  +  ( l_a**a )  +  ( l_ru**ru )  +  ( l_r**r + l_r2**r2 )  +  ( ks**r2**s/v )  +  ( ks_u**r2s );

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
	if( w < K1/K ){
		*p10 -= 1;
		*p00 += 1;
		*a += 1;
	}
	/* a dissociation from p11 */
	else if( w >= K1/K && w < K2/K ){
		*p11 -= 1;
		*p01 += 1;
		*a += 1;
	}
	/* a dissociation from p12 */
	else if( w >= K2/K && w < K3/K ){
		*p12 -= 1;
		*p02 += 1;
		*a += 1;
	}
	/*a association to p00 */
	else if( w >= K3/K && w < K4/K ){
		*p00 -= 1;
		*p10 += 1;
		*a -= 1;
	}
	/*a association to p01 */
	else if( w >= K4/K && w < K5/K ){
		*p01 -= 1;
		*p11 += 1;
		*a -= 1;
	}
	/*a association to p02 */
	else if( w >= K5/K && w < K6/K ){
		*p02 -= 1;
		*p12 += 1;
		*a -= 1;
	}
	/*r2 dissociation from p01 */
	else if( w >= K6/K && w < K7/K ){
		*p01 -= 1;
		*p00 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p02 */
	else if( w >= K7/K && w < K8/K ){
		*p02 -= 1;
		*p01 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p11 */
	else if( w >= K8/K && w < K9/K ){
		*p11 -= 1;
		*p10 += 1;
		*r2 += 1;
	}
	/*r2 dissociation from p12 */
	else if( w >= K9/K && w < K10/K ){
		*p12 -= 1;
		*p11 += 1;
		*r2 += 1;
	}
	/*r2 association to p00 */
	else if( w >= K10/K && w < K11/K ){
		*p00 -= 1;
		*p01 += 1;
		*r2 -= 1;
	}
	/*r2 association to p01 */
	else if( w >= K11/K && w < K12/K ){
		*p01 -= 1;
		*p02 += 1;
		*r2 -= 1;
	}
	/*r2 association to p10 */
	else if( w >= K12/K && w < K13/K ){
		*p10 -= 1;
		*p11 += 1;
		*r2 -= 1;
	}
	/*r2 association to p11 */
	else if( w >= K13/K && w < K14/K ){
		*p11 -= 1;
		*p12 += 1;
		*r2 -= 1;
	}
	/*transcription with p00 */
	else if( w >= K14/K && w < K15/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p01 */
	else if( w >= K15/K && w < K16/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p10 */
	else if( w >= K16/K && w < K17/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p11 */
	else if( w >= K17/K && w < K18/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p12 */
	else if( w >= K18/K && w < K19/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*transcription with p02 */
	else if( w >= K19/K && w < K20/K ){
		*m+=1;
		*ru+=1;
		*au+=1;
	}
	/*m degradation */
	else if( w >= K20/K && w < K21/K ){
		*m -= 1;
	}
	/*a translation */
	else if( w >= K21/K && w < K22/K ){
		*au += 1;
	}
	/*r translation */
	else if( w >= K22/K && w < K23/K ){
		*ru += 1;
	}
	/*a folding */
	else if( w >= K23/K && w < K24/K ){
		*au -= 1;
		*a += 1;
	}
	/*r folding */
	else if( w >= K24/K && w < K25/K ){
		*ru -= 1;
		*r += 1;
	}
	/*r dimerization */
	else if( w >= K25/K && w < K26/K ){
		*r2 += 1;
		*r -= 2;
	}
	/*r2 dimer disruption */
	else if( w >= K26/K && w < K27/K ){
		*r2 -= 1;
		*r += 2;
	}
	/*au degradation */
	else if( w >= K27/K && w < K28/K ){
		*au -= 1;
	}
	/*a degradation */
	else if( w >= K28/K && w < K29/K ){
		*a -= 1;
	}
	/*ru degradation */
	else if( w >= K29/K && w < K30/K ){
		*ru -= 1;
	}
	/*r degradation */
	else if( w >= K30/K && w < K31/K ){
 		*r -= 1;
	}
	/*r2 degradation */
	else if( w >= K31/K && w < K32/K ){
		*r2 -= 1;
	}
	/*r2 and s association */
	else if( w >= K32/K && w < K33/K ){
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

	/*SI LOS HAGO INT LA EMBARRO, NO ES MEJOR POR EFICIENCIA?*/
	double t=0.0;
	double p00 = 1;
	double p01 = 0;
	double p02 = 0;
	double p10 = 0;
	double p11 = 0;
	double p12 = 0;
	double m = 0;
	double au = 0;
	double a = 0;
	double ru = 0;
	double r = 0;
	double r2 = 0;
	double s = 0;
	double r2s = 0;

	FILE *out = fopen("one_cell.dat", "w");

	fprintf(out, "%f %f %f %f %f\n",t,a,r,r2,r2s);
	while(t<t_total){
		step(&t, &p00, &p01, &p02, &p10, &p11, &p12, &m, &au, &a, &ru, &r, &r2, &s, &r2s);
		fprintf(out, "%f %f %f %f %f\n",t,a,r,r2,r2s);
	}
	fclose(out);
}
