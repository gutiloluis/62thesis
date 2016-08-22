#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;
int n_cells = 10000;
int n_elements = 8;
int t_total;

const gsl_rng_type *T;
gsl_rng *R;

double k_r = 0.01;
double tau_r = 120.0;
double tau_p = 3600.0;
double g_r;
double g_p;
double b = 20.0;
double k_p;

void step(double *t, double *r, double *p);
void cell(double t_total, double *r_f, double *p_f);
void cells(double n_cells, double t_total, char *filename);
void param_B();
void param_K_r();
void param_Tau_p();
int main(int argc, char **argv){

  *&t_total = 10.0*tau_p;

  *&g_r = log(2)/tau_r;
  *&g_p = log(2)/tau_p;
  *&k_p = b*g_r;

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));

  param_B();
  /*param_K_r();
  param_Tau_p();*/

  return 0;
}

void step(double *t, double *r, double *p){

  double K = k_r + g_r**r + k_p**r + g_p**p;

  double Dt = -log(gsl_rng_uniform(R))/K;
  double which = gsl_rng_uniform(R);

  *t+=Dt;

  if(which < k_r/K){
    *r+=1;
  } else if(which >= k_r/K && which < (k_r+g_r**r)/K){
    *r-=1;
  } else if(which >= (k_r+g_r**r)/K && which < (k_r+g_r**r+k_p**r)/K){
    *p+=1;
  } else {
    *p-=1;
  }
}

void cell(double t_total, double *r_f, double *p_f){

  double t = 0;
  double r = 0;
  double p = 0;

  while(t<t_total){
    step(&t, &r, &p);
  }
  *r_f = r;
  *p_f = p;
}

void cells(double n_cells, double t_total, char *filename){
  int i;
  double ave_p = 0;
  double ave_p2 = 0;

  double r_f;
  double p_f;

  FILE *out = fopen(filename, "a");

  for(i=0;i<n_cells;i++){
    cell(t_total, &r_f, &p_f);
    ave_p += p_f;
    ave_p2 += p_f*p_f;
  }

  ave_p/=n_cells;
  ave_p2/=n_cells;

  fprintf(out, "%f %f\n", ave_p, (ave_p2 - ave_p*ave_p)/ave_p);
  fclose(out);
}

void param_B(){

  char filename[10] = "b.dat";
  for(i=0;i<n_elements;i++){
    *&b = 5.0*(i+1);
    *&k_p = b*g_r;
    cells(n_cells, t_total, filename);
  }
  *&b = 20.0;
  *&k_p = b*g_r;
}

void param_K_r(){

  char filename[10] = "kr.dat";
  for(i=0;i<n_elements;i++){
    *&k_r = 0.0025*(i+1);
    cells(n_cells, t_total, filename);
  }
  *&k_r = 0.01;
}

void param_Tau_p(){

  char filename[10] = "taup.dat";
  for(i=0;i<n_elements;i++){
    *&tau_p = 900.0*(i+1);
    *&g_p = log(2)/tau_p;
    cells(n_cells, t_total, filename);
  }
  *&tau_p = 3600.0;
  *&g_p = log(2)/tau_p;
}
