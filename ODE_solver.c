#include <stdio.h>
#include <math.h>

#define PATH "euler.dat"
#define sizeof_array(array) sizeof array / sizeof array[0]
typedef float function(int *i, float *x);
int DIMENSION;
float DELTA_t = 0.015625;


float lorenz(int *i, float x[])
{
  //uses splot
  float sigma = 10;
  float r = 28;//reynolds
  float b = 2.667;
  switch(*i)
  {
    case 0: return sigma * (-x[0] + x[1]);
    case 1: return r * x[0] - x[1] -x[0] * x[2];
    case 2: return x[0] * x[1] - b * x[2];
  }
  return -1;
}

float harmonic_osc(int *i, float x[])
{
  float K = 10;
  switch(*i)
  {
    case 0: return x[1];
    case 1: return -K * x[0];
  }
  return -1;
}

float damped_harmonic_osc(int *i, float x[])
{
  float K = 10;
  float P = 0.5;
  switch(*i)
  {
    case 0: return x[1];
    case 1: return -K * x[0] -P * x[1];
  }
  return -1;
}

float poincare_osc(int *i, float x[])
{
  switch(*i)
  {
    case 0: return 5;
    case 1: return x[0] * (1 - x[0]);
  }
  return -1;
}
float vanderPol(int *i, float x[])
{
  float K = 10;
  float Mu = 4;
  switch(*i)
  {
    case 0: return x[1];
    case 1: return -K * x[0] + Mu * (1 - x[0] * x[0]) * x[1];
  }
  return -1;
}
float pendulum(int *i, float x[])
{
  switch(*i)
  {
    case 0: return x[1];
    case 1: return -9.81/0.01 * sin(x[0]);
  }
  return -1;
}

float duffing(int *i, float x[])
{
	float beta = 10;
	switch(*i)
  {
    case 0: return x[1];
    case 1: return x[0] - beta * x[0]*x[0]*x[0];
  }
  return -1;
}

void euler(function f, float t, float t_f, float x[])
{
	FILE *fp = fopen(PATH, "w");
	while(t <= t_f)
  {
    fprintf(fp, "%f", t);
    for (int i = 0; i < DIMENSION; ++i)
      fprintf(fp, " %f", x[i]);
    fprintf(fp, "\n");

    t += DELTA_t;
    for (int i = 0; i < DIMENSION; ++i)
     x[i] += f(&i, x) * DELTA_t;
 }
 fclose(fp);
}

int main()
{
  //initial conditions:
	float x_0[] = {1,1,1};
  //defining dimension:
  DIMENSION = sizeof_array(x_0);
  //Euler method:
  euler(lorenz, 0, 100, x_0);

  //Plotting with gnuplot
  char * commandsForGnuplot[] = {"splot 'euler.dat' u 2:3:4 w lines"};
  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
  for (int i=0; i < sizeof_array(commandsForGnuplot); i++)
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
  return 0;
}