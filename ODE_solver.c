#include <stdio.h>
#include <math.h>

#define PATH "euler.dat"
#define sizeof_array(array) sizeof array / sizeof array[0]
const double PI = 4*atan(1);
typedef double function(int *i, double *q);
int DIMENSION;
//double DELTA_t = 0.015625; //1/2^6
//double DELTA_t = 0.001953125; //1/2^9
double DELTA_t = 0.00006103515625; //1/2^14

double two_body_problem(int *i, double q[])
{
  /*
  Initial conditions at perihelion:
  x(0) = .98 AU, y(0) = 0 AU
  v_x(0) = 0 AU/yr
  v_y(0) = .98 * 2π AU/yr
  */
  double G = 4*PI*PI;
  double r = sqrt(q[0]*q[0]+q[2]*q[2]);
  switch(*i)
  {
    //x:
    case 0: return q[1];
    case 1: return -G*q[0]/r*r*r;
    //y:
    case 2: return q[3];
    case 3: return -G*q[2]/r*r*r;
    default: return 0;
  }
}

double lorenz(int *i, double q[])
{
  //uses splot
  double sigma = 10;
  double r = 28;//reynolds
  double b = 8/3.0;
  switch(*i)
  {
    case 0: return sigma * (-q[0] + q[1]);
    case 1: return r * q[0] - q[1] -q[0] * q[2];
    case 2: return q[0] * q[1] - b * q[2];
    default: return 0;
  }
}

double harmonic_osc(int *i, double q[])
{
  double omega_squared = 10; //angular frequency
  switch(*i)
  {
    case 0: return q[1];
    case 1: return -omega_squared * q[0];
    default: return 0;
  }
}

double damped_harmonic_osc(int *i, double q[])
{
  double K = 10;
  double P = 0.5;
  switch(*i)
  {
    case 0: return q[1];
    case 1: return -K * q[0] -P * q[1];
    default: return 0;
  }
}

double pendulum(int *i, double q[])
{
  double g = 9.81; //acceleration due to gravity
  double L = 10; //String's length
  switch(*i)
  {
    case 0: return q[1];
    case 1: return -g/L * sin(q[0]);
    default: return 0;
  }
}

void rk4(function f, double t, double t_f, double q[])
{
	FILE *fp = fopen(PATH, "w");
	while(t <= t_f)
  {
    fprintf(fp, "%f", t);
    for (int i = 0; i < DIMENSION; ++i)
      fprintf(fp, " %f", q[i]);
    fprintf(fp, "\n");

    t += DELTA_t;
    double k1[DIMENSION];
    double k2[DIMENSION];
    double k3[DIMENSION];
    int i;
    for (i = 0; i < DIMENSION; ++i)
      k1[i] = q[i] + f(&i,q) * DELTA_t/2;
    for (i = 0; i < DIMENSION; ++i)
      k2[i] = q[i] + f(&i,k1) * DELTA_t/2;
    for (i = 0; i < DIMENSION; ++i)
      k3[i] = q[i] + f(&i,k2) * DELTA_t;
    for (i = 0; i < DIMENSION; ++i)
      q[i] += ( f(&i, q)
        + 2 * ( f(&i,k1) + f(&i,k2) )
        + f(&i,k3) ) * DELTA_t/6;
  }
  fclose(fp);
}

int main()
{
  //initial conditions:
 double q_0[] = {1,0,0,2*PI};
  //defining dimension:
 DIMENSION = sizeof_array(q_0);
  //Runge-Kutta 4th Order:
 rk4(two_body_problem, 0, 1, q_0);

  //Plotting with gnuplot
 char *gnuplot_command = {"plot 'euler.dat' u 2:4 w lines"};
 FILE *gnuplot_pipe = popen ("gnuplot -persistent", "w");
 fprintf(gnuplot_pipe, "%s\n", gnuplot_command);
 pclose(gnuplot_pipe);
 return 0;
}