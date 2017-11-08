#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef void function(double, double[], double[]);

const double sigma = 10;
const double r = 28;
const double b = 8/3.0;

void dydt(double t, double y[], double dydt[])
{
    dydt[0] = sigma * (-y[0] + y[1]);
    dydt[1] = r * y[0] - y[1] -y[0] * y[2];
    dydt[2] = y[0] * y[1] - b * y[2];
}
void rk4(int dof, double dx, double x, double y[dof], function dydx)
{
    double hh = dx*0.5;
    double h6 = dx/6.0;
    double xh = x + hh;

    double yt[dof];
    double k1[dof];
    double k2[dof];
    double k3[dof];
    dydx(x,y,k1);
    for (int i = 0; i < dof; ++i)
        yt[i] = y[i] + hh*k1[i];
    dydx(xh,yt,k2);
    for (int i = 0; i < dof; ++i)
        yt[i] = y[i] + hh*k2[i];
    dydx(xh,yt,k3);
    for (int i = 0; i < dof; ++i){
        yt[i] = y[i] + dx*k3[i];
        k3[i] += k2[i];
    }
    dydx(x+dx,yt,k2);
    for (int i = 0; i < dof; ++i)
        y[i] = y[i] + h6*(k1[i]+k2[i]+2.0*k3[i]);
}

int main()
{
    const int dof = 3;
    const double dt = .001;
    double y[dof] = {0,10,1};
    FILE *gp = popen("gnuplot","w");
    char *cmd = "set terminal epslatex standalone lw 3\n\
                 set output 'gptemp.tex'\n\
                 set colorsequence podo\n\
                 set sample 300\n\
                 set xlabel '$x$'\n\
                 set ylabel '$y$'\n\
                 plot '-' w l not\n";
    fprintf(gp,"%s",cmd);
    for (double t=0; t < 40; t+=dt) {
        fprintf(gp, "%f %f\n", y[0], y[1]);
        rk4(dof,dt,t,y,dydt);
    }
    fprintf(gp, "e\n");
    pclose(gp);
    system("pdflatex -interaction batchmode gptemp.tex && mv gptemp.pdf output.pdf && rm -f gptemp*");
  return 0;
}
