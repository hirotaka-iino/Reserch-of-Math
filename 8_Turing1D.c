#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define      PI     (3.14159265358979)
#define      N      (200)

double a = 3.0;//3.0;
double b = 2.0;//2.0;
double c = 0.0;//0.0;
double d = 1.0;//1.0;

double f(double u, double v){ return d * (u * (1 - u * u) - v);}
double g(double u, double v){ return a * (u - c) - b * v;}
double rand_range(double min,double max){ return (max - min) * (rand()) / RAND_MAX + min;}

int main(void)
{
    double  Lx = 100.0;
    double  dx = Lx / N;
    double  dt = 0.0005;
    double  du = 1.0;
    double  dv = 1.0;
    double  lambda_u = du * dt / (dx * dx);
    double  lambda_v = dv * dt / (dx * dx);
    
	FILE *gp;
	gp = popen("gnuplot -persist","w");
	fprintf(gp,"set terminal x11\n");
    //fprintf(gp,"set terminal aqua\n");
	fprintf(gp,"set xrange[0:%f]\n",Lx);
	fprintf(gp,"set yrange[-0.5:0.5]\n");
	double	u0[N + 2];
	double	u1[N + 2];
    double  v0[N + 2];
    double  v1[N + 2];

    int INTV = 100;

	for(int i = 1;i <= N;i ++)
	{
        double Max =  0.5; //Maximum of random number
        double Min = -0.5; //Minimum of random number
		u0[i] = rand_range(Min,Max);
        v0[i] = rand_range(Min,Max);
	}
	for(int i_time = 0;;i_time ++)
	{
      
		//B.C.
		u0[0] = u0[1];
		u0[N + 1] = u0[N];
        v0[0] = v0[1];
        v0[N + 1] = v0[N];

		//plot
        if(i_time % INTV == 0)
        {
            fprintf(gp,"plot '-' with lines\n");
            for(int i = 1;i <= N;i ++)
            {
                double x = (i - 0.5) * dx;
                double y = v0[i];
                //double m = 19;
                //y = 0.2 * cos(m * PI * x / Lx);
                fprintf(gp,"%f %f\n",x,y);
            }
            fprintf(gp,"e\n");
            fflush(gp);
        }
		//Calc Eq
		for(int i = 1;i <= N;i ++)
		{
			u1[i] = u0[i]
            + lambda_u * (u0[i - 1] - 2 * u0[i] + u0[i + 1])
            + dt * f(u0[i], v0[i]);
            v1[i] = v0[i] + lambda_v * (v0[i - 1] - 2 * v0[i] + v0[i + 1])
            + dt * g(u0[i], v0[i]);

		}

		//Subs
		for(int i = 1;i <= N;i ++)
		{
			u0[i] = u1[i];
            v0[i] = v1[i];
		}
	}
	fclose(gp);
	return 0;
}

