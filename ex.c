#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

const int N = 1000;
const double G = 1;
//const double G = 6.673*pow(10, -11);
const int k = 4;
const double timeStop=0.000001;

typedef struct body {
double mass;
double p_x;
double p_y;
double v_x;
double v_y;
double a_x;
double a_y;
} body;


double findDist(body i, body j){
	return pow(pow(j.p_x-i.p_x,2)+pow(j.p_y-i.p_y, 2), .5);
}

double findDistTrap(body i, body j){
	printf("\nBodies with distance 0:\nbody at (%f, %f)\nbody at (%f, %f)\n", i.p_x, i.p_y, j.p_x, j.p_y);
	double x = pow(j.p_x-i.p_x,2);
	double y = pow(j.p_y-i.p_y, 2);
	printf("difference of x squared: %.17g\ndifference of y squared: %.17g\n", x,y);
	double res = pow(pow(j.p_x-i.p_x,2)+pow(j.p_y-i.p_y, 2), .5);
	double res2 = pow(x+y, .5);
	printf("res1: %f\nres2: %f\n", res,res2);
	return pow(pow(j.p_x-i.p_x,2)+pow(j.p_y-i.p_y, 2), .5);
}
/*  
written by Brennan Proudfoot
I certify that no unauthorized assistance has been received or
given in the completion of this work
*/  
int main() {
	body bodies[N];
	int i, j, t;
	for(i=0; i<N; i++){
    	body b = {1, i, i, 0, 0, 0, 0};
		bodies[i]=b;
        //printf("\nThe mass of body #%d is %f\nThe x position is %f\nThe y position is %f\n, the x velocity is %f\n, the y velocity is %f\n\n",i, bodies[i].mass, bodies[i].p_x, bodies[i].p_y, bodies[i].v_x, bodies[i].v_y);    
    }
    double start = omp_get_wtime() ;  
	for(t=0; t<k; t++){
		double forces[N][2];
		for(i=0; i<N; i++){
			for(j=0; j<N; j++){
				if(j<=i){continue;}
				double position = findDist(bodies[i], bodies[j]);
				double force_x = (G*bodies[i].mass*bodies[j].mass*(bodies[j].p_x-bodies[i].p_x))/(pow(position, 3));	
				double force_y = (G*bodies[i].mass*bodies[j].mass*(bodies[j].p_y-bodies[i].p_y))/(pow(position, 3));
				
				forces[i][0]+=force_x;
				forces[i][1]+=force_y;

				forces[j][0]-=force_x;
				forces[j][1]-=force_y;
				// bodies[i].a_x+=(force_x/bodies[i].mass);
				// bodies[i].a_y+=(force_y/bodies[i].mass);

				// bodies[j].a_x-=(force_x/bodies[j].mass);
				// bodies[j].a_y-=(force_y/bodies[j].mass);
			}
			bodies[i].a_x=(forces[i][0]/bodies[i].mass);
			bodies[i].a_y=(forces[i][1]/bodies[i].mass);
		}

		for(i=0; i<N; i++){
        	bodies[i].v_x=bodies[i].v_x+timeStop*bodies[i].a_x;
			bodies[i].v_y=bodies[i].v_y+timeStop*bodies[i].a_y;
			bodies[i].p_x=bodies[i].p_x+timeStop*bodies[i].v_x;
			bodies[i].p_y=bodies[i].p_y+timeStop*bodies[i].v_y;
			bodies[i].a_x=0;
			bodies[i].a_y=0;
		}

	}

	double end = omp_get_wtime();
	double elapsed = end - start;
	printf("Interactions per second = %.17g\n\n", (k*N*N)/elapsed);

	// for(i=0; i<N; i++){
 //        	printf("\nThe mass of body #%d is %f\nThe x position is %.17g\nThe y position is %.17g\n, the x velocity is %.17g\n, the y velocity is %.17g\n",i, bodies[i].mass, bodies[i].p_x, bodies[i].p_y, bodies[i].v_x, bodies[i].v_y);
	// }
	return 0;
}
