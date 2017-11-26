/*
 * COMP 633:  sequential implementation for n bodies (pa1a)
 *            all-pair interactions
 *
 * compile using icc or gcc:
 *   icc -O3 -fopenmp pa1-ref.c -lm -o pa1-ref
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// maximum value of n
#define  NMAX  20000 

/*
 * factored array of structs representation
 */

// body state required for computation of forces
struct Body {
  double x, y;
  double m;
};

// force on body
struct Force {
  double fx, fy;
};

// body state only used in integration
struct Velocity {
  double vx, vy;
};

static struct Body     P[NMAX];
static struct Force    F[NMAX];
static struct Velocity V[NMAX];

/* 
 * nbody simulation parameters
 */
double G      = 6.673E-11;
double DeltaT = 0.001;
int    nts    = 4; 

double fullPairs(int n){
  int i, j, k;
  double t1,t2;

  /* 
   * arbitrary initial configuration:
   *    n bodies evenly spaced along perimeter of a circle with radius n.  
   *    Each body has velocity 0 and unit mass.
   */
  #pragma omp parallel for schedule(static)
  for (i = 0; i < n; i++) {
    double theta = i * (2.0 * M_PI / n);
    P[i].x = n * cos(theta);
    P[i].y = n * sin(theta);
    P[i].m = 1.0;
    V[i].vx = V[i].vy = 0.0;
  }

  /*  
   *  time n-body simulation for nts steps
   */
  t1 = omp_get_wtime();
  #pragma omp parallel private(calcs, i, j)
  {
    for (k = 0; k < nts; k++) {

      #pragma omp for schedule(static) collapse(1)
      for (i = 0; i < n; i++) {

        double Fx = 0.0, Fy = 0.0;

        for (j = 0 ; j < n;  j++) {

            if (i != j) { 
              double rx, ry, d2, d, c;
              
              rx = P[j].x - P[i].x;   
              ry = P[j].y - P[i].y;   
              d2 = (rx * rx) + (ry * ry);
              c  = P[j].m / (d2 * sqrt(d2));           
              Fx += c * rx;              
              Fy += c * ry;               
            } /* i != j */

        } /* j */

        F[i].fx =  P[i].m * Fx;
        F[i].fy =  P[i].m * Fy;

      } /* i */

      /*
       * advance positions and velocities
       */
      for (i = 0; i < n; i++) {
        P[i].x  += V[i].vx * DeltaT;
        P[i].y  += V[i].vy * DeltaT;
        V[i].vx += ((G * F[i].fx) / P[i].m) * DeltaT; 
        V[i].vy += ((G * F[i].fy) / P[i].m) * DeltaT; 
      }

    } /* nts */
  }

  t2 = omp_get_wtime();
  double interactions = nts*n*n;
  double dt = t2 - t1;
  return dt;
}

double halfPairs(int n){
  int i, j, k;
  double t1,t2;

  /* 
   * arbitrary initial configuration:
   *    n bodies evenly spaced along perimeter of a circle with radius n.  
   *    Each body has velocity 0 and unit mass.
   */
  //#pragma omp parallel for schedule(static)
  for (i = 0; i < n; i++) {
    double theta = i * (2.0 * M_PI / n);
    P[i].x = n * cos(theta);
    P[i].y = n * sin(theta);
    P[i].m = 1.0;
    V[i].vx = V[i].vy = 0.0;
  }

  /*  
   *  time n-body simulation for nts steps
   */
  t1 = omp_get_wtime();
  //#pragma omp parallel private(calcs, i, j)
  {
    for (k = 0; k < nts; k++) {

     //#pragma omp for schedule(static) collapse(1)
      for (i = 0; i < n; i++) {


        for (j = i+1 ; j < n;  j++) {
          double rx, ry, d2, d, c, Fx, Fy;
          
          rx = P[j].x - P[i].x;   
          ry = P[j].y - P[i].y;   
          
          d2 = (rx * rx) + (ry * ry);
          c  = P[j].m / (d2 * sqrt(d2));           
          
          Fx = c * rx;              
          Fy = c * ry;   

          F[i].fx =  P[i].m * Fx;
          F[i].fy =  P[i].m * Fy;
          F[j].fx =  P[j].m * -Fx;
          F[j].fy =  P[j].m * -Fy;
        } /* j */
      } /* i */

      /*
       * advance positions and velocities
       */
      for (i = 0; i < n; i++) {
        P[i].x  += V[i].vx * DeltaT;
        P[i].y  += V[i].vy * DeltaT;
        V[i].vx += ((G * F[i].fx) / P[i].m) * DeltaT; 
        V[i].vy += ((G * F[i].fy) / P[i].m) * DeltaT; 
      }

    } /* nts */
  }

  t2 = omp_get_wtime();
  double interactions = nts*n*n;
  double dt = t2 - t1;
  return dt;
}  

/*
 * reference implementation
 * 
 * pa1-ref <n>
 *  
 * reports sequential performance in units of millions of 
 * interactions per second for given argument n
 */
int main(int argc, char * argv[])
{
  FILE* fp = fopen("simTimes.csv","w+");
  int n[11] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000,20000};
  int i;
  for(i = 0; i<11; i++){
    double fullPairs = fullPairs(n[i]);//(RUNS*pow(n[i],2))/sim(n[i]);
    double halfPairs = halfPairs(n[i]);
    fprintf(fp, "%d,%.10lg,%.10lg\n",n[i],fullPairs,halfPairs);
    printf("%d, %.10lg,%.10lg\n",n[i], fullPairs,halfPairs);
    // double slowTime = simSlow(n[i]);//(RUNS*pow(n[i],2))/simSlow(n[i]);
    // fprintf(fp, "%d,%.10lg,%.10lg;\n",n[i],slowTime, time);
    // printf("%d, %lg, %lg\n",n[i], slowTime, time);
  }
  fclose(fp);
} 


