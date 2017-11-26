#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
const double G = 1.000000000000;
const double TIMESTEP = 0.000001;
const double BODY_MASS = 1.0;
const double TOTAL_TIME = 6.325914;
const int RUNS = 4;
struct Body{
		double mass;
		double a[2]; //accel vector
		double v[2]; //velocity vector
		double r[2]; //pos vector
};

void printSystem(struct Body* b, int size){
	int i;
	for(i=0; i<size; i++){
		printf("%.8lg, %.8lg\n", b[i].r[0], b[i].r[1]);
	}
}

double sim(int NUM_BODIES){
	int i,j;
	struct Body bodies[NUM_BODIES];
	for(i = 0; i<NUM_BODIES; i++){
		double nextPos = -1.0 + i * 2.0/(NUM_BODIES-1);
		struct Body body = {BODY_MASS, {0, 0}, {0, 0}, {nextPos, nextPos}};
		bodies[i] = body;
	}

	int calcs;
	double start = omp_get_wtime();
	for(calcs = 0; calcs < RUNS; calcs++){
		for(i = 0; i<NUM_BODIES; i++){
			for(j = i+1; j<NUM_BODIES; j++){
					double rij[2] = {bodies[j].r[0] - bodies[i].r[0], bodies[j].r[1] - bodies[i].r[1]};
					double sharedCalc = (G * bodies[i].mass * bodies[j].mass) / 
						pow(pow(rij[0], 2) + pow(rij[1], 2),1.5);
					double currForce[2] = {rij[0] * sharedCalc, rij[1] * sharedCalc};
					bodies[i].a[0] += currForce[0] / bodies[i].mass;
					bodies[i].a[1] += currForce[1] / bodies[i].mass;
					bodies[j].a[0] += -currForce[0] / bodies[j].mass;
					bodies[j].a[1] += -currForce[1] / bodies[j].mass; 
			}
		}

		for(i = 0; i<NUM_BODIES; i++){
			bodies[i].v[0] += bodies[i].a[0] * TIMESTEP;
			bodies[i].v[1] += bodies[i].a[1] * TIMESTEP;
			bodies[i].r[0] += bodies[i].v[0] * TIMESTEP;
			bodies[i].r[1] += bodies[i].v[1] * TIMESTEP;
			bodies[i].a[0] = 0;
			bodies[i].a[1] = 0;
		}
	}
	return omp_get_wtime() - start;
}

double simSlow(int NUM_BODIES){
	int i,j;
	struct Body bodies[NUM_BODIES];
	for(i = 0; i<NUM_BODIES; i++){
		double nextPos = -1.0 + i * 2.0/(NUM_BODIES-1);
		struct Body body = {BODY_MASS, {0, 0}, {0, 0}, {nextPos, nextPos}};
		bodies[i] = body;
	}
	int calcs;
	double start = omp_get_wtime();
	for(calcs = 0; calcs <RUNS; calcs++){
		for(i = 0; i<NUM_BODIES; i++){
			double totalForce[2] = {0,0};
			for(j = 0; j<NUM_BODIES; j++){
				if(i!=j){
					double rij[2] = {bodies[j].r[0] - bodies[i].r[0], bodies[j].r[1] - bodies[i].r[1]};
					double sharedCalc = (G * bodies[i].mass * bodies[j].mass) / 
						pow(pow(rij[0], 2) + pow(rij[1], 2),1.5);
					totalForce[0] += rij[0] * sharedCalc;
					totalForce[1] += rij[1] * sharedCalc;
				}
			}
			bodies[i].a[0] = totalForce[0] / bodies[i].mass;
			bodies[i].a[1] = totalForce[1] / bodies[i].mass; 
		}

		for(i = 0; i<NUM_BODIES; i++){
			bodies[i].v[0] += bodies[i].a[0] * TIMESTEP;
			bodies[i].v[1] += bodies[i].a[1] * TIMESTEP;
			bodies[i].r[0] += bodies[i].v[0] * TIMESTEP;
			bodies[i].r[1] += bodies[i].v[1] * TIMESTEP;
		}
	}
	return omp_get_wtime() - start;
}

int main(){
	FILE* fp = fopen("simTimes.csv","w+");
	int n[10] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000};
	int i,j;
	for(j=0; j< 5; j++){
	for(i = 0; i<10; i++){
		double time = (RUNS*pow(n[i],2))/sim(n[i]);
		double slowTime = (RUNS*pow(n[i],2))/simSlow(n[i]);
		fprintf(fp, "%d,%.17g,%.17g\n",n[i],slowTime, time);
		printf("%d, %.17g, %.17g\n",n[i], slowTime, time);
	}
	}
	fclose(fp);
}
