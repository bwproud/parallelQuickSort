#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>
#include <cilk/cilk.h>

#define  NMAX  75000000 

static double   N[NMAX];
static int lt[NMAX];
static int gt[NMAX];
static double local[NMAX];

void printArray(int n){
    int j;
    printf("[");
    int t =0;
    for(j = 0; j<n; j++){
          if(t){
            printf(", %f", N[j]);
          }else{
            t=1;
            printf("%f", N[j]);
          }
    }
    printf("]\n");
}

double drand ( double low, double high )
{   
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

void fillArrayRandom(int n){
    int j;
    for(j = 0; j<n; j++){

        double r = drand(0,1000);
        N[j]=r;
    }
}

void fillSorted(int n){
    int j;
    N[0]=0.5;
    for(j = 1; j<n; j++){
        double r = N[j-1]+1;
        N[j]=r;
    }
}

void fillReverseSorted(int n){
    int j;
    N[0]=n+0.5;
    for(j = 1; j<n; j++){
        double r = N[j-1]-1;
        N[j]=r;
    }
}

void fillSame(int n){
    int j;
    double r = drand(0,1000);
    for(j = 0; j<n; j++){
        N[j]=r;
    }
}

int partition(int p, int r){
  double key=N[r];
  int i=p-1;
  int j;
  double temp;
  for(j=p; j<r; j++){
    if(N[j]<=key){
      i+=1;
      temp = N[i];
      N[i]=N[j];
      N[j]=temp;
    }  
  } 
  temp = N[i+1];   
  N[i+1]=N[r];
  N[r]=temp;
  return i+1;
}

void quickSortHelper(int p, int r){
    if(p<r){
        int q=partition(p,r);
   	quickSortHelper(p,q-1);
	quickSortHelper(q+1,r);   
    }
}

double sequentialQuickSort(int n){
    double t1, t2;
    t1 = omp_get_wtime();
    
    quickSortHelper(0, n-1);
    
    t2 = omp_get_wtime();
    
    return t2-t1;
}


void insertionSortHelper(int p, int r){
    double key;
    int j, i;

    for (i = p+1; i<r+1 ; i++){
      key = N[i];
      j = i-1;

      while (j >= p && N[j] > key){
         N[j+1] = N[j];
         j--;
      }
      N[j+1] = key;
    }
}
void prefixSum(int arr[], int p, int r){
    int i;
    for(i=p+1;i<r+1;i++){
        arr[i]+=arr[i-1];
    }
}

int log_2(int n){
  int i=0;
  while(n >>= 1) {++i;}
  return i;
}

void parallelPrefixSum(int p, int r){
        int len = r-p+1;
        int shift, j, h;
        int k = log_2(len);
        for(h=1; h<k+1;h++){
                shift = 1<<h;
                cilk_for (j=1; j<(len/shift)+1;j++){
                        lt[p+j*shift-1]+=lt[p+j*shift-(shift/2)-1];
                        gt[p+j*shift-1]+=gt[p+j*shift-(shift/2)-1];
                }
        }

        for(h=k; h>-1;h--){
                shift = 1<<h;
                cilk_for (j=2; j<(len/shift)+1;j++){
                        if(j%2==1){
                                lt[p+j*shift-1]+=lt[p+j*shift-shift-1];
                                gt[p+j*shift-1]+=gt[p+j*shift-shift-1];
                        }
                }
        }
}

int parallelPartition(int p, int r){
  double key=N[r];
  int i,j;
  double temp;

    cilk_for (i=p; i<r+1; i++){
      lt[i]=0;
      gt[i]=0;
      local[i]=N[i];
    }

    cilk_for (i = p; i <r; i++){
        if(N[i]<key){
            lt[i]=1;
            gt[i]=0;
        }else{
            lt[i]=0;
            gt[i]=1;
        }
    }
    
    parallelPrefixSum(p,r);
  
    int pivot = lt[r];
    N[pivot+p]=key;

    cilk_for (i=p; i<r; i++){
        if(local[i]<key){
            int index = p+lt[i]-1;
            N[index]=local[i];
        }else{
            int index = p+pivot+gt[i];
            N[index]=local[i];
        }    
    }
    

  return pivot+p;
}

void psqHelper(int p, int r, int size){
    if(p<r){
        if(r-p<=50){
            insertionSortHelper(p,r);
        }else{
	int q;
	if(r-p < 0.5*size){
            q = partition(p,r);
        }else{
            q=parallelPartition(p,r);
        }
            cilk_spawn psqHelper(p,q-1, size);
            psqHelper(q+1,r, size);
        }  
    }    
}

double parallelQuickSort(int n){
    double t1, t2;
    #pragma omp master
    t1 = omp_get_wtime();
    
      psqHelper(0, n-1, n);
    
    #pragma omp master
    t2 = omp_get_wtime();
    
    return t2-t1;
}

int checkArray(int n){
    int j;
    for(j = 0; j<n-1; j++){
        if(N[j]>N[j+1]){
          return -1;
        }
    }
    return 0;
}

int main(int argc, char * argv[]){
  __cilkrts_set_param("nworkers", argv[1]);
  FILE* fp = fopen("simTimes.csv","a+");
  int len=16;
  int n[] = {3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147, 531441, 1594323, 4782969, 14348947, 43046721};
  int i;
  srand(getpid());
  if (atoi(argv[1]) == 1){
	  for(i=0; i<len; i++){
	  	fillArrayRandom(n[i]);
	  	double t = sequentialQuickSort(n[i]);
	  	fprintf(fp,"1,%d,%f\n",n[i],t);
	  }
   }
  else{
	  for(i = 0; i<len; i++){
	    fillArrayRandom(n[i]);
	    double t = parallelQuickSort(n[i]);
	    int numworkers = __cilkrts_get_nworkers();
	    printf("%d elements sorted in %f time with %d workers\n", n[i], t, numworkers);
	    fprintf(fp,"%d,%d,%f\n",numworkers,n[i],t);
	    if(checkArray(n[i])==-1){
	      printf("SORT FAILED\n");
	    }else{
	      printf("SUCCESSFUL SORT\n");
	    }
	  }
	}
  fclose(fp);
} 


