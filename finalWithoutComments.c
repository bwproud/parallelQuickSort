#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

// maximum value of n
#define  NMAX  133500000 
//#define  NMAX  200
#define CHUNKSIZE 20

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

int cmpfunc (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b)
    return 1;
  else if (*(double*)a < *(double*)b)
    return -1;
  else
    return 0;  
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
    double t1;
    t1 = omp_get_wtime();
    
    quickSortHelper(0, n-1);
    
    double t2;
    t2 = omp_get_wtime();
    
    return (double)(t2-t1)/ CLOCKS_PER_SEC;
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
//                #pragma omp parallel for schedule(static) private(j)
                for(j=1; j<(len/shift)+1;j++){
                        lt[p+j*shift-1]+=lt[p+j*shift-(shift/2)-1];
                        gt[p+j*shift-1]+=gt[p+j*shift-(shift/2)-1];
                }
        }

        for(h=k; h>-1;h--){
                shift = 1<<h;
  //              #pragma omp parallel for schedule(static) private(j)
                for(j=2; j<(len/shift)+1;j++){
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

//  #pragma omp parallel
//  {
  //  #pragma omp for schedule(static) private(i)
    for(i=p; i<r+1; i++){
      lt[i]=0;
      gt[i]=0;
      local[i]=N[i];
    }

  //  #pragma omp for schedule(static) private(i)
    for(i = p; i <r; i++){
        if(N[i]<key){
            lt[i]=1;
            gt[i]=0;
        }else{
            lt[i]=0;
            gt[i]=1;
        }
    }
 // }


  //parallelPrefixSum(p,r);
  prefixSum(lt, p,r);
  prefixSum(gt,p,r);

  int pivot = lt[r];
  N[pivot+p]=key;

//  #pragma omp parallel
//  {
  //  #pragma omp for schedule(static) private(i)
    for(i=p; i<r; i++){
        if(local[i]<key){
            int index = p+lt[i]-1;
            N[index]=local[i];
        }else{
            int index = p+pivot+gt[i];
            N[index]=local[i];
        }    
    }
 // }  

  return pivot+p;
}

void psqHelper(int p, int r){
    if(p<r){
        if(r-p<=50){
            insertionSortHelper(p,r);
        }else{
          int q=parallelPartition(p,r);
          
//         #pragma omp parallel
  //       {
    //        #pragma omp task
            psqHelper(p,q-1);

      //      #pragma omp task
	    psqHelper(q+1,r);

	//    #pragma omp taskwait
        //  }
        }  
    }    
}

double parallelQuickSort(int n){
    time_t t1;
    t1 = omp_get_wtime();
    
    psqHelper(0, n-1);
    
    time_t t2;
    t2 = omp_get_wtime();
    
    return (double)(t2-t1)/ CLOCKS_PER_SEC;
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

void tester(int n){
    srand(getpid());
    fillArrayRandom(n);
    printArray(n);
    double t = parallelQuickSort(n);
    printArray(n);
}
int main(int argc, char * argv[]){
  FILE* fp = fopen("simTimes.csv","w+");
  int len=15;
  int n[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000,20000,200000,2000000,20000000,133500000};
  int i;
//  srand(getpid());
  for(i = 0; i<len; i++){
    fillArrayRandom(n[i]);
    double t = parallelQuickSort(n[i]);
    printf("%d elements sorted in %f time\n", n[i], t);
    if(checkArray(n[i])==-1){
      printf("SORT FAILED\n");
    }else{
      printf("SUCCESSFUL SORT\n");
    }
  }
  fclose(fp);
} 


