#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <omp.h>
#include <cilk/cilk.h>

// maximum value of n
#define  NMAX  75000000 
//#define  NMAX  200
#define CHUNKSIZE 20

typedef struct{
    int lt;
    int gt;
} flag;

flag *flags;
double *local;
int size;
// void printArray(int n){
//     int j;
//     printf("[");
//     int t =0;
//     for(j = 0; j<n; j++){
//         //if(j<10){
//           if(t){
//             printf(", %f", N[j]);
//           }else{
//             t=1;
//             printf("%f", N[j]);
//           }
          
//         //}
//     }
//     printf("]\n");
// }

double drand ( double low, double high )
{   
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

void fillArrayRandom(double *arr, int n){
    int j;
    cilk_for(j = 0; j<n; j++){

        double r = drand(0,1000);
        arr[j]=r;
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

void swap(double *xp, double *yp){
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
// void merge(int l, int m, int r){
//     int i, j, k;
//     int n1 = m - l + 1;
//     int n2 =  r - m;
 
//     /* create temp arrays */
//     int L[n1], R[n2];
 
//     /* Copy data to temp arrays L[] and R[] */
//     for (i = 0; i < n1; i++)
//         L[i] = N[l + i];
//     for (j = 0; j < n2; j++)
//         R[j] = N[m + 1+ j];
 
//     /* Merge the temp arrays back into arr[l..r]*/
//     i = 0; // Initial index of first subarray
//     j = 0; // Initial index of second subarray
//     k = l; // Initial index of merged subarray
//     while (i < n1 && j < n2){
//         if (L[i] <= R[j]){
//             N[k] = L[i];
//             i++;
//         }else{
//             N[k] = R[j];
//             j++;
//         }
//         k++;
//     }
 
//     /* Copy the remaining elements of L[], if there
//        are any */
//     while (i < n1){
//         N[k] = L[i];
//         i++;
//         k++;
//     }
 
//     /* Copy the remaining elements of R[], if there
//        are any */
//     while (j < n2){
//         N[k] = R[j];
//         j++;
//         k++;
//     }
// }

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
// void mergeSortHelper(int l, int r){
//     if (l < r){
//         int m = l+(r-l)/2;
 
//         // Sort first and second halves
//         mergeSortHelper(l, m);
//         mergeSortHelper(m+1, r);
 
//         merge(l, m, r);
//     }
// }

// double mergeSort(int n){
//     double t1, t2;
//     #pragma omp master
//     t1 = omp_get_wtime();
    
//     mergeSortHelper(0,n-1);
    
//     #pragma omp master
//     t2 = omp_get_wtime();
//     return t2-t1;
// }

// double insertionSort(int n){
//     int key, j, i;
//     double t1, t2;
//     #pragma omp master
//     t1 = omp_get_wtime();

//     for (i = 1; i < n; i++){
//       key = N[i];
//       j = i-1;

//       while (j >= 0 && N[j] > key){
//          N[j+1] = N[j];
//          j--;
//       }
//       N[j+1] = key;
//     }
//     #pragma omp master
//     t2 = omp_get_wtime();
    
//     return t2-t1;
// }

int partition(double *arr, int p, int r){
  double key=arr[r];
  int i=p-1;
  int j;
  double temp;
  for(j=p; j<r; j++){
    if(arr[j]<=key){
      i+=1;
      temp = arr[i];
      arr[i]=arr[j];
      arr[j]=temp;
    }  
  } 
  temp = arr[i+1];   
  arr[i+1]=arr[r];
  arr[r]=temp;
  return i+1;
}

// void quickSortHelper(int p, int r){
//     if(p<r){
//         int q=partition(p,r);
	

// 		#pragma omp task
//     {
//     			quickSortHelper(p,q-1);
// 		}	

// 	quickSortHelper(q+1,r);   
// }
// }

// double sequentialQuickSort(int n){
//     double t1, t2;
//     #pragma omp master
//     t1 = omp_get_wtime();
//     #pragma omp parallel
//   {
//      We only want our master thread to be executed once, thus we use the singel construct here.
//       nowait is used becuse we have no need for synchronization at the end of the region 
//     #pragma omp single nowait
//     {
//       quickSortHelper(0, n-1);
//     }
//   }
    
    
//     #pragma omp master
//     t2 = omp_get_wtime();
    
//     return t2-t1;
// }


void insertionSortHelper(double *arr, int p, int r){
    double key;
    int j, i;

    for (i = p+1; i<r+1 ; i++){
      key = arr[i];
      j = i-1;

      while (j >= p && arr[j] > key){
         arr[j+1] = arr[j];
         j--;
      }
      arr[j+1] = key;
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
                // #pragma omp parallel for schedule(static) private(j)
                cilk_for (j=1; j<(len/shift)+1;j++){
                        flags[p+j*shift-1].lt+=flags[p+j*shift-(shift/2)-1].lt;
                        flags[p+j*shift-1].gt+=flags[p+j*shift-(shift/2)-1].gt;
                }
        }

        for(h=k; h>-1;h--){
                shift = 1<<h;
                // #pragma omp parallel for schedule(static) private(j)
                cilk_for (j=2; j<(len/shift)+1;j++){
                        if(j%2==1){
                                flags[p+j*shift-1].lt+=flags[p+j*shift-shift-1].lt;
                                flags[p+j*shift-1].gt+=flags[p+j*shift-shift-1].gt;
                        }
                }
        }
}

int parallelPartition(double *arr, int p, int r){
//  printf("p %d r %d\n", p, r);
  double key=arr[r];
  int i,j;
  double temp;

//  printf("%d: before first parallel region\n",omp_get_thread_num());
    //#pragma omp for schedule(static) private(i)
    cilk_for (i=p; i<r+1; i++){
      flags[i].lt=0;
      flags[i].gt=0;
      local[i]=arr[i];
    }

    //#pragma omp for schedule(static) private(i)
    cilk_for (i = p; i <r; i++){
        if(arr[i]<key){
            flags[i].lt=1;
            flags[i].gt=0;
        }else{
            flags[i].lt=0;
            flags[i].gt=1;
        }
    }

//  printf("before less than: [");
  // for(i=p; i<r+1; i++){
  //    printf("%d, ", lt[i]);
  // }
  // printf("]\n");
  // printf("before greater than: [");
  // for(i=p; i<r+1; i++){
  //    printf("%d, ", gt[i]);
  // }
  // printf("]\n");
  
  
 // printf("%d: before prefix sum\n",omp_get_thread_num());

  parallelPrefixSum(p,r);
  // prefixSum(lt, p,r);
  // prefixSum(gt,p,r);

 // printf("%d: after prefix sum\n",omp_get_thread_num());
  //prefixSum(lt, gt, p, r);

  // printf("after less than: [");
  // for(i=p; i<r+1; i++){
  //    printf("%d, ", lt[i]);
  // }
  // printf("]\n");
  // printf("after greater than: [");
  // for(i=p; i<r+1; i++){
  //    printf("%d, ", gt[i]);
  // }
  // printf("]\n");

  int pivot = flags[r].lt;
//  printf("pivot point is %d\n",pivot);
  arr[pivot+p]=key;

    //#pragma omp for schedule(static) private(i)
    cilk_for (i=p; i<r; i++){
        if(local[i]<key){
            int index = p+flags[i].lt-1;
            arr[index]=local[i];
        }else{
            int index = p+pivot+flags[i].gt;
            arr[index]=local[i];
        }    
    }
    

 // printf("%d: after second parallel\n", omp_get_thread_num());
  return pivot+p;
}

void psqHelper(double *arr, int p, int r){
    int q;
    if(p<r){
        if(r-p<=50){
            insertionSortHelper(arr, p,r);
        }else{
          if(r-p < 0.5*size){
            q = partition(arr,p,r);
          }else{
            q=parallelPartition(arr,p,r);
	  }
//	  printf("left p: %d left r: %d\n", p, q-1);          
//          printf("right p: %d right r: %d\n", q+1, r); 
            cilk_spawn psqHelper(arr,p,q-1);
            psqHelper(arr,q+1,r);
        }  
    }    
}

double parallelQuickSort(double *arr, int n){
    double t1, t2;
    flags = malloc(sizeof(flag)*n);
    local = malloc(sizeof(double)*n);
    int i;
    cilk_for(i=0; i<n; i++){
      local[i]=0;
      flags[i].lt=0;
      flags[i].gt=0;
    }
    #pragma omp master
    t1 = omp_get_wtime();
    
  //   #pragma omp parallel
  // {
  //    We only want our master thread to be executed once, thus we use the singel construct here.
  //     nowait is used becuse we have no need for synchronization at the end of the region 
  //   #pragma omp single nowait
  //   {
      psqHelper(arr, 0, n-1);
  //   }
  // }
    
    #pragma omp master
    t2 = omp_get_wtime();
    
    return t2-t1;
}

// double selectionSort(int n){
//     int j, min_idx,i;
//     double t1,t2;
//     #pragma omp master
//     t1 = omp_get_wtime();

//     // One by one move boundary of unsorted subarray
//     for (i = 0; i < n-1; i++){
        
//         // Find the minimum element in unsorted array
//         min_idx = i;

//         for (j = i+1; j < n; j++){
//             if (N[j] < N[min_idx]){
//                 min_idx = j;
//             }  
//         }
//         double temp = N[i];
//         N[i] = N[min_idx];
//         N[min_idx]=temp;
//     }
//     #pragma omp master
//     t2 = omp_get_wtime();
    
//     return t2-t1;
// }

int checkArray(double *arr, int n){
    int j;
    for(j = 0; j<n-1; j++){
        if(arr[j]>arr[j+1]){
          return -1;
        }
    }
    return 0;
}

// void tester(int n){
//     srand(getpid());
//     fillArrayRandom(n);
//     printArray(n);
//     double t = parallelQuickSort(n);
//     printArray(n);
// }
int main(int argc, char * argv[]){
  FILE* fp = fopen("simTimes.csv","w+");
  int len=15;
  //int len=5;
  int n[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000,20000,200000,2000000,20000000,75000000};
  int i;
  srand(getpid());
  //tester(10);
  for(i = 0; i<len; i++){
    size = n[i];
    double *arr = malloc(sizeof(double)*n[i]);
    fillArrayRandom(arr, n[i]);
    //printf("Before for array of size %d:\n", n[i]);
    //printArray(n[i]);
    double t = parallelQuickSort(arr, n[i]);
    printf("%d elements sorted in %f time\n", n[i], t);
    if(checkArray(arr, n[i])==-1){
      printf("SORT FAILED\n");
    }else{
      printf("SUCCESSFUL SORT\n");
    }
    //printf("after for array of size %d:", n[i]);
    //printArray(n[i]);
  }
  fclose(fp);
} 


