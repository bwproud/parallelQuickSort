
/* multithreaded quickSort
usage with gcc (version 4.2 or higher required):
     gcc -O -fopenmp -o quickSort-openmp quickSort-openmp.c 
     ./quickSort-openmp size numWorkers
*/

#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#define MAXSIZE 50000000  /* maximum array size */
#define MAXWORKERS 10 /* maximum number of workers */

void quickSort(int *inputArray, int size);
void swap(int *inputArray, int leftIndex, int rightIndex);
int partition(int*N, int p, int r){
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

void quickSortHelper(int* N, int p, int r){
    if(p<r){
        int q=partition(N,p,r);
        quickSortHelper(N,p,q-1);
        quickSortHelper(N,q+1,r);
    }    
}

double sequentialQuickSort(int* N, int n){
    double t1, t2;
    #pragma omp master
    t1 = omp_get_wtime();
    
    quickSortHelper(N,0, n-1);
    
    #pragma omp master
    t2 = omp_get_wtime();
    
    return t2-t1;
}

double start_time, end_time, s_start_time, s_end_time; /* start and end times */
long size; /* array size */
/* ---------------------------------------------------------------------------- */
/* read command line, initialize, and create threads */
int main(int argc, char *argv[]) {
  int i, numWorkers;
  
  /* read command line args if any */
  size = (argc > 1)? atoi(argv[1]) : MAXSIZE;
  if (size > MAXSIZE){
    size = MAXSIZE;
  }

  if (size < 2){
    printf("Array is only one element!\n");
    return 0; /* Invalid input array*/
  }

  numWorkers = (argc > 2)? atoi(argv[2]) : MAXWORKERS;
  if (numWorkers > MAXWORKERS) numWorkers = MAXWORKERS;
  omp_set_num_threads(numWorkers); /* Set number of threads */

  int *inputArray; /* testDataArray */
  inputArray = malloc(sizeof(int) * size); /* Allocate in memory instead */
  srand(time(NULL));
  /* Create testData array */
  for (i = 0; i < size; i++) {
    inputArray[i] = rand()%99999;
  }

  #ifdef DEBUG
    printf("array size: %ld \n", size);
    printf("numWorkers: %d \n", numWorkers);
    printf("[ ");
    for (i = 0; i < size; i++) {
      printf(" %d", inputArray[i]);
    }
    printf(" ]\n");
  #endif

  start_time = omp_get_wtime();
  /* Call the quickSort function to sort the list */
  #pragma omp parallel
  {
    /* We only want our master thread to be executed once, thus we use the singel construct here.
      nowait is used becuse we have no need for synchronization at the end of the region */
    #pragma omp single nowait
    {
      quickSort(inputArray, size);
    }
  }

  /* get end time */
  end_time = omp_get_wtime();
  /*int myid = omp_get_thread_num(); */
  /* print results */
  for (i = 0; i < size; i++) {
    inputArray[i] = rand()%99999;
  }
  s_start_time = omp_get_wtime();
  sequentialQuickSort(inputArray, size);
  s_end_time = omp_get_wtime();
  #ifdef PRINT
  printf("[ ");
  for (i = 0; i < size; i++) {
    printf(" %d", inputArray[i]);
  }
  printf(" ]\n");
  #endif
  free(inputArray);
  printf("The execution time is %g sec\n Seq execution in %g sec\n", end_time - start_time, s_end_time - s_start_time);
  return 0;
}

void quickSort(int *inputArray, int size){
  int pivot, leftIndex, rightIndex;
  /* End of reccursion */
  if (size <= 1) { return; }
  /* Set pivot */
  pivot = inputArray[size/2];
  for(leftIndex = 0, rightIndex = size -1;; leftIndex++, rightIndex--) {
    while(inputArray[leftIndex] < pivot){
      leftIndex++;
    }
    while(pivot < inputArray[rightIndex]){
      rightIndex--;
    }
    if(rightIndex <= leftIndex){
      break;
    }
    swap(inputArray, leftIndex, rightIndex);    
  }
  #pragma omp task if(rightIndex-leftIndex > 1000)
  {
    quickSort(inputArray, leftIndex); /* Sort lower */
  }
  //#pragma omp task
  //{
    quickSort(inputArray + rightIndex + 1, size - rightIndex -1); /* Sort upper */
  //}
}
/* ---------------------------------------------------------------------------- */
/* Swaps two elements */
void swap(int *inputArray, int leftIndex, int rightIndex){
  int temp;
  temp = inputArray[leftIndex];
  inputArray[leftIndex] = inputArray[rightIndex];
  inputArray[rightIndex] = temp;
}