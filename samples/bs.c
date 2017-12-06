/* Standard includes */
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

/* Custom includes */
#include "quicksort_common.h"

/* Parameters for degrading to sequential work. */
#define SERIAL_INSERTION_NSIZE 32       // Insertion sort when subproblem size = 32
#define SERIAL_PARTITION_N_FACTOR 0.7   // Sequential partition when size' = 0.7*size

/* Structure for comparison flags used in parallel partition. */
typedef struct{
    int lte;
    int gt;
} lte_gt;

/* Global arrays: flags and copy array used for swapping in parallel. */
lte_gt *flags;
long *copyArray;    

/* Some constants that will be referenced */
int workers, problem_size;

/* Inclusive, in-place parallel prefix sum.
   Computes the prefix sum of a flag array within a range S[left...left+n]
   n = 2^k
*/
void parallel_prefix_sum(lte_gt *S, int left, int n, int k) {

    int i, h;
    for(h = 1 ; h <= k ; h++) {
        cilk_for (i = 1; i <= (n >> h); i++) {
            S[left+i * (1 << h) - 1].lte += S[left+(1 << h) * i - (1 << (h-1)) - 1].lte;
            S[left+i * (1 << h) - 1].gt += S[left+(1 << h) * i - (1 << (h-1)) - 1].gt;
        }
    }

    for (h = k ; h >= 1; h--) {
        cilk_for(i = 2; i <= (n >> (h-1)); i++){
            if (i % 2) {
                S[left+i * (1 << (h-1)) -1].lte += S[left+i * (1 << (h-1)) - (1 << (h-1)) - 1].lte;
                S[left+i * (1 << (h-1)) -1].gt += S[left+i * (1 << (h-1)) - (1 << (h-1)) - 1].gt;
            }
        }
    }

}

/* An efficient insertion sort used for small subproblems */
void sequential_insertionSort(long *array, int left, int right) {

    int i, j, val;
    for(i = left; i <= right; i++) {
        val = array[i];
        j = i - 1;
        while( j >= 0 && array[j] > val) {
            array[j+1] = array[j];
            j--;
        }
        array[j+1] = val;
    }
}

/* Parallel partition that makes use of prefix sums. */
int partition(long *array, int left, int right){

    // Compute n, k for use by prefix sum function.
    int n = (right - left + 1),
        k = (int) log2(n),
        i;

    // Use right as pivot since we already randomized the array.
    long pivot = array[right];

    // Set flags in comparison flag arrays
    cilk_for (i = left; i <= right; i++) {

        // Go ahead and copy array[i] for swap.
        copyArray[i] = array[i];

        // Strictly less than pivot
        if (array[i] < pivot) {
            flags[i].lte = 1;
            flags[i].gt = 0;
        }
        // Strictly greater than pivot
        else if(array[i]>pivot) {
            flags[i].lte = 0;
            flags[i].gt = 1;
        }
        // Equal to pivot: shove into less than range.
        else{
            if(i==right){
                flags[i].lte=0;
            }
            else {
                flags[i].lte=1;
            }
            flags[i].gt=0;
        }
    }

    // Compute the prefix sum of the flag array (in 2 dimensions).
    parallel_prefix_sum(flags, left, n,k);

    // The pivot belongs in this location.
    int pivotIndex = left+flags[right].lte;
    array[pivotIndex] = pivot;

    // Now use these mappings to swap in parallel
    // Note, we don't look at i = right, since it is the pivot
    cilk_for (i = left; i < right; i++){
        // Value was LTE pivot
        if(copyArray[i]<=pivot){
            array[left+flags[i].lte-1] = copyArray[i];
        }
        // Value was GT pivot
        else if (copyArray[i]>pivot){
            array[pivotIndex+flags[i].gt] = copyArray[i];
        }
    }

    // Return the pivot index that ensures the invariant of quicksort!
    return pivotIndex;
}

/* Recursive helper that manages degradation of quicksort */
void quicksort_recursive(long *array,int left,int right){

    // Base case.
    if(left >= right) {
        return;
    }

    //First see if we've degraded to serial insertion sort.
    if (right - left + 1 < SERIAL_INSERTION_NSIZE) {
        return sequential_insertionSort(array, left, right);
    }
    else {
        int splitPoint;

        // Degrade to sequential partition?
        if (right - left + 1 <= SERIAL_PARTITION_N_FACTOR * problem_size) {
            splitPoint = sequential_partition(array,left, right);
        } else{
            splitPoint = partition(array,left, right);
        }

        // Spawn a new task for one sub problem, compute another, implicit sync.
        cilk_spawn quicksort_recursive(array,left,splitPoint-1);
        quicksort_recursive(array,splitPoint+1,right);
    }
}

/* Top-level function for quicksort. Prepares helper arrays and warms caches
 * before triggering recursion. 
 */
void quicksort(long *array, int size) {
    // Allocate helper arrays.
    copyArray = (long *) malloc (sizeof(long) * size);
    flags = (lte_gt *) malloc (sizeof(lte_gt) * size);

    int i;
    // Warm the caches via initialization
    cilk_for(i = 0 ; i < size ; i++) {
        copyArray[i] = 0;
        flags[i].lte = 0;
        flags[i].gt = 0;
    }

    // Recurse.
    quicksort_recursive(array, 0, size-1);
}

int main(int argc, char **argv) {

    // Time keeping related variables
    double start, stop, time_elapsed;

    // Validate run-time arguments
    if(argc < 2){
        dbg_printf("Incorrect number of arguments, expected 1 argument\n");
        return EXIT_FAILURE;
    }
    // Optional N_workers command line argument
    else if (argc == 3) {
        workers = (atoi(argv[2]) < __cilkrts_get_nworkers()) ? atoi(argv[2]): __cilkrts_get_nworkers();
    } else { // Unspecified ==> use default # available
        workers =  __cilkrts_get_nworkers(); 
    }

    // Set number of workers.
     __cilkrts_set_param("nworkers", argv[2]);
    dbg_printf("Using %d available workers.\n", workers);

    // Set up input array.
    problem_size = atoi(argv[1]);
    long* array = (long *) malloc(problem_size * sizeof(long));

    // Initialize array with uniformly random values. Do in cilk_for to warm cache.
    srand(time(NULL));
    int i;
    cilk_for(i = 0; i < problem_size; i++){
        array[i] = random_int(0, problem_size);
    }

    // If we are in PRINTMODE, print the array before sorting.
    #ifdef PRINTMODE
    printf("Unsorted Array\n");
    printArray(array,0, problem_size-1);
    #endif

    // Start the clock, run quicksort, stop the clock.
    start = wctime();
    quicksort(array, problem_size);
    stop = wctime();
    time_elapsed = (double) (stop - start);

    // If we are in PRINTMODE, print the array after sorting.
    #ifdef PRINTMODE
    printf("Sorted Array\n");
    printArray(array,0, problem_size-1);
    #endif

    // Always print a CSV-ready row for data summary:
    // (algorithm type, problem size, number of workers, time elapsed (sec))
    printf("parallel, %d, %d, %f\n", problem_size, workers, time_elapsed);
    return 0;
}