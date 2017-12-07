# final633
final 633 project

# Description
Our basic goal was to take an efficient sequential quicksort and parallelize both the recursion and partitioning steps in an effort to reduce the runtime of the algorithm. In order to accomplish this goal we used cilk with the icc compiler. While our parallel recursive step is quite simple (though this is by far the biggest contributor to the runtime reduction), our partitioning step runs in expected O(lgn) step complexity through the use of a scan based algorithm. We also added a number of features that look to prevent degenerate cases from sinking our runtime.
