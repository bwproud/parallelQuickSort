#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#define N 1000
#define CHUNKSIZE 4

int * scan(int * num, int len);
int main(int argc, char *argv[]) {
int i=0;
int read;
int k = 0;
printf("Enter the number k to set the length of the array(len=2^k):\n");
scanf("%d", &read);
k=read;
printf("Ok, enter %d values:\n", (int)pow(2, k));
int len=(int)pow(2, k);
int num[len];
while(i<len){
        scanf("%d", &read);
        num[i++]=read;
}

int * s = scan(num, len);
for(int i = 0; i<len; i++){
	printf("s[%d]=%d\n", i, s[i]);
}
return 0; 
}

int * scan(int * num, int len){
if(len==1){
	return num;
}
int chunk = CHUNKSIZE;
int i=0;
int * s = malloc(len * sizeof(int));
int y[(int)(len/2)];
int z[(int)(len/2)];
for(i=0; i<len; i++){
        s[i]=num[i];
}
 #pragma omp parallel shared(y,num) private(i)
 {
   #pragma omp for schedule(static,chunk)
   for (i=0; i < (int)(len/2); i++){
   	y[i]=num[2*i]+num[(2*i)+1];
   }
}
   
   int * ret=scan(y, len/2);

#pragma omp parallel shared(s,z,num) private(i)
 {
   #pragma omp for schedule(static,chunk)
   for (i=0; i < (int)(len/2); i++){
        z[i]=ret[i];
   }   
   #pragma omp for schedule(static,chunk)
   for (i=0; i < len; i++){
        if(i%2==1){
		s[i]=z[(int)i/2];
	}else if(i==0){
		s[i]=num[i];
	}else{
		s[i]=z[(int)((i-1)/2)]+num[i];
	}
   }

}   /* end of parallel region */
return s;
}
