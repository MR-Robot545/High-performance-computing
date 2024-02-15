#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
// to run
// gcc -fopenmp main.c -o main -lm
// export OMP_NUM_THREADS=1
// ./main

int n, localN, rem;
double arr[131072+20], localS[100000+5], sum=0, localSum=0,localSquaredSum=0, mean=0, variance=0, stdDev=0, totalSum = 0, totalSquaredSum=0;
    int numOfThreads = 1;
   double start_time, end_time, total_time;
int main(int argc, char **argv) {
    
    start_time = omp_get_wtime();
    
    printf("Enter array size:\n");
    scanf("%d", &n);
    printf("Enter array elements:\n");
    for (int i=0; i<n; i++) {
        arr[i]+=(i+1)*1.0;
    }

    #pragma omp parallel
    {
        #pragma omp master
        {
            numOfThreads = omp_get_num_threads();
        }
    }

    localN = n / numOfThreads;
    rem = n % numOfThreads;

    #pragma omp parallel private(localSum) reduction(+:totalSum)
    {
        int tid = omp_get_thread_num();
        #pragma omp critical
        {
            int start = tid * localN;
            int end = start + localN;
            if ((tid == numOfThreads-1) && (rem)) {
                end += rem;
            }
            localSum = 0.0;
            for (int i=start; i<end; i++) {
                localS += arr[i];
            }
            totalSum += localSum;
        }
    }

    mean = totalSum / n;



    totalSquaredSum = 0.0;
    #pragma omp parallel private(localSquaredSum) reduction(+:totalSquaredSum)
    {
        int tid = omp_get_thread_num();
        #pragma omp critical
        {
            int start = tid * localN;
            int end = start + localN;
            if ((tid == numOfThreads-1) && (rem)) {
                end += rem;
            }
            for (int i=start; i<end; i++) {
                localData[i-start] = arr[i];
            }
            localSquaredSum = 0.0;
            for (int i=start; i<end; i++) {
                localSquaredSum += pow(localData[i-start] - mean, 2);
            }
            totalSquaredSum += localSquaredSum;
        }
    }
    variance = totalSquaredSum / n;
    stdDev = sqrt(variance);
    printf("Mean = %.1f\nVariance = %.4f\nStandard deviation = %.4f\n", mean, variance, stdDev);

   
     end_time = omp_get_wtime();
      total_time = end_time - start_time;
      total_time/=1000.0;
    printf("Total time: %f seconds\n", total_time);

    
    total_time = end_time - start_time;

    return 0;
}
