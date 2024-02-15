#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <math.h>
// to run 
// mpicc -o task3.out task3.c-lm
// mpiexec -np 3 ./task3.out
int rank,key,cnt;
int range;
int siz;
MPI_Status status;
float arr[100000+5],final_sum=0;
float variance=0.0,final_sqrtsum=0,standardDeviation=0.0;
float mean=0.0;
int n,range;
float summation(){
    float sum=0;
    int low=(rank-1)*range;
    int high=low+range;
    if(rank==siz-1){
        high=n;
    }
    for (int i = low; i <high ; ++i) {
        sum+=arr[i];
    }
    return sum;
}

float squaredsum(){
    float sqrtsum=0;
    int low=(rank-1)*range;
    int high=low+range;
    if(rank==siz-1){
        high=n;
    }

    for (int i = low; i <high ; ++i) {
        float diff=(arr[i]-mean);
        sqrtsum+=(diff*diff);
    }
    return sqrtsum;
}
void mpi(){
    if(siz==1){
        printf("enter the sizeof array \n");
        scanf("%d",&n);

        for (int i = 0; i < n; ++i) {
            arr[i]=i+1;
        }
        for (int i = 0; i <n ; ++i) {
            final_sum+=arr[i];
        }
        mean=(double)(final_sum*1.0/n*1.0);
        printf("mean= %f\n",mean);
        for (int i = 0; i <n ; ++i) {
            float diff=(arr[i]-mean);
            final_sqrtsum+=(diff*diff);
        }
        variance=(double)((final_sqrtsum*1.0)/(n*1.0));
        printf("variance= %f\n",variance);
        standardDeviation=(float) sqrt((double)variance);
        printf("standard deviation= %f\n",standardDeviation);
    }else{
        if(!rank){
            printf("enter the sizeof array \n");
            scanf("%d",&n);

            for (int i = 0; i < n; ++i) {
                arr[i]=i+1;
            }
            range=(n/(siz-1));
        }
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&arr,n,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(&range,1,MPI_INT,0,MPI_COMM_WORLD);

        float sum=0;
        if(rank){
            sum=summation();
        }
        MPI_Reduce(&sum,&final_sum,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
        if(!rank){
            mean=(double)((final_sum*1.0)/(n*1.0));
            printf("mean= %f\n",mean);
        }
        MPI_Bcast(&mean,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        float sqrtsum=0;
        if(rank){
            sqrtsum=squaredsum();
        }
        MPI_Reduce(&sqrtsum,&final_sqrtsum,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
        if(!rank){
            variance=(double)((final_sqrtsum*1.0)/(n*1.0));
            printf("variance= %f\n",variance);
            standardDeviation=(float) sqrt((double)variance);
            printf("standard deviation= %f\n",standardDeviation);
        }
    }

}
int main(){
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &siz);
     double time=0.0;
     MPI_Barrier(MPI_COMM_WORLD);
     time-= MPI_Wtime();
     mpi();
     time+=MPI_Wtime();
     if(!rank){
     	 printf("total time= %f\n",time/1000);
     }
    MPI_Finalize();
    return 0;
}
