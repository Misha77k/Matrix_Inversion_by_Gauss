#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<errno.h>
#include<pthread.h>
#include<sys/time.h>
#include"Parallel_Gauss.h"
#define EPS 10e-50
int main(){
    double *A;
    double *E;
    int *mas;
    int N;
    Args *args;
    pthread_t *threads;
    int nthreads;
    printf("Введи число потоков\n");
    scanf("%d",&nthreads);
    if (!(threads = (pthread_t*)
          malloc (nthreads * sizeof (pthread_t))))
    {
        printf ("Not enough memory!\n");
        return -1;
    }
    if (!(args = (Args*) malloc (nthreads * sizeof (Args))))
    {
        printf ("Not enough memory!\n");
        return -1;
    }
    printf("Введи размер матрицы\n");
    scanf("%d",&N);
    A=(double *)malloc(N*N*sizeof(double));
    if(A==NULL)
        return -1;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            A[i*N+j]=abs(i - j);
    E=(double *)malloc(N*N*sizeof(double));
    if(E==NULL)
        return -1;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            if(i==j)
                E[i*N+j]=1;
            else E[i*N+j]=0;
        }
    mas=(int *)malloc(N*sizeof(int));
    if(mas==NULL)
        return -1;
    for(int i=0;i<N;i++)
        mas[i]=i;
    int error = 0;
    for (int i = 0; i < nthreads; i++)
    {
        args[i].A= A;
        args[i].E = E;
        args[i].mas = mas;
        args[i].N=N;
        args[i].thread_num = i;
        args[i].total_threads = nthreads;
        args[i].residual = 0.0;
        args[i].error = error;
    }
    //double t_full = get_full_time ();
    /* Запускаем задачи */
    for (int i = 1; i < nthreads; i++)
    {
        if (pthread_create (threads + i, 0,Gauss_method_threaded,args + i)){
            printf ("cannot create thread #%d!\n",i);
            return 10;
        }
    }
    Gauss_method_threaded(args + 0); // для главного потока
    /* Ожидаем окончания задач */
    
    printf("time is %lf\n", args[0].time);
    if(args[0].error == -1){
        printf("error\n");
        return -1;
    }
    printf("residual is %.17le\n",args[0].residual);
    /* очень быстрый компьютер... */
    /*for(int i=0;i<N;i++){
     for(int j=0;j<N;j++)
     printf("%lf ",E[i*N+j]);
     printf("\n");
     }*/
    free (threads);
    free (args);
    free (A);
    free (E);
    free (mas);
    return 0;
}

