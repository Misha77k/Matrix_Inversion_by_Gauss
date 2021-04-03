#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<errno.h>
#include<pthread.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<sys/sysctl.h>
#define EPS 10e-50
void synchronize (int total_threads);
void synchronize (int total_threads)
{
    static pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in=PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out=PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    pthread_mutex_lock (&mutex);
    threads_in++;
    if (threads_in >= total_threads)
    {
        threads_out = 0;
        pthread_cond_broadcast (&condvar_in);
    }
    else
    {
        while (threads_in < total_threads)
        {
            pthread_cond_wait (&condvar_in, &mutex);
        }
    }
    threads_out++;
    if (threads_out >= total_threads)
    {
        threads_in = 0;
        pthread_cond_broadcast (&condvar_out);
    }
    else
    {
        while (threads_out < total_threads)
        {
            pthread_cond_wait (&condvar_out, &mutex);
        }
    }
    pthread_mutex_unlock (&mutex);
}
double norma(double *matr,int N);
double norma(double *matr,int N){
    double max=0.0;
    double cur;
    for(int i=0;i<N;i++){
        cur = 0.0;
        for(int j=0;j<N;j++)
            cur+=fabs(matr[i*N+j]);
        if(cur>max)
            max=cur;
    }
    return max;
}
double nevyazka(double *A,double *B,int N);
double nevyazka(double *A,double *B,int N){
    double max=0.0;
    double cur;
    double p;
    for(int i=0;i<N;i++){
        cur=0.0;
        for(int j=0;j<N;j++){
            p=0.0;
            for(int k=0;k<N;k++)
                p=p+A[i*N+k]*B[k*N+j];
            if(j==i)
                p=p-1.0;
            cur=cur+fabs(p);
        }
        if(cur>max)
            max=cur;
    }
    return max;
}
double get_full_time();
double get_full_time(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec+(t.tv_usec)/1000000.0);
}
int max_line(double *mas,int N,int y, int l);
int max_line(double *mas,int N,int y, int l){
    double max=0.0;
    int k=0;
    for(int i=y*N + l;i<(y+1)*N;i++){
        if((fabs(mas[i])>fabs(max))){
            max=mas[i];
            k=i;
        }
    }
    if(fabs(max)<EPS){
        return -1;
    }
    return (k-y*N);
}
int maxLineInThreads(double max,int k);
int maxLineInThreads(double max,int k){
    static double max_in_row=0.0;
    static int k_in_row=0;
    static int count=0;
    static pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&mutex);
    if(count==0){
        max_in_row=fabs(max);
        k_in_row=k;
        count++;
    }else{
        if(fabs(max_in_row)<fabs(max)){
            max_in_row=max;
            k_in_row=k;
        }
    }
    if(fabs(max_in_row<EPS)){
        pthread_mutex_unlock(&mutex);
        return -1;
    }
    pthread_mutex_unlock(&mutex);
    return k_in_row;
    
}
int maxLine(double *mas,int N,int y,int l,int thread_num,int total_threads);
int maxLine(double *mas,int N,int y,int l,int thread_num,int total_threads){
    double max=0.0;
    int k=0;
    int current=thread_num;
    while(current<l){
        current=current+total_threads;
    }
    while(current<N){
        if(fabs(mas[y*N+current])>fabs(max)){
            max=mas[y*N+current];
            k=current;
        }
        current=current+total_threads;
    }
    k=maxLineInThreads(max,k);
    synchronize(total_threads);
    return k;
}
void change(double *A,int k,int j,int N,int *mas);
void change(double *A,int k,int j,int N,int *mas){
    double c;
    int d;
    for(int i=0;i<N;i++){
        c=A[i*N+k];
        A[i*N+k]=A[i*N+j];
        A[i*N+j]=c;
    }
    d=mas[k];
    mas[k]=mas[j];
    mas[j]=d;
}
void swapColumns(double *A,int k,int j,int N,int *mas,int thread_num,int total_threads);
void swapColumns(double *A,int k,int j,int N,int *mas,int thread_num,int total_threads){
    int current=thread_num;
    double c;
    int d;
    while(current<N){
        c=A[current*N+j];
        A[current*N+j]=A[current*N+k];
        A[current*N+k]=c;
        current=current+total_threads;
    }
    synchronize(total_threads);
    if(thread_num==0){
        d=mas[k];
        mas[k]=mas[j];
        mas[j]=d;
    }
}
void Gauss_method_straight(double *A,double *E,int N,int k,int thread_num,int total_threads);
void Gauss_method_straight(double *A,double *E,int N,int k,int thread_num,int total_threads){
    
    int current = thread_num;
    while(current < N){
        if(k < current){
            A[k*N + current] /= A[k*N + k];
        }
        E[k*N + current] /= A[k*N + k];
        current += total_threads;
    }
    synchronize(total_threads);
    for(int i = k+1; i < N; i++){
        current = thread_num;
        while (current < N) {
            if(current > k){
                A[i*N + current] -= A[k*N + current]*A[i*N + k];
            }
            E[i*N + current] -= E[k*N + current]*A[i*N + k];
            current += total_threads;
        }
    }
    
}


void straight(double *A,double *E,int N,int k,int thread_num,int total_threads);
void straight(double *A,double *E,int N,int k,int thread_num,int total_threads){
    int current=thread_num;
    double c;
    while(current<N){
        
        if(current>k){
            A[k*N+current]=A[k*N+current]/A[k*N+k];
        }
        E[k*N+current]=E[k*N+current]/A[k*N+k];
        current=current+total_threads;
    }
    synchronize(total_threads);
    current=thread_num;
    while(current<k)
        current=current+total_threads;
    if(current==k){
        current=current+total_threads;
    }
    while(current<N){
        c=A[current*N+k];
        for(int i=0;i<N;i++){
            if(i>k){
                A[current*N+i]=A[current*N+i]-c*A[k*N+i];
            }
            E[current*N+i]=E[current*N+i]-c*E[k*N+i];
        }
        current=current+total_threads;
    }
}
void Gauss_method_reverse(double *A,double *E,int N,int k,int thread_num,int total_threads);
void Gauss_method_reverse(double *A,double *E,int N,int k,int thread_num,int total_threads){
    int current=0;
    for(int i=k;i>0;i--){
        current=thread_num;
        while (current<N) {
            E[(i-1)*N + current]=E[(i-1)*N + current]-A[(i-1)*N + k]*E[k*N + current];
            current=current+total_threads;
        }
    }
}

void reverse(double *A,double *E,int N,int k,int thread_num,int total_threads);
void reverse(double *A,double *E,int N,int k,int thread_num,int total_threads){
    int current=thread_num;
    while(current<k)
        current=current+total_threads;
    if(current>=k){
        current=current-total_threads;
    }
    while(current>=0){
        for(int j=0;j<N;j++){
            E[current*N+j]=E[current*N+j]-A[current*N+k]*E[k*N+j];
        }
        current=current-total_threads;
    }
}
void swap_matr(double *A,int i,int j,int N);
void swap_matr(double *A,int i,int j,int N){
    double c;
    for(int k=0;k<N;k++){
        c=A[i*N+k];
        A[i*N+k]=A[j*N+k];
        A[j*N+k]=c;
    }
}
void swap_arr(int *mas,int j,int i);
void swap_arr(int *mas,int j,int i){
    int s;
    s=mas[i];
    mas[i]=mas[j];
    mas[j]=s;
}
int Gauss_method(double *A,double *E,int *mas,int N,int thread_num,int total_threads,int &error);


int Gauss_method(double *A,double *E,int *mas,int N,int thread_num,int total_threads,int &error){
    int k;
    static int q = 0;
    for(int j=0;j<N;j++){
        if(thread_num == 0){
            k=max_line(A,N,j,j);
            if(k == -1){
                error = -1;
                q = -1;
                return 0;
            }
            change(A,k,j,N,mas);
        }
        synchronize(total_threads);
        straight(A, E, N, j, thread_num, total_threads);
        synchronize(total_threads);
    }
    for(int i=N-1;i>=0;i--){
        Gauss_method_reverse(A,E,N,i, thread_num, total_threads);
    }
    synchronize(total_threads);
    if(thread_num == 0){
        for(int i=0;i<N;){
            if(mas[i]!=i){
                swap_matr(E,mas[i],i,N);
                swap_arr(mas,mas[i],i);
            }
            else i++;
        }
    }
    return 0;
}
struct Args{
    double *A;
    double *E;
    int *mas;
    int N;
    int thread_num;
    int total_threads;
    double time;
    int error;
};
void *Gauss_method_threaded(void *pa);
void * Gauss_method_threaded (void *pa)
{
    Args *pargs = (Args*)pa;
    int k=0;
    pargs->time = get_full_time ();
    k=Gauss_method(pargs->A,pargs->E,pargs->mas,pargs->N,pargs->thread_num,pargs->total_threads,pargs->error);
    pargs->time = get_full_time () - pargs->time;
    return 0;
}
int main(){
    double *A;
    double *E;
    double *B;
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
    B=(double *)malloc(N*N*sizeof(double));
    if(B==NULL)
        return -1;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            B[i*N+j]=A[i*N+j];
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
        args[i].error = error;
    }
    for (int i = 0; i < nthreads; i++)
    {
        if (pthread_create (threads + i, 0,Gauss_method_threaded,args + i)){
            printf ("cannot create thread #%d!\n",i);
            return 10;
        }
    }
    //Gauss_method_threaded(args + 0);
    
    for (int i = 0; i < nthreads; i++)
    {
        if (pthread_join (threads[i], 0))
            fprintf (stderr, "cannot wait thread #%d!\n", i);
    }
    
    printf("time is %lf\n", args[0].time);
    if(args[0].error == -1){
        printf("error 1\n");
        return -1;
    }
    double nev;
    nev=nevyazka(B,E,N);
    printf("Невязка равна %.17le\n",nev);
    free (threads);
    free (args);
    free (A);
    free (E);
    free (mas);
    return 0;
}
