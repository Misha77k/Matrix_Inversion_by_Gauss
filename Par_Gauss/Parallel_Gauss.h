#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<errno.h>
#include<pthread.h>
#include<sys/time.h>
#define EPS 10e-50
void synchronize (int total_threads);
double func(int i,int j);
int ReadfromFile(char *fp,double *A,int N);
double norma(double *matr,int N);
double p_nevyazka(double max);
double nevyazka(double *A,double *B,int N, int thread_num, int total_threads);
double get_full_time();
int max_line(double *mas,int N,int y, int l);
void change(double *A,int k,int j,int N,int *mas);
void Gauss_method_straight(double *A,double *E,int N,int k,int thread_num,int total_threads);
void Gauss_method_reverse(double *A,double *E,int N,int k,int thread_num,int total_threads);
void swap_matr(double *A,int i,int j,int N);
void swap_arr(int *mas,int j,int i);
int Gauss_method(double *A,double *E,int *mas,int N,int thread_num,int total_threads,double &residual,int &error);
struct Args{
    double *A;
    double *E;
    int *mas;
    int N;
    int thread_num;
    int total_threads;
    double residual;
    double time;
    int error;
};
void *Gauss_method_threaded(void *pa);
void print(double *A,int N);
void print_part(double *A,int N,int K);
