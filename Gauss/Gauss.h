#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cstdlib>
#include<time.h>
#define EPS 10e-20
double func(int i,int j);
int ReadfromFile(char *fp,double *A,int N);
int max_line(double *mas,int N,int y);
void change(double *A,int k,int j,int N,int *mas);
void Gauss_method_straight(double *A,double *E,int N,int k);
void Gauss_method_reverse(double *A,double *E,int N,int k);
int Gauss(double *A,double *E,int N,int *mas);
double nevyazka(double *A,double *B,int N);
void print(double *A,int N);
void print_part(double *A,int N,int K);
