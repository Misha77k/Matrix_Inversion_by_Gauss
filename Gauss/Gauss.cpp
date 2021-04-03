#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cstdlib>
#include<time.h>
#include<string.h>
#include"Gauss.h"
#define EPS 10e-20
double func(int i,int j){
    return abs(i-j);
}
int ReadfromFile(char*fp,double *A,int N){
    FILE *f;
    int K=0;
    // int M=0;
    char a;
    f=fopen(fp,"r");
    if(f==NULL){
        printf("файл не открылся\n");
        return -1;
    }
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            if(fscanf(f,"%lf",&A[i*N+j])!=1){
                printf("error\n");
                return -1;
            }
            K++;
        }
    int q=fscanf(f,"%c",&a);
    if (fscanf(f,"%c",&a)==1){
        printf("Несовпадение размерностей\n");
        fclose(f);
        return -1;
    }
    if(K!=N*N){
        printf("Несовпадение размерностей\n");
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}
int max_line(double *mas,int N,int y){
    double max=0.0;
    int k=0;
    for(int i=y*N;i<(y+1)*N;i++){
        // printf("mas[%d]=%lf\n",i,mas[i]);
        if((fabs(mas[i])>fabs(max))){
            //  printf(" %d:",y);
            // printf("   %lf",fabs(mas[i]));
            max=mas[i];
            
            k=i;
        }
    }
    //    printf("max=%.17lf\n",max);
    if(fabs(max)<EPS){
        printf("матрица не имеет обратной\n");
        return -1;
    }
    //  printf("%d, ",k-y*N);
    return (k-y*N);
}
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
void Gauss_method_straight(double *A,double *E,int N,int k){
    double c=A[k*N+k];
    double d;
    for(int i=0;i<N;i++){
        A[k*N+i]=A[k*N+i]/c;
        E[k*N+i]=E[k*N+i]/c;
    }
    for(int i=k+1;i<N;i++){
        d=A[i*N+k];
        for(int j=0;j<N;j++){
            A[i*N+j]=A[i*N+j]-d*A[k*N+j];
            E[i*N+j]=E[i*N+j]-d*E[k*N+j];
        }
    }
}
void Gauss_method_reverse(double *A,double *E,int N,int k){
    double c=A[k*N+k];
    double d;
    for(int i=0;i<=k;i++){
        A[k*N+i]=A[k*N+i]/c;
        E[k*N+i]=E[k*N+i]/c;
    }
    for(int i=k-1;i>=0;i--){
        d=A[i*N+k];
        for(int j=N-1;j>=0;j--){
            A[i*N+j]=A[i*N+j]-d*A[k*N+j];
            E[i*N+j]=E[i*N+j]-d*E[k*N+j];
        }
    }
}
void swap_matr(double *A,int i,int j,int N){
    double c;
    for(int k=0;k<N;k++){
        c=A[i*N+k];
        A[i*N+k]=A[j*N+k];
        A[j*N+k]=c;
    }
}
void swap_arr(int *mas,int j,int i){
    int s;
    s=mas[i];
    mas[i]=mas[j];
    mas[j]=s;
}
int Gauss(double *A,double *E,int N,int *mas){
    int k;
    double d;
    for(int j=0;j<N;j++){
        k=max_line(A,N,j);
        if(k==-1)
            return -1;
        change(A,k,j,N,mas);
        Gauss_method_straight(A,E,N,j);
    }
    for(int i=N-1;i>=0;i--){
        Gauss_method_reverse(A,E,N,i);
    }
    for(int i=0;i<N;){
        if(mas[i]!=i){
            swap_matr(E,mas[i],i,N);
            swap_arr(mas,mas[i],i);
        }
        else i++;
    }
    return 0;
}
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
                p=p-1;
            cur=cur+fabs(p);
        }
        if(cur>max)
            max=cur;
    }
    return max;
}
void print(double *A,int N){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++)
            printf("%lf  ",A[i*N+j]);
        printf("\n");
    }
}
void print_part(double *A,int N,int K){
    if(K<N){
        for(int i=0;i<K;i++){
            for(int j=0;j<K;j++)
                printf("%lf ",A[i*N+j]);
            printf("\n");
        }
    }
    if(K>=N)
        print(A,N);
}
