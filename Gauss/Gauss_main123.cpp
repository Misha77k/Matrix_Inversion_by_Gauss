#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cstdlib>
#include<time.h>
#include"Gauss.h"
int main(int argc,char** argv){
    double *A,*B;
    double *E;
    int N;
    int M;
    clock_t t;
 //   printf("Введи размерность матрицы\n");
    //scanf("%d",&N);
    if(argc!=2 && argc!=3){
        printf("error\n");
        return -1;
    }
    else{
        N=atoi(argv[1]);
    }
    if(N<=0){
        printf("Неверная размерность\n");
        return -1;
    }
    A=(double*)malloc(N*N*sizeof(double));
    if(argc==2){
       // N=atoi(argv[1]);
        //A=(double*)malloc(N*N*sizeof(double));
        /*for(int i=0;i<N*N;i++)
         if(scanf("%lf",&A[i])!=1)
         return -1;*/
        for(int i=0;i<N;i++)
            for(int j=0;j<N;j++)
                A[i*N+j]=func(i,j);
    }
    if(argc==3){
       int k = ReadfromFile(argv[2],A,N);
       if(k==-1){
           free(A);
           return -1;
       }
    }
   // A=(double*)malloc(N*N*sizeof(double));
    /*for(int i=0;i<N*N;i++)
     if(scanf("%lf",&A[i])!=1)
     return -1;*/
    /*for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            A[i*N+j]=func(i,j);*/
    E=(double*)malloc(N*N*sizeof(double));
    for(int i=0;i<N*N;i++)
        if (i % (N+1)==0)
            E[i]=1;
        else E[i]=0;
    int *mas;
    mas=(int*)malloc(N*sizeof(int));
    for(int i=0;i<N;i++)
        mas[i]=i;
  //  double d;
    B=(double*)malloc(N*N*sizeof(double));
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            B[i*N+j]=A[i*N+j];
    t=clock();
    int r;
    r=Gauss(A,E, N,mas);
    if(r==-1){
        free(mas);
	free(A);
	free(E);
	free(B);
        printf("матрица не имеет обратной\n");
        t=clock()-t;
        printf("Время=%lf\n",((double)t)/CLOCKS_PER_SEC);
        return -1;
    }
   /* for(int i=0;i<N;i++)
        if(mas[i]!=i){
            for(int j=0;j<N;j++){
                d=E[i*N+j];
                E[i*N+j]=E[mas[i]*N+j];
                E[mas[i]*N+j]=d;
            }
            mas[mas[i]]=mas[i];
        }*/
    t=clock()-t;
    printf("Время=%lf\n",((double)t)/CLOCKS_PER_SEC);
   // print(E,M);
    double nev;
    nev=nevyazka(B, E, N);
    printf("невязка=%le\n",nev);
    printf("Введите размерность главного минора\n");
    scanf("%d",&M);
    if(M>=N)
        print(E,N);
    else print_part(E,N,M);
    free(mas);
    free(A);
    free(B);
    free(E);
    return 0;
}

