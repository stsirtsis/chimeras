#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

const char* getfield(char* line, int num){
  const char* tok;
  for (tok = strtok(line, ",");
  tok && *tok;
  tok = strtok(NULL, ",\n"))
  {
    if (!--num)
    return tok;
  }
  return NULL;
}

__global__ void calcPot(double *u, double *unext, int *currRefracIter, double *currTime, int *d_it, double *w, double *lastTime, double *d_refracTime,
int *d_maxRefracIter, int *d_N, int *d_R, double *d_uth, int *d_minMPVIter, double *d_dt, double *d_mi, double *d_sumCoeff){

  double pi=3.14159265359;
  int it = *d_it;
  double refracTime = *d_refracTime;
  int maxRefracIter = *d_maxRefracIter;
  int N = *d_N;
  int R = *d_R;
  double uth = *d_uth;
  int minMPVIter = *d_minMPVIter;
  double dt = *d_dt;
  double mi = *d_mi;
  double sumCoeff = *d_sumCoeff;
  int myId = blockDim.x * blockIdx.x + threadIdx.x;

  /*******Refractory Period*******/
  if (*(u+myId)==0 && *(currRefracIter+myId)<maxRefracIter){
    (*(currRefracIter+myId))++;
    return;
  }
  else{
    *(currRefracIter+myId)=0;
  }

  /*******Sum Calculation*******/
  double sumVar=0.0;
  int k,l;
  int iLeftCorner=N+myId/N-R;
  int jLeftCorner=N+myId%N-R;
  for  (k=iLeftCorner; k<iLeftCorner+2*R+1; k++){
    for (l=jLeftCorner; l<jLeftCorner+2*R+1; l++){
      sumVar+=*(u+myId)-*(u+(k%N)*N+l%N);
    }
  }

  *(unext+myId)=*(u+myId)+dt*(mi-*(u+myId)+sumCoeff*sumVar);
  *(currTime+myId)+=dt;

  if(*(unext+myId)>=uth){   //Threshold crossed
    *(unext+myId)=0.0;
    if (it>=minMPVIter){
      *(w+myId)=((*(w+myId))*(*(lastTime+myId))+2*pi)/((*(lastTime+myId))+(*(currTime+myId))+refracTime);
      *(lastTime+myId)+=(*(currTime+myId))+refracTime;
    }
    *(currTime+myId)=0.0;
  }

  return;
}

int main(int argc, char** argv){

  FILE *file1;
  char filename[100];

  /*******Parameter Declarations*******/
  int N=100;  //Grid dimension
  double dt=0.001; //0.001
  int totalTime=10000;  //Simulation time
  int it=0;
  int totalIter=totalTime/dt; //Total iterations
  int R=22; //Square radius
  double sigma=0.7; //Coupling strength
  double sumCoeff=sigma/((2*R+1)*(2*R+1)-1);  //Potential sum coefficient
  double mi=1.0; //Integrator floor
  double uth=0.98;

  double Ts=log(mi/(mi-uth));
  double refracTime=0.22*Ts;  //Refractory period time
  int maxRefracIter=(int)ceil(refracTime/dt); //Refractory period iterations
  int i,j;

  double u[N*N];
  double unext[N*N];

  int currRefracIter[N*N];  //Current iterations already in refractory period
  int maxMPVIter=30000;
  int minMPVIter=2000000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.

  double currTime[N*N];
  double lastTime[N*N];
  double w[N*N];
  double t=0.0;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      (*(unext+i*N+j))=0.0;
      (*(currTime+i*N+j))=0.0;
      (*(lastTime+i*N+j))=0.0;
      (*(currRefracIter+i*N+j))=0.0;
    }
  }
  file1=fopen(argv[1],"r"); //argv[1]
  char line[2048];
  i=0;
  while(fgets(line, 2048, file1)){
    for(j=1;j<=N;j++){
      char* tmp = strdup(line);
      (*(u+N*i+j-1))=atof(getfield(tmp,j));
      free(tmp);
    }
    i++;
  }
  fclose(file1);

  double *d_u, *d_unext, *d_currTime, *d_w, *d_lastTime, *d_refracTime, *d_uth, *d_dt, *d_mi, *d_sumCoeff;
  int *d_currRefracIter, *d_it, *d_maxRefracIter, *d_N, *d_R, *d_minMPVIter;

  cudaMalloc(&d_u, N*N*sizeof(double));
  cudaMalloc(&d_unext, N*N*sizeof(double));
  cudaMalloc(&d_currRefracIter, N*N*sizeof(int));
  cudaMalloc(&d_currTime, N*N*sizeof(double));
  cudaMalloc(&d_it, sizeof(int));
  cudaMalloc(&d_w, N*N*sizeof(double));
  cudaMalloc(&d_lastTime, N*N*sizeof(double));

  cudaMalloc(&d_refracTime, sizeof(double));
  cudaMalloc(&d_maxRefracIter, sizeof(int));
  cudaMalloc(&d_N, sizeof(int));
  cudaMalloc(&d_R, sizeof(int));
  cudaMalloc(&d_uth, sizeof(double));
  cudaMalloc(&d_minMPVIter, sizeof(int));
  cudaMalloc(&d_dt, sizeof(double));
  cudaMalloc(&d_mi, sizeof(double));
  cudaMalloc(&d_sumCoeff, sizeof(double));

  cudaMemcpy(d_refracTime, &refracTime, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_maxRefracIter, &maxRefracIter, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_N, &N, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_R, &R, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_uth, &uth, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_minMPVIter, &minMPVIter, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dt, &dt, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_mi, &mi, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_sumCoeff, &sumCoeff, sizeof(double), cudaMemcpyHostToDevice);

  time_t benchBegin = time(NULL);
  /*******Simulation*******/
  while (it<totalIter){

    if (it%10000==0) printf("Iteration %d of %d\n", it, totalIter);

    cudaMemcpy(d_u, u, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_unext, unext, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_currRefracIter, currRefracIter, N*N*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_currTime, currTime, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_it, &it, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_w, w, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lastTime, lastTime, N*N*sizeof(double), cudaMemcpyHostToDevice);

    //printf("STARTING\n");
    calcPot<<<100,100>>>(d_u, d_unext, d_currRefracIter, d_currTime, d_it, d_w, d_lastTime, d_refracTime, d_maxRefracIter, d_N, d_R, d_uth, d_minMPVIter, d_dt, d_mi, d_sumCoeff);
    cudaThreadSynchronize();

    cudaMemcpy(unext, d_unext, N*N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(currRefracIter, d_currRefracIter, N*N*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(currTime, d_currTime, N*N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(w, d_w, N*N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(lastTime, d_lastTime, N*N*sizeof(double), cudaMemcpyDeviceToHost);

    // for(i=0; i<N; i++){
    //   for(j=0; j<N; j++){
    //     printf("%lf ", *(unext+N*i+j));
    //   }
    //   printf("\n");
    // }
    // printf("FINISHED\n");

    if(it%10000==0){
      sprintf(filename, "ResultsCUDA%s/Results_POT_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",argv[2],sigma,R,t);
      file1=fopen(filename,"w");
      for(i=0;i<N;i++){
        for(j=0;j<N;j++){
          fprintf(file1, "%lf,",*(unext+N*i+j));
        }
        fprintf(file1,"\n");
      }
      fclose(file1);
    }
    if (it>minMPVIter){
      if ((it-minMPVIter)%maxMPVIter==0){
        sprintf(filename, "ResultsCUDA%s/Results_MPV_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",argv[2],sigma,R,t);
        file1=fopen(filename,"w");
        for(i=0;i<N;i++){
          for(j=0;j<N;j++){
            fprintf(file1,"%lf,",*(w+N*i+j));
          }
          fprintf(file1,"\n");
        }
        fclose(file1);
      }
    }
    if (it == 2000000){
      time_t benchEnd = time(NULL);
      sprintf(filename, "ResultsCUDA%s/execTime.dat",argv[2]);
      file1=fopen(filename,"w");
      fprintf(file1,"Execution time for 2000 time units: %ld seconds\n",benchEnd-benchBegin);
      fclose(file1);
    }
    for (i=0; i<N; i++){
      for (j=0; j<N; j++){
        (*(u+N*i+j))=*(unext+N*i+j);
      }
    }
    t+=dt;
    it++;
  } //edw kleinei h while.

  cudaFree(d_u);
  cudaFree(d_unext);
  cudaFree(d_currRefracIter);
  cudaFree(d_currTime);
  cudaFree(d_it);
  cudaFree(d_w);
  cudaFree(d_lastTime);
  cudaFree(d_refracTime);
  cudaFree(d_maxRefracIter);
  cudaFree(d_N);
  cudaFree(d_R);
  cudaFree(d_uth);
  return(0);
}
