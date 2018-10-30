#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

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

int main(int argc, char** argv){

  FILE *file1;
  char filename[100];
  double pi=3.14159265359 ;

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
  int i,j,k,l;
  int iLeftCorner, jLeftCorner;

  double u[N][N];
  double unext[N][N];

  double sumVar=0.0;  //Potential sum coefficient

  int currRefracIter[N][N];  //Current iterations already in refractory period
  int maxMPVIter=30000;
  int minMPVIter=2000000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.

  double currTime[N][N];
  double lastTime[N][N];
  double w[N][N];
  double t=0.0;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      unext[i][j]=0.0;
      currTime[i][j]=0.0;
      lastTime[i][j]=0.0;
      currRefracIter[i][j]=0.0;
    }
  }

  file1=fopen(argv[1],"r"); //argv[1]
  char line[2048];
  i=0;
  while(fgets(line, 2048, file1)){
    for(j=1;j<=N;j++){
      char* tmp = strdup(line);
      u[i][j-1]=atof(getfield(tmp,j));
      free(tmp);
    }
    i++;
  }
  fclose(file1);

  time_t benchBegin = time(NULL);
  /*******Simulation*******/
  while (it<totalIter){

    if (it%10000==0) printf("Iteration %d of %d\n", it, totalIter);
    #pragma omp parallel num_threads(5)
    {
      #pragma omp for private(sumVar,iLeftCorner,jLeftCorner,k,l) schedule(dynamic) collapse(2)
      for (i=0; i<N; i++){
        for (j=0; j<N; j++){

          /*******Refractory Period*******/
          if (u[i][j]==0 && currRefracIter[i][j]<maxRefracIter){
            currRefracIter[i][j]++;
            continue;
          }
          else{
            currRefracIter[i][j]=0;
          }

          /*******Sum Calculation*******/
          sumVar=0.0;
          iLeftCorner=N+i-R;
          jLeftCorner=N+j-R;
          for  (k=iLeftCorner; k<iLeftCorner+2*R+1; k++){
            for (l=jLeftCorner; l<jLeftCorner+2*R+1; l++){
              sumVar+=(u[i][j]-u[k%N][l%N]);
            }
          }

          unext[i][j]=u[i][j]+dt*(mi-u[i][j]+sumCoeff*sumVar);
          currTime[i][j]+=dt;

          if(unext[i][j]>=uth){   //Threshold crossed
            unext[i][j]=0.0;
            if (it>=minMPVIter){
              w[i][j]=(w[i][j]*lastTime[i][j]+2*pi)/(lastTime[i][j]+currTime[i][j]+refracTime);
              lastTime[i][j]+=(currTime[i][j]+refracTime);
            }
            currTime[i][j]=0.0;
          }
        }//edw kleinei h for j.
      }//edw kleinei h for i.
    }
    if(it%10000==0){
      sprintf(filename, "ResultsOpenMP%s/Results_POT_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",argv[2],sigma,R,t);
      file1=fopen(filename,"w");
      for(i=0;i<N;i++){
        for(j=0;j<N;j++){
          fprintf(file1, "%lf,",unext[i][j]);
        }
        fprintf(file1,"\n");
      }
      fclose(file1);
    }
    if (it>minMPVIter){
      if ((it-minMPVIter)%maxMPVIter==0){
        sprintf(filename, "ResultsOpenMP%s/Results_MPV_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",argv[2],sigma,R,t);
        file1=fopen(filename,"w");
        for(i=0;i<N;i++){
          for(j=0;j<N;j++){
            fprintf(file1,"%lf,",w[i][j]);
          }
          fprintf(file1,"\n");
        }
        fclose(file1);
      }
    }
    if (it == 2000000){
      time_t benchEnd = time(NULL);
      sprintf(filename, "ResultsOpenMP%s/execTime.dat",argv[2]);
      file1=fopen(filename,"w");
      fprintf(file1,"Execution time for 2000 time units: %ld seconds\n",benchEnd-benchBegin);
      fclose(file1);
    }
    for (i=0; i<N; i++){
      for (j=0; j<N; j++){
        u[i][j]=unext[i][j];
      }
    }
    t+=dt;
    it++;
  } //edw kleinei h while.

  return(0);
}
