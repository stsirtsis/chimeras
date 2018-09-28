#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

int main(){

  FILE *file1;
  char filename[80];
  double pi=3.14159265359 ;
  int seed=823093819;   //seed1
  //int seed= 395849566; //seed2 πλέγμα χιμαιρών
  //int seed= 95679999;
  // int seed= 65243455;
  // int seed= 97243562;
  // int seed= 943746563;

  srand(seed);

  /*******Parameter Declarations*******/
  int N=100;  //Grid dimension
  double dt=0.01; //0.001
  int totalTime=10000;  //Simulation time
  int it=0;
  int totalIter=totalTime/dt; //Total iterations
  int R=22; //Square radius
  double sigma=0.7; //Coupling strength
  double sumCoeff=sigma/((2.0*R+1.0)*(2.0*R+1.0)-1.0);  //Potential sum coefficient
  double mi=1.0; //Integrator floor
  double uth=0.98;

  double Ts=log(mi/(mi-uth));
  double refracTime=0.22*Ts;  //Refractory period time
  int maxRefracIter=round(refracTime/dt); //Refractory period iterations
  int i,j,k,l;
  int iLeftCorner, jLeftCorner;

  double u[N][N];
  double unext[N][N];

  double sumVar=0.0;  //Potential sum coefficient

  double currRefracIter[N][N];  //Current iterations already in refractory period
  int currThresCrossings[N][N];
  int maxMPVIter=20000;
  int minMPVIter=50000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.
  double w;
  double t=0.0;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      unext[i][j]=0.0;
      currRefracIter[i][j]=0.0;
      currThresCrossings[i][j]=0;
    }
  }

  file1=fopen("initial.csv","r");
  char line[2048];
  i=0;
  while(fgets(line, 2048, file1)){
    for(j=1;j<=N;j++){
      char* tmp = strdup(line);
      u[i][j-1]=atof(getfield(tmp,j));
      //printf("Field %d,%d is %s\n", i, j-1, getfield(tmp,j));
      free(tmp);
    }
    i++;
  }
  fclose(file1);

  /*******Simulation*******/
  while (it<totalIter){

    if (it%1000==0) printf("Iteration %d of %d\n", it, totalIter);
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

        if(unext[i][j]>=uth){   //Threshold crossed
          unext[i][j]=0.0;
          if (it>=minMPVIter){
            currThresCrossings[i][j]++;
          }
        }

        u[i][j]=unext[i][j];
      }//edw kleinei h for j.
    }//edw kleinei h for i.

    if (it>=minMPVIter){
      if ((it-minMPVIter)%maxMPVIter==0){
        sprintf(filename, "Results_MPV_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",sigma,R,t);
        file1=fopen(filename,"w");
        for(i=0;i<N;i++){
          for(j=0;j<N;j++){
            w=2.0*pi*currThresCrossings[i][j]/((it-minMPVIter)*dt);
            //currThresCrossings[i][j]=0;
            fprintf(file1,"%lf,",w);
          }
          fprintf(file1,"\n");
        }
        fclose(file1);
      }
    }
    if(it%1000==0){
      sprintf(filename, "Results_POT_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",sigma,R,t);
      file1=fopen(filename,"w");
      for(i=0;i<N;i++){
        for(j=0;j<N;j++){
          fprintf(file1, "%lf,",u[i][j]);
        }
        fprintf(file1,"\n");
      }
      fclose(file1);
    }
    t+=dt;
    it++;
  } //edw kleinei h while.

  return(0);
}
