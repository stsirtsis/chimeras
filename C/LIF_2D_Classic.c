#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Returns a real number in the interval [0,1]
double myrand(){
  return((double)rand()/(double)RAND_MAX);
}

int main(){

  FILE *file1;
  file1=fopen("Potential.dat","w");
  fprintf(file1, "u[x][y]\n");
  fclose(file1);

  FILE *file2;
  char filename[60];
  double pi=3.14159265359 ;
  //int seed=823093819;   //seed1
  // int seed= 395849566; //seed2 πλέγμα χιμαιρών
  int seed= 95679999;
  // int seed= 65243455;
  // int seed= 97243562;
  // int seed= 65345981;

  srand(seed);

  /*******Parameter Declarations*******/
  int N=100;  //Grid dimension
  double dt=0.001;
  int totalTime=10000;  //Simulation time
  int it=0;
  int totalIter=totalTime/dt; //Total iterations
  int R=25; //Square radius
  double sigma=0.2; //Coupling strength
  double sumCoeff=sigma/((2.0*R+1.0)*(2.0*R+1.0)-1.0);  //Potential sum coefficient
  double mi=1.0; //Integrator floor

  double refracTime;  //Refractory period time
  int maxRefracIter=500; //Refractory period iterations
  refracTime=maxRefracIter*dt;

  //Neuron coordinates for diagram purposes
  int x=5;
  int y=78;

  int i,ki,kki,ll;
  int j,kj,kkj,pp;

  double uth=0.98;
  double u[N][N];
  double unext[N][N];

  double sumVar=0.0;  //Potential sum coefficient

  double currRefracIter[N][N];  //Current iterations already in refractory period
  int currThresCrossings[N][N];
  int currMPVIter=0;
  int maxMPVIter=500;
  int minMPVIter=10;//2000000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.
  double w;
  double t=0.0;
  //arxikopoihsh tou pinaka dunamikwn, ekf , currRefracIter.
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      unext[i][j]=0.0;
      currRefracIter[i][j]=0.0;
      currThresCrossings[i][j]=0;
      do{
        u[i][j]=myrand();
      }while(u[i][j]>=uth);
    }
  }


  /*******Simulation*******/

  while (it<totalIter){

    printf("Iteration %d of %d\n", currMPVIter, maxMPVIter);
    for (i=0; i<N; i++){
      for (j=0; j<N; j++){

        /*******Refractory Period*******/
        if (u[i][j]==0 && currRefracIter[i][j]<maxRefracIter){
          currRefracIter[i][j]++;
          break;
        }
        else{
          currRefracIter[i][j]=0;
        }

        /*******Sum Calculation*******FIX!!!!!!!!!!!!*/
        sumVar=0.0;
        for  (ki=0; ki<N; ki++){
          for ( kj=0; kj<N; kj++){
            sumVar=sumVar+(u[i][j]-u[ki][kj]);
          }
        }

        unext[i][j]=u[i][j]+dt*(mi-u[i][j]+sumCoeff*sumVar) ;

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
      if (currMPVIter==maxMPVIter){
        for(i=0;i<N;i++){
          for(j=0;j<N;j++){
            w=2.0*pi*currThresCrossings[i][j]/(maxMPVIter*dt);
            currThresCrossings[i][j]=0;
            if (i==0 && j==0){
              //snprintf (filename, "LIF_2D_Classic_sigma_time_%d.dat",it);
              sprintf(filename, "LIF_2D_Classic_sigma.dat");
              file2=fopen(filename,"w");
            }
            //fprintf(file2,"%d, %d, %lf\n",i,j,w);
            fprintf(file2,"%lf,",w);
            if(j==N-1) fprintf(file2,"\n");
            if (i==N-1 && j==N-1){
              fclose(file2);
            }
          }
        }
      }
      if (currMPVIter==maxMPVIter) currMPVIter=0;
      else currMPVIter++;
    }
    t+=dt;
    it++;
    printf("%lf\n", u[x][y]);
    file1=fopen("Potential.dat","a");
    fprintf(file1,"%lf\n",u[x][y]);
    fclose(file1);
  } //edw kleinei h while.

  //printf(" ****** Telos!! ******\n");
  return(0);
}
