#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

const size_t NUMTHREADS = 4;
int done = 0;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
pthread_cond_t finCond = PTHREAD_COND_INITIALIZER;

typedef struct {
  int *N;
  double *u;
  double *unext;
  int *currRefracIter;
  int *maxRefracIter;
  int *currThresCrossings;
  int *R;
  double *dt;
  double *mi;
  double *sumCoeff;
  double *uth;
  int *minMPVIter;
  int *it;
  int *totalIter;
} threadArgs;

typedef struct {
  threadArgs *ta;
  int id;
} allThreadArgs;

void* threadJob(void *args){//Needs thread id, N, u, currRefracIter, maxRefracIter, R, unext, dt, mi, sumCoeff, uth, minMPVIter

  allThreadArgs *allArgs = args;
  threadArgs *myArgs = allArgs->ta;
  int threadId = allArgs->id;
  int N = *(myArgs->N);
  int maxRefracIter = *(myArgs->maxRefracIter);
  int R = *(myArgs->R);
  double dt = *(myArgs->dt);
  double mi = *(myArgs->mi);
  double sumCoeff = *(myArgs->sumCoeff);
  double uth = *(myArgs->uth);
  int minMPVIter = *(myArgs->minMPVIter);
  int it = *(myArgs->it);
  int totalIter = *(myArgs->totalIter);
  //printf("TEST %d\n", *(myArgs->it));
  int *currThresCrossings = myArgs->currThresCrossings;
  double *u = myArgs->u;
  double *unext = myArgs->unext;
  int *currRefracIter = myArgs->currRefracIter;

  while(it<totalIter){
    //printf("Thread %d beginning iteration %d\n", threadId, it);
    for (int i=threadId*N/((int)NUMTHREADS); i<(threadId+1)*N/((int)NUMTHREADS); i++){
      for (int j=0; j<N; j++){

        /*******Refractory Period*******/
        if ((*(u+N*i+j))==0 && (*(currRefracIter+N*i+j))<maxRefracIter){
          (*(currRefracIter+N*i+j))++;
          continue;
        }
        else{
          (*(currRefracIter+N*i+j))=0;
        }

        /*******Sum Calculation*******/
        double sumVar=0.0;
        int iLeftCorner=N+i-R;
        int jLeftCorner=N+j-R;
        for  (int k=iLeftCorner; k<iLeftCorner+2*R+1; k++){
          for (int l=jLeftCorner; l<jLeftCorner+2*R+1; l++){
            sumVar+=(*(u+N*i+j)-*(u+N*(k%N)+l%N));
          }
        }

        (*(unext+N*i+j))=(*(u+N*i+j)+dt*(mi-(*(u+N*i+j))+sumCoeff*sumVar));

        if((*(unext+N*i+j))>=uth){   //Threshold crossed
          (*(unext+N*i+j))=0.0;
          if (it>=minMPVIter){
            (*(currThresCrossings+N*i+j))++;
          }
        }
        (*(u+N*i+j))=(*(unext+N*i+j));
      }//edw kleinei h for j.
    }//edw kleinei h for i.
    // we're going to manipulate done and use the cond, so we need the mutex
    //printf("Thread %d finished iteration %d\n", threadId, it);
    pthread_mutex_lock( &mutex );

    // increase the count of threads that have finished their work.
    done++;

    // wait up the main thread (if it is sleeping) to test the value of done
    pthread_cond_signal( &cond );
    pthread_cond_wait(&finCond, &mutex);
    //printf("Thread %d waiting for next iteration\n", threadId);
    pthread_mutex_unlock( & mutex );

    it++;
    //printf("Thread %d stepping to iteration %d\n", threadId, it);
  }
  //free(allArgs);
  //free(myArgs);
  return NULL;
}

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

  pthread_t threads[NUMTHREADS];
  allThreadArgs* allArgs[NUMTHREADS];
  threadArgs *args = malloc(sizeof *args);
  FILE *file1;
  char filename[80];
  double pi=3.14159265359 ;
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
  int i,j;

  double u[N][N];
  double unext[N][N];

  double currRefracIter[N][N];  //Current iterations already in refractory period
  int currThresCrossings[N][N];
  int maxMPVIter=20000;
  int minMPVIter=100000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.
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
      free(tmp);
    }
    i++;
  }
  fclose(file1);

  args->N=&N;
  args->R=&R;
  args->dt=&dt;
  args->mi=&mi;
  args->it=&it;
  args->maxRefracIter=&maxRefracIter;
  args->uth=&uth;
  args->sumCoeff=&sumCoeff;
  args->minMPVIter=&minMPVIter;
  args->totalIter=&totalIter;

  args->u=(double *)u;
  args->unext=(double *)unext;
  args->currRefracIter=(int *)currRefracIter;
  args->currThresCrossings=(int *)currThresCrossings;
  for( int t=0; t<(int)NUMTHREADS; t++ ){
    allArgs[t]=malloc(sizeof(allThreadArgs));
    allArgs[t]->ta=args;
    allArgs[t]->id=t;
    pthread_create(&threads[t], NULL, threadJob, allArgs[t]);
  }

  /*******Simulation*******/
  while (it<totalIter){

    if (it%1000==0) printf("Iteration %d of %d\n", it, totalIter);


    pthread_mutex_lock( &mutex );
    while(done < (int)NUMTHREADS ){
      //printf( "MAIN THREAD %d threads are done which is < %d so waiting on cond\n", done, (int)NUMTHREADS );

      /* block this thread until another thread signals cond. While
      blocked, the mutex is released, then re-aquired before this
      thread is woken up and the call returns. */
      pthread_cond_wait(&cond, &mutex);
      //puts( "MAIN THREAD condition was signalled!" );

      /* we go around the loop with the lock held */
    }
    done=0;
    pthread_mutex_unlock( & mutex );


    if (it>minMPVIter){
      if ((it-minMPVIter)%maxMPVIter==0){
        sprintf(filename, "ResultsPar_MPV_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",sigma,R,t);
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
      sprintf(filename, "ResultsPar_POT_LIF_2D_Classic_sigma_%lf_R_%d_time_%lf_.dat",sigma,R,t);
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
    //printf("MAIN THREAD writing is done, increasing it!\n");

    //pthread_mutex_lock( &mutex );
    it++;
    pthread_cond_broadcast(&finCond);
    //pthread_mutex_unlock( & mutex );

  } //edw kleinei h while.

  return(0);
}
