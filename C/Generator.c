#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Returns a real number in the interval [0,1]
double myrand(){
  return((double)rand()/(double)RAND_MAX);
}

int main() {

  FILE *file1;

  int N=100;
  double uth=0.98;
  int i,j;
  double tmp;

int seed = time(NULL);
  srand(seed);

  file1=fopen("initial.csv","w");
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      do{
        tmp=myrand();
      }while(tmp>=uth);
      fprintf(file1,"%lf,",tmp);
    }
    fprintf(file1,"\n");
  }
  fclose(file1);

  return 0;
}
