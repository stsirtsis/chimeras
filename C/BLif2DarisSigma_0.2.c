#include <stdio.h>
#include <math.h> 
#include <stdlib.h>


double myrand()
{
 return((double)rand()/(double)RAND_MAX);
}

int main()

{
//FILE *Lif_2_chim1_2D_0=fopen("Lif_2_chim1_2D_0.dat","w");

//FILE *Lif_2_chim1_2D_pm3d=fopen("Lif_2_chim1_2D_pm3d.dat","w");
//FILE *Lif_2_2D_freq=fopen("Lif_2_2D_freq.dat","w");
FILE *file1;
FILE *file2;
FILE *Lif_2_chim2_2D=fopen("BLif_2_2D_sigma_0.1.dat","w");
FILE *Lif_2_2D= fopen("BLif2-u[l][t]_2D_sigma_0.1.dat","w");//Gia dhmiourgia arxeiou gia diagramma enos u[i] ws pros xrono.
FILE *Lif_2_chim_tomh=fopen("BLif_2_chim_tomh_sigma_0.1.dat","w");
//int seed=823093819;   //seed1
 double pi=3.14159265359 ;
 // int seed= 395849566; //seed2 πλέγμα χιμαιρών
 int seed= 95679999;
// int seed= 65243455;
// int seed= 97243562;
// int seed= 65345981;

srand(seed);

//dhlwsh metablhtwn

int N=100;
//printf(" Posous neurwnes 8elete ana diastash; px an balete 10 o sunolikos # neurwnwn 8a einai 10*10=100 PROSOXH!! mono akeraious!!\n");
//scanf("%d",&N);

double dt=0.001; // bhma xronou
//printf("Poio 8a einai to bhma tou xronou?\n");
//scanf("%lf",&dt);
int TIME=10000;
//printf("Gia poso xrono na tre3ei to programma? (Mono fusikous ari8mous!!)\n");
//scanf("%d",&TIME);  
int tn=TIME/dt; //sunolika bhmata xronou
//printf("8a ginoun %d bhmata xronou",tn);


int a=25;
//printf(" H pleura tou tetragwnou al/shs einai 2*a+1 , poio 8a einai to a;\n a=");
//scanf("%d",&a);

double sigmaconst=0.2;
//printf(" Poia 8a einai h sta8era suzeu3hs?\n sigmaconst=");
//scanf("%lf",&sigmaconst);


double m=1.0; // sta8era m
double t=0.0; // o xronos

double tabs_arxiko ; 
int n=500;
//printf("Posa bhmata xronou 8a diarkei h efhsyxash;\n");
//scanf("%d",&n);
tabs_arxiko=n*dt;

int x=5;
int y=78;
//printf("\n Gia poion neurwna na kanw diagramma; Prepei na einai apo (0,0) mexri (%d,%d)\t",N-1,N-1);
//scanf("%d %d",&x,&y); 

//printf("\n check1 \n"); // Gia epi8eorhsh rohs .

double sigma[N][N][N][N] ;

//printf(" check2 \n"); // Gia epi8eorhsh rohs .

int i,ki,kki,ll;
int j,kj,kkj,pp;







//arxikopoihsh tou pinaka suzeu3hs

for ( i=0 ; i<N ; i++)
   {
    
     for ( j=0 ; j<N ; j++)
       {
        
         for ( ki=0 ; ki<N ; ki++)
           {
           
              for ( kj=0 ; kj<N ; kj++)
               {
                 sigma[i][j][ki][kj]=0.0;
                
                 //printf(" i=%d j=%d ki=%d kj=%d \n",i,j,ki,kj);
               }
           }
        }
   }

//printf(" \n check3 \n"); // Gia epi8eorhsh rohs .


for ( i=0 ; i<N ;i++)
   {
    for ( j=0 ; j<N ; j++)
       {
        for ( ki=i-a ; ki<=i+a ; ki++)
           {
            for ( kj=j-a ; kj<=j+a ; kj++)
               {
                kki=ki;
                kkj=kj;
                 if (ki>N-1)
                   {
                     kki=ki-N;
                   }               
                 else if (ki<0)
                   {
                     kki=N+ki ;
                   }
                 
                 if (kj>N-1)
                   {
                     kkj=kj-N;
                   }               
                 else if (kj<0)
                   {
                     kkj=N+kj ;
                   }
                 
                 sigma[i][j][kki][kkj]=sigmaconst;
                 //printf("%lf \n",sigma[i][j][kki][kkj]); 
  


               }
           }
        }
     }
//printf(" check4 \n"); // Gia epi8eorhsh rohs .



double uth=0.98;
//printf("Poio 8a einai to Dunamiko katwfliou; \n uth= ");
//scanf("%lf",&uth);
double u1[N][N];


                 
//***********************************************************************************************//


double Sum=0.0 ;
double tabs[N][N];// xronos euhshxashs
double u[N][N];
double w[N][N];
double tfreq[N][N];
double tfreq0[N][N];

int k=0; // sta8era pou anebainei kata ena me to bhma tou xronou wste na thn xrhsimopoihsw ston pinaka u[i][t] anti gia t (pou den pairnei aparaithta fusikes times) vazw to k.

int ekf[N][N] ; //ekf apo ekforthsh. otan einai 1 sto epomeno loop mhdenizetai to dunamiko kai arxizei na metra antistrofa o xronos efhsuxashs.




//arxikopoihsh tou pinaka dunamikwn, ekf , tabs.
for ( i=0; i<N ; i++)
   {
  for ( j=0 ;j<N ;j++) 
     {
      u[i][j]=0.0;
      ekf[i][j]=0;
      tabs[i][j]=0.0;
      w[i][j]=0.0;  //gia upologismo suxnothtas N/T apo N o ari8mos ekfortisewn.
      tfreq[i][j]=0.0;
      tfreq0[i][j]=0.0;
      
      do      
        {
             u1[i][j]=myrand(); 
        }       
      while(u1[i][j]>=uth);
      
     }
   }
//printf(" check 5 \n"); // Gia epi8eorhsh rohs .

int deikths=0;
int kfreq=2000000; //bhma meta to opoio me endiaferei na arxisw na upologizw thn syxnothta twn neurwnwn.

//Edw 3ekina h roh tou xronou.


while (t<=TIME)
  { 
  //  printf("\n 3ekinhse o xronos kai twra eimaste sto %lf",t);
  for ( i=0 ; i<N ; i++)
      {
        for ( j=0 ;j<N ;j++)
          {
           Sum=0.0;
             for  ( ki=0 ; ki<N ; ki++)
                 {
                 for ( kj=0 ; kj<N ; kj++)
                     {
                      Sum=Sum+ sigma[i][j][ki][kj]/((2*a+1)*(2*a+1))*(u1[i][j]-u1[ki][kj]);
                     }                                                   
                  }
        //   printf("\n to a8roisma einai %lf",Sum);
                           
       //   fprintf(Lif_2_sum,"%d \t %lf \n",i,Sum);
            
        
//entoles pou ulopoioun ekforthsh kai efhsuxash.
        if (ekf[i][j]==1)
            {
             
        //    printf(" \n exw ekfwrtish tou %d \t neurwna thn %d xronikh stigmh.",i,k);
             tabs[i][j]=tabs_arxiko;
             u[i][j]=0.0;
             
            }
        else if (tabs[i][j]>0.0) 
            {
                tabs[i][j]=tabs[i][j]-dt;
                u[i][j]=0.0;
                
        //       printf(" \n exw efhsuxash tou %d \t neurwna thn %d xronikh stigmh.",i,k);
                
            }
         
        else 
            {
             
               
               u[i][j]=u1[i][j]+dt*(m-u1[i][j] + Sum) ;
              if (k>=kfreq)
                {
                 tfreq[i][j]=tfreq[i][j]+dt ;
                }     
            }
            



         if(u[i][j]>=uth ) 
            {
             w[i][j]=(w[i][j]*tfreq0[i][j]+1)/(tfreq0[i][j]+tfreq[i][j]+tabs_arxiko);
             if (k>=kfreq)
               {
                tfreq0[i][j]=tfreq0[i][j]+tfreq[i][j]+tabs_arxiko ;
                tfreq[i][j]=0 ;
               }
             u[i][j]=0.0;
             ekf[i][j]=1 ;
             
             
             }
         else
          {
          ekf[i][j]=0 ;
          }
          
      
          if (k%(30000)==0)
           {

            
            
              if (i==0 && j==0)
                { 
                 if (k>=kfreq)
                  {
                    char filefreq[30];   
                    sprintf (filefreq, "BLif_2D_sigma%lf_freq%d.dat" ,sigmaconst, deikths);
                    file2=fopen(filefreq,"w");
                  }
                  char filename[30];
                  sprintf (filename, "BLif_2D_sigma%lf_pote%d.dat",sigmaconst, deikths);  
                  file1=fopen(filename ,"w");
                  
                }  
              
              fprintf(file1,"%d \t %d \t %lf \n",i,j,u[i][j]); 
            
            if (k>=kfreq)
               {
                fprintf(file2," %d \t %d \t %lf \n",i,j,w[i][j]);
               }
            
            if (i==N-1 && j==N-1)
              {
               fclose(file1);
               if (k>=kfreq)
                 {
                  fclose(file2);
                 }
              }
           }
          
          // if (k==tn-1 &&  j<N-1) 
            // {
            //   fprintf(Lif_2_chim1_2D_pm3d,"%d \t %d \t %lf \n ",i,j,u[i][j]);             
            // }
         //  else if (k==tn-1)
           //  {
              // fprintf(Lif_2_chim1_2D_pm3d,"%d \t %d \t %lf \n \n ",i,j,u[i][j]);    
            // }
          //if (k%(tn/1000)==0 &i<N-1)
          //  {
           //  fprintf(Lif_2_chim2,"%d \t %lf \n",i,u[i]);
     //       }
       //   else if (k%(tn/1000)==0)
         //   {
           //  fprintf(Lif_2_chim2,"%d \t %lf \n\n\n",i,u[i]);
           // }
           if (j==1)
             {
              if(k%(10000)==0 && i<N-1 )
                 {       
                  fprintf(Lif_2_chim_tomh,"%d \t %lf \n",i,u[i][j]);
                 }
              else if (k%(10000)==0)
                 {
                  fprintf(Lif_2_chim_tomh,"%d \t %lf \n \n ",i,u[i][j]);
                 }
             }
           
              
             
      }//edw kleinei h for j.
     }//edw kleinei h for i.     
     
     fprintf(Lif_2_2D,"%lf \t %lf \n",t,u[x][y]);
           for ( ll=0 ;  ll<N ; ll++)
             {
              for ( pp=0 ; pp<N ; pp++)
                {
                  u1[ll][pp]=u[ll][pp];
                }
             }
     k=k+1;
     if (k%(30000)==0)
       {
        deikths=deikths+1;
       }
     t=t+dt;      
     
       
  } //edw kleinei h while.

//printf(" ****** Telos!! ******\n");

return(0);
}
