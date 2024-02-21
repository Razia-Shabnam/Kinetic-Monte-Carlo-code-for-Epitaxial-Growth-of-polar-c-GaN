#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<time.h>
#include "variables.h"
#include "Height_Calculation.c"


int main()
{

 int I, J, K;
 outNfiles = 867;
 rank = 40 ;
 char outfile[256];
 FILE *out;
 int ON,OR;

 for(ON=867;ON<=outNfiles;ON++){
     for(OR=0;OR<rank;OR++){
          sprintf(outfile,"outputGaN%dx%d",ON,OR);
          out = fopen(outfile,"r");
          if(out==NULL){ printf("error opening outGaN files\n"); exit(1); }

          for(K=(Nz-10);K<Nmax;K++){
                if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
                     for(I=0;I<Nx;I++){
                             for(J=0;J<Ny;J++){
                                 fscanf(out,"%d ",&Box[I][J][K]);}  fscanf(out,"\n");}  fscanf(out,"\n"); }
                                   }

 Height_Calculation();

 int Ga, N, AdGa, X, Y, Z ;
 Ga =0, N=0,AdGa =0;

 for(Z=Nz; Z<Nmax; Z++){
     if(((Z%8)==0) || (((Z+5)%8)==0) || (((Z+4)%8)==0) || (((Z+1)%8)==0)){
          for(X=0; X<Nx; X++){
              for(Y=0; Y<Ny; Y++){
                  if     (Box[X][Y][Z]==1) { Ga++;   }
                  else if(Box[X][Y][Z]==2) { N++;    }
                  else if(Box[X][Y][Z]==3) { AdGa++; } }}}}

 printf( "ON %d OR %d Ga %d N%d AdGa%d\n", ON, OR, Ga, N, AdGa);}}
 return 0;
}  // End of main of PP code;

