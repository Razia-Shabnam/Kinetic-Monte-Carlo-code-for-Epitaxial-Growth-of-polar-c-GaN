//Subroutine to initialize an array for a hexagonal lattice having AaBb... stacking.
 

#include "variables.h"     

void Initial_Conditions()
{
 int I,J,K;

 // Initialization
 for(I=0; I<Nx; I++){
    for(J=0; J<Ny; J++){
       for(K=0; K<Nmax; K++){ 
	   Box[I][J][K] = 0; }}}

 // Ga and N Sublattices //
 for(K=0; K<Nz; K++){
    for(I=0; I<Nx; I++){
       for(J=0; J<Ny; J++){
          if((((I%4)==0) && (J%2)==0) || (((I+2)%4)==0 && ((J+1)%2==0))){
            if((K%8)==0) { Box[I][J][K] = 1; }
              else if((K+5)%8==0) { Box[I][J][K] = 2;   }}}}}

 for(K=0; K<Nz; K++){
    for(I=0; I<Nx; I++){
       for(J=0; J<Ny; J++){
          if(((((I+3)%4)==0) && ((J+1)%2)==0) || ((((I+1)%4)==0) && ((J%2)==0))){
            if(((K+4)%8)==0) { Box[I][J][K] = 1; }
              else if(((K+1)%8==0)) { Box[I][J][K] = 2; }}}}}



// Storing initial configuration in a file
 
 FILE *inp;
 inp = fopen("input","w");
 if(inp==NULL){ printf("error in opening input file to write\n"); exit(1);}

 for(K=(Nz-10); K<Nz; K++){
    if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
      for(I=0; I<Nx; I++){
         for(J=0; J<Ny; J++){
             H[I][J] = K;  
             fprintf(inp,"%d ",Box[I][J][K]); }
             fprintf(inp,"\n"); }
             fprintf(inp,"\n"); }}
 fclose(inp);
}   
