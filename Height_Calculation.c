//Subroutine to calculate Height of an any atom positioned at certain array index (i,j).


#include "variables.h"

void Height_Calculation()
{
 int I,J,K;

 for(I=0;I<Nx;I++){
    for(J=0;J<Ny;J++){      
       if((((I%4)==0) && (J%2)==0) || (((I+2)%4)==0 && ((J+1)%2==0)))
         {
          K=0;
          while(K<Nmax){
               if((K%8)==0 || ((K+5)%8)==0){
                 if(Box[I][J][K]!=0){
                    H[I][J]=K;}
                   }                    
                 K++;
	   }
	 }          
 
       else if(((((I+3)%4)==0) && ((J+1)%2)==0) || ((((I+1)%4)==0) && ((J%2)==0)))
              {
               K=0;
               while(K<Nmax){  
                     if((((K+4)%8)==0) || (((K+1)%8)==0)){
                        if(Box[I][J][K]!=0){
                           H[I][J]=K;
			  }}                  
                        K++;
                  }}}}
}
