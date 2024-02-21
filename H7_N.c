//Subroutine for diffusion of N atom at  height multiple of 7, 15, 23,....This subroutine checks all the 18 diffuison sites. If no site found diffusion fails.


#include "variables.h"

int H7_N(int X,int Y)           
{

 int i, j, k, l, a, b, temp, site, rej_sit;
 
 int Xp1 = (X+1)%Nx,     Xm1 = (X-1+Nx)%Nx,    Yp1 = (Y+1)%Ny,     Ym1 = (Y-1+Ny)%Ny,    Xp2 = (X+2)%Nx,     Xm2 = (X-2+Nx)%Nx,
     Yp2 = (Y+2)%Ny,     Ym2 = (Y-2+Ny)%Ny,    Xp3 = (X+3)%Nx,     Xm3 = (X-3+Nx)%Nx,    Yp3 = (Y+3)%Ny,     Ym3 = (Y-3+Ny)%Ny;  
     
     Neighb = 0, rej_sit = 0, temp = 0;     
                  
     for (a=0; a<tot_site; a++) { NX_H4[a] = 0;  NY_H4[a] = 0 ; }
     for (b=0; b<tot_site; b++) { reject_site[b] ; }

     NX_H4[0]  = Xp1 ;   NY_H4[0]  = Y   ;       NX_H4[1]  = Xm1 ;   NY_H4[1]  = Ym1 ;       NX_H4[2]  = Xm1 ;   NY_H4[2]  = Yp1 ;

     NX_H4[3]  = Xp2 ;   NY_H4[3]  = Ym1 ;       NX_H4[4]  = Xp2 ;   NY_H4[4]  = Yp1 ;       NX_H4[5]  = X   ;   NY_H4[5]  = Ym2 ;
     NX_H4[6]  = Xm2 ;   NY_H4[6]  = Ym1 ;       NX_H4[7]  = X   ;   NY_H4[7]  = Yp2 ;       NX_H4[8]  = Xm2 ;   NY_H4[8]  = Yp1 ;

     NX_H4[9]  = Xp3 ;   NY_H4[9]  = Ym1 ;       NX_H4[10] = Xp1 ;   NY_H4[10] = Ym2 ;       NX_H4[11] = Xp3 ;   NY_H4[11] = Yp1 ;
     NX_H4[12] = Xp1 ;   NY_H4[12] = Yp2 ;       NX_H4[13] = Xm1 ;   NY_H4[13] = Ym3 ;       NX_H4[14] = Xm3 ;   NY_H4[14] = Ym2 ;
     NX_H4[15] = Xm3 ;   NY_H4[15] = Y   ;       NX_H4[16] = Xm1 ;   NY_H4[16] = Yp3 ;       NX_H4[17] = Xm3 ;   NY_H4[17] = Yp2 ;
           
     site = tot_site;

     while(site > 0){
                     double Rdir = gsl_rng_uniform_pos(rdir);
                     double randiff = 1.0/ site;

                     for (i=0; i<site; i++){
                               if(((i*randiff) < Rdir) && (Rdir < ((i+1)*randiff))){

                                    if(rej_sit >= 1){
                                                     for(j=0; j<rej_sit; j++){
                                                              if(i >= reject_site[j]) { i = i+1; }}}				    
  
	                            if((Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]]) ==3) {
                                              Neighb = countNeighb(NX_H4[i],NY_H4[i],H[NX_H4[i]][NY_H4[i]]);    
                                              if(Neighb!=0){
							    H[NX_H4[i]][NY_H4[i]] = H[NX_H4[i]][NY_H4[i]]    + 5 ;
                                                            Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]]   = 1 ; 
						            Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]-5] = 2 ;								  
                                                            Box[X][Y][H[X][Y]] = 0;        H[X][Y] = H[X][Y] - 3 ;
                                                            N_move[cntGa]++ ;  break; }}
                       	          
			            else if((Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]]) ==1) {
                                                   Neighb = countNeighb(NX_H4[i],NY_H4[i],H[NX_H4[i]][NY_H4[i]]);
                                                   if(Neighb!=0){
								 H[NX_H4[i]][NY_H4[i]] = H[NX_H4[i]][NY_H4[i]]  + 3 ;
							         Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]] = 2 ;
							         Box[X][Y][H[X][Y]] = 0;      H[X][Y] = H[X][Y] - 3 ;
							         N_move[cntGa]++;  stable_N(NX_H4[i], NY_H4[i]); break; }} 			

			            break; }}		     

                     if(Box[X][Y][H[X][Y]] == 1){ break; }		     
                     
                     else  if(site!=0){
                                       site--;
                                       reject_site[rej_sit] = i;  rej_sit++;
                                       for(k=1; k<rej_sit; k++){
                                                                temp = reject_site[k]; l = k-1;
                                                                while(l>=0 && temp <= reject_site[l]) { reject_site[l+1] = reject_site[l]; l = l-1; }
                                                                reject_site[l+1] = temp; }
                                      }

                     else  if(site==0) { NoNsites++; }
                    }
     return 0;
}  
