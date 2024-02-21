// Subroutine for diffusion of sub surface N atom at  height multiple of 3, 11, 19,.....


#include "variables.h"

int H3_sbN(int X,int Y)           
{
 int i, j, k, l, a, b, c, temp, site, rej_sit, NN;
 
 int Xp1 = (X+1)%Nx,     Xm1 = (X-1+Nx)%Nx,    Yp1 = (Y+1)%Ny,     Ym1 = (Y-1+Ny)%Ny,    Xp2 = (X+2)%Nx,     Xm2 = (X-2+Nx)%Nx,
     Yp2 = (Y+2)%Ny,     Ym2 = (Y-2+Ny)%Ny,    Xp3 = (X+3)%Nx,     Xm3 = (X-3+Nx)%Nx,    Yp3 = (Y+3)%Ny,     Ym3 = (Y-3+Ny)%Ny;  
     
     Neighb = 0, rej_sit = 0, temp = 0;  NN =0;    
             
     for (a=0; a<tot_site; a++) { NX_H8[a] = 0;  NY_H8[a] = 0 ; }
     for (b=0; b<tot_site; b++) { reject_site[b] ; }
     
     NX_H8[0]  = Xm1 ;   NY_H8[0]  = Y   ;       NX_H8[1]  = Xp1 ;   NY_H8[1]  = Ym1 ;       NX_H8[2]  = Xp1 ;   NY_H8[2]  = Yp1 ;

     NX_H8[3]  = Xm2 ;   NY_H8[3]  = Ym1 ;       NX_H8[4]  = Xm2 ;   NY_H8[4]  = Yp1 ;       NX_H8[5]  = X   ;   NY_H8[5]  = Ym2 ;
     NX_H8[6]  = Xp2 ;   NY_H8[6]  = Ym1 ;       NX_H8[7]  = X   ;   NY_H8[7]  = Yp2 ;       NX_H8[8]  = Xp2 ;   NY_H8[8]  = Yp1 ;

     NX_H8[9]  = Xm3 ;   NY_H8[9]  = Ym1 ;       NX_H8[10] = Xm1 ;   NY_H8[10] = Ym2 ;       NX_H8[11] = Xm3 ;   NY_H8[11] = Yp1 ;
     NX_H8[12] = Xm1 ;   NY_H8[12] = Yp2 ;       NX_H8[13] = Xp1 ;   NY_H8[13] = Ym3 ;       NX_H8[14] = Xp3 ;   NY_H8[14] = Ym2 ;
     NX_H8[15] = Xp3 ;   NY_H8[15] = Y   ;       NX_H8[16] = Xp1 ;   NY_H8[16] = Yp3 ;       NX_H8[17] = Xp3 ;   NY_H8[17] = Yp2 ;
     
            
     site = tot_site;

     while(site > 0){
                     double Rdir = gsl_rng_uniform_pos(rdir);
                     double randiff = 1.0/ site;

                     for (i=0; i<site; i++){
                               if(((i*randiff) < Rdir) && (Rdir < ((i+1)*randiff))){

                                    if(rej_sit >= 1){
                                                     for(j=0; j<rej_sit; j++){
                                                              if(i >= reject_site[j]) { i = i+1; }}}				    
                                      
	                            if((Box[NX_H8[i]][NY_H8[i]][H[NX_H8[i]][NY_H8[i]]]) ==3) {
                                            Neighb = countNeighb(NX_H8[i],NY_H8[i],H[NX_H8[i]][NY_H8[i]]);    
                                            if(Neighb!=0){
							  H[NX_H8[i]][NY_H8[i]] = H[NX_H8[i]][NY_H8[i]]    + 5 ;
                                                          Box[NX_H8[i]][NY_H8[i]][H[NX_H8[i]][NY_H8[i]]]   = 1 ; 
						          Box[NX_H8[i]][NY_H8[i]][H[NX_H8[i]][NY_H8[i]]-5] = 2 ;								  
                                                          Box[X][Y][H[X][Y]] = 0;        H[X][Y] = H[X][Y] - 5 ;
							  Box[X][Y][H[X][Y]] = 3;  sbN_move++;   break; }}				                           	          				    
			            break; }}		     
		     
                     if(Box[X][Y][H[X][Y]] == 3){			     
                                                 for(c=0; c<3; c++){
                                                                    if((Box[NX_H8[c]][NY_H8[c]][H[NX_H8[c]][NY_H8[c]]]) == 1){	
								        NN = countNeighb(NX_H8[c],NY_H8[c],(H[NX_H8[c]][NY_H8[c]]-5));}								    
								    else if((Box[NX_H8[c]][NY_H8[c]][H[NX_H8[c]][NY_H8[c]]]) == 3){
									     NN = countNeighb(NX_H8[c],NY_H8[c],(H[NX_H8[c]][NY_H8[c]]));}									    
								    if(NN == 0){ 									   
				                                                H[X][Y] = H[X][Y]    + 5;
                                                                                Box[X][Y][H[X][Y]]   = 1;
                                                                                Box[X][Y][H[X][Y]-5] = 2;
                                                                                Box[NX_H8[i]][NY_H8[i]][H[NX_H8[i]][NY_H8[i]]] = 0 ;
                                                                                H[NX_H8[i]][NY_H8[i]] = H[NX_H8[i]][NY_H8[i]]  - 5 ;
                                                                                Box[NX_H8[i]][NY_H8[i]][H[NX_H8[i]][NY_H8[i]]] = 3 ;
                                                                                rejected_sb_move++;}}  						 
                                                break; }		     
 
                     else  if(site!=0){
                                       site--;
                                       reject_site[rej_sit] = i;  rej_sit++;
                                       for(k=1; k<rej_sit; k++){
                                                                temp = reject_site[k]; l = k-1;
                                                                while(l>=0 && temp <= reject_site[l]) { reject_site[l+1] = reject_site[l]; l = l-1; }
                                                                reject_site[l+1] = temp; }
                                      }

                     else  if(site==0) {  NosbNsites++; } 
                    }
 return 0;
}  
