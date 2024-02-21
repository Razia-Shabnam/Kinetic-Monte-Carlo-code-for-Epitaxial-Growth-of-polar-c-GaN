//Subroutine for diffusion of Ga atom at  height multiple of 4, 12, 20... This subroutine checks all the 18 diffuison sites. If no site found diffusion fails.


#include "variables.h"

int H4_Ga(int X,int Y)           
{	
 int i, j, k, l, a, b, c, d, temp, site, rej_sit, NN;
 
 int Xp1 = (X+1)%Nx,     Xm1 = (X-1+Nx)%Nx,    Yp1 = (Y+1)%Ny,     Ym1 = (Y-1+Ny)%Ny,    Xp2 = (X+2)%Nx,     Xm2 = (X-2+Nx)%Nx,
     Yp2 = (Y+2)%Ny,     Ym2 = (Y-2+Ny)%Ny,    Xp3 = (X+3)%Nx,     Xm3 = (X-3+Nx)%Nx,    Yp3 = (Y+3)%Ny,     Ym3 = (Y-3+Ny)%Ny;  
     
     Neighb = 0;   rej_sit = 0;   temp = 0;   NN = 0;
             
     for (a=0; a<tot_site; a++) { NX_H4[a] = 0;  NY_H4[a] = 0 ; }
     for (b=0; b<tot_site; b++) { reject_site[b] ; }

/*Sites for diffusion*/

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

		 		    if((Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]]) ==2) {
                                            Neighb = countNeighb(NX_H4[i],NY_H4[i],H[NX_H4[i]][NY_H4[i]]);
                                            if(Neighb!=0){
                                                          H[NX_H4[i]][NY_H4[i]] = H[NX_H4[i]][NY_H4[i]]  + 5 ;
                                                          Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]] = 1 ;
                                                          Box[X][Y][H[X][Y]] = 0;      H[X][Y] = H[X][Y] - 5 ;
                                                          Ga_move[cntN][cntAdGa]++ ; break; }}

                                    else if((Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]]) ==1) {
					     if(i>2){ Neighb = countNeighb(NX_H4[i],NY_H4[i],H[NX_H4[i]][NY_H4[i]]);						     
						      if(Neighb!=0){
                                                                    stable_AdGa(NX_H4[i],NY_H4[i],H[NX_H4[i]][NY_H4[i]]) ;
                                                                    if(Ga_Ga == 3){
                                                                                   H[NX_H4[i]][NY_H4[i]] = H[NX_H4[i]][NY_H4[i]]  + 3 ;
                                                                                   Box[NX_H4[i]][NY_H4[i]][H[NX_H4[i]][NY_H4[i]]] = 3 ;
                                                                                   Box[X][Y][H[X][Y]] = 0;      H[X][Y] = H[X][Y] - 5 ;
                                                                                   Ga_move[cntN][cntAdGa]++ ; break; } } } }
				    break;} }
		     

                     if(Box[X][Y][H[X][Y]] == 2){
                                                 for(c=0; c<3; c++){
                                                     if((Box[NX_H4[c]][NY_H4[c]][H[NX_H4[c]][NY_H4[c]]]) == 1){
                                                         if((H[NX_H4[c]][NY_H4[c]] - H[X][Y])>0){
                                                               NN = countNeighb(NX_H4[c],NY_H4[c],(H[NX_H4[c]][NY_H4[c]]-5));
                                                               if(NN == 0){
                                                                           H[X][Y] = H[X][Y]  + 5;
                                                                           Box[X][Y][H[X][Y]] = 1;
                                                                           Box[NX_H4[c]][NY_H4[c]][H[NX_H4[c]][NY_H4[c]]] = 0;
                                                                           H[NX_H4[c]][NY_H4[c]] = H[NX_H4[c]][NY_H4[c]]  - 5;
                                                                           side_move_N++;  stable_N(NX_H4[c], NY_H4[c]); break; }}}}
                                                 if(Box[X][Y][H[X][Y]] == 2) {
                                                                              stable_N(X, Y);
                                                                              if(Box[X][Y][H[X][Y]] == 1) {
										 for(d=0; d<3; d++){     
                                                                                     if((Box[NX_H4[d]][NY_H4[d]][H[NX_H4[d]][NY_H4[d]]]) == 2){
                                                                                             if((H[NX_H4[d]][NY_H4[d]] - H[X][Y])>0){
                                                                                                 NN = countNeighb(NX_H4[d],NY_H4[d],(H[NX_H4[d]][NY_H4[d]]-3));
                                                                                                 if(NN == 0){
                                                                                                             H[X][Y] = H[X][Y]  + 3;
                                                                                                             Box[X][Y][H[X][Y]] = 2;
                                                                                                             Box[NX_H4[d]][NY_H4[d]][H[NX_H4[d]][NY_H4[d]]] = 0;
                                                                                                             H[NX_H4[d]][NY_H4[d]] = H[NX_H4[d]][NY_H4[d]]  - 3;
                                                                                                             side_move_Ga++; stable_N(X, Y); break; }}}}}}
                                                                               break; }

                     else  if(site!=0){
                                       site--;
                                       reject_site[rej_sit] = i;  rej_sit++;
                                       for(k=1; k<rej_sit; k++){
                                                                temp = reject_site[k]; l = k-1;
                                                                while(l>=0 && temp <= reject_site[l]) { reject_site[l+1] = reject_site[l]; l = l-1; }
                                                                reject_site[l+1] = temp; }
                                      }

		     else  if(site==0) { NoGasites++; move_sbN( X, Y); }
                    }                                             
     return 0;
}  

