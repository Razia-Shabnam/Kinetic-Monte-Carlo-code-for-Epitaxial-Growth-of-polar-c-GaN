// Subroutine is written to get the stable position for both N and AdGa. If diffused or deposited N encounter any AdGa in nearest neighbor of it, then it will dispalce AdGa n occupy its place.  


#include "variables.h"

int stable_N(int X, int Y)  
{
 
  int Xp1 = (X+1)%Nx,   Xm1 = (X-1+Nx)%Nx,   Yp1 = (Y+1)%Ny,   Ym1 = (Y-1+Ny)%Ny;
      

      if((H[X][Y]+1)%8==0){	      		  
                           if((Box[Xp1][Y][H[Xp1][Y]]) ==3){
                                                 if((H[X][Y] - H[Xp1][Y]) > 0){
                                                                               H[Xp1][Y] = H[Xp1][Y] + 5; Box[Xp1][Y][H[Xp1][Y]] = 1; Box[Xp1][Y][H[Xp1][Y]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; }
                                                 else if((H[X][Y] - H[Xp1][Y]) < 0){
                                                                                    H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]] = 1; Box[Xp1][Y][H[Xp1][Y]] = 0; H[Xp1][Y] = H[Xp1][Y] - 3 ; AdGa_N++ ; }}
                                                                                           
                           else if((Box[Xm1][Ym1][H[Xm1][Ym1]]) ==3){
                                                 if((H[X][Y] - H[Xm1][Ym1]) > 0){
                                                                                 H[Xm1][Ym1] = H[Xm1][Ym1] + 5; Box[Xm1][Ym1][H[Xm1][Ym1]] = 1; Box[Xm1][Ym1][H[Xm1][Ym1]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; }
						 else if((H[X][Y] - H[Xm1][Ym1]) < 0){
                                                                                      H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]] = 1; Box[Xm1][Ym1][H[Xm1][Ym1]] = 0; H[Xm1][Ym1] = H[Xm1][Ym1] - 3; AdGa_N++ ; }}       
                                                                                      
                           else if((Box[Xm1][Yp1][H[Xm1][Yp1]]) ==3){
                                                 if((H[X][Y] - H[Xm1][Yp1]) > 0){
                                                                                 H[Xm1][Yp1] = H[Xm1][Yp1] + 5; Box[Xm1][Yp1][H[Xm1][Yp1]] = 1; Box[Xm1][Yp1][H[Xm1][Yp1]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; }
                                                 else if((H[X][Y] - H[Xm1][Yp1]) < 0){
						          	                      H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]]  = 1; Box[Xm1][Yp1][H[Xm1][Yp1]] = 0; H[Xm1][Yp1] = H[Xm1][Yp1] - 3 ;  AdGa_N++ ; }}			   
                          }                   
			       
      else if((H[X][Y]+1)%8!=0){		  
                           if((Box[Xm1][Y][H[Xm1][Y]]) ==3){
                                                 if((H[X][Y] - H[Xm1][Y]) > 0){
							                       H[Xm1][Y] = H[Xm1][Y] + 5; Box[Xm1][Y][H[Xm1][Y]] = 1; Box[Xm1][Y][H[Xm1][Y]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; }
                                                 else if((H[X][Y] - H[Xm1][Y]) < 0){
							                            H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]] = 1; Box[Xm1][Y][H[Xm1][Y]] = 0; H[Xm1][Y] = H[Xm1][Y] - 3 ; AdGa_N++ ; }}			   
                                                                 
                           else if((Box[Xp1][Ym1][H[Xp1][Ym1]]) ==3) {
                                                 if((H[X][Y] - H[Xp1][Ym1]) > 0){
							                         H[Xp1][Ym1] = H[Xp1][Ym1] + 5; Box[Xp1][Ym1][H[Xp1][Ym1]] = 1; Box[Xp1][Ym1][H[Xp1][Ym1]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; } 
                                                 else if((H[X][Y] - H[Xp1][Ym1]) < 0){
							                              H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]] = 1; Box[Xp1][Ym1][H[Xp1][Ym1]] = 0; H[Xp1][Ym1] = H[Xp1][Ym1] - 3; AdGa_N++; }} 			   
                                  
                           else if((Box[Xp1][Yp1][H[Xp1][Yp1]]) ==3) {
                                                 if((H[X][Y] - H[Xp1][Yp1]) > 0){
                                                                                 H[Xp1][Yp1] = H[Xp1][Yp1] + 5; Box[Xp1][Yp1][H[Xp1][Yp1]] = 1; Box[Xp1][Yp1][H[Xp1][Yp1]-5] = 2; Box[X][Y][H[X][Y]] = 0; H[X][Y] = H[X][Y] - 3; N_AdGa++; }
                                                 else if((H[X][Y] - H[Xp1][Yp1]) < 0){
							                              H[X][Y] = H[X][Y] + 5; Box[X][Y][H[X][Y]] = 1; Box[Xp1][Yp1][H[Xp1][Yp1]] = 0; H[Xp1][Yp1] = H[Xp1][Yp1] - 3; AdGa_N++ ; }}   			   
                              }
return 0;
}

