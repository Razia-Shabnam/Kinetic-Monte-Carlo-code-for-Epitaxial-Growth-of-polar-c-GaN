//Subroutine for diffusion of AdGa atom. AdGa atom can diffuse over N atom having atleast one Ga atom as neighbor. 
//Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached.



#include"variables.h"

int move_AdGa(int X, int Y)            
{ 
  int Xp1 = (X+1)%Nx,   Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny ;	
    
      cntGa = 0, Prob = 0.0;
          
      if(((H[X][Y]+1)%8) ==0){
                              if((Box[Xp1][Y]  [H[X][Y]+1]) == 1)       { cntGa++;  }
                              if((Box[Xm1][Ym1][H[X][Y]+1]) == 1)       { cntGa++;  }
                              if((Box[Xm1][Yp1][H[X][Y]+1]) == 1)       { cntGa++;  }
                             }


      else if(((H[X][Y]+1)% 8)!=0){
                                  if((Box[Xm1][Y]  [H[X][Y]+1]) == 1)  { cntGa++;  }
                                  if((Box[Xp1][Ym1][H[X][Y]+1]) == 1)  { cntGa++;  }
                                  if((Box[Xp1][Yp1][H[X][Y]+1]) == 1)  { cntGa++;  }
                                 }


      if(cntGa==0) { /* Prob = 1.0 */ if(((H[X][Y]+1) % 8) ==0) { H7_AdGa( X, Y) ; }  else if (((H[X][Y]+1) % 8)!=0 ) { H3_AdGa( X, Y) ; } }     
	 			      	
      else {
            if     (cntGa==1) { Prob = Prob_AdGa1Ga; }	      
	    else if(cntGa==2) { Prob = Prob_AdGa2Ga; }
            else if(cntGa==3) { Prob = Prob_AdGa3Ga; } 
           
            double Racc = gsl_rng_uniform_pos(racc);    
	    if(Racc < Prob) {  
		             if(((H[X][Y]+1) % 8) ==0)        { H7_AdGa( X, Y) ; } 
 			     else if (((H[X][Y]+1) % 8) !=0 ) { H3_AdGa( X, Y) ; } }

	    else { /* failed diffusion */  NoAdGadiff++; }
           }
 return 0;
}                          
