//Subroutine for diffusion of Ga atom.  Following process have been included in deposition of atom.
//Ga adatom diffuses over N atom having atleast one Ga atom as neighbor. Ga adatom can also get diffuse over surface Ga atoms having Ga/AdGa as nearest neighbor, thus forming  AdGa atoms.
//Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached. 
//Formation of AdGa will take place only for next to next neighbr and next to next to next neighbor as dffusion to nearest neighbor will give rise N atom and that is not a stable configuration right now in our code.
//Subsurfce diffusion of N atom will be attempted if surface Ga atom moves fails due to no site avalaible or the Probabilty too high.



#include"variables.h"

int move_Ga (int X,int Y)           
{   
 int Xp1 = (X+1)%Nx,  Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny; 
     
     cntN = 0;   cntAdGa = 0;   Prob = 0.0;      

     /*counting bonds to calculate probability*/

     if((H[X][Y] % 8) == 0){ 
                            if((Box[Xm1][Y]  [H[X][Y]-1]) == 2)       { cntN++;    }
                            if((Box[Xp1][Ym1][H[X][Y]-1]) == 2)       { cntN++;    }
                            if((Box[Xp1][Yp1][H[X][Y]-1]) == 2)       { cntN++;    }
 	  	            if((Box[Xm1][Y]  [H[X][Y]-1]) == 3)       { cntAdGa++; }
                            if((Box[Xp1][Ym1][H[X][Y]-1]) == 3)       { cntAdGa++; }
                            if((Box[Xp1][Yp1][H[X][Y]-1]) == 3)       { cntAdGa++; }
                           }

     else if((H[X][Y] % 8)!= 0){ 
                                if((Box[Xp1][Y]  [H[X][Y]-1]) == 2)  { cntN++;    }
                                if((Box[Xm1][Ym1][H[X][Y]-1]) == 2)  { cntN++;    }
                                if((Box[Xm1][Yp1][H[X][Y]-1]) == 2)  { cntN++;    }
                                if((Box[Xp1][Y]  [H[X][Y]-1]) == 3)  { cntAdGa++; }
                                if((Box[Xm1][Ym1][H[X][Y]-1]) == 3)  { cntAdGa++; }
                                if((Box[Xm1][Yp1][H[X][Y]-1]) == 3)  { cntAdGa++; }   			   
                               } 

    
     if     ((cntN==0) && (cntAdGa==0))             { Prob = Prob_Gaad         ; }
     else if((cntN==0) && (cntAdGa==1))             { Prob = Prob_Ga1AdGa      ; }
     else if((cntN==0) && (cntAdGa==2))             { Prob = Prob_Ga2AdGa      ; }
     else if((cntN==0) && (cntAdGa==3))             { Prob = Prob_Ga3AdGa      ; }
     else if((cntN==1) && (cntAdGa==1))             { Prob = Prob_Ga1N_1AdGa   ; }
     else if((cntN==1) && (cntAdGa==2))             { Prob = Prob_Ga1N_2AdGa   ; }
     else if((cntN==2) && (cntAdGa==1))             { Prob = Prob_Ga2N_1AdGa   ; }
     else if((cntN==1) && (cntAdGa==0))             { Prob = Prob_Ga1N         ; }
     else if((cntN==2) && (cntAdGa==0))             { Prob = Prob_Ga2N         ; }
     else if((cntN==3) && (cntAdGa==0))             { NoGadiff++;      return 0; }
             
     double Racc = gsl_rng_uniform_pos(racc);                     
	  
     if(Racc < Prob) {
		      if((H[X][Y] % 8) == 0)      { H8_Ga( X, Y) ; } 
		      else if((H[X][Y] % 8) !=0 ) { H4_Ga( X, Y) ; }  
	                   }
     
     else { NoGadiff++; move_sbN( X, Y); }
	     	                   
 return 0;
}                  
