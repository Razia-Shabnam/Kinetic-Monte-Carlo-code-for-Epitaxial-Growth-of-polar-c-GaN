//Subrotine is written to check if Ga atom has zero direct neighbor. In that case N get deposited over Ga or  excite AdGa atom (w.r.t to incoming site) and get depoisted in its place.


#include "variables.h"

int stable_N_OGa(int X, int Y, int Z)  
{
	
 int Xp1 = (X+1)%Nx,
     Xm1 = (X-1+Nx)%Nx,
     Yp1 = (Y+1)%Ny,
     Ym1 = (Y-1+Ny)%Ny;
     

     if(Z%8==0){ 
                if     ((Box[Xm1][Y][H[Xm1][Y]]) == 3)     { H[Xm1][Y]   = H[Xm1][Y] + 5;   Box[Xm1][Y][H[Xm1][Y]] = 1;     Box[Xm1][Y][H[Xm1][Y]-5] = 2;     Natoms++; N_AdGa++; dcnt++; }
		else if((Box[Xm1][Y][H[Xm1][Y]]) == 1)     { H[Xm1][Y]   = H[Xm1][Y] + 3;   Box[Xm1][Y][H[Xm1][Y]] = 2;     Natoms++; dcnt++; stable_N(Xm1,Y);   }

		else if((Box[Xp1][Ym1][H[Xp1][Ym1]]) == 3) { H[Xp1][Ym1] = H[Xp1][Ym1] + 5; Box[Xp1][Ym1][H[Xp1][Ym1]] = 1; Box[Xp1][Ym1][H[Xp1][Ym1]-5] = 2; Natoms++; N_AdGa++; dcnt++; }
		else if((Box[Xp1][Ym1][H[Xp1][Ym1]]) == 1) { H[Xp1][Ym1] = H[Xp1][Ym1] + 3; Box[Xp1][Ym1][H[Xp1][Ym1]] = 2; Natoms++; dcnt++; stable_N(Xp1,Ym1); }

	        else if((Box[Xp1][Yp1][H[Xp1][Yp1]]) == 3) { H[Xp1][Yp1] = H[Xp1][Yp1] + 5; Box[Xp1][Yp1][H[Xp1][Yp1]] = 1; Box[Xp1][Yp1][H[Xp1][Yp1]-5] = 2; Natoms++; N_AdGa++; dcnt++; }
                else if((Box[Xp1][Yp1][H[Xp1][Yp1]]) == 1) { H[Xp1][Yp1] = H[Xp1][Yp1] + 3; Box[Xp1][Yp1][H[Xp1][Yp1]] = 2; Natoms++; dcnt++; stable_N(Xp1,Yp1); }
               }

     else if(Z%8!=0){
                     if     ((Box[Xp1][Y][H[Xp1][Y]]) == 3)     { H[Xp1][Y]   = H[Xp1][Y] + 5;   Box[Xp1][Y][H[Xp1][Y]] = 1;     Box[Xp1][Y][H[Xp1][Y]-5] = 2;     Natoms++; N_AdGa++; dcnt++; }
                     else if((Box[Xp1][Y][H[Xp1][Y]]) == 1)     { H[Xp1][Y]   = H[Xp1][Y] + 3;   Box[Xp1][Y][H[Xp1][Y]] = 2;     Natoms++; dcnt++; stable_N(Xp1,Y);   }

                     else if((Box[Xm1][Ym1][H[Xm1][Ym1]]) == 3) { H[Xm1][Ym1] = H[Xm1][Ym1] + 5; Box[Xm1][Ym1][H[Xm1][Ym1]] = 1; Box[Xm1][Ym1][H[Xm1][Ym1]-5] = 2; Natoms++; N_AdGa++; dcnt++; }
                     else if((Box[Xm1][Ym1][H[Xm1][Ym1]]) == 1) { H[Xm1][Ym1] = H[Xm1][Ym1] + 3; Box[Xm1][Ym1][H[Xm1][Ym1]] = 2; Natoms++; dcnt++; stable_N(Xm1,Ym1); }

                     else if((Box[Xm1][Yp1][H[Xm1][Yp1]]) == 3) { H[Xm1][Yp1] = H[Xm1][Yp1] + 5; Box[Xm1][Yp1][H[Xm1][Yp1]] = 1; Box[Xm1][Yp1][H[Xm1][Yp1]-5] = 2; Natoms++; N_AdGa++; dcnt++; }
                     else if((Box[Xm1][Yp1][H[Xm1][Yp1]]) == 1) { H[Xm1][Yp1] = H[Xm1][Yp1] + 3; Box[Xm1][Yp1][H[Xm1][Yp1]] = 2; Natoms++; dcnt++; stable_N(Xm1,Yp1); }
                    }
return 0;

}

        


