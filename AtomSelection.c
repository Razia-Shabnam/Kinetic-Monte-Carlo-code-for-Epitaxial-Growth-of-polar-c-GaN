// Subroutine for selection of type of atom (Ga, AdGa or N) randomly for diffusion.


#include "variables.h"

void AtomSelection()                                  
{                                                                      
    int Xc, Yc, X, Y, randomsite;

    randomsite = (int) (double)((nx*ny)*gsl_rng_uniform_pos(rsel));

    Xc = randomsite / ny;                                
    Yc = randomsite % ny;

    Y = Yc;

    if((Yc%2) == 0){
        if((Xc%2) == 0){X = (2*Xc);}
        else if((Xc%2) != 0){X = (2*Xc+1);}}
                 
    else if((Yc%2) != 0){
             if((Xc%2) == 0){X = (2*Xc+1);}
             else if((Xc%2) != 0){X = (2*Xc);}}
                                                     
     
    if((Box[X][Y][H[X][Y]]) == 1) 
      {
       if(((H[X][Y]) == 0) || ((H[X][Y]) == 4)) { printf("Bulk layers, No diffsion (%d,%d) H[X][Y] %d \n",X,Y,H[X][Y]);  exit(0)        ; }
       if((Box[X][Y][(H[X][Y])-5]) == 0)        { printf("overhang Ga atom, atom below Ga = %d\n",Box[X][Y][H[X][Y]-5]); exit(0)        ; }
       else if((Box[X][Y][(H[X][Y])-5])==2)     {/*printf("Ga atom selected for diffusion = %d\n",Box[X][Y][H[X][Y]]);*/ move_Ga( X, Y) ; }
      } 

    else if((Box[X][Y][H[X][Y]]) == 2)
           { 
            if((Box[X][Y][(H[X][Y])-3]) == 0)      { printf("overhang N atom, Height = %d, atom below N = %d\n",H[X][Y],Box[X][Y][H[X][Y]-3]); exit(0)  ; }
            else if((Box[X][Y][(H[X][Y])-3]) == 1) {/*printf("N atom selected  for diffusion = %d\n",Box[X][Y][H[X][Y]]);*/  move_N( X, Y) ; }
           }

    else if((Box[X][Y][H[X][Y]]) == 3)
           { 
            if((Box[X][Y][(H[X][Y])-3]) == 0)      { printf("overhang AdGa atom, Height = %d, atom below AdGa = %d\n",H[X][Y],Box[X][Y][H[X][Y]-3]); exit(0) ; }
            else if((Box[X][Y][(H[X][Y])-3]) == 1) {/*printf("AdGa atom selected  for diffusion = %d\n",Box[X][Y][H[X][Y]]);*/  move_AdGa( X, Y) ; }
           }   
}    
