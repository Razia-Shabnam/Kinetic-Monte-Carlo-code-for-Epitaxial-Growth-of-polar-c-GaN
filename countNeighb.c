// Sept 28, 2020: Program to count neighbours for X,Y coordinates. This sub routine was written just to make sure that no ad atom (Ga, N or AdGa ) get deposited over isolated atom. 


#include "variables.h"

int countNeighb(int X, int Y, int Z)  
{
   int neighb = 0;     
 
   int Xp1 = (X+1)%Nx,  Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny;
      
       if((Box[X][Y][Z]) ==1)
          {
           if(Z%8==0){		 
                      if((Box[Xm1][Y]  [Z-1]) ==2)              { neighb++; }
                      if((Box[Xp1][Ym1][Z-1]) ==2)              { neighb++; }
	              if((Box[Xp1][Yp1][Z-1]) ==2)              { neighb++; }		      
                     }

 	   else if(Z%8!=0){
 		           if((Box[Xp1][Y]  [Z-1]) ==2)         { neighb++; }
                           if((Box[Xm1][Ym1][Z-1]) ==2)         { neighb++; }
                           if((Box[Xm1][Yp1][Z-1]) ==2)         { neighb++; }			   
 	                  } 
	  }   
        
       else if(((Box[X][Y][Z]) ==2) || ((Box[X][Y][Z]) ==3))
                {  
   	         if((Z+1)%8==0){	      
                                if((Box[Xp1][Y]  [Z+1]) ==1)          { neighb++; }
                                if((Box[Xm1][Ym1][Z+1]) ==1)          { neighb++; }
                                if((Box[Xm1][Yp1][Z+1]) ==1)          { neighb++; }     				
	                       }          

                 else if((Z+1)%8!=0){
                                     if((Box[Xm1][Y]  [Z+1]) ==1)     { neighb++; }
                                     if((Box[Xp1][Ym1][Z+1]) ==1)     { neighb++; }
                                     if((Box[Xp1][Yp1][Z+1]) ==1)     { neighb++; }     				     
                                    }  	  
                }

   return (neighb); 
}


