//Subroutine for deposition atom. Following process have been included in deposition of atom.
//Ga adatom can deposit over N atom having atleast one Ga atom as neighbor. But in case of an isolated N atom, Ga adatom will get deposit over N present as neighbor w.r.t to slected N atom site. Ga adatom can also get deposited over surface Ga atoms having Ga/AdGa as nearest neighbor, thus formatiing  Adlayer Ga atoms.
//N adatom deposits over Ga atom having atleast one N atom as neighbor. But in case of an isolated N atom, N adatom will get deposit over Ga present as neighbor w.r.t selected Ga atom site. If incoming N adatom encounters Adlayer Ga atom at the deposition site, then incoming N adatom pushes Adlayer Ga atom over itself thus occupying adlayer Ga atom site. 
// 1 is used for Ga, 2 for N and 3 for adlayer Gaatom (AdGa). Atoms get deposited at height following lattice ratio for GaN(0001) system. Our code doesn't allow overhangs while deposition.



#include "variables.h"

void External_Flux()
{                          
  int Neighb, Xc, Yc, AdatomX, AdatomY;
      
      Rdepsite = (int) ((double)(nx*ny)*gsl_rng_uniform_pos(rdep));                  
       
      Xc = Rdepsite / ny;             
      Yc = Rdepsite % ny;
   
      /* finding actual coordinates from randomly generated deposition site */
      AdatomY = Yc;
      if((Yc%2)==0){
        if((Xc%2)==0){AdatomX = (2*Xc);}
        else if((Xc%2)!=0){AdatomX = (2*Xc+1);}}
            
      else if((Yc%2)!=0){
             if((Xc%2)==0){AdatomX = (2*Xc+1);}
             else if((Xc%2)!=0){AdatomX = (2*Xc);}}


/* Periodic boundary condition used in code. */

  int Xp1 = (AdatomX+1)%Nx,
      Xm1 = (AdatomX-1+Nx)%Nx,
      Yp1 = (AdatomY+1)%Ny,
      Ym1 = (AdatomY-1+Ny)%Ny;


  double Rtype;       
  Rtype = gsl_rng_uniform_pos(rflux);  


//*******************************************************  Ga DEPOSITION *****************************************************************
  
  if(Rtype < Flux_Ratio)
    {   
     Attempt_Ga++;	

     if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 2){                         
       if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-3]) == 1){                 
             Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
             if(Neighb!=0){
                           H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 5;
                           Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 1;  
                           Gaatoms++;     dcnt++;   }
	     else if(Neighb==0){ stable_Ga(AdatomX,AdatomY,H[AdatomX][AdatomY]); }}}                                                                                           

     else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 1)                   
             {                                                   
              if(((H[AdatomX][AdatomY])%8) == 0)
                {               
                 if((Box[Xm1][AdatomY][H[Xm1][AdatomY]])==2){
                    if((Box[Xm1][AdatomY][H[Xm1][AdatomY]-3])==1){ 
                        Neighb = countNeighb(Xm1,AdatomY,H[Xm1][AdatomY]);
                        if(Neighb!=0){
                                      H[Xm1][AdatomY] = H[Xm1][AdatomY]  + 5;
                                      Box[Xm1][AdatomY][H[Xm1][AdatomY]] = 1;
                                      Gaatoms++;     dcnt++; }
			else if(Neighb==0){ stable_Ga(Xm1,AdatomY,H[Xm1][AdatomY]); }}}		 
                                                                                                                                                     
                 else if((Box[Xp1][Ym1][H[Xp1][Ym1]])==2){
                        if((Box[Xp1][Ym1][H[Xp1][Ym1]-3])==1){
                            Neighb = countNeighb(Xp1,Ym1,H[Xp1][Ym1]);
                            if(Neighb!=0){
                                          H[Xp1][Ym1] = H[Xp1][Ym1]  + 5;
                                          Box[Xp1][Ym1][H[Xp1][Ym1]] = 1;
                                          Gaatoms++;     dcnt++;   } 
                            else if(Neighb==0){ stable_Ga(Xp1,Ym1,H[Xp1][Ym1]); }}}			    
			   
                 else if((Box[Xp1][Yp1][H[Xp1][Yp1]])==2){
                        if((Box[Xp1][Yp1][H[Xp1][Yp1]-3])==1){
                            Neighb = countNeighb(Xp1,Yp1,H[Xp1][Yp1]);
                            if(Neighb!=0){
                                          H[Xp1][Yp1] = H[Xp1][Yp1]  + 5;
                                          Box[Xp1][Yp1][H[Xp1][Yp1]] = 1;
                                          Gaatoms++;     dcnt++;   }
			    else if(Neighb==0){ stable_Ga(Xp1,Yp1,H[Xp1][Yp1]); }}}
                                                  
                 else{
                      Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                      if(Neighb!=0){
                                    stable_AdGa(AdatomX,AdatomY,H[AdatomX][AdatomY]) ;
                                    if(Ga_Ga == 3){
                                                   H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                                                   Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 3;
                                                   AdGaatoms++;     dcnt++;  }}}			      
		}                  
        
              else if(((H[AdatomX][AdatomY])%8) != 0)
                     {                           
                      if((Box[Xp1][AdatomY][H[Xp1][AdatomY]])==2){
                        if((Box[Xp1][AdatomY][H[Xp1][AdatomY]-3])==1){
                            Neighb = countNeighb(Xp1,AdatomY,H[Xp1][AdatomY]);
                            if(Neighb!=0){
                                          H[Xp1][AdatomY] = H[Xp1][AdatomY]  + 5;
                                          Box[Xp1][AdatomY][H[Xp1][AdatomY]] = 1;
                                          Gaatoms++;  dcnt++; }
			    else if(Neighb==0){ stable_Ga(Xp1,AdatomY,H[Xp1][AdatomY]); }}}
		      
                      else if((Box[Xm1][Ym1][H[Xm1][Ym1]])==2){
                             if((Box[Xm1][Ym1][H[Xm1][Ym1]-3])==1){
                                 Neighb = countNeighb(Xm1,Ym1,H[Xm1][Ym1]);
                                 if(Neighb!=0){
                                               H[Xm1][Ym1] = H[Xm1][Ym1]  + 5;
                                               Box[Xm1][Ym1][H[Xm1][Ym1]] = 1;
                                               Gaatoms++; dcnt++; }
	                         else if(Neighb==0){ stable_Ga(Xm1,Ym1,H[Xm1][Ym1]);}}}			 

                      else if((Box[Xm1][Yp1][H[Xm1][Yp1]])==2){
                             if((Box[Xm1][Yp1][H[Xm1][Yp1]-3])==1){
                                 Neighb = countNeighb(Xm1,Yp1,H[Xm1][Yp1]);
                                 if(Neighb!=0){
                                               H[Xm1][Yp1] = H[Xm1][Yp1]  + 5;
                                               Box[Xm1][Yp1][H[Xm1][Yp1]] = 1;
                                               Gaatoms++; dcnt++; }
		                 else if(Neighb==0){ stable_Ga(Xm1,Yp1,H[Xm1][Yp1]);}}}	     
	             
		      else{
                           Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                           if(Neighb!=0){
                                         stable_AdGa(AdatomX,AdatomY,H[AdatomX][AdatomY]) ;
                                         if(Ga_Ga == 3){
                                                        H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                                                        Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 3;
                                                        AdGaatoms++;     dcnt++;  }}}	
		     }  
             }   
    }

// **********************************************************  N DEPOSITION  *****************************************************************

  else if((Flux_Ratio < Rtype) && (Rtype < 1.00))
          { 
           Attempt_N++;		  

           if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 1){                     
             if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-5]) == 2){
                 Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                 if(Neighb!=0){			       
                               H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                               Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 2;
			       Natoms++;    dcnt++;   stable_N(AdatomX,AdatomY); }
	         else if(Neighb==0){ stable_N_OGa(AdatomX,AdatomY,H[AdatomX][AdatomY]); }}}

           else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]])==2)                        
                  {
                   if(((H[AdatomX][AdatomY]+1)%8)==0)  
                     {        	 	 
                      if((Box[Xp1][AdatomY][H[Xp1][AdatomY]])==1){                         
                        if((Box[Xp1][AdatomY][H[Xp1][AdatomY]-5])==2){
		            Neighb = countNeighb(Xp1,AdatomY,H[Xp1][AdatomY]);
                            if(Neighb!=0){ 
                                          H[Xp1][AdatomY] = H[Xp1][AdatomY]  + 3;
                                          Box[Xp1][AdatomY][H[Xp1][AdatomY]] = 2;
					  Natoms++;   dcnt++;   stable_N(Xp1,AdatomY); } 
			    else if(Neighb==0){ stable_N_OGa(Xp1,AdatomY,H[Xp1][AdatomY]); }}} 		      
		                           
                      else if((Box[Xm1][Ym1][H[Xm1][Ym1]])==1){
                             if((Box[Xm1][Ym1][H[Xm1][Ym1]-5])==2){
                                 Neighb = countNeighb(Xm1,Ym1,H[Xm1][Ym1]);
                                 if(Neighb!=0){
                                               H[Xm1][Ym1] = H[Xm1][Ym1]  + 3;
                                               Box[Xm1][Ym1][H[Xm1][Ym1]] = 2;
					       Natoms++;   dcnt++;   stable_N(Xm1,Ym1); } 
				 else if(Neighb==0){ stable_N_OGa(Xm1,Ym1,H[Xm1][Ym1]); }}}

                      else if((Box[Xm1][Yp1][H[Xm1][Yp1]])==1){
                             if((Box[Xm1][Yp1][H[Xm1][Yp1]-5])==2){
                                 Neighb = countNeighb(Xm1,Yp1,H[Xm1][Yp1]);
                                 if(Neighb!=0){
                                               H[Xm1][Yp1] = H[Xm1][Yp1]  + 3;
                                               Box[Xm1][Yp1][H[Xm1][Yp1]] = 2;
					       Natoms++;   dcnt++;   stable_N(Xm1,Yp1); } 
				 else if(Neighb==0){ stable_N_OGa(Xm1,Yp1,H[Xm1][Yp1]); }}}				 
		     } 

		   else if(((H[AdatomX][AdatomY]+1)%8)!=0)
                          { 
                           if((Box[Xm1][AdatomY][H[Xm1][AdatomY]])==1){                  
                             if((Box[Xm1][AdatomY][H[Xm1][AdatomY]-5])==2){  
                                 Neighb = countNeighb(Xm1,AdatomY,H[Xm1][AdatomY]);
                                 if(Neighb!=0){
                                               H[Xm1][AdatomY] = H[Xm1][AdatomY]  + 3;
                                               Box[Xm1][AdatomY][H[Xm1][AdatomY]] = 2;
					       Natoms++;   dcnt++;   stable_N(Xm1,AdatomY); } 
				 else if(Neighb==0){stable_N_OGa(Xm1,AdatomY,H[Xm1][AdatomY]); }}}
                    
                           else if((Box[Xp1][Ym1][H[Xp1][Ym1]])==1){
                                  if((Box[Xp1][Ym1][H[Xp1][Ym1]-5])==2){
                                      Neighb = countNeighb(Xp1,Ym1,H[Xp1][Ym1]);
                                      if(Neighb!=0){
                                                    H[Xp1][Ym1] = H[Xp1][Ym1]  + 3;
                                                    Box[Xp1][Ym1][H[Xp1][Ym1]] = 2;
                                                    Natoms++;   dcnt++;   stable_N(Xp1,Ym1); } 
				      else if(Neighb==0){ stable_N_OGa(Xp1,Ym1,H[Xp1][Ym1]); }}}   			   

                           else if((Box[Xp1][Yp1][H[Xp1][Yp1]])==1){
                                  if((Box[Xp1][Yp1][H[Xp1][Yp1]-5])==2){
                                      Neighb = countNeighb(Xp1,Yp1,H[Xp1][Yp1]);
                                      if(Neighb!=0){
                                                    H[Xp1][Yp1] = H[Xp1][Yp1]  + 3;
                                                    Box[Xp1][Yp1][H[Xp1][Yp1]] = 2;
						    Natoms++;   dcnt++;   stable_N(Xp1,Yp1); }
				      else if(Neighb==0){ stable_N_OGa(Xp1,Yp1,H[Xp1][Yp1]); }}}				      
                          } 
                  }		   

	   else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]])==3) {
		  if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-3])==1){
                      H[AdatomX][AdatomY] = H[AdatomX][AdatomY]    + 5 ;
		      Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]   = 1 ;
		      Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-5] = 2 ;
	              Natoms++; N_AdGa++;  dcnt++;   }}             	      
       } 
}    
