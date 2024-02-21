//Nov23, 2021: Main program to run code.

#include "headerfiles.h"
#include "variables.h"
#include <mpi.h>

int main( int argc, char *argv[])
{

       MPI_Init(&argc,&argv);
       MPI_Comm_size( MPI_COMM_WORLD, &size );
       MPI_Comm_rank( MPI_COMM_WORLD, &rank );

       int r, m;
          
       for (r=0; r<4; r++) { for (m=0; m<4; m++)   { Ga_move[r][m] = 0 ; }  AdGa_move[r] = 0;  N_move[r] = 0; }

       nNrot = 0, nNhop = 0, AdGa_N = 0, N_AdGa = 0, NoAdGadiff = 0, NoGadiff = 0, NoNdiff = 0, NosbNdiff = 0;
       NoGasites = 0, NoNsites = 0, NoAdGasites = 0, NosbNsites = 0, side_move_AdGa= 0, side_move_Ga = 0, side_move_N =0, rejected_sb_move = 0;
       Attempt_N = 0, Attempt_Ga = 0, sbN_Island = 0, sbN_move = 0, Ga3_sbN = 0 ;
       
       FILE *accpmoves;
       FILE *accpdep;

       int I, J;

       clock_t start, end;
       double cpu_time_used;
       start = clock();

       unsigned long int seed;
       seed = time(NULL);
       
       if(rank==0){printf("Today's Date and Time is = %s\n",ctime(&seed));}

       const gsl_rng_type *T;
       T = gsl_rng_mt19937;

                           gsl_rng_env_setup();
                           gsl_rng_default_seed=abs(seed)+rank;
                           rsel = gsl_rng_alloc(T);              // random no. for selecting an atomic site
                           if(rsel==NULL){printf("Problem in rsel, rng_instance, in main\n"); exit(1);}
                           seed+=1;

                           gsl_rng_env_setup();
                           gsl_rng_default_seed=abs(seed)+rank;
                           rdir = gsl_rng_alloc(T);             // random no. to select a direction to move the selected atom 
                           if(rdir==NULL){printf("Problem in rdir, rng_instance, in main\n"); exit(1);}   
                           seed+=1;
 
                           gsl_rng_env_setup();
                           gsl_rng_default_seed=abs(seed)+rank;
                           racc = gsl_rng_alloc(T);             // random no. to compare the Probability and random no. fo acceptance/rejection
                           if(racc==NULL){printf("Problem in racc, rng_instance, in main\n"); exit(1);}
                           seed+=1;

                           gsl_rng_env_setup();
                           gsl_rng_default_seed=abs(seed)+rank;
                           rdep = gsl_rng_alloc(T);            // random no. to select a site to deposit an atom
                           if(rdep==NULL){printf("Problem in rdep, rng_instance, in main\n"); exit(1);}
                           seed+=1;
       
                           gsl_rng_env_setup();
                           gsl_rng_default_seed=abs(seed)+rank;
                           rflux = gsl_rng_alloc(T);           // random no. to choose a type of atom (Ga/N)
                           if(rflux==NULL){printf("Problem in rflux, rng_instance, in main\n"); exit(1);}
                           seed+=1;
                        
     
       outNfiles = 0 ;

       if   (outNfiles == 0) { Initial_Conditions(); ReadInput();}
       else { ReadOutput();}

    /* Epitaxial Ga Probabilty calculation */
       Prob_Gaad       = exp(-(Energy_Gaad        - Energy_fast)/(kB*Temp));     
       Prob_Ga1AdGa    = exp(-(Energy_Ga1AdGa     - Energy_fast)/(kB*Temp)); 
       Prob_Ga2AdGa    = exp(-(Energy_Ga2AdGa     - Energy_fast)/(kB*Temp)); 
       Prob_Ga3AdGa    = exp(-(Energy_Ga3AdGa     - Energy_fast)/(kB*Temp)); 
       Prob_Ga1N_1AdGa = exp(-(Energy_Ga1N_1AdGa  - Energy_fast)/(kB*Temp)); 
       Prob_Ga1N_2AdGa = exp(-(Energy_Ga1N_2AdGa  - Energy_fast)/(kB*Temp)); 
       Prob_Ga2N_1AdGa = exp(-(Energy_Ga2N_1AdGa  - Energy_fast)/(kB*Temp));  
       Prob_Ga1N       = exp(-(Energy_Ga1N        - Energy_fast)/(kB*Temp)); 
       Prob_Ga2N       = exp(-(Energy_Ga2N        - Energy_fast)/(kB*Temp)); 
       Prob_Ga3N       = exp(-(Energy_Ga3N        - Energy_fast)/(kB*Temp)); 
                                    

    /* Epitaxial AdGa Probabilty calculation */
       Prob_AdGa0Ga    = exp(-(Energy_AdGa0Ga     - Energy_fast)/(kB*Temp)); 
       Prob_AdGa1Ga    = exp(-(Energy_AdGa1Ga     - Energy_fast)/(kB*Temp)); 
       Prob_AdGa2Ga    = exp(-(Energy_AdGa2Ga     - Energy_fast)/(kB*Temp)); 
       Prob_AdGa3Ga    = exp(-(Energy_AdGa3Ga     - Energy_fast)/(kB*Temp));
			 
    /* On surface N Probabilty calculation */
       Prob_N0Ga1N     = exp(-(Energy_N0Ga1N      - Energy_fast)/(kB*Temp));
       Prob_N0Ga2N     = exp(-(Energy_N0Ga2N      - Energy_fast)/(kB*Temp)); 
       Prob_N0Ga3N     = exp(-(Energy_N0Ga3N      - Energy_fast)/(kB*Temp));
       Prob_NRot       = exp(-(Energy_N1GaRot     - Energy_fast)/(kB*Temp)); 
       Prob_NHop       = exp(-(Energy_N1GaHop     - Energy_fast)/(kB*Temp));                             
       Prob_N2Ga       = exp(-(Energy_N2Ga        - Energy_fast)/(kB*Temp));       
       Prob_N3Ga       = exp(-(Energy_N3Ga        - Energy_fast)/(kB*Temp));

   /* sub surface N Probabilty calculation */
       Prob_NGa_Ga    = exp(-(Energy_NGa_Ga    - Energy_fast)/(kB*Temp));
       
 
       start_point = ((outNfiles*50)+1) ; // Restart point

       for(I = start_point; I<=  Deposition ; I++){ 
  	       for(J=1; J<=Nsteps;J++){
                        if(((J-1)%Nsteps)==0){
         	          	              External_Flux(); 			        
                                              cnt++;     
                                              
                                                   char accpdepfile[256]; 
                                                   sprintf(accpdepfile,"20ML_accpdep%d",1,rank);   
                                                   accpdep = fopen(accpdepfile,"a");  
                                                   fprintf(accpdep,"ExtFlux called= %d\t  Atoms Deposited= %d\t  GaDep= %d\t  NDep= %d\t  AdGaDep= %d\t Attempt_Ga= %d\t Attempt_N= %d\t seed= %d\t rdep= %d\t rflux= %d\t racc = %d\t rdir= %d\t rsel= %d\n" ,cnt,dcnt,Gaatoms,Natoms,AdGaatoms,Attempt_Ga,Attempt_N, seed,rdep, rflux, racc, rdir, rsel);
                                                   fclose(accpdep);
                                                  }
                          AtomSelection();       
                                          } 
       if((I%(50))==0){     		  
                       outNfiles++;
                       OutGaNfiles();
                       
	                    char accpmovefile[256];
                            sprintf(accpmovefile,"20ML_accpmoves%d",1,rank);
                            accpmoves = fopen(accpmovefile,"a");                          
		
		            fprintf(accpmoves,"Gaad= %ld\t Ga1N= %ld\t Ga2N= %ld\t Ga3N= %ld\t Ga2N_1AdGa= %ld\t Ga1N_2AdGa= %ld\t Ga1N_1AdGa= %ld\t  Ga3AdGa= %ld\t Ga2AdGa= %ld\t Ga1AdGa= %ld\t NoGadiff= %ld\n",Ga_move[0][0], Ga_move[1][0], Ga_move[2][0], Ga_move[3][0], Ga_move[2][1], Ga_move[1][2], Ga_move[1][1], Ga_move[0][3], Ga_move[0][2],Ga_move[0][1], NoGadiff);	
                            fprintf(accpmoves,"Nad= %ld\t N1Ga= %ld\t Nrot= %ld\t Nhop= %ld\t N2Ga= %ld\t N3Ga= %ld\t NoNdiff= %ld\n",N_move[0],N_move[1], nNrot, nNhop, N_move[2], N_move[3], NoNdiff);			  
	                    fprintf(accpmoves,"AdGa0Ga= %ld\t AdGa1Ga= %ld\t AdGa2Ga= %ld\t AdGa3Ga= %ld\t NoAdGadiff= %ld\t sbN_move= %ld\t NosbNdiff= %ld\n",AdGa_move[0], AdGa_move[1], AdGa_move[2], AdGa_move[3], NoAdGadiff, sbN_move, NosbNdiff);		       
                            fprintf(accpmoves,"NoGasites= %ld\t NoNsites= %ld\t NosbNsites= %ld\t NoAdGasites= %ld\t AdGa_N= %ld\t N_AdGa= %ld\t side_move_N= %ld\t side_move_Ga= %ld\t rejected_sb_move=%ld\t sbN_Island=%ld\t Ga3_sbN =%ld\n", NoGasites, NoNsites, NoAdGasites, NosbNsites, AdGa_N, N_AdGa, side_move_N, side_move_Ga, rejected_sb_move, sbN_Island, Ga3_sbN);			
	                    fclose(accpmoves);
                           }
                      }

       end = clock();
       cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
       double time = (double) (cpu_time_used/3600.00);  // will be time in hours
       if(rank==0){ printf("\nno. of times extflux is called = %d\n",cnt);
                    printf("program exceutition time = %lf hours\n",time);
                    printf("**** End of Program ****\n");}

      gsl_rng_free(rsel);
      gsl_rng_free(rdir);
      gsl_rng_free(rdep);
      gsl_rng_free(racc);
      gsl_rng_free(rflux); 
//    ReadOutput();
      MPI_Finalize();
      return 0;
}




