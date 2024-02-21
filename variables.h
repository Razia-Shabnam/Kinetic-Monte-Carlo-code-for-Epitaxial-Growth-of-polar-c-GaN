/* Epi Ga atoms */


#define Energy_Gaad                     0.50      //when Epi Ga has Neighbor = 0
#define Energy_Ga1N                     0.66      //when Epi Ga has Neighbor = 1 (N)
#define Energy_Ga1AdGa                  0.55      //when Epi Ga has Neighbor = 1 (AdGa)
#define Energy_Ga2N                     1.45      //when Epi Ga has Neighbor = 2 (N)
#define Energy_Ga2AdGa                  0.60      //when Epi Ga has Neighbor = 2 (AdGa)       
#define Energy_Ga1N_1AdGa               0.76      //when Epi Ga has Neighbor = 2 (1N and 1 AdGa) 
#define Energy_Ga3N                     2.20      //when Epi Ga has Neighbor = 3 (N)
#define Energy_Ga3AdGa                  0.80      //when Epi Ga has Neighbor = 3 (AdGa)       
#define Energy_Ga2N_1AdGa               1.55      //when Epi Ga has Neighbor = 3 (2N and 1 AdGa)
#define Energy_Ga1N_2AdGa               0.86      //when Epi Ga has Neighbor = 3 (1N and 2 AdGa) 
                       

/*For Ad Ga atoms */

#define Energy_AdGa0Ga                  0.40      //when Ad Ga has Neighbor = 0   
#define Energy_AdGa1Ga                  0.50      //when Ad Ga has Neighbor = 1 (Ga)  
#define Energy_AdGa2Ga                  0.60      //when Ad Ga has Neighbor = 2 (Ga)
#define Energy_AdGa3Ga                  0.70      //when Ad Ga has Neighbor = 3 (Ga)

/*For on surface N atoms */
#define Energy_N0Ga1N                   0.65      //when N has Neighbor = 0
#define Energy_N0Ga2N                   0.90      //when N has Neighbor = 0
#define Energy_N0Ga3N                   1.16      //when N has Neighbor = 0
#define Energy_N1GaRot                  0.79      //when N has Neighbor = 1 (Ga - Rot) 
#define Energy_N1GaHop                  1.22      //when N has Neighbor = 1 (Ga - hop)   
#define Energy_N2Ga                     1.90      //when N has Neighbor = 2 (Ga)
#define Energy_N3Ga                     2.60      //when N has Neighbor = 3 (Ga)

/* For sub surafce N atoms */

#define Energy_NGa_Ga                  0.60      //for all sub surafce N atom combination.

#define Energy_fast 0.4

//********************************************************************************************************************************************

#define Nmax 72                                                                     // maximum height given to system
#define Nx 100
#define Ny 100
#define Nz 10                                                                       // out of Nmax, atoms are filled only upto Nz
#define nx Nx/2                                                                     // nx and ny are for randomly choosing the atoms
#define ny Ny                                                                       // nx and ny are for randomly choosing the atoms

#define Flux_Ratio 0.29                                                              //A variable Flux_Ratio to control Ga & N incoming atoms.
#define ML 20.0                                                                     //no. of ML required to de 1 ML = ((nx*ny)/2)

#define flux_rate 0.7                                                               //(ML/sim_time)
#define Deposition        (int) ((ML*((double)nx)*((double)ny))/2)                  // (nx*ny) represents 1 BL. Fog 1 ML, have to divide by 2

#define Sweeps  205167198L                                                          //10^11*exp(-Energy_fast/kBT)
#define Total_timesteps   (unsigned long int) ((Sweeps*Deposition)/flux_rate)       // total no. of steps to do the particular simulation
#define Nsteps            (unsigned long int) (Total_timesteps/Deposition)          // 1 atom will be deposited after every Nsteps.

#define kB 1.0                                                                      // (8.6173×10^−5) bond-strength b/w N-Ga 
#define Temp 0.06462975                                                             // 750 K 

#define tot_site 18


//*******************************************************************************************************************************************

int H[Nx][Ny], Box[Nx][Ny][Nmax];                         

double Prob, Prob_Gaad, Prob_Ga1N, Prob_Ga2N, Prob_Ga3N, Prob_Ga1N_1AdGa, Prob_Ga1N_2AdGa, Prob_Ga2N_1AdGa, Prob_Ga1AdGa, Prob_Ga2AdGa,Prob_Ga3AdGa;
double Prob_N0Ga3N, Prob_N0Ga2N, Prob_N0Ga1N, Prob_NRot, Prob_NHop, Prob_N2Ga, Prob_N3Ga;
double Prob_AdGa0Ga, Prob_AdGa1Ga, Prob_AdGa2Ga, Prob_AdGa3Ga, Prob_NGa_Ga;

//********************************************************************************************************************************************
int  cntAdGa, cntGa, cntN, Ga_Ga, Neighb, av_s, NN, NN_sb, N4N, start_point, Rdepsite ;
int NX_H8[tot_site], NY_H8[tot_site], NX_H4[tot_site], NY_H4[tot_site], av_site[tot_site], av_NX[tot_site], av_NY[tot_site], GaX[3], GaY[3], atom[tot_site], reject_site[tot_site]; 
long int N_AdGa, AdGa_N, NoGasites, NoNsites, NoAdGasites, NosbNsites, side_move_AdGa, side_move_Ga, side_move_N, rejected_sb_move, sbN_Island;
long int Ga_move[4][4], AdGa_move[4], N_move[4], sbN_move, NosbNdiff, Ga3_sbN;
int cnt, dcnt, Gaatoms, Natoms, AdGaatoms,  Attempt_N, Attempt_Ga;   
long int nNad , nNhop, nNrot, nN2Ga, nN3Ga, NoNdiff;
long int nGaad, nGa1AdGa, nGa2AdGa, nGa3AdGa, nGa1N1AdGa, nGa1N2AdGa, nGa2N1AdGa, nGa1N, nGa2N, nGa3N, NoGadiff;
long int nAdGa0Ga, nAdGa1Ga, nAdGa2Ga, nAdGa3Ga, NoAdGadiff;

int Nfiles, outNfiles;    // Nfiles are the no. of files to be created for output and FC.xyz files.
int size, rank;
//FILE *accpmoves, *cntout, *dev, *adout, *dep_atom, *unoccup;
//FILE *FCptr, *Iscnt, *IsSize;   // FCptr & Iscnt (island counting), IsSize (island size) would have to be defined here for the PP need. 
//int N;  // N is for applying loop over the total outNfiles for Post-Processing.
//int Label[Nx][Ny];  // for PP.

