#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "variables.h"

#include <gsl/gsl_rng.h>

gsl_rng *rsel, *rdir, *racc, *rdep, *rflux;     // rsel : random no. for selecting a site for move;  rdep : random no. to select a site for deposition; rflux : to get either type of atom
                                                // rdir : random no. to select a direction to move the selected atom;  racc : random no. to decide whether to accept or reject a move.
                                                //
#include "allsubfunctions.h"
                                               
