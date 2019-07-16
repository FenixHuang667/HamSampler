//
//  PF.h
//  Inv_seq
//
//  Created by Fenix Huang on 8/19/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#ifndef __Inv_seq__PF__
#define __Inv_seq__PF__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"
#include "stat.h"

#endif /* defined(__Inv_seq__PF__) */



extern double *Q[7];
extern paramT *P;

extern char Outputfilename[500];
extern char *Outfile;
extern int Interval;
extern int refold;
extern int compute_ratio;
extern int show_sequence;
extern int show_fold_structure;
extern int show_MI;
extern int show_energy;
extern int eval;


extern long double ** PFunction(int *ptable, int *order, int *patten);
//Compute the PF on a given structrue [ptabel]

extern char * Sample_seq (int *ptable, int *order, long double **Q, int *pat);
//Sample a sequence on a structure
//In condition that PF has been computed

extern long double EXP(int G);

extern void Sampler (int mul, int *ptable, int *pat);

extern void helix_e(void);

extern void Hull(int *ptable, int *pat);

extern void debug();

extern patten *seq_sampler_ham (char *structure, char *ref_sequence, char *sequence_patten, int distance, int mul, FILE *out);

extern patten *sample_compatible_seq(char *structure, int mul);

extern int* BPair_order (int *ptable);
extern int energy_eval(char *seq, int *ptable, int *order); 
