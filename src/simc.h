//
//  simc.h
//  Inv_seq
//
//  Created by Fenix Huang on 4/28/17.
//  Copyright (c) 2017 Fenix Huang. All rights reserved.
//

#ifndef __Inv_seq__simc__
#define __Inv_seq__simc__

#include <stdio.h>


extern void Overlap(char *seqA, char *seqB, char *ref_sequence, int *struc_pat, FILE *out);

extern void Compare_Boltzmann(char *seqA, char *seqB, int *struc_pat, FILE *out);

extern void Inverse_fold_rate(char *struc, char *ref_seq, char *sequence_patten, int distance, int mul, FILE *out);

#endif /* defined(__Inv_seq__simc__) */
