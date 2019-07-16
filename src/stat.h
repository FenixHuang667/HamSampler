//
//  stat.h
//  Inv_seq
//
//  Created by Fenix Huang on 8/28/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#ifndef __Inv_seq__stat__
#define __Inv_seq__stat__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

typedef struct patten patten;

struct patten
{
    char *pat;
    unsigned long cnt;
    patten *next;
};



#endif /* defined(__Inv_seq__stat__) */

extern void ratio (char *seq, int *rat);
//nucleotides ratio


extern patten *patten_detector(char *seq, int i, int j, int frequency, patten *List);

extern float clustering (char **set, int start, int end, int n);

extern int *SSampler(int l);

extern void sortList (patten *List);
