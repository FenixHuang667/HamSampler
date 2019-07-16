//
//  pfunc.h
//  Bolt
//
//  Created by Fenix Huang on 5/27/16.
//  Copyright (c) 2016 Fenix Huang. All rights reserved.
//

#ifndef Bolt_pfunc_h

#include "stat.h"

#define MAXENG 100000

#define R GASCONST

extern int fold_sec (char *seq, int * ptable);

extern patten* struc_sampler(char *seq, int mul, int *pat);

#endif
