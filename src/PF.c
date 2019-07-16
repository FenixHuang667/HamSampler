//
//  PF.c
//  Inv_seq
//
//  Created by Fenix Huang on 8/19/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include <string.h>
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "energy_par.h"
#include "utils.h"
#include "stat.h"
#include "pfunc.h"
#include "fold.h"

#include "PF.h"


//#define EXP(G) exp(-1*(10*(G))/(R*(P->temperature+K0)))

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

#define PUBLIC
#define PRIVATE static

#define R GASCONST


int compute_ratio;
char Outputfilename[500];
char *Outfile;
int Interval, refold, show_sequence, show_fold_structure, show_energy, show_MI, eval;

FILE *fp;
PUBLIC paramT *P; //energy parameters

long double EXP(int G)
{
    return (long double) exp(-1*(10*(G))/(R*(P->temperature+K0)));
}


int* BPair_order (int *ptable)
{
    int i, n, j;
    int *order;
    
    n= ptable[0];
    order = (int *) space(sizeof(int)*n);
    
    j=1;
    for (i=1; i<=n; i++) {
        if (ptable[i] != 0 && i > ptable [i]) {
            order[j] = ptable[i];
            j++;
        }
    }
    order[0] = j-1;
    return order;
}
//Order the computational order of base pairs


int nest_bp (int *ptable, int i)
{
    int j, n;
    
    n=0;
    j=i+1;
    while (j<ptable[i]) {
        if (ptable[j]!=0) {
            j=ptable[j]+1;
            n++;
        } else {
            j++;
        }
    }
    return n;
}

int nest_up (int *ptable, int i)
{
    int j, n;
    
    n=0;
    j=i+1;
    while (j<ptable[i]) {
        if (ptable[j] == 0) {
            n++;
            j++;
        } else {
            j=ptable[j]+1;
        }
    }
    return n;
}

void type_to_bp(int type, int i, int j, int *s)
{
    switch (type) {
        case 1:
            s[i] = 2;
            s[j] = 3;
            break;
        case 2:
            s[i] = 3;
            s[j] = 2;
            break;
        case 3:
            s[i] = 3;
            s[j] = 4;
            break;
        case 4:
            s[i] = 4;
            s[j] = 3;
            break;
        case 5:
            s[i] = 1;
            s[j] = 4;
            break;
        case 6:
            s[i] = 4;
            s[j] = 1;
            break;
        default:
            s[i] = 0;
            s[j] = 0;
            break;
    }
}
void make_sequence(int *seq, char *res, int i, int j)
{
    int r;
    for (r=i; r<=j; r++) {
        switch (seq[r]) {
            case 1:
                res[r]='A';
                break;
            case 2:
                res[r]='C';
                break;
            case 3:
                res[r]='G';
                break;
            case 4:
                res[r]='U';
                break;
            default :
                res[r]='_';
                break;
        }
    }
}
//map the coded sequnce to char
//1:A 2:C 3:G 4:U

int undertermined(int *patten, int i, int j)
{
    int r, ud=0;
    for (r=i; r<=j; r++) {
        if (patten[r] == 0) ud++;
    }
    return ud;
}

int n_choose_k (int n, int k)
{
    int result = 1, i;
    if (n==0 && k==0) return 1;
    else if (n==0 && k>0) return 0;
    
    for (i=1; i<=k; i++) {
        result *= n - (k-i);
        result /= i;
    }
    return result;
}

int random_seq_hamming(int i, int j, int d, int *patten, int *ref_seq)
{
    int r, diff = 0, ud=0, m;
    if (j-i < 0) {
        if (d == 0) return 1;
        else return 0;
    }
    
    for (r = i; r<=j; r++) {
        if (patten[r] != 0 && patten[r] != ref_seq[r]) diff++;
    }
    
    if (d - diff < 0) return 0;
    
    for (r=i; r<=j; r++) {
        if (patten[r] == 0) ud++;
    }
    m =  n_choose_k(ud, d-diff);
    return m * pow(3, d-diff);
}

int random_seq_hamming_gap(int i, int j, int p, int q, int d, int *patten, int *ref_seq) //[i..j] [p..q]
{
    int r, diff = 0, ud = 0, m;
    if (j-i < 0 && q-p < 0) {
        if (d==0) return 1;
        else return 0;
    }

    for (r = i; r<=j; r++) {
        if (patten[r] != 0 && patten[r] != ref_seq[r]) diff++;
    }
    for (r = p; r<=q; r++) {
        if (patten[r] != 0 && patten[r] != ref_seq[r]) diff++;
    }
    if (d - diff < 0) return 0;
    
    for (r=i; r<=j; r++) {
        if (patten[r] == 0) ud++;
    }
    
    for (r=p; r<=q; r++) {
        if (patten[r] == 0) ud++;
    }
    m = n_choose_k(ud, d-diff); 
    return n_choose_k(ud, d-diff) * pow(3, d-diff);
}


long double **PFunction(int *ptable, int *order, int * patten)
{
    int bp=order[0], nbp, ubp, ud;
    int i, j, energy, r, p, q, l, next_bp;
    int type, type2;
    int s[9];
    char str[12];
    double long Qtemp;
    double long **Q;
    
    Q = (long double **) space(sizeof(long double *)*7);
    
    str[0]='_';
    for (i=1; i<=6; i++) {
        Q[i] = (long double *) space(sizeof(long double)*(bp+5));
        for (j=0; j<=bp; j++)
            Q[i][j]=0;
    }

    
    // initial the Q arrays.
    
    for (r=1; r<=bp; r++) {
        i=order[r];
        nbp = nest_bp(ptable, i);
        j=ptable[i];
        
        if (nbp == 0) { // hairpin loop
            if (j-i-1==4) { //tetra loop
                s[0] = 6;
                str[7]=0;
                
                for (type = 1; type <= 6; type++) {
                    type_to_bp(type, 1, 6, s);
                    if (patten[i]!=0 && patten[i]!= s[1]) continue;
                    if (patten[j]!=0 && patten[j]!= s[6]) continue;
                    
                    for (s[2]=1; s[2]<=4; s[2]++) {
                        for (s[3]=1; s[3]<=4; s[3]++) {
                            for (s[4]=1; s[4]<=4; s[4]++) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                                    if (patten[i+2]!=0 && patten[i+2]!= s[3]) continue;
                                    if (patten[i+3]!=0 && patten[i+3]!= s[4]) continue;
                                    if (patten[i+4]!=0 && patten[i+4]!= s[5]) continue;
                                    
                                    make_sequence(s, str, 1, 6);
                                    energy = HairpinE(j-i-1, type, s[2], s[5], str+1);
                                    Q[type][r]+=EXP(energy);
                                }
                            }
                        }
                    }
                }
            } else if (j-i-1==3) {
                s[0] = 5;
                str[6] = 0;
                        
                for (type = 1; type <= 6; type++) {
                    type_to_bp(type, 1, 5, s);
                    if (patten[i]!=0 && patten[i]!= s[1]) continue;
                    if (patten[j]!=0 && patten[j]!= s[5]) continue;
                    
                    for (s[2]=1; s[2]<=4; s[2]++) {
                        for (s[3]=1; s[3]<=4; s[3]++) {
                            for (s[4]=1; s[4]<=4; s[4]++) {
                                if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                                if (patten[i+2]!=0 && patten[i+2]!= s[3]) continue;
                                if (patten[i+3]!=0 && patten[i+3]!= s[4]) continue;
                                
                                make_sequence(s, str, 1, 5);
                                energy = HairpinE(j-i-1, type, s[2], s[4], str+1);
                                Q[type][r]+=EXP(energy);
                            }
                        }
                    }
                }
            } else if (j-i-1>4){ // normal hairpins
                s[0]=4;
                str[5] = 0;
                
                for (type = 1; type <= 6; type++) {
                    type_to_bp(type, 1, 4, s);
                    if (patten[i]!=0 && patten[i]!= s[1]) continue;
                    if (patten[j]!=0 && patten[j]!= s[4]) continue;
                    
                    for (s[2]=1; s[2]<=4; s[2]++) {
                        for (s[3]=1; s[3]<=4; s[3]++) {
                            
                            if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                            if (patten[j-1]!=0 && patten[j-1]!= s[3]) continue;
                            
                            ud = undertermined(patten, i+2, j-2);
                            
                            make_sequence(s, str, 1, 4);
                            energy = HairpinE(j-i-1, type, s[2], s[3], str+1);
                            Q[type][r]+=pow(4, ud) * EXP(energy);
                        }
                    }
                }
            }
            
            //for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
            
        } else if (nbp == 1) { // interior loop
            p=i+1;
            while (ptable[p]==0) p++;
            for (l=1; l<=order[0]; l++) {
                if (order[l] == p) {
                    next_bp = l;
                    break;
                }
            }
            q = ptable[p];
            
//          printf("%d %d\n", i, l);
            for (type = 1; type<=6; type++) {
                type_to_bp(type, 1, 4, s);
                if (patten[i]!=0 && patten[i]!= s[1]) continue;
                if (patten[j]!=0 && patten[j]!= s[4]) continue;
                
                for (s[2]=1; s[2]<=4; s[2]++) {
                    for (s[3]=1; s[3]<=4; s[3]++) {
                        if (patten[p]!=0 && patten[p]!= s[2]) continue;
                        if (patten[q]!=0 && patten[q]!= s[3]) continue;
                        
                        type2 = BP_pair[s[2]][s[3]];
                        if (!type2) continue;
                        type2 = rtype[type2];
                            
                        if (p-i-1 == 0 && j-q-1 == 0) {
                            energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[3], s[1], s[4]);
                            Q[type][r] += EXP(energy) * Q[rtype[type2]][next_bp];
                                
                        } else if (p-i-1 == 1 && j-q-1 == 0) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                            
                                energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[3], s[5], s[4]);
                                Q[type][r] += EXP(energy) * Q[rtype[type2]][next_bp];
                            }
                        } else if (p-i-1 >= 2 && j-q-1 == 0) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                    if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                    ud = undertermined(patten, i+2, p-2);
                                    
                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[3], s[6], s[4]);
                                    Q[type][r] += pow(4, ud) * EXP(energy) * Q[rtype[type2]][next_bp];
                                }
                            }
                        } else if (p-i-1 == 0 && j-q-1 == 1) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                if (patten[j-1]!=0 && patten[j-1]!= s[5]) continue;
                                
                                energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[5], s[1], s[5]);
                                Q[type][r] += EXP(energy) * Q[rtype[type2]][next_bp];
                            }
                                
                        } else if (p-i-1 == 1 && j-q-1 == 1) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                    if (patten[j-1]!=0 && patten[j-1]!= s[6]) continue;
                                    
                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[6], s[5], s[6]);
                                    Q[type][r] += EXP(energy) * Q[rtype[type2]][next_bp];
                                }
                            }
                                
                        } else if (p-i-1 >= 2 && j-q-1 == 1) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    for (s[7]=1; s[7]<=4; s[7]++) {
                                        if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                        if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                        if (patten[j-1]!=0 && patten[j-1]!= s[7]) continue;
                                        
                                        energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[7], s[6], s[7]);
                                        ud = undertermined(patten, i+2, p-2);
                                        
                                        Q[type][r] += pow(4,ud) * EXP(energy) * Q[rtype[type2]][next_bp];
                                    }
                                }
                            }
                        } else if (p-i-1 == 0 && j-q-1 >= 2) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    if (patten[q+1]!=0 && patten[q+1]!= s[5]) continue;
                                    if (patten[j-1]!=0 && patten[j-1]!= s[6]) continue;
                                    
                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[6], s[1], s[5]);
                                    ud = undertermined(patten, q+2, j-2);
                                    
                                    Q[type][r] += pow(4,ud) * EXP(energy) * Q[rtype[type2]][next_bp];
                                }
                            }
                                
                        } else if (p-i-1 == 1 && j-q-1 >= 2) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    for (s[7]=1; s[7]<=4; s[7]++) {
                                        if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                        if (patten[q+1]!=0 && patten[q+1]!= s[6]) continue;
                                        if (patten[j-1]!=0 && patten[j-1]!= s[7]) continue;
                                        
                                        energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[7], s[5], s[6]);
                                        ud = undertermined(patten, q+2, j-2);
                                        
                                        Q[type][r] += pow(4,ud) * EXP(energy) * Q[rtype[type2]][next_bp];
                                    }
                                }
                            }
                        } else if (p-i-1 >= 2 && j-q-1 >= 2) {
                            for (s[5]=1; s[5]<=4; s[5]++) {
                                for (s[6]=1; s[6]<=4; s[6]++) {
                                    for (s[7]=1; s[7]<=4; s[7]++) {
                                        for (s[8]=1; s[8]<=4; s[8]++) {
                                            if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                            if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                            if (patten[q+1]!=0 && patten[q+1]!= s[7]) continue;
                                            if (patten[j-1]!=0 && patten[j-1]!= s[8]) continue;
                                            
                                            energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[8], s[6], s[7]);
                                            ud = undertermined(patten, i+2, p-2) + undertermined(patten, q+2, j-2);
                                            
                                            Q[type][r] += pow(4, ud) * EXP(energy) * Q[rtype[type2]][next_bp];
                                                
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
    
       // for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
            
            
        } else if (nbp > 1) { //multi loop
            for (type = 1; type<=6; type++) {
                type_to_bp(type, 1, 2, s);
                
                if (patten[i]!=0 && patten[i]!= s[1]) continue;
                if (patten[j]!=0 && patten[j]!= s[2]) continue;
                
                energy = P->MLclosing + P->MLintern[type];
                
                ud=0;
                p=i+1;
                while (p<j) {
                    if (ptable[p]==0) {
                        if (patten[p] == 0) ud++;
                        p++;
                    } else {
                        p = ptable[p]+1;
                    }
                }
                
                p=i+1;
                Q[type][r] = 1;
                while (p<j) {
                    while (ptable[p]==0 && p<j) p++;
                    q=ptable[p];
                    Qtemp = 0;
                    
                    for (l=1; l<=order[0]; l++) {
                        if (order[l] == p) {
                            next_bp = l;
                            break;
                        }
                    }
                    
                    for (s[3]=1; s[3]<=4; s[3]++){
                        for (s[4]=1; s[4]<=4; s[4]++) {
                            type2 = BP_pair[s[3]][s[4]];
                            if (patten[p]!=0 && patten[p]!= s[3]) continue;
                            if (patten[q]!=0 && patten[q]!= s[4]) continue;
                            
                            if (!type2) continue;
        
                            Qtemp += EXP (P->MLintern[type2]) * Q[type2][next_bp];
                        }
                    }
                        
                    Q[type][r] *= Qtemp;
                        
                    p=q+1;
                    while (ptable[p]==0 && p<j) p++;
                        
                }
                
                
                Q[type][r] *= pow(4, ud) * EXP(energy);
            }
            
            //for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
        }
    }
    return Q;
}

double long PH_hamming_multiloop(int i, int j, int d_remain, double long ***Q_ham, int *ptable, int *order, int *patten, int *ref_seq, int is_mul)
{
    int uh, index, d_unpaired, d_sub, G, r, s, type, next_bp, l, d_c, max_unpaired, max_sub;
    double long Q = 0, Qtemp;
    int base[12];
    
    if (d_remain < 0) {
        return 0;
    }
    if (i > j && d_remain ==0) {
        return 1;
    }
    else if (i>j && d_remain >0) {
        return 0;
    }
    
    index = i;
    while (ptable[index]==0) index++;
    if (index > j) {
        uh = random_seq_hamming(i, j, d_remain, patten, ref_seq);
        if (is_mul) G = (j-i+1) * P->MLbase; else G = 0;
        Q = uh *EXP(G);
        
        // printf("%Le\n", Q);
        return Q;
    }
    
    next_bp = 0;
    r = index; s = ptable[r];
    for (l=1; l<=order[0]; l++) {
        if (order[l] == r) {
            next_bp = l;
            break;
        }
    }
    
    max_unpaired = index-i;
    if (d_remain < max_unpaired) max_unpaired = d_remain;
    for (d_unpaired = 0; d_unpaired <= max_unpaired; d_unpaired++) {
        uh = random_seq_hamming(i, index-1, d_unpaired, patten, ref_seq);
        G = (r-i) * P->MLbase;
        max_sub = s-r+1;
        if (d_remain-d_unpaired < max_sub) max_sub = d_remain-d_unpaired;
        for (d_sub = 0; d_sub <= max_sub; d_sub++) {
            if (d_remain-d_sub-d_unpaired>0 && s+1>j) continue;
            Q_ham[0][next_bp][d_remain-d_sub-d_unpaired] = PH_hamming_multiloop(s+1, j, d_remain-d_sub-d_unpaired, Q_ham, ptable, order, patten, ref_seq, is_mul);
            for (type = 1; type <=6; type++) {
                type_to_bp(type, 1, 2, base);
                
                if (patten[r]!=0 && patten[r]!= base[1]) continue;
                if (patten[s]!=0 && patten[s]!= base[2]) continue;
                
                d_c = 0;
                
                if (base[1] != ref_seq[r]) d_c++;
                if (base[2] != ref_seq[s]) d_c++;
                
                if (d_sub < d_c) continue;
                Qtemp = uh;
                if (is_mul) Qtemp *= EXP(G) * EXP (P->MLintern[type]);
                Qtemp *= Q_ham[type][next_bp][d_sub];
                Qtemp *= Q_ham[0][next_bp][d_remain-d_sub-d_unpaired];
                
                Q += Qtemp;
                //if (i==9 && j==66) printf("PF %d: %Le\n", d_remain, Q);
            }
        }
    }
    return Q;
}

double long PH_hamming_multiloop_Z(int i, int j, int d_remain, double long ***Q_ham, double long ***Z_ham, int *ptable, int *order, int *patten, int *ref_seq, int is_mul)
{
    int uh, index, d_unpaired, d_sub, G, r, s, type, next_bp, l, d_c, max_unpaired, max_sub;
    double long Q = 0, Qtemp, Q1, Q2, Q3, Z1, Z2, Z3;
    int base[12];
    
    if (d_remain < 0) {
        return 0;
    }
    if (i > j && d_remain ==0) {
        return 0;
    }
    else if (i>j && d_remain >0) {
        return 0;
    }
    
    index = i;
    while (ptable[index]==0) index++;
    if (index > j) {
        uh = random_seq_hamming(i, j, d_remain, patten, ref_seq);
        if (is_mul) G = (j-i+1) * P->MLbase; else G=0;
        Q = uh * (((float) G)/100) * EXP(G);
        
        // printf("%Le\n", Q);
        return Q;
    }
    
    next_bp = 0;
    r = index; s = ptable[r];
    for (l=1; l<=order[0]; l++) {
        if (order[l] == r) {
            next_bp = l;
            break;
        }
    }
    
    max_unpaired = index-i;
    if (d_remain < max_unpaired) max_unpaired = d_remain;
    for (d_unpaired = 0; d_unpaired <= max_unpaired; d_unpaired++) {
        uh = random_seq_hamming(i, index-1, d_unpaired, patten, ref_seq);
        G = (r-i) * P->MLbase;
        max_sub = s-r+1;
        if (d_remain-d_unpaired < max_sub) max_sub = d_remain-d_unpaired;
        for (d_sub = 0; d_sub <= max_sub; d_sub++) {
            if (d_remain-d_sub-d_unpaired>0 && s+1>j) continue;
            // Q_ham[0][next_bp][d_remain-d_sub-d_unpaired] = PH_hamming_multiloop(s+1, j, d_remain-d_sub-d_unpaired, Q_ham, ptable, order, patten, ref_seq, is_mul);
            Z_ham[0][next_bp][d_remain-d_sub-d_unpaired] = PH_hamming_multiloop_Z(s+1, j, d_remain-d_sub-d_unpaired, Q_ham, Z_ham, ptable, order, patten, ref_seq, is_mul);
            for (type = 1; type <=6; type++) {
                type_to_bp(type, 1, 2, base);
                
                if (patten[r]!=0 && patten[r]!= base[1]) continue;
                if (patten[s]!=0 && patten[s]!= base[2]) continue;
                
                d_c = 0;
                
                if (base[1] != ref_seq[r]) d_c++;
                if (base[2] != ref_seq[s]) d_c++;
                
                if (d_sub < d_c) continue;
                
                if (is_mul) Q1 = EXP(G); else Q1 = 1;
                if (is_mul) Q2 = Q_ham[type][next_bp][d_sub] * EXP(P->MLintern[type]); else Q2 = Q_ham[type][next_bp][d_sub];
                Q3 = Q_ham[0][next_bp][d_remain-d_sub-d_unpaired];
                
                if (is_mul) Z1 = (((float) G) / 100) * EXP(G); else Z1 = 0;
                if (is_mul) Z2 = Q_ham[type][next_bp][d_sub] * EXP(P->MLintern[type]) * (((float) P->MLintern[type]) / 100) + EXP(P->MLintern[type]) * Z_ham[type][next_bp][d_sub];
                else Z2 = Z_ham[type][next_bp][d_sub];
                Z3 = Z_ham[0][next_bp][d_remain-d_sub-d_unpaired];
                
                Q += uh * (Z1 * Q2 * Q3 + Q1 * Z2 * Q3 + Q1 * Q2 * Z3);
                
                //if (i==9 && j==66) printf("PF %d: %Le\n", d_remain, Q);
            }
        }
    }
    return Q;
}

void PFunction_hamming(int *ptable, int *order, int * patten, int *ref_seq, int distance, long double ***Q_ham, long double ***Z_ham, long double ***X_ham)
{
    int bp=order[0], nbp, ubp, ud;
    int i, j, energy, n, r, p, q, l, d, next_bp, ref_i, h, uh, h_interior, d_c;
    int type, type2;
    int s[50];
    int ref_segment[50];
    
    char str[50], segment[50];
    double long Qtemp;
    
    n = l = ref_seq[0];
    
    str[0]='_';
    // initial the Q arrays.
    
    for (d = 0; d <= distance; d++) {
        for (r=1; r<=bp; r++) {
            i=order[r];
            nbp = nest_bp(ptable, i);
            j=ptable[i];
            
            if (nbp == 0) { // hairpin loop
                if (j-i-1==4) { //tetra loop
                    s[0] = 6;
                    str[7]=0;
                    
                    ref_segment[0] = 6; segment[7] = 0;
                    for (ref_i = 1; ref_i <= j-i+1; ref_i++) ref_segment[ref_i] = ref_seq[i+ref_i-1];
                    
                    for (type = 1; type <= 6; type++) {
                        type_to_bp(type, 1, 6, s);
                        if (patten[i]!=0 && patten[i]!= s[1]) continue;
                        if (patten[j]!=0 && patten[j]!= s[6]) continue;
                        
                        for (s[2]=1; s[2]<=4; s[2]++) {
                            for (s[3]=1; s[3]<=4; s[3]++) {
                                for (s[4]=1; s[4]<=4; s[4]++) {
                                    for (s[5]=1; s[5]<=4; s[5]++) {
                                        if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                                        if (patten[i+2]!=0 && patten[i+2]!= s[3]) continue;
                                        if (patten[i+3]!=0 && patten[i+3]!= s[4]) continue;
                                        if (patten[i+4]!=0 && patten[i+4]!= s[5]) continue;
                                        
                                        make_sequence(s, str, 1, 6);
                                        make_sequence(ref_segment, segment, 1, 6);
                                        h = hamming(segment+1, str+1);
                                        
                                        if (h == d) {
                                            energy = HairpinE(j-i-1, type, s[2], s[5], str+1);
                                            Q_ham[type][r][h] += EXP(energy);
                                            Z_ham[type][r][h] += (((float) energy)/100) * EXP(energy);
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (j-i-1==3) {
                    s[0] = 5;
                    str[6] = 0;
                    
                    ref_segment[0] = 5; segment[6] = 0;
                    for (ref_i = 1; ref_i <= j-i+1; ref_i++) ref_segment[ref_i] = ref_seq[i-1+ref_i];
                    
                    for (type = 1; type <= 6; type++) {
                        type_to_bp(type, 1, 5, s);
                        if (patten[i]!=0 && patten[i]!= s[1]) continue;
                        if (patten[j]!=0 && patten[j]!= s[5]) continue;
                        
                        for (s[2]=1; s[2]<=4; s[2]++) {
                            for (s[3]=1; s[3]<=4; s[3]++) {
                                for (s[4]=1; s[4]<=4; s[4]++) {
                                    if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                                    if (patten[i+2]!=0 && patten[i+2]!= s[3]) continue;
                                    if (patten[i+3]!=0 && patten[i+3]!= s[4]) continue;
                                    
                                    make_sequence(s, str, 1, 5);
                                    make_sequence(ref_segment, segment, 1, 5);
                                    h = hamming(segment+1, str+1);
                                    if (h == d) {
                                        energy = HairpinE(j-i-1, type, s[2], s[4], str+1);
                                        Q_ham[type][r][h]+=EXP(energy);
                                        Z_ham[type][r][h] += (((float) energy)/100) * EXP(energy);
                                    }
                                }
                            }
                        }
                    }
                } else if (j-i-1>4){ // normal hairpins
                    s[0]=4;
                    str[5] = 0;
                    
                    ref_segment[0] = 4; segment[5] = 0;
                    ref_segment[1] = ref_seq[i];
                    ref_segment[2] = ref_seq[i+1];
                    ref_segment[3] = ref_seq[j-1];
                    ref_segment[4] = ref_seq[j];
                    
                    
                    for (type = 1; type <= 6; type++) {
                        type_to_bp(type, 1, 4, s);
                        if (patten[i]!=0 && patten[i]!= s[1]) continue;
                        if (patten[j]!=0 && patten[j]!= s[4]) continue;
                        
                        for (s[2]=1; s[2]<=4; s[2]++) {
                            for (s[3]=1; s[3]<=4; s[3]++) {
                                
                                if (patten[i+1]!=0 && patten[i+1]!= s[2]) continue;
                                if (patten[j-1]!=0 && patten[j-1]!= s[3]) continue;
                                
                                make_sequence(s, str, 1, 4);
                                make_sequence(ref_segment, segment, 1, 4);
                                h = hamming(segment+1, str+1);
                                
                                if (h <= d) {
                                    uh = random_seq_hamming(i+2, j-2, d-h, patten, ref_seq);
                                    energy = HairpinE(j-i-1, type, s[2], s[3], str+1);
                                    Q_ham[type][r][d] += uh * EXP(energy);
                                    Z_ham[type][r][d] += uh * (((float) energy)/100) * EXP(energy);
                                }
                            }
                        }
                    }
                }
                
                //for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
                
            } else if (nbp == 1) { // interior loop
                p=i+1;
                while (ptable[p]==0) p++;
                for (l=1; l<=order[0]; l++) {
                    if (order[l] == p) {
                        next_bp = l;
                        break;
                    }
                }
                q = ptable[p];
                
                ref_segment[1] = ref_seq[i];
                ref_segment[2] = ref_seq[p];
                ref_segment[3] = ref_seq[q];
                ref_segment[4] = ref_seq[j];
                
                //          printf("%d %d\n", i, l);
                for (type = 1; type<=6; type++) {
                    type_to_bp(type, 1, 4, s);
                    if (patten[i]!=0 && patten[i]!= s[1]) continue;
                    if (patten[j]!=0 && patten[j]!= s[4]) continue;
                    
                    for (s[2]=1; s[2]<=4; s[2]++) {
                        for (s[3]=1; s[3]<=4; s[3]++) {
                            if (patten[p]!=0 && patten[p]!= s[2]) continue;
                            if (patten[q]!=0 && patten[q]!= s[3]) continue;
                            
                            type2 = BP_pair[s[2]][s[3]];
                            if (!type2) continue;
                            type2 = rtype[type2];
                            
                            d_c = 0;
                            if (s[2] != ref_segment[2]) d_c++;
                            if (s[3] != ref_segment[3]) d_c++;
                            
                            if (p-i-1 == 0 && j-q-1 == 0) {
                                s[0] = 4; ref_segment[0] = 4;
                                str[5] = 0; segment[5] = 0;
                                make_sequence(s, str, 1, 4);
                                make_sequence(ref_segment, segment, 1, 4);
                                h = hamming(segment+1, str+1)-d_c;
                                if (h <= d) {
                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[3], s[1], s[4]);
                                    Q_ham[type][r][d] += EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h];
                                    
                                    Z_ham[type][r][d] += EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h] + (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h];
                                }
                            } else if (p-i-1 == 1 && j-q-1 == 0) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    
                                    if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                    ref_segment[5] = ref_seq[i+1];
                                    
                                    s[0] = 5; ref_segment[0] = 5;
                                    str[6] = 0; segment[6] = 0;
                                    make_sequence(s, str, 1, 5);
                                    make_sequence(ref_segment, segment, 1, 5);
                                    h = hamming(segment+1, str+1);
                                    
                                    if (h <= d) {
                                        energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[3], s[5], s[4]);
                                        Q_ham[type][r][d] += EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                        
                                        Z_ham[type][r][d] += EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h+d_c] + (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                    }
                                }
                            } else if (p-i-1 >= 2 && j-q-1 == 0) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                        if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                        
                                        ref_segment[5] = ref_seq[i+1];
                                        ref_segment[6] = ref_seq[p-1];
                                        
                                        s[0] = 6; ref_segment[0] = 6;
                                        str[7] = 0; segment[7] = 0;
                                        make_sequence(s, str, 1, 6);
                                        make_sequence(ref_segment, segment, 1, 6);
                                        h = hamming(segment+1, str+1);
                                        
                                        if (h <= d) {
                                            for (h_interior = 0; h_interior <= d-h; h_interior++) {
                                                uh = random_seq_hamming(i+2, p-2, h_interior, patten, ref_seq);
                                        
                                                energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[3], s[6], s[4]);
                                                Q_ham[type][r][d] += uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                
                                                Z_ham[type][r][d] += uh * EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h-h_interior+d_c] + uh * (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                            }
                                        }
                                    }
                                }
                            } else if (p-i-1 == 0 && j-q-1 == 1) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    if (patten[j-1]!=0 && patten[j-1]!= s[5]) continue;
                                    ref_segment[5] = ref_seq[j-1];
                                    
                                    s[0] = 5; ref_segment[0] = 5;
                                    str[6] = 0; segment[6] = 0;
                                    make_sequence(s, str, 1, 5);
                                    make_sequence(ref_segment, segment, 1, 5);
                                    h = hamming(segment+1, str+1);
                                    
                                    if (h <= d) {
                                        energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[5], s[1], s[5]);
                                        Q_ham[type][r][d] += EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                        
                                        Z_ham[type][r][d] += EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h+d_c] + (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                    }
                                }
                                
                            } else if (p-i-1 == 1 && j-q-1 == 1) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                        if (patten[j-1]!=0 && patten[j-1]!= s[6]) continue;
                                        
                                        ref_segment[5] = ref_seq[i+1];
                                        ref_segment[6] = ref_seq[j-1];
                                        
                                        s[0] = 6; ref_segment[0] = 6;
                                        str[7] = 0; segment[7] = 0;
                                        make_sequence(s, str, 1, 6);
                                        make_sequence(ref_segment, segment, 1, 6);
                                        h = hamming(segment+1, str+1);
                                        
                                        if (h <= d) {
                                            energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[6], s[5], s[6]);
                                            Q_ham[type][r][d] += EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                            
                                            Z_ham[type][r][d] += EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h+d_c] + (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h+d_c];
                                        }
                                    }
                                }
                                
                            } else if (p-i-1 >= 2 && j-q-1 == 1) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        for (s[7]=1; s[7]<=4; s[7]++) {
                                            if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                            if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                            if (patten[j-1]!=0 && patten[j-1]!= s[7]) continue;
                                            
                                            ref_segment[5] = ref_seq[i+1];
                                            ref_segment[6] = ref_seq[p-1];
                                            ref_segment[7] = ref_seq[j-1];
                                            
                                            s[0] = 7; ref_segment[0] = 7;
                                            str[8] = 0; segment[8] = 0;
                                            make_sequence(s, str, 1, 7);
                                            make_sequence(ref_segment, segment, 1, 7);
                                            h = hamming(segment+1, str+1);
                                            
                                            if (h <= d) {
                                                for (h_interior = 0; h_interior <= d-h; h_interior++) {
                                                    uh = random_seq_hamming(i+2, p-2, h_interior, patten, ref_seq);
                                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[7], s[6], s[7]);
                                            
                                                    Q_ham[type][r][d] += uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                    
                                                    Z_ham[type][r][d] += uh * EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h-h_interior+d_c] + (((float) energy)/100) * uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                }
                                            }
                                        }
                                    }
                                }
                            } else if (p-i-1 == 0 && j-q-1 >= 2) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        if (patten[q+1]!=0 && patten[q+1]!= s[5]) continue;
                                        if (patten[j-1]!=0 && patten[j-1]!= s[6]) continue;
                                        
                                        ref_segment[5] = ref_seq[q+1];
                                        ref_segment[6] = ref_seq[j-1];
                                        
                                        s[0] = 6; ref_segment[0] = 6;
                                        str[7] = 0; segment[7] = 0;
                                        make_sequence(s, str, 1, 6);
                                        make_sequence(ref_segment, segment, 1, 6);
                                        h = hamming(segment+1, str+1);
                                        
                                        if (h <= d) {
                                            for (h_interior = 0; h_interior <= d-h; h_interior++) {
                                                uh = random_seq_hamming(q+2, j-2, h_interior, patten, ref_seq);
                                                energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[2], s[6], s[1], s[5]);
                                                Q_ham[type][r][d] += uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                
                                                Z_ham[type][r][d] += uh * EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h-h_interior+d_c] + uh * (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                            }
                                        }
                                    }
                                }
                                
                            } else if (p-i-1 == 1 && j-q-1 >= 2) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        for (s[7]=1; s[7]<=4; s[7]++) {
                                            if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                            if (patten[q+1]!=0 && patten[q+1]!= s[6]) continue;
                                            if (patten[j-1]!=0 && patten[j-1]!= s[7]) continue;
                                            
                                            ref_segment[5] = ref_seq[i+1];
                                            ref_segment[6] = ref_seq[q+1];
                                            ref_segment[7] = ref_seq[j-1];
                                            
                                            s[0] = 7; ref_segment[0] = 7;
                                            str[8] = 0; segment[8] = 0;
                                            make_sequence(s, str, 1, 7);
                                            make_sequence(ref_segment, segment, 1, 7);
                                            h = hamming(segment+1, str+1);
                                            
                                            if (h <= d) {
                                                for (h_interior = 0; h_interior <= d-h; h_interior++) {
                                                    uh = random_seq_hamming(q+2, j-2, h_interior, patten, ref_seq);
                                                    energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[7], s[5], s[6]);
                                                    Q_ham[type][r][d] += uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                    
                                                    Z_ham[type][r][d] += uh * EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h-h_interior+d_c] + uh * (((float) energy)/100) *EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                }
                                            }
                                        }
                                    }
                                }
                            } else if (p-i-1 >= 2 && j-q-1 >= 2) {
                                for (s[5]=1; s[5]<=4; s[5]++) {
                                    for (s[6]=1; s[6]<=4; s[6]++) {
                                        for (s[7]=1; s[7]<=4; s[7]++) {
                                            for (s[8]=1; s[8]<=4; s[8]++) {
                                                if (patten[i+1]!=0 && patten[i+1]!= s[5]) continue;
                                                if (patten[p-1]!=0 && patten[p-1]!= s[6]) continue;
                                                if (patten[q+1]!=0 && patten[q+1]!= s[7]) continue;
                                                if (patten[j-1]!=0 && patten[j-1]!= s[8]) continue;
                                                
                                                ref_segment[5] = ref_seq[i+1];
                                                ref_segment[6] = ref_seq[p-1];
                                                ref_segment[7] = ref_seq[q+1];
                                                ref_segment[8] = ref_seq[j-1];
                                                
                                                s[0] = 8; ref_segment[0] = 8;
                                                str[9] = 0; segment[9] = 0;
                                                make_sequence(s, str, 1, 8);
                                                make_sequence(ref_segment, segment, 1, 8);
                                                h = hamming(segment+1, str+1);

                                                if (h <= d) {
                                                    for (h_interior = 0; h_interior <= d-h; h_interior++) {
                                                        uh = random_seq_hamming_gap(i+2, p-2, q+2, j-2, h_interior, patten, ref_seq);
                                                        energy = LoopEnergy(p-i-1, j-q-1, type, type2, s[5], s[8], s[6], s[7]);
                                                        
                                                        Q_ham[type][r][d] += uh * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                        
                                                        Z_ham[type][r][d] += uh * EXP(energy) * Z_ham[rtype[type2]][next_bp][d-h-h_interior+d_c] + uh * (((float) energy)/100) * EXP(energy) * Q_ham[rtype[type2]][next_bp][d-h-h_interior+d_c];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                // for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
                
                
            } else if (nbp > 1) { //multi loop
                for (type = 1; type<=6; type++) {
                    type_to_bp(type, 1, 2, s);
                    
                    if (patten[i]!=0 && patten[i]!= s[1]) continue;
                    if (patten[j]!=0 && patten[j]!= s[2]) continue;
                    
                    d_c = 0;
                    if (s[1] != ref_seq[i]) d_c ++;
                    if (s[2] != ref_seq[j]) d_c ++;
                    
                    energy = P->MLclosing + P->MLintern[type];
                    Q_ham[0][r][d-d_c] = PH_hamming_multiloop(i+1, j-1, d-d_c, Q_ham, ptable, order, patten, ref_seq, 1);
                    Z_ham[0][r][d-d_c] = PH_hamming_multiloop_Z(i+1, j-1, d-d_c, Q_ham, Z_ham, ptable, order, patten, ref_seq, 1);
                    Q_ham[type][r][d] = EXP(energy) * Q_ham[0][r][d-d_c];
                    Z_ham[type][r][d] = EXP(energy) * Z_ham[0][r][d-d_c] + (((float) energy)/100) * EXP(energy) * Q_ham[0][r][d-d_c];
                }
                
                //for (l=1; l<=6; l++) printf("%d: %f\n", i, Q[l][r]);
            }
            
            //if (1) {
            //    printf("%d: \n", d);
            //    for (l=1; l<=6; l++) printf("%d: %Le\n", r, Q_ham[l][r][d]);
            //    printf("\n");
            //}
        }
        
        //exterior loop
        
        Q_ham[0][bp+1][d] = PH_hamming_multiloop(1, n, d, Q_ham, ptable, order, patten, ref_seq, 0);
        Z_ham[0][bp+1][d] = PH_hamming_multiloop_Z(1, n, d, Q_ham, Z_ham, ptable, order, patten, ref_seq, 0);
        //printf("%d: %Le\n", d, Q_ham[1][bp+1][d]);
    }
}


int random_nucleotide (void)
{
    float seed;
    seed=(float) ((drand48()) * 4);
    if (seed < 1) {
        return 1;
    } else if (seed < 2 && seed >=1) {
        return 2;
    } else if (seed < 3 && seed >=2) {
        return 3;
    } else if (seed < 4 && seed >=3) {
        return 4;
    }
    return 0;
}

int random_nucleotide_hamming (int avoid)
{
    float seed;
    int res;
    seed=(float) ((drand48()) * 3);
    if (seed < 1) {
        res = 1;
    } else if (seed < 2 && seed >=1) {
        res = 2;
    } else if (seed < 3 && seed >=2) {
        res = 3;
    } else {
        res = 0;
    }
    if (res >= avoid) res++;
    return res;
}

char * Sample_seq (int *ptable, int *order, long double **Q, int *pat) // In condition that the PF has been computed
{
    char * seq;
    int bp, i, j, p, q, r, l, nbp, ubp, n, next_bp, ud, k;
    float seed;
    int type, type2;
    double long Qtemp, Qv, Qsum, *Qmul, Qdetermined;
    int energy;
    short done_sample;
    int *s;
    int y[9];
    char str[12];
    
    n=ptable[0];
    seq = (char *) calloc (n+3, sizeof(char));
    s = (int *) calloc (n+3, sizeof(int));

    bp = order[0];
    for (l=1; l<=n; l++) s[l] = pat[l];
    s[0] = n;
    
    
    p=1;
    while (p<=n) {
        if (ptable[p] == 0) {
            if (pat[p]==0) s[p] = random_nucleotide();
            p++;
        } else {
            q = ptable[p];
            Qsum = 0;
            for (l=1; l<=order[0]; l++) {
                if (order[l] == p) {
                    next_bp = l;
                    break;
                }
            }
            for (l=1; l<=6; l++) {
                Qsum+=Q[l][next_bp];
            }
            seed = (float) (1-drand48());
            Qv = seed * Qsum;
            Qtemp = 0;
            
            for (type=1; type<=6; type++) {
                type_to_bp(type, 1, 2, y);
                if (pat[p]!=0 && pat[p]!=y[1]) continue;
                if (pat[q]!=0 && pat[q]!=y[2]) continue;
                
                Qtemp += Q[type][next_bp];

                if (Qv <=Qtemp) {
                    s[p] = y[1];
                    s[q] = y[2];
                    break;
                }
            }
            p = ptable[p]+1;
        }
    }
    //sample exterior loop
    
    for (r=bp; r>=1; r--) {
        i = order[r];
        j = ptable[i];
        type = BP_pair[s[i]][s[j]];
        
        nbp = nest_bp(ptable, i);
        seed = (float) (1-drand48());
        Qv = seed * Q[type][r];
        
        if (nbp == 0) {
            if (j-i-1 == 4) {
                done_sample = 0;
                Qtemp = 0;
                y[1] = s[i];
                y[6] = s[j];
                str[7] = 0;
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                                if (pat[i+2]!=0 && pat[i+2]!= y[3]) continue;
                                if (pat[i+3]!=0 && pat[i+3]!= y[4]) continue;
                                if (pat[i+4]!=0 && pat[i+4]!= y[5]) continue;
                                
                                make_sequence(y, str, 1,6);
                                energy = HairpinE(j-i-1, type, y[2], y[5], str+1);
                                Qtemp += EXP(energy);
                                if (Qv <=Qtemp) {
                                    for (l=2; l<=5; l++) s[i+l-1] = y[l];
                                    done_sample = 1;
                                }
                            }
                        }
                    }
                }
                
            } else if (j-i-1 == 3) {
                done_sample = 0;
                Qtemp = 0;
                y[1] = s[i];
                y[5] = s[j];
                str[6]=0;
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                            if (pat[i+2]!=0 && pat[i+2]!= y[3]) continue;
                            if (pat[i+3]!=0 && pat[i+3]!= y[4]) continue;
                            
                            make_sequence(y, str, 1, 5);
                            energy = HairpinE(j-i-1, type, y[2], y[4], str+1);
                            Qtemp += EXP(energy);
                            
                            if (Qv <=Qtemp) {
                                for (l=2; l<=4; l++) s[i+l-1] = y[l];
                                done_sample = 1;
                            }
                        }
                    }
                }
            } else if (j-i-1 > 4) {
                done_sample = 0;
                Qtemp = 0;
                y[1] = s[i];
                y[4] = s[j];
                str[5]=0;
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                        if (pat[j-1]!=0 && pat[j-1]!= y[3]) continue;
                        
                        make_sequence(y, str, 1, 4);
                        energy = HairpinE(j-i-1, type, y[2], y[3], str+1);
                        
                        ud = undertermined(pat, i+2, j-2);
                        
                        Qtemp += pow(4, ud) * EXP(energy);
                        
                        if (Qv <=Qtemp) {
                            done_sample = 1;
                            s[i+1] = y[2];
                            s[j-1] = y[3];
                            for (l=i+2; l<=j-2 && pat[l]==0; l++) {
                                s[l] = random_nucleotide();
                            }
                        }
                        
                    }
                }
            }
        } else if (nbp == 1) {
            p=i+1;
            while (ptable[p]==0) p++;
            for (l=1; l<=order[0]; l++) {
                if (order[l] == p) {
                    next_bp = l;
                    break;
                }
            }
            q = ptable[p];
            done_sample = 0;
            Qtemp = 0;
            
            for (type2=1; type2<=6 && !done_sample; type2++) {
                type_to_bp(type2, 1, 2, y);
                if (pat[p]!=0 && pat[p]!=y[1]) continue;
                if (pat[q]!=0 && pat[q]!=y[2]) continue;
                
                if (p-i-1 == 0 && j-q-1 == 0) {
                    energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[1], y[2], s[i], s[j]);
                    Qtemp += EXP(energy) * Q[type2][next_bp];
                    
                    if (Qv <=Qtemp) {
                        done_sample = 1;
                        s[p] = y[1];
                        s[q] = y[2];
                    }
                } else if (p-i-1 == 1 && j-q-1 == 0) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                            
                        energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[2], y[3], s[j]);
                        Qtemp += EXP(energy) * Q[type2][next_bp];
                            
                        if (Qv <=Qtemp) {
                            done_sample = 1;
                            s[p] = y[1];
                            s[q] = y[2];
                            s[i+1] = y[3];
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 == 0) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                            if (pat[p-1]!=0 && pat[p-1]!= y[4]) continue;
                        
                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[2], y[4], s[j]);
                            ud = undertermined(pat, i+2, p-2);
                            Qtemp += pow(4,ud) * EXP(energy) * Q[type2][next_bp];
                        
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[p] = y[1];
                                s[q] = y[2];
                                s[i+1] = y[3];
                                s[p-1] = y[4];
                                for (l=i+2; l<=p-2 && pat[l]==0; l++) {
                                    if (pat[l]==0) s[l] = random_nucleotide();
                                }
                            }
                        }
                    }
                } else if (p-i-1 == 0 && j-q-1 == 1) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        if (pat[j-1]!=0 && pat[j-1]!= y[3]) continue;
                        
                        energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[1], y[3], s[i], y[3]);
                        Qtemp += EXP(energy) * Q[type2][next_bp];
                        if (Qv <=Qtemp) {
                            done_sample = 1;
                            s[p] = y[1];
                            s[q] = y[2];
                            s[j-1] = y[3];
                        }
                    }
                } else if (p-i-1 == 1 && j-q-1 == 1) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                            if (pat[j-1]!=0 && pat[j-1]!= y[4]) continue;
                            
                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[4], y[3], y[4]);
                            Qtemp += EXP(energy) * Q[type2][next_bp];
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[p] = y[1];
                                s[q] = y[2];
                                s[i+1] = y[3];
                                s[j-1] = y[4];
                            }
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 == 1) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                                if (pat[p-1]!=0 && pat[p-1]!= y[4]) continue;
                                if (pat[j-1]!=0 && pat[j-1]!= y[5]) continue;
                                
                                energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[5], y[4], y[5]);
                                ud = undertermined(pat, i+2, p-2);
                                Qtemp += pow(4,ud) * EXP(energy) * Q[type2][next_bp];
                                if (Qv <=Qtemp) {
                                    done_sample = 1;
                                    s[p] = y[1];
                                    s[q] = y[2];
                                    s[i+1] = y[3];
                                    s[p-1] = y[4];
                                    s[j-1] = y[5];
                                    for (l=i+2; l<=p-2 && pat[l]==0; l++) {
                                        if (pat[l]==0) s[l] = random_nucleotide();
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 == 0 && j-q-1 >= 2) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            if (pat[q+1]!=0 && pat[q+1]!= y[3]) continue;
                            if (pat[j-1]!=0 && pat[j-1]!= y[4]) continue;
                            
                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[1], y[4], s[i], y[3]);
                            ud = undertermined(pat, q+2, j-2);
                            Qtemp += pow(4,ud) * EXP(energy) * Q[type2][next_bp];
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[p] = y[1];
                                s[q] = y[2];
                                s[q+1] = y[3];
                                s[j-1] = y[4];
                                for (l=q+2; l<=j-2 && pat[l]==0; l++) {
                                    if (pat[l]==0) s[l] = random_nucleotide();
                                }
                            }
                            
                        }
                    }
                } else if (p-i-1 == 1 && j-q-1 >= 2) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                                if (pat[q+1]!=0 && pat[q+1]!= y[4]) continue;
                                if (pat[j-1]!=0 && pat[j-1]!= y[5]) continue;
                                
                                energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[5], y[3], y[4]);
                                ud = undertermined(pat, q+2, j-2);
                                Qtemp += pow(4,ud) * EXP(energy) * Q[type2][next_bp];
                                if (Qv <=Qtemp) {
                                    done_sample = 1;
                                    s[p] = y[1];
                                    s[q] = y[2];
                                    s[i+1] = y[3];
                                    s[q+1] = y[4];
                                    s[j-1] = y[5];
                                    for (l=q+2; l<=j-2 && pat[l]==0; l++) {
                                        if (pat[l]==0) s[l] = random_nucleotide();
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 >= 2) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                                for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                                    if (pat[i+1]!=0 && pat[i+1]!= y[3]) continue;
                                    if (pat[p-1]!=0 && pat[p-1]!= y[4]) continue;
                                    if (pat[q+1]!=0 && pat[q+1]!= y[5]) continue;
                                    if (pat[j-1]!=0 && pat[j-1]!= y[6]) continue;
                                    
                                    ud = undertermined(pat, i+2, p-2) + undertermined(pat, q+2, j-2);
                                    energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[3], y[6], y[4], y[5]);
                                    Qtemp += pow(4, ud) * EXP(energy) * Q[type2][next_bp];
                                    
                                    if (Qv <=Qtemp) {
                                        s[p] = y[1];
                                        s[q] = y[2];
                                        s[i+1] = y[3];
                                        s[p-1] = y[4];
                                        s[q+1] = y[5];
                                        s[j-1] = y[6];
                                        for (l=i+2; l<=p-2 && pat[l]==0; l++) {
                                            if (pat[l]==0) s[l] = random_nucleotide();
                                        }
                                        for (l=q+2; l<=j-2 && pat[l]==0; l++) {
                                            if (pat[l]==0) s[l] = random_nucleotide();
                                        }
                                        done_sample = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            if (s[p] == 0) {
                printf("Error\n");
            }
        } else if (nbp > 1) {
            ud=0;
            p=i+1;
            while (p<j) {
                if (ptable[p]==0) {
                    if (pat[p] == 0) {
                        ud++;
                        s[p] = random_nucleotide();
                    }
                    p++;
                } else {
                    p = ptable[p]+1;
                }
            }
            
            energy = P->MLclosing + P->MLintern[type] + P->MLbase * ud;
            
            Qmul = (long double*) space(sizeof(long double) *(nbp+2));
            for (l=0;l<=nbp;l++) Qmul[l] = 0;
            p=i+1;
            l=1;
            while (p<j) {
                if (ptable[p]==0) {
                    p++;
                } else {
                    q=ptable[p];
                    for (k=1; k<=order[0]; k++) {
                        if (order[k]==p) {
                            next_bp=k;
                            break;
                        }
                    }
                    for (type2=1; type2<=6; type2++) {
                        Qmul[l] += EXP(P->MLintern[type2]) * Q[type2][next_bp];
                    }
                    p=q+1;
                    l++;
                }
            }
            
            p=i+1;
            l=1;
            Qdetermined = 1;
            while (p<j) {
                if (ptable[p] == 0){
                    p++;
                } else {
                    q = ptable[p];
                    for (k=1; k<=order[0]; k++) {
                        if (order[k]==p) {
                            next_bp=k;
                            break;
                        }
                    }
                    
                    Qsum =1;
                    for (k=l+1; k<=nbp; k++) Qsum *= Qmul[k];
                    Qtemp = 0;
                    for (type2 = 1; type2 <= 6; type2++) {
                        type_to_bp(type2, 1, 2, y);
                        if (pat[p]!=0 && pat[p]!=y[1]) continue;
                        if (pat[q]!=0 && pat[q]!=y[1]) continue;
                        Qtemp += pow(4,ud) * Qdetermined * EXP(energy) * EXP(P->MLintern[type2]) * Q[type2][next_bp] * Qsum;
                        
                        if (Qv <=Qtemp) {
                            s[p]=y[1];
                            s[q]=y[2];
                            Qdetermined *= EXP(P->MLintern[type2]) * Q[type2][next_bp];
                            seed = (float) (1-drand48());
                            Qv = seed * Qdetermined * pow(4,ud) * EXP(energy) * Qsum;
                            break;
                        }
                    }
                    p=q+1;
                    l++;
                }
                
            }
        }
    }
//    for (l=1; l<=n; l++) printf("%d ", s[l]);
//    printf("\n");
    
    make_sequence(s, seq, 0, n);
    seq[n+1] = 0;
    
    free(s);
    

    return seq;
}

void sample_random_seq_segment_hamming(int i, int j, int d_remain, int *sample_seq, int *pat, int *ref_seq)
{
    int k, diff, d, uh, uh1, uh0;
    float seed, Q;
    
    if (j<i) return;
    diff = 0;
    for (k = i; k<=j; k++) {
        if (pat[k]!=0 && ref_seq[k] != pat[k]) diff++;
    }
    if (diff > d_remain) {
        printf("Error, too many differences between reference sequence and patten seqeunce! \n");
        exit(0);
    }
    d = d_remain - diff;
    for (k = i; k<=j; k++) {
        if (pat[k] != 0) {
            sample_seq[k] = pat[k];
        } else {
            uh = random_seq_hamming(k, j, d, pat, ref_seq);
            uh0 = random_seq_hamming(k+1, j, d, pat, ref_seq);
            uh1 = 3*random_seq_hamming(k+1, j, d-1, pat, ref_seq);
            
            if (uh0+uh1!=uh) {
                printf("Error! %d %d %d\n");
            }
            
            seed = (float) (1-drand48());
            Q = seed * uh;
            if (Q <= uh0) {
                sample_seq[k] = ref_seq[k];
            } else {
                sample_seq[k] = random_nucleotide_hamming(ref_seq[k]);
                d--;
            }
        }
    }
    
}

void sample_random_seq_segment_gap_hamming(int i, int j, int r, int s, int d_remain, int *sample_seq, int *pat, int *ref_seq)
{
    int ref_seq_segment[80], segment[80], patten_segment[80], k;
    if (j<i && s<r) return;
    for (k=1; k<=j-i+1; k++) {
        ref_seq_segment[k] = ref_seq[i+k-1];
        patten_segment[k] = pat[i+k-1];
    }
    for (k=j-i+2; k<= j-i+1 + s-r+1; k++) {
        ref_seq_segment[k] = ref_seq[r+k-1 - (j-i+1)];
        patten_segment[k] = pat[r+k-1 - (j-i+1)];
    }
    sample_random_seq_segment_hamming(1, j-i+1 + s-r+1, d_remain, segment, patten_segment, ref_seq_segment);
    
    for (k=1; k<=j-i+1; k++) {
        sample_seq[i+k-1] = segment[k];
    }
    
    for (k=j-i+2; k<= j-i+1 + s-r+1; k++) {
        sample_seq[r+k-1 - (j-i+1)] = segment[k];
    }
}


void Sample_multiloop_hamming(int i, int j, int d_remain, int current_bp, int *remain_distance_sub, int *sample_seq, int *ptable, int *order, long double ***Q_ham, int *pat, int *ref_seq, int is_mul)
{
    double long Qsum = 0, Qtemp, Qpro, Q_multi, Qdebug;
    float seed;
    int energy, uh, index, d_unpaired, d_sub, G, r, s, type, next_bp, l, d_c, done_sample, max_unpaired, max_sub;
    int base[12];
    
    if (d_remain < 0) {
        printf("Error, negative hamming distance! [%d, %d]:%d\n", i, j ,d_remain);
        return;
    }
    if (i > j) {
        if (d_remain!=0) printf("Error!");
        return;
    }
    
    index = i;
    while (ptable[index]==0) index++;
    if (index > j) {
        sample_random_seq_segment_hamming(i, j, d_remain, sample_seq, pat, ref_seq);
        return;
    }
    
    next_bp = 0;
    r = index; s = ptable[r];
    for (l=1; l<=order[0]; l++) {
        if (order[l] == r) {
            next_bp = l;
            break;
        }
    }
    
    Qtemp = Q_ham[0][current_bp][d_remain];
    // Qdebug = PH_hamming_multiloop(i, j, d_remain, Q_ham, ptable, order, pat, ref_seq, is_mul);

    
    //if (!is_mul) {
    //    Qtemp = Q_ham[0][current_bp][d_remain];
    //} else {
        // Qtemp = PH_hamming_multiloop(i, j, d_remain, Q_ham, ptable, order, pat, ref_seq, 1);
    //    Qtemp = Q_ham[0][current_bp][d_remain];
    //}
    seed = (float) (1-drand48());
    Qpro = Qtemp *seed;
    // printf("Total: %Le %Le\n", Qtemp, Qpro);
    
    done_sample = 0;
    max_unpaired = index-i;
    if (d_remain < max_unpaired) max_unpaired = d_remain;
    for (d_unpaired = 0; d_unpaired <= max_unpaired; d_unpaired++) {
        uh = random_seq_hamming(i, index-1, d_unpaired, pat, ref_seq);
        G = (r-i) * P->MLbase;
        max_sub = s-r+1;
        if (d_remain-d_unpaired<max_sub) max_sub = d_remain-d_unpaired;
        for (d_sub = 0; d_sub <= max_sub; d_sub++) {
            if (d_remain-d_sub-d_unpaired >0 && s+1>j) continue;
            Q_multi = Q_ham[0][next_bp][d_remain-d_sub-d_unpaired];
            for (type = 1; type <=6; type++) {
                type_to_bp(type, 1, 2, base);
                
                if (pat[r]!=0 && pat[r]!= base[1]) continue;
                if (pat[s]!=0 && pat[s]!= base[2]) continue;
                
                d_c = 0;
                
                if (base[1] != ref_seq[r]) d_c++;
                if (base[2] != ref_seq[s]) d_c++;
                
                if (d_sub < d_c) continue;
                
                Qtemp = uh;
                if (is_mul) Qtemp *= EXP(G) * EXP (P->MLintern[type]);
                Qtemp *= Q_ham[type][next_bp][d_sub];
                // Qtemp *= PH_hamming_multiloop(s+1, j, d_remain-d_sub-d_unpaired, Q_ham, ptable, order, pat, ref_seq, is_mul);
                Qtemp *= Q_multi;
                
                Qsum += Qtemp;
                
                // if (i==9 && j==66) printf("Sample: %d %Le\n", d_remain, Qsum);
                
                if (Qpro <= Qsum && !done_sample) {
                    sample_random_seq_segment_hamming(i, index-1, d_unpaired, sample_seq, pat, ref_seq);
                    sample_seq[r] = base[1];
                    sample_seq[s] = base[2];
                    
                    remain_distance_sub[next_bp] = d_sub;
                    
                    Sample_multiloop_hamming(s+1, j, d_remain-d_sub-d_unpaired, next_bp, remain_distance_sub, sample_seq, ptable, order, Q_ham, pat, ref_seq, is_mul);
                    
                    //if (is_mul) printf("%d %d %d\n", d_unpaired, d_sub, d_remain-d_sub-d_unpaired);
                    done_sample = 1;
                    return;
                }
            }
            //printf("%d %d %d %Le %Le\n",d_unpaired, d_sub, d_remain-d_sub-d_unpaired, Qpro, Qsum);
        }
    }
    printf("Multiloop sample error!");
}

char * Sample_seq_ham (int *ptable, int *ref_seq, int *order, int *remain_distance_sub, long double ***Q_ham, int *pat, int distance) // In condition that the PF_ham has been computed
{
    char * seq;
    int bp, i, j, p, q, r, l, nbp, n, next_bp, ud, k, d, h, uh, d_c, h_interior, debug_i;
    float seed;
    int type, type2;
    double long Qtemp, Qv, Qsum, Qdetermined, Qtotal;
    int energy;
    short done_sample;
    int *s;
    int y[20];
    char str[50], segment[50];
    int ref_segment[50], gap[50];
    
    n=ptable[0];
    seq = (char *) calloc (n+3, sizeof(char));
    s = (int *) calloc (n+3, sizeof(int));
    
    bp = order[0];
    for (l=1; l<=n; l++) s[l] = pat[l];
    s[0] = n;
    
    
    Qtotal = Q_ham[1][order[0]+1][distance];
    
    Sample_multiloop_hamming(1, n, distance, order[0]+1, remain_distance_sub, s, ptable, order, Q_ham, pat, ref_seq, 0);
    // s[2] = 3; s[73] = 2; remain_distance_sub[23]=5;

    p=1;
    for (r=bp; r>=1; r--) {
        
        //----------
        
        /*
        make_sequence(s, seq, 0, n);
        seq[n+1] = 0;
        printf("%s\n", seq+1);
        for (debug_i=1; debug_i<=order[0]; debug_i++) printf("%d ", remain_distance_sub[debug_i]);
        printf("\n");
        */
        
        //----------
        
        i = order[r];
        j = ptable[i];
        type = BP_pair[s[i]][s[j]];
        
        nbp = nest_bp(ptable, i);
        d = remain_distance_sub[r];
        seed = (float) (1-drand48());
        Qv = seed * Q_ham[type][r][d];
        // printf("Probability: %f, %Le, %Le\n", seed, Qv, Q_ham[type][r][d]);
        
        if (nbp == 0) {
            if (j-i-1 == 4) {
                done_sample = 0;
                Qtemp = 0;
                y[0] = 6;
                y[1] = s[i];
                y[6] = s[j];
                str[7] = 0;
                
                ref_segment[0] = 6; segment[7] = 0;
                ref_segment[1] = ref_seq[i];
                ref_segment[2] = ref_seq[i+1];
                ref_segment[3] = ref_seq[i+2];
                ref_segment[4] = ref_seq[i+3];
                ref_segment[5] = ref_seq[i+4];
                ref_segment[6] = ref_seq[i+5];
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                                if (pat[i+2]!=0 && pat[i+2]!= y[3]) continue;
                                if (pat[i+3]!=0 && pat[i+3]!= y[4]) continue;
                                if (pat[i+4]!=0 && pat[i+4]!= y[5]) continue;
                                
                                make_sequence(y, str, 1, 6);
                                make_sequence(ref_segment, segment, 1, 6);
                                h = hamming(segment+1, str+1);
                                
                                if (h == d) {
                                    energy = HairpinE(j-i-1, type, y[2], y[5], str+1);
                                    Qtemp += EXP(energy);
                                    if (Qv <=Qtemp) {
                                        for (l=2; l<=5; l++) s[i+l-1] = y[l];
                                        done_sample = 1;
                                    }
                                }
                            }
                        }
                    }
                }
                
            } else if (j-i-1 == 3) {
                done_sample = 0;
                Qtemp = 0;
                y[0] = 5;
                y[1] = s[i];
                y[5] = s[j];
                str[6]=0;
                
                ref_segment[0] = 5; segment[6] = 0;
                ref_segment[1] = ref_seq[i];
                ref_segment[2] = ref_seq[i+1];
                ref_segment[3] = ref_seq[i+2];
                ref_segment[4] = ref_seq[i+3];
                ref_segment[5] = ref_seq[i+4];
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        for (y[4]=1; y[4]<=4 && !done_sample; y[4]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                            if (pat[i+2]!=0 && pat[i+2]!= y[3]) continue;
                            if (pat[i+3]!=0 && pat[i+3]!= y[4]) continue;
                            
                            make_sequence(y, str, 1, 5);
                            make_sequence(ref_segment, segment, 1, 5);
                            h = hamming(segment+1, str+1);
                            if (h == d) {
                                energy = HairpinE(j-i-1, type, y[2], y[4], str+1);
                                Qtemp += EXP(energy);
                            
                                if (Qv <=Qtemp) {
                                    for (l=2; l<=4; l++) s[i+l-1] = y[l];
                                    done_sample = 1;
                                }
                            }
                        }
                    }
                }
            } else if (j-i-1 > 4) {
                done_sample = 0;
                Qtemp = 0;
                y[0] = 4;
                y[1] = s[i];
                y[4] = s[j];
                str[5]=0;
                
                ref_segment[0] = 4; segment[5] = 0;
                ref_segment[1] = ref_seq[i];
                ref_segment[2] = ref_seq[i+1];
                ref_segment[3] = ref_seq[j-1];
                ref_segment[4] = ref_seq[j];
                
                for (y[2]=1; y[2]<=4 && !done_sample; y[2]++) {
                    for (y[3]=1; y[3]<=4 && !done_sample; y[3]++) {
                        if (pat[i+1]!=0 && pat[i+1]!= y[2]) continue;
                        if (pat[j-1]!=0 && pat[j-1]!= y[3]) continue;
                        
                        make_sequence(y, str, 1, 4);
                        make_sequence(ref_segment, segment, 1, 4);
                        h = hamming(segment+1, str+1);
                        if (h <= d) {
                            energy = HairpinE(j-i-1, type, y[2], y[3], str+1);
                        
                            uh = random_seq_hamming(i+2, j-2, d-h, pat, ref_seq);
                        
                            Qtemp += uh * EXP(energy);
                        
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[i+1] = y[2];
                                s[j-1] = y[3];
                                
                                sample_random_seq_segment_hamming(i+2, j-2, d-h, s, pat, ref_seq);
                                
                            }
                        }
                        
                    }
                }
            }
        } else if (nbp == 1) {
            p=i+1;
            next_bp = 0;
            while (ptable[p]==0) p++;
            for (l=1; l<=order[0]; l++) {
                if (order[l] == p) {
                    next_bp = l;
                    break;
                }
            }
            q = ptable[p];
            done_sample = 0;
            Qtemp = 0;
            
            y[1] = s[i];
            y[4] = s[j];
            ref_segment[1] = ref_seq[i];
            ref_segment[2] = ref_seq[p];
            ref_segment[3] = ref_seq[q];
            ref_segment[4] = ref_seq[j];
            
            for (type2=1; type2<=6 && !done_sample; type2++) {
                type_to_bp(type2, 2, 3, y);
                if (pat[p]!=0 && pat[p]!=y[2]) continue;
                if (pat[q]!=0 && pat[q]!=y[3]) continue;
                
                d_c = 0;
                if (y[2] != ref_segment[2]) d_c ++;
                if (y[3] != ref_segment[3]) d_c ++;
                
                if (p-i-1 == 0 && j-q-1 == 0) {
                    y[0] = 4; ref_segment[0] = 4;
                    str[5] = 0; segment[5] = 0;
                    make_sequence(y, str, 1, 4);
                    make_sequence(ref_segment, segment, 1, 4);
                    h = hamming(segment+1, str+1);
                    if (h <= d) {
                        energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[2], y[3], y[1], y[4]);
                        Qtemp += EXP(energy) * Q_ham[type2][next_bp][d-h+d_c];
                        
                        if (Qv <= Qtemp) {
                            done_sample = 1;
                            s[p] = y[2];
                            s[q] = y[3];
                            remain_distance_sub[next_bp] = d - h + d_c;
                        }
                    }
                } else if (p-i-1 == 1 && j-q-1 == 0) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                        
                        y[0] = 5; ref_segment[0] = 5;
                        str[6] = 0; segment[6] = 0;
                        
                        ref_segment[5] = ref_seq[i+1];
                        
                        make_sequence(y, str, 1, 5);
                        make_sequence(ref_segment, segment, 1, 5);
                        h = hamming(segment+1, str+1);
                        if (h <= d) {
                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[3], y[5], y[4]);
                            Qtemp += EXP(energy) * Q_ham[type2][next_bp][d-h+d_c];
                        
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[p] = y[2];
                                s[q] = y[3];
                                s[i+1] = y[5];
                                remain_distance_sub[next_bp] = d - h + d_c;
                            }
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 == 0) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                            if (pat[p-1]!=0 && pat[p-1]!= y[6]) continue;
                            
                            ref_segment[5] = ref_seq[i+1];
                            ref_segment[6] = ref_seq[p-1];
                            
                            y[0] = 6; ref_segment[0] = 6;
                            str[7] = 0; segment[7] = 0;
                            make_sequence(y, str, 1, 6);
                            make_sequence(ref_segment, segment, 1, 6);
                            h = hamming(segment+1, str+1);
                            
                            if (h<=d) {
                                for (h_interior = 0; h_interior <= d-h && !done_sample; h_interior++) {
                                    uh = random_seq_hamming(i+2, p-2, h_interior, pat, ref_seq);
                                    energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[3], y[6], y[4]);
                                    Qtemp += uh * EXP(energy) * Q_ham[type2][next_bp][d-h-h_interior+d_c];
                            
                                    if (Qv <=Qtemp) {
                                        done_sample = 1;
                                        s[p] = y[2];
                                        s[q] = y[3];
                                        s[i+1] = y[5];
                                        s[p-1] = y[6];
                                        sample_random_seq_segment_hamming(i+2, p-2, h_interior, s, pat, ref_seq);
                                        remain_distance_sub[next_bp] = d-h-h_interior+d_c;
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 == 0 && j-q-1 == 1) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        if (pat[j-1]!=0 && pat[j-1]!= y[5]) continue;
                        ref_segment[5] = ref_seq[j-1];
                        
                        y[0] = 5; ref_segment[0] = 5;
                        str[6] = 0; segment[6] = 0;
                        
                        ref_segment[5] = ref_seq[j-1];
                        
                        make_sequence(y, str, 1, 5);
                        make_sequence(ref_segment, segment, 1, 5);
                        
                        h = hamming(segment+1, str+1);
                        
                        if (h <= d) {
                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[2], y[5], y[1], y[5]);
                            Qtemp += EXP(energy) * Q_ham[type2][next_bp][d-h+d_c];
                            if (Qv <=Qtemp) {
                                done_sample = 1;
                                s[p] = y[2];
                                s[q] = y[3];
                                s[j-1] = y[5];
                                remain_distance_sub[next_bp] = d-h+d_c;
                            }
                        }
                    }
                } else if (p-i-1 == 1 && j-q-1 == 1) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                            if (pat[j-1]!=0 && pat[j-1]!= y[6]) continue;
                            
                            ref_segment[5] = ref_seq[i+1];
                            ref_segment[6] = ref_seq[j-1];
                            
                            y[0] = 6; ref_segment[0] = 6;
                            str[7] = 0; segment[7] = 0;
                            make_sequence(y, str, 1, 6);
                            make_sequence(ref_segment, segment, 1, 6);
                            h = hamming(segment+1, str+1);
                            
                            if (h <= d) {
                                energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[6], y[5], y[6]);
                                Qtemp += EXP(energy) * Q_ham[type2][next_bp][d-h+d_c];
                                if (Qv <=Qtemp) {
                                    done_sample = 1;
                                    s[p] = y[2];
                                    s[q] = y[3];
                                    s[i+1] = y[5];
                                    s[j-1] = y[6];
                                    remain_distance_sub[next_bp] = d-h+d_c;
                                }
                            }
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 == 1) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            for (y[7]=1; y[7]<=4 && !done_sample; y[7]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                                if (pat[p-1]!=0 && pat[p-1]!= y[6]) continue;
                                if (pat[j-1]!=0 && pat[j-1]!= y[7]) continue;
                                
                                ref_segment[5] = ref_seq[i+1];
                                ref_segment[6] = ref_seq[p-1];
                                ref_segment[7] = ref_seq[j-1];
                                
                                y[0] = 7; ref_segment[0] = 7;
                                str[8] = 0; segment[8] = 0;
                                make_sequence(y, str, 1, 7);
                                make_sequence(ref_segment, segment, 1, 7);
                                h = hamming(segment+1, str+1);
                                
                                if (h <= d) {
                                    for (h_interior = 0; h_interior <= d-h && !done_sample; h_interior++) {
                                        uh = random_seq_hamming(i+2, p-2, h_interior, pat, ref_seq);
                                        
                                        energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[7], y[6], y[7]);
    
                                        Qtemp += uh * EXP(energy) * Q_ham[type2][next_bp][d-h-h_interior+d_c];
                                        if (Qv <=Qtemp) {
                                            done_sample = 1;
                                            s[p] = y[2];
                                            s[q] = y[3];
                                            s[i+1] = y[5];
                                            s[p-1] = y[6];
                                            s[j-1] = y[7];
                                            sample_random_seq_segment_hamming(i+2, p-2, h_interior, s, pat, ref_seq);
                                            remain_distance_sub[next_bp] = d-h-h_interior+d_c;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 == 0 && j-q-1 >= 2) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            if (pat[q+1]!=0 && pat[q+1]!= y[5]) continue;
                            if (pat[j-1]!=0 && pat[j-1]!= y[6]) continue;
                            
                            ref_segment[5] = ref_seq[q+1];
                            ref_segment[6] = ref_seq[j-1];
                            
                            y[0] = 6; ref_segment[0] = 6;
                            str[7] = 0; segment[7] = 0;
                            make_sequence(y, str, 1, 6);
                            make_sequence(ref_segment, segment, 1, 6);
                            h = hamming(segment+1, str+1);
                            
                            if (h <= d) {
                                for (h_interior = 0; h_interior <= d-h && !done_sample; h_interior++) {
                                    uh = random_seq_hamming(q+2, j-2, h_interior, pat, ref_seq);
                                
                                    energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[2], y[6], y[1], y[5]);
                                    Qtemp += uh * EXP(energy) * Q_ham[type2][next_bp][d-h-h_interior+d_c];
                                    if (Qv <=Qtemp) {
                                        done_sample = 1;
                                        s[p] = y[2];
                                        s[q] = y[3];
                                        s[q+1] = y[5];
                                        s[j-1] = y[6];
                                        sample_random_seq_segment_hamming(q+2, j-2, h_interior, s, pat, ref_seq);
                                        remain_distance_sub[next_bp] = d-h-h_interior+d_c;
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 == 1 && j-q-1 >= 2) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            for (y[7]=1; y[7]<=4 && !done_sample; y[7]++) {
                                if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                                if (pat[q+1]!=0 && pat[q+1]!= y[6]) continue;
                                if (pat[j-1]!=0 && pat[j-1]!= y[7]) continue;
                                
                                y[0] = 7; ref_segment[0] = 7;
                                str[8] = 0; segment[8] = 0;
                                
                                ref_segment[5] = ref_seq[i+1];
                                ref_segment[6] = ref_seq[q+1];
                                ref_segment[7] = ref_seq[j-1];
                                
                                make_sequence(y, str, 1, 7);
                                make_sequence(ref_segment, segment, 1, 7);
                                h = hamming(segment+1, str+1);
                                
                                if (h <= d) {
                                    for (h_interior = 0; h_interior <= d-h && !done_sample; h_interior++) {
                                        uh = random_seq_hamming(q+2, j-2, h_interior, pat, ref_seq);
                                
                                        energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[7], y[5], y[6]);
                                        Qtemp += uh * EXP(energy) * Q_ham[type2][next_bp][d-h-h_interior+d_c];
                                        
                                        if (Qv <=Qtemp) {
                                            done_sample = 1;
                                            s[p] = y[2];
                                            s[q] = y[3];
                                            s[i+1] = y[5];
                                            s[q+1] = y[6];
                                            s[j-1] = y[7];
                                            sample_random_seq_segment_hamming(q+2, j-2, h_interior, s, pat, ref_seq);
                                            remain_distance_sub[next_bp] = d-h-h_interior+d_c;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (p-i-1 >= 2 && j-q-1 >= 2) {
                    for (y[5]=1; y[5]<=4 && !done_sample; y[5]++) {
                        for (y[6]=1; y[6]<=4 && !done_sample; y[6]++) {
                            for (y[7]=1; y[7]<=4 && !done_sample; y[7]++) {
                                for (y[8]=1; y[8]<=4 && !done_sample; y[8]++) {
                                    if (pat[i+1]!=0 && pat[i+1]!= y[5]) continue;
                                    if (pat[p-1]!=0 && pat[p-1]!= y[6]) continue;
                                    if (pat[q+1]!=0 && pat[q+1]!= y[7]) continue;
                                    if (pat[j-1]!=0 && pat[j-1]!= y[8]) continue;
                                    
                                    ref_segment[5] = ref_seq[i+1];
                                    ref_segment[6] = ref_seq[p-1];
                                    ref_segment[7] = ref_seq[q+1];
                                    ref_segment[8] = ref_seq[j-1];
                                    
                                    y[0] = 8; ref_segment[0] = 8;
                                    str[9] = 0; segment[9] = 0;
                                    make_sequence(y, str, 1, 8);
                                    make_sequence(ref_segment, segment, 1, 8);
                                    h = hamming(segment+1, str+1);
                                    
                                    if (h <= d) {
                                        for (h_interior = 0; h_interior <= d-h && !done_sample; h_interior++) {
                                            uh = random_seq_hamming_gap(i+2, p-2, q+2, j-2, h_interior, pat, ref_seq);
                                            energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type2], y[5], y[8], y[6], y[7]);
                                            Qtemp += uh * EXP(energy) * Q_ham[type2][next_bp][d-h-h_interior+d_c];
                                            
                                            if (Qv <=Qtemp) {
                                                s[p] = y[2];
                                                s[q] = y[3];
                                                s[i+1] = y[5];
                                                s[p-1] = y[6];
                                                s[q+1] = y[7];
                                                s[j-1] = y[8];
                                                
                                                sample_random_seq_segment_gap_hamming(i+2, p-2, q+2, j-2,h_interior, s, pat, ref_seq);
                                                remain_distance_sub[next_bp] = d-h-h_interior+d_c;
                                                done_sample = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            if (s[p] == 0) {
                printf("Error\n");
            }
        } else if (nbp > 1) {
            d = remain_distance_sub[r];
            d_c = 0;
            if (s[i] != ref_seq[i]) d_c ++;
            if (s[j] != ref_seq[j]) d_c ++;
            
            Sample_multiloop_hamming(i+1, j-1, d-d_c, r, remain_distance_sub, s, ptable, order, Q_ham, pat, ref_seq, 1);
        }
    }
    //    for (l=1; l<=n; l++) printf("%d ", s[l]);
    //    printf("\n");
    
    make_sequence(s, seq, 0, n);
    seq[n+1] = 0;
    
    free(s);
    return seq;
}


int energy_eval(char *seq, int *ptable, int *order)
{
    int r, i, j, p, q, energy=0, nbp, ubp, n=ptable[0], type, type2;
    int *s;
    
    s = (int *) space (sizeof(int)*(n+5));
    for (r=0; r<n; r++) {
        switch(seq[r]) {
            case 'A':
                s[r+1]=1;
                break;
            case 'C':
                s[r+1]=2;
                break;
            case 'G':
                s[r+1]=3;
                break;
            case 'U':
                s[r+1]=4;
                break;
            default:
                s[r+1]=0;
                break;
        }
    }
    
    for (r=1; r<=order[0]; r++){
        i = order[r];
        j = ptable[i];
        type = BP_pair[s[i]][s[j]];
        
        nbp = nest_bp(ptable, i);
        if (nbp == 0) {
            energy += HairpinE(j-i-1, type, s[i+1], s[j-1], seq+i-1);
        } else if (nbp == 1) {
            p=i+1;
            while (ptable[p]==0 && p<j) p++;
            q=ptable[p];
            type2 = BP_pair[s[q]][s[p]];
            energy += LoopEnergy(p-i-1, j-q-1, type, type2, s[i+1], s[j-1], s[p-1], s[q+1]);
            
        } else if (nbp >1) {
            ubp = nest_up(ptable, i);
            energy += P->MLclosing+P->MLintern[type]+ubp*P->MLbase;
            p=i+1;
            while (p<j) {
                while (ptable[p]==0 && p<j) p++;
                q=ptable[p];
                if (p<q && q<j) {
                    type2 = BP_pair[s[p]][s[q]];
                    energy += P->MLintern[type2];
                    p=q+1;
                }
            }
        }
    }
    
    free(s);
    
    return energy;
}
float Ppatten(int *ptable, char *patten, int i, int j, int *order, long double ** gen_pf)// In condition that the PF has been computed
{
    float prob=0;
    int r, n=ptable[0], l, closure;
    int *pat;
    long double ** con_pf, Q1, Q2;
    
    pat = (int *) space (sizeof(int)*(n+3));
    
    for (r=0; r<=n;r++) {
        if (r>=i && r<=j) {
            switch (patten[r-i]){
                case 'A':
                    pat[r] = 1;
                    break;
                case 'C':
                    pat[r] = 2;
                    break;
                case 'G':
                    pat[r] = 3;
                    break;
                case 'U':
                    pat[r] = 4;
                    break;
                default:
                    pat[r] = 0;
                    break;
            }
        } else pat[r] = 0;
    }
    
    con_pf = PFunction(ptable, order, pat);
    
    Q1 = 0; Q2 = 0;
    for (r=1; r<=6; r++) {
        Q1+=con_pf[r][order[0]];
        Q2+=gen_pf[r][order[0]];
    }
    prob = (float) (Q1/ Q2);

    free(pat);
    return prob;
}

double Info(int n, int *fre)
{
    double res = 0;
    int i;
    for (i = 1; i<=n; i++) {
        if (fre[i]){
            res += (-1) * ((double) fre[i] / (double) n) * (log((double) fre[i] / (double) n));
        }
    }
    return res;
}



int H_dis (int *a, int *b){
    int i, h=0;
    for (i=1; i<=a[0]; i++) {
        if (a[i]!=b[i]) h++;
    }
    return h;
}

float random_seq_clustering(int l, int n)
{
    char **seq_set;
    int i, j;
    float C;
    
    seq_set = (char **) space (sizeof(char*)*(n+2));
    for (i=1; i<=n; i++) {
        seq_set[i]=(char *) space(sizeof(char)*(l+2));
        for (j=0;j<l;j++) {
            seq_set[i][j] = random_nucleotide();
        }
        seq_set[i][l] = 0;
    }
    
    C = clustering(seq_set, 1, l, n);
    
    for (i=1; i<=n; i++) {
        free(seq_set[i]);
    }
    free(seq_set);
    
    return C;
}

double random_seq_patten (int n, int l)
{
    char * seq;
    int j;
    patten *List=NULL, *index;
    int *fre;
    double res;
    
    seq = (char *) space (sizeof(char)*(l+2));
    fre = (int *) space (sizeof(int)*n);
    
    for (j=0; j<n; j++) fre[j]=0;
    
    for (j = 1; j<=n; j++) {
        randomseq(l, seq);
        List = patten_detector(seq, 0, l-1, 1, List);
    }
    
    index = List;
    while (index!=NULL) {
        // printf("%s  %lu\n", index->pat, index->cnt);
        res += (-1) * ((double) index->cnt / (double) n) * (log((double) index->cnt / (double) n));
        //fre[index->cnt]++;
        index = index->next;
    }
    
    free(seq);
    free(fre);
    return res/log(4);
}


double MI(int *ptable, char *seq)
{
    long double PF_struc, PF_seq, **PF, Qtemp;
    int *order, pat[500], ud;
    int i,r,bp,energy;
    double mi;
    
    for (i=0;i<500;i++) pat[i]=0;
    
    update_fold_params();
    P = scale_parameters();
    
    order = BPair_order(ptable);
    
    PF = PFunction(ptable, order, pat);
    PF_seq = pfunc(seq);
    ud = 0;
    i=1;
    PF_struc = 1;
    while (i<=ptable[0]) {
        if (ptable[i]==0) {
            i++;
            if (pat[i]==0) ud++;
        } else {
            Qtemp=0;
            for (r=1;r<=order[0];r++) {
                if (order[r]==i) {
                    bp=r;
                    break;
                }
            }
            
            for (r=1;r<=6;r++) {
                Qtemp+=PF[r][bp];
            }
            PF_struc *= Qtemp;
            i = ptable[i]+1;
        }
    }
    PF_struc *= pow(4,ud);
    
    energy = energy_eval(seq, ptable, order);
    
    //printf("%Le, %Le\n", PF_seq, PF_struc);
    
    mi = (double) (-1) * EXP (energy) * log ((EXP(energy)/ (PF_seq * PF_struc)));
    
    free(order);
    for (r=1;i<=6;r++) free(PF[r]);
    free(PF);
    
    return log(mi);
}


void Sampler (int mul, int *ptable, int *pat)
{
    int *order, *sam_p;
    char *str, *sam_str;
    long double **PF;
    int i,sum;
    char *sample, *ref_sample;
    int correct, total, correct100;
    int rat[5];
    int energy, energy_fold, h = 0;
    double entropy,mi1, mi2;
    int rank[500], rank_overall[500],delta_rank[500];
    
    
    update_fold_params();
    P = scale_parameters();
    
    order = BPair_order(ptable);
    str = pair2structure(ptable);
    
    PF = PFunction(ptable, order, pat);
    // printf("%s\n", str);
    
    for (i=0; i<=4; i++) rat[i]=0;
    correct100=0;
    for (i=0;i<500; i++) {
        rank[i]=0;
        rank_overall[i]=0;
        delta_rank[i]=0;
    }
    
    fp = fopen(Outputfilename, "w");
    if (fp==NULL) printf("Can't open file \"%s\".\n", Outputfilename);
    
    
    
    for (i=1; i<=mul;i++) {
        sample = Sample_seq(ptable, order, PF, pat);
        
        if (i == 1) {
            ref_sample = (char *) calloc (strlen(sample)+5, sizeof(char));
            strcpy(ref_sample, sample+1);
        }
        
        h += hamming(sample+1, ref_sample);
        
        if (compute_ratio) ratio(sample+1, rat);
        if (show_sequence) {
            fprintf(fp,">%d\n",i);
            fprintf(fp,"%s\n", sample+1);
        }
        
        
        if (refold) {
            //refold
            sam_p = (int *) calloc (ptable[0]+5, sizeof(int));
            energy_fold = fold_sec(sample+1, sam_p);
            sam_str = pair2structure(sam_p);
            
            if (show_fold_structure) fprintf(fp, "%s\n", sam_str);
        
            correct=compare(sam_p, ptable);
            total=number(sam_p);
        
            energy = energy_eval(sample+1, ptable, order);
        
        
            //mi1 = MI(ptable, sample+1);
            //mi2 = MI(sam_p, sample+1);
        
            if (show_energy) fprintf(fp, ">%d Energy: %d  Energy_fold %d \n", i, energy, energy_fold);
            
            if (show_MI) fprintf(fp, "MI: %lf  MI_fold:  %lf\n", mi1, mi2);
        
        //printf("Energy: %d\n", energy);
            energy = energy/100; 
            rank_overall[-energy]++;
            delta_rank[(abs(energy-energy_fold))/100]++;
            if (correct == order[0] && correct == total) {
               if (energy < 0) rank[-energy]++;
               correct100++;
            }
        //refold end
            
            free(sam_str);
            free(sam_p);
        }
        free(sample);
    }
    if (compute_ratio) {
        sum=0;
        for (i=1; i<=4; i++) sum+=rat[i];
    
        fprintf(fp, "A: %f ", ((float) rat[1]) / ((float) sum));
        fprintf(fp, "C: %f ", ((float) rat[2]) / ((float) sum));
        fprintf(fp, "G: %f ", ((float) rat[3]) / ((float) sum));
        fprintf(fp, "U: %f ", ((float) rat[4]) / ((float) sum));
        
        fprintf(fp, "\n");
        fprintf(fp, "A/U: %f ", ((float) (rat[1]+rat[4])) / ((float) sum));
        fprintf(fp, "C/G: %f ", ((float) (rat[2]+rat[3])) / ((float) sum));
        fprintf(fp, "\n");
    }
    
    if (refold) {
        fprintf(fp, "100 Correct: %d\n", correct100);
    
    
        fprintf(fp, "Energy Distributon of Inverse fold solution: \n");
        for (i=0; i<=50; i++) {
            fprintf(fp, "%d,", rank[i]);
        }
        fprintf(fp, "\n");
    
        fprintf(fp, "Energy Distribution: \n");
        for (i=0; i<=50; i++) {
            fprintf(fp, "%d,", rank_overall[i]);
        }
    
        fprintf(fp, "\n");
        fprintf(fp, "Delta Eenergy: \n");
        for (i=0; i<=100; i++) {
            fprintf(fp, "%d,", delta_rank[i]);
        }
    
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "Ave hamming ditance: %f \n", ((float) h / mul));
    
    fclose(fp);
    free(order);
}

void mutaion_rate(patten *List, char *ref_sequence, FILE *out, float rat[1000][5][5])
{
    int i, n = strlen(ref_sequence), so, sn, cnt = 0, r, s, n_ref;
    patten *index;
    char no, nn;
    index = List;
    while (index!=NULL) {
        for (i=0; i<n; i++) {
            so = nucleotide2code(ref_sequence[i]);
            sn = nucleotide2code(index->pat[i]);
            rat[i][so][sn] += index->cnt;
        }
        cnt += index->cnt;
        index = index->next;
    }
    
    for (i=0; i<n; i++) {
        for (r=0; r<=4; r++) {
            for (s=0; s<=4; s++) {
                rat[i][r][s] /= cnt;
            }
        }
    }
    
    if (out) {
        for (i=0; i<n; i++) {
            n_ref = nucleotide2code(ref_sequence[i]);
            for (r=1; r<=4; r++) {
                for (s=1; s<=4; s++) {
                    no = code2nucleotide(r);
                    nn = code2nucleotide(s);
                    
                    if (rat[i][r][s]==0 && r==n_ref) fprintf(out, "%d: %c->%c: %f\n",i, no, nn, rat[i][r][s]);
                }
            }
        }
    }
}



patten *sample_compatible_seq(char *structure, int mul)
{
    int *ptable, *order;
    int i, j, type;
    int seq[1000];
    char *sequence;
    patten *List = NULL;
    float seed;
    
    ptable = structure2pair(structure);
    order = BPair_order(ptable);
    
    for (i=0; i<mul; i++) {
        for (j = 1; j<=ptable[0]; j++) {
            if (ptable[j] == 0) {
                seq[j] = random_nucleotide();
            } else if (ptable[j] != 0 && j < ptable[j]) {
                seed=(float) ((drand48()) * 6);
                type = ((int) seed) + 1;
                type_to_bp(type, j, ptable[j], seq);
            }
        }
        sequence = (char *) calloc (ptable[0]+5, sizeof(char));
        sequence[0] = '_';
        make_sequence(seq, sequence, 1, ptable[0]);
        List = patten_detector(sequence, 1, ptable[0], 1, List);
    }
    
    
    return List;
}


patten *sample_compatible_seq_ham(char *structure, char *ref_sequence, int distance, int mul)
{
    return NULL;
}


patten *seq_sampler_ham (char *structure, char *ref_sequence, char *sequence_patten, int distance, int mul, FILE *out)
{
    int *order, *sam_p, *remain_d_sub, *ptable, *seq_pat, *ref_seq, energy, ave_energy=0, fold_struc[1000], h, n_seq=0;
    char *str, *sam_str;
    long double ***PF_ham=NULL, ***PF_ham_Z=NULL, ***PF_ham_X=NULL;
    int i,j,r;
    char *sample, *fold_structure;
    patten * List = NULL, *index;
    float rat[1000][5][5];
    int bp, l;
    int temp_str[1000];
    
    update_fold_params();
    P = scale_parameters();
    
    ptable = structure2pair(structure);
    seq_pat = seq2code(sequence_patten);
    ref_seq = seq2code(ref_sequence);
    
    order = BPair_order(ptable);
    
    PF_ham = (long double ***) space(sizeof(long double **)*7);
    PF_ham_Z = (long double ***) space(sizeof(long double **)*7);
    PF_ham_X = (long double ***) space(sizeof(long double **)*7);
    
    bp = order[0];
    l = ptable[0];
    for (i=0; i<=6; i++) {
        PF_ham[i] = (long double **) calloc (bp+5, sizeof(long double));
        PF_ham_Z[i] = (long double **) calloc (bp+5, sizeof(long double));
        for (j=0; j<=bp+2; j++) {
            PF_ham[i][j] = (long double *) calloc (l+5, sizeof(long double));
            PF_ham_Z[i][j] = (long double *) calloc (l+5, sizeof(long double));
            for (r=0; r<=l; r++) {
                PF_ham[i][j][r] = 0;
                PF_ham_Z[i][j][r] = 0;
            }
        }
    }

    PFunction_hamming(ptable, order, seq_pat, ref_seq, distance, PF_ham, PF_ham_Z, PF_ham_X);
    
    
    energy = energy_eval(ref_sequence, ptable, order);
    //printf("%d\n", energy);
    //printf("%Le\n", PF_ham[0][order[0]+1][distance]);
    //printf("%Le\n", PF_ham_Z[0][order[0]+1][distance]);
    if (out) fprintf(out, "Expected Energy: %Lf\n", PF_ham_Z[0][order[0]+1][distance]/PF_ham[0][order[0]+1][distance]);
    // printf("%Le\n", EXP(energy));
    
    remain_d_sub = (int *) calloc (ptable[0]+5, sizeof(int));
    for (i=1; i<=mul; i++) {
        sample = Sample_seq_ham(ptable, ref_seq, order, remain_d_sub, PF_ham, seq_pat, distance);
        List = patten_detector(sample, 1, ptable[0], 1, List);
    }
    
    //sortList(List);
    
    /*
    for (i=0;i<1000;i++) {
        for (j=0;j<5; j++) {
            for (r=0;r<5;r++) {
                rat[i][j][r] = 0;
            }
        }
    }
    mutaion_rate(List, ref_sequence, out, rat);
    */
    
    if (out) {
        index = List;
        while (index!=NULL) {
            fprintf(out, "%s\n", index->pat);
            fprintf(out, "%d\n", index->cnt);
            // h = hamming(index->pat, ref_sequence);
            // fprintf(out, "D: %d\n", h);
            if (eval) {
                energy = energy_eval(index->pat, ptable, order);
                fprintf(out, "%d\n", energy);
                ave_energy += energy * index->cnt;
            }
            if (refold) {
                energy = fold_sec(index->pat, fold_struc);
                fold_structure = pair2structure(fold_struc);
                fprintf(out, "%s\n", fold_structure);
                free(fold_structure);
            }
            n_seq++;
            index = index->next;
        }
        
    
        fprintf(out, "Num_seq: %d\n", n_seq);
        fprintf(out, "Sampled average energy: %f\n", (((float) ave_energy)/ (mul*100)));
    }
    
    
    free(remain_d_sub);
    free(order);
    
    free(ptable);
    free(seq_pat);
    free(ref_seq);
    
    return List;
}



void multi_sample(int mul, int *ptable, int *pat, int flag)
{
    int *order;
    char *str, *sam_str;
    int i, j,k, l;
    char *sample;
    patten * List, *indx, *current;
    long double **PF;
    int rat[5], sum;
    int *sam_p;
    int correct,total, correct100;
    int *sort_energy, *sort_id, *energy_fold;
    int energy;
    int rank[500], rank_overall[500],delta_rank[500];
    char **sample_set;
    float C[100][100];
    int fre[500];
    double entropy,mi1, mi2;
    
    
    update_fold_params();
    P = scale_parameters();
    
    order = BPair_order(ptable);
    str = pair2structure(ptable);
    
    PF = PFunction(ptable, order, pat);
    printf("%s\n", str);
    
    
    for (i=0; i<=4; i++) rat[i]=0;
    List = NULL;
    correct100=0;
    for (i=0;i<500; i++) {
        rank[i]=0;
        rank_overall[i]=0;
        delta_rank[i]=0;
    }
    energy_fold = (int *) space(sizeof(int));
    
    if (flag==5) {
        printf("%f\n", Ppatten(ptable, "GAAACCCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GGGACCCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GGAAGGGC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GAAAGGCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GUGACCCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GAGACCCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GAAAGGGC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GUGAGGCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GCAAGGCC", 38, 45, order, PF));
        printf("%f\n", Ppatten(ptable, "GGAAGGCC", 38, 45, order, PF));
        
        return;
    }
    
    sample_set = (char **) space (sizeof(char*)*(mul+2));
    
    

    printf("Native MI %lf\n", MI(ptable, "GCAGCAGGGAACUCACGCUUGCGUAGAGGCUAAGUGCUUCGGCACAGCACAAGCCCGCUGCG"));  //2N3R
    //printf("Native MI %le\n", MI(ptable, "GGAGGUAGUAGGUCGAAAGACCAUUCUGCCUCC"));  //2JXV
    
    if (flag==0) {
        sample = Sample_seq(ptable, order, PF, pat);
        printf("%s\n", sample+1);
    }
    
    for (i=1; i<=mul; i++) {
        
        //sample
        sample = Sample_seq(ptable, order, PF, pat);


        if (flag ==1) {
            ratio(sample+1, rat);
            //printf(">%d\n",i);
            //printf("%s\n", sample+1);
        
            //refold
            energy_fold = fold_sec(sample+1,sam_p);
            sam_str = pair2structure(sam_p);
        
            correct=compare(sam_p, ptable);
            total=number(sam_p);
        
            energy = energy_eval(sample+1, ptable, order);
            
            
            mi1 = MI(ptable, sample+1);
            mi2 = MI(sam_p, sample+1);
            
            if (50 < mi1 && mi1 < 51) {
                
                printf("%s  %lf\n", sample+1, mi1);
            }
            
            //printf(">%d Energy: %d  Energy_fold %d   MI: %lf  MI_fold:  %lf\n", i, energy, *energy_fold, mi1, mi2);
            
            //printf("Energy: %d\n", energy);
            rank_overall[-energy]++;
            delta_rank[(abs(energy-*energy_fold))/100]++;
            if (correct == order[0] && correct == total) {
                if (energy < 0) rank[-energy]++;
                correct100++;
            }
            //refold end
        }
        
        if (flag==2 || flag == 3) {
            sample_set[i] = Sample_seq(ptable, order, PF, pat);
        }
        
        if (flag ==4) List = patten_detector(sample+1, 38-1, 45-1, 1, List);
        
        if (flag==1) {
            free(sam_str);
            free(sam_p);
            free(sample);
        }
    }
    
    if (flag == 2) {
        for (j=1; j<=10; j++) {
            for (k=1; j+k<=ptable[0]; k++) {
                C[k][j+k] = clustering(sample_set, k,j+k,mul);
                printf("(%d,%d): %f\n", k, j+k, C[k][j+k]);
            }
        }
    }
    
    if (flag == 1) {
        sum=0;
        for (i=1; i<=4; i++) sum+=rat[i];
    
        for (i=1; i<=4; i++) {
            printf("%f ", ((float) rat[i]) / ((float) sum));
        }
        printf("\n");
        printf("A/U: %f ", ((float) (rat[1]+rat[4])) / ((float) sum));
        printf("C/G: %f ", ((float) (rat[2]+rat[3])) / ((float) sum));
        printf("\n");
    
    
        printf("100 Correct: %d\n", correct100);
        
        
        printf("ED of Inverse fold solution: \n");
        for (i=0; i<=50; i++) {
            printf("%d,", rank[i]);
        }
        printf("\n");
        
        printf("ED: \n");
        for (i=0; i<=50; i++) {
            printf("%d,", rank_overall[i]);
        }
        
        printf("\n");
        printf("delta ED: \n");
        for (i=0; i<=100; i++) {
            printf("%d,", delta_rank[i]);
        }
        
        printf("\n");
    }
    
    
    if (flag==3) {
        for (l=0; l<=7; l++) {
            printf("L: %d\n", l);
            for (i=1; i+l<=ptable[0]; i++) {
                List = NULL;
                entropy = 0;
                //printf("[%d %d]: ", i, i+l);
            //for (j=0; j<500; j++) fre[j]=0;
                for (j=1; j<=mul; j++) {
                    List = patten_detector(sample_set[j]+1, i-1, i+l-1,1, List);
                }
                indx = List;
                while (indx != NULL) {
                    current = indx;
                    //printf("%s  %lu\n", current->pat, current->cnt);
                    entropy += (-1) * ((double) indx->cnt / (double) mul) * (log((double) indx->cnt / (double) mul));
                    //fre[current->cnt]++;
                    indx = indx->next;
                    free(current->pat);
                    free(current);
                }
                //for (j=1; j<100; j++) printf("%d ", fre[j]);
                //printf("\n");
                printf("%f, ", 1-((entropy/log(4))/(l+1)));
            }
            printf("\n");
        }
    }
    if (flag==4) {
        indx =List;
        while (indx!=NULL){
            if (indx->cnt>100) printf("%s  %d\n", indx->pat, indx->cnt);
            indx=indx->next;
        }
    }

    free(energy_fold);
    free(order);
}

void Entropy (int i, int j, int *pat, int l, long double **PF_gen, int *ptable, int *order, float *E) {

    int r,ud=0, udp=0, k;
    long double ** PF_con,Q1=0,Q2=0;
    
    if (l == j-i+1) {
        PF_con = PFunction(ptable, order, pat);
        r = 1;
        while (r<=ptable[0]) {
            if (ptable[r] == 0) {
                ud++;
                if (pat[r]==0) udp++;
                r++;
            } else {
                for (k=1; k<=6; k++) {
                    Q1 += PF_con[k][order[0]];
                    Q2 += PF_gen[k][order[0]];
                    //only work for a uniqe exteiror arc
                }
                r = ptable[r]+1;
            }
        }
        Q1 *= pow(4,udp);
        Q2 *= pow(4,ud);
        *E += (float) (-1)* ((Q1/Q2) * log(Q1/Q2) / log(4));
        
        for (r=1; r<=6; r++) {
            free(PF_con[r]);
        }
        free(PF_con);
    } else {
        for (r=1; r<=4; r++) {
            pat[i+l] = r;
            Entropy(i,j,pat,l+1,PF_gen,ptable, order,E);
            pat[i+l]=0;
        }
    }
}

void HP(int *ptable)
{
    int *order,r,l,i;
    float *E;
    int pat[500];
    long double ** PF;
    
    update_fold_params();
    P = scale_parameters();
    
    order = BPair_order(ptable);
    E = (float *) space (sizeof(float));
    
    for (r=1; r<500;r++) pat[r]=0;

    PF = PFunction(ptable, order, pat);
    
    
    for (l=0; l<=8; l++) {
        printf("L: %d\n", l);
        for (i=1; i+l<=ptable[0]; i++) {
            *E=0;
            Entropy (i,i+l,pat, 0,PF,ptable,order,E);
            printf("%f, ", 1-(*E/(l+1)));
        }
        printf("\n");
    }
}


void Hull(int *ptable, int *pat)
{
    int *order, *sam_p, i;
    char *str, *sam_str;
    long double **PF;
    char *sample;
    int *energy_fold;
    patten *List = NULL, *index;
    
    
    update_fold_params();
    P = scale_parameters();
    
    order = BPair_order(ptable);
    str = pair2structure(ptable);
    
    PF = PFunction(ptable, order, pat);
    // printf("%s\n", str);
    
    energy_fold = (int *) space(sizeof(int));
    
    fp = fopen(Outputfilename, "w");
    if (fp==NULL) printf("Can't open file \"%s\".\n", Outputfilename);
    
    
    for (i=1; i<10000; i++) {
        sample = Sample_seq(ptable, order, PF, pat);
        //printf("%s\n", sample);
        sam_p = (int *) calloc (strlen(sample)+2, sizeof(int));
        
        energy_fold = fold_sec(sample+1,sam_p);
        sam_str = pair2structure(sam_p);
        List = patten_detector(sam_str, 1, ptable[0], 1, List);
        
        free(sam_p);
    }
    
    index = List;
    while (index!=NULL) {
        if (index->cnt > 50) printf("%s %d\n", index->pat, index->cnt);
        index = index->next;
    }
    
}






void debug()
{
    

    
}


void helix_e(void)
{
    int p1, p2, p3, p4;
    int energy;
    int seq[10];
    char string[10];
    
    seq[0]=8;
    string[0] = '_'; string[9]=0;
    
    for (p1=1;p1<=6; p1++) {
        for (p2=1; p2<=6; p2++) {
            for (p3=1;p3<=6; p3++) {
                for (p4=1;p4<=6; p4++) {
                    type_to_bp(p1, 1, 8, seq);
                    type_to_bp(p2, 2, 7, seq);
                    type_to_bp(p3, 3, 6, seq);
                    type_to_bp(p4, 4, 5, seq);
                    
                    energy = 0;
                    energy+=LoopEnergy(0, 0, p1, rtype[p2], seq[2], seq[7], seq[1], seq[8]);
                    energy+=LoopEnergy(0, 0, p2, rtype[p3], seq[3], seq[6], seq[2], seq[7]);
                    energy+=LoopEnergy(0, 0, p3, rtype[p4], seq[4], seq[5], seq[3], seq[6]);
                    
                    make_sequence(seq, string, 1, 8);
                    printf("%c%c%c%c-%c%c%c%c: %d\n", string[1],string[2],string[3],string[4],string[5],string[6],string[7],string[8], energy);
                }
            }
        }
    }
}


