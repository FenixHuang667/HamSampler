#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fold_vars.h"
#include "pair_mat.h"
#include "energy_par.h"
#include "utils.h"
#include "stat.h"
#include "energy_const.h"
#include "fold.h"

#include "PF.h"
#include "pfunc.h"



#define R GASCONST
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))


typedef struct block block;
struct block
{
    int i;
    int j;
    int type;
    // 10:c
    // 11:f5 12:fML
    // 14:f51 15:fML1
    
    block *Next;
};


double *Qc,*QfML,*QfML1;
double *Qf5,*Qf51;

int *c,*fML, *fPL,*indx2, *f5, *sym;
int *DMLi, *DMLi1;

int *indx2;
block *List;




void pfunc_initial (int n)
{
    int i;
    
    indx2 = (int *) calloc ((n+5), sizeof(int));
    for (i = 1; i <= n; i++)
        indx2[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */
    
    Qc     = (double *) calloc(((n*(n+1))/2+2), sizeof(double));
    Qf5    = (double *) calloc(((n*(n+1))/2+2), sizeof(double));
    Qf51   = (double *) calloc(((n*(n+1))/2+2), sizeof(double));
    QfML   = (double *) calloc(((n*(n+1))/2+2), sizeof(double));
    QfML1  = (double *) calloc(((n*(n+1))/2+2), sizeof(double));
}

void pfunc_freevar()
{
    free(Qc); free(QfML); free(QfML1);
    free(Qf5); free(Qf51);
    free(indx2);
}

int * encode_seq(char *seq) {
    int i,l, *s;
    
    l = strlen(seq);
    s = (int *) calloc (l+2, sizeof(int));
    s[0] = l;
    
    for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
        s[i]= (short) encode_char(toupper(seq[i-1]));
    }
    return s;
}

block * newblock()
{
    block * b;
    b = (block *) calloc (1, sizeof(block));
    b->i = 0;
    b->j = 0;
    b->type = 0;
    b->Next = NULL;
    return b;
}

void pushblock (int i, int j, int type)
{
    block * new;
    new = newblock();
    new->i = i;
    new->j = j;
    new->type = type;
    new->Next = List;
    List = new;
}


void initial_sec (int n)
{
    int i,k;
    
    indx2 = (int *) space(sizeof(int)*(n+5));
    for (i = 1; i <= n; i++)
        indx2[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */
    
    c     = (int *) calloc (((n*(n+1))/2+2), sizeof(int));
    f5    = (int *) calloc (((n*(n+1))/2+2), sizeof(int));
    DMLi    = (int *) calloc (n+3, sizeof(int));
    DMLi1    = (int *) calloc (n+3, sizeof(int));
    sym    = (int *) calloc (((n*(n+1))/2+2), sizeof(int));
    fML   = (int *) calloc (((n*(n+1))/2+2), sizeof(int));
    //fPL   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
}

void free_sec (void)
{
    free(indx2);
    free(c);
    free(f5);
    free(DMLi);
    free(DMLi1);
    free(sym);
    free(fML);
}

void free_pf (void)
{
    free (Qc);
    free (Qf5);
    free (Qf51);
    free (QfML);
    free (QfML1);
}

void fill_array_sec(const char *seq, int *s)
{
    int i,j,k1,k2,p,q,ij;
    int length,type2,G,type;
    int new_c, new_ML, new_P, new_f5;
    int *FF;
    
    length=strlen(seq);
    s = encode_seq(seq);
    
    if (length<=TURN) return;
    
    for (i=1;i<=length;i++)
        for (j=i;j<=i+TURN+1;j++) {
            ij=indx2[j]+i;
            c[ij]=MAXENG;
            f5[ij]=0;
            fML[ij]=MAXENG;
        }
    for (j=1;j<=length;j++) {
        DMLi[j]=MAXENG;
        DMLi1[j]=MAXENG;
    }
    
    for (i=length-TURN-1; i>=1; i--) {
        for (j=i+TURN+1; j<=length; j++) {
            
            type = BP_pair[s[i]][s[j]];
            ij=indx2[j]+i;
            
            if (type) {
                
                new_c = HairpinE(j-i-1, type, s[i+1], s[j-1], seq+i-1); //Hairpin loop
                
                for (p = i+1; p <=  MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
                    int minq = j-i+p-MAXLOOP-2;
                    if (minq<p+1+TURN) minq = p+1+TURN;
                    for (q = minq; q < j; q++) {
                        type2 = BP_pair[s[p]][s[q]];
                        
                        if (type2==0) continue;
                        type2 = rtype[type2];
                        
                        G = LoopEnergy(p-i-1, j-q-1, type, type2, s[i+1], s[j-1], s[p-1], s[q+1]);
                        new_c = MIN2(new_c, G+c[indx2[q]+p]); // Interior loop
                    }
                }
                
                G=P->MLintern[type]+P->MLclosing+DMLi1[j-1];
                if (G < new_c) {
                    new_c=G;  // Multi-loop
                }
                
                c[ij]=new_c;
                
            } else
                c[ij]=MAXENG; // Tight sec struct
            
            
            sym[ij]=0;
            new_f5=0;
            
            if (type>2) G=P->TerminalAU; else G=0;
            if (G+c[ij] < new_f5) {
                new_f5=G+c[ij];
                sym[ij]=1;
            }
            if (f5[ij+1] < new_f5) {
                new_f5=f5[ij+1];
                sym[ij]=0;
            }
            
            for (k1=i+TURN+1;k1<j;k1++) {
                //if (!sym[indx2[k1]+i]) continue;
                type2 = BP_pair[s[i]][s[k1]];
                if (type2>2) G=P->TerminalAU; else G=0;
                if (type2 && c[indx2[k1]+i]+f5[indx2[j]+k1+1]+G < new_f5) {
                    new_f5=c[indx2[k1]+i]+f5[indx2[j]+k1+1]+G;
                    sym[ij]=0;
                }
            }
            f5[ij]=new_f5; // f5 structure
            
            new_ML=MAXENG;
            if (fML[ij+1]+P->MLbase < new_ML) new_ML=fML[ij+1]+P->MLbase;
            if (fML[indx2[j-1]+i]+P->MLbase < new_ML) new_ML=fML[indx2[j-1]+i]+P->MLbase;
            if (type && c[ij]+P->MLintern[type]< new_ML)
                new_ML=c[ij]+P->MLintern[type];
            
            DMLi[j]=DMLi[j-1]+P->MLbase;
            for (k1=i+TURN+1;k1<j-TURN-1;k1++) {
                type2=BP_pair[s[k1+1]][s[j]];
                if (type2) {
                    G=P->MLintern[type2];
                    if (c[indx2[j]+k1+1]+fML[indx2[k1]+i]+G < DMLi[j]) {
                        DMLi[j]=c[indx2[j]+k1+1]+fML[indx2[k1]+i]+G;
                    }
                }
            }
            
            new_ML = MIN2(new_ML, DMLi[j]);
            fML[ij]=new_ML; //fML
        }
        
        FF=DMLi1; DMLi1=DMLi; DMLi=FF;
        for (j=1;j<=length;j++) {
            DMLi[j]=MAXENG;
        }
    }
}
// secondary structure S done!


void de_I(block *T, const char *seq, int *s, int *ptable)
{
    int p,q,k1,k2, minq,G,i,j, ij;
    int typeij, type_l;
    
    i=T->i;
    j=T->j;
    ij=indx2[j]+i;
    
    switch (T->type) {
        case 10: //decompose c
            typeij=BP_pair[s[i]][s[j]];
            if (typeij) {
                ptable[i]=j;
                ptable[j]=i;
                
                if (c[indx2[j]+i]==HairpinE(j-i-1, typeij, s[i+1], s[j-1], seq+i-1)) {
                    return;
                }
                
                for (p=i+1;p<=MIN2(j-2-TURN,i+MAXLOOP+1);p++) {
                    minq = j-i+p-MAXLOOP-2;
                    if (minq<p+1+TURN) minq = p+1+TURN;
                    for (q = minq; q<j; q++) {
                        type_l=BP_pair[s[q]][s[p]];
                        if (type_l==0) continue;
                        G=LoopEnergy(p-i-1, j-q-1, typeij, type_l, s[i+1], s[j-1], s[p-1], s[q+1]);
                        if (c[indx2[j]+i]==G+c[indx2[q]+p]) {
                            pushblock(p,q,10);
                            return;
                        }
                    }
                }
                
                for (k1=i+1; k1<j; k1++) {
                    for (k2=k1+TURN; k2<j; k2++) {
                        type_l = BP_pair[s[k1]][s[k2]];
                        if (type_l == 0) continue;
                        G = P->MLbase * (k1-i-1) + P->MLintern[type_l] + c[indx2[k2]+k1] +fML[indx2[j-1]+k2+1];
                        if (c[ij] == P->MLintern[typeij]+P->MLclosing +G) {
                            pushblock(k1, k2, 10);
                            pushblock(k2+1, j-1, 12);
                            return;
                        }
                    }
                }
            } else {
                printf("Error! Do not match the energy, please check your program carefully!\n");
                return;
            }
            break;
            
        case 11: //decompose f5
            if (j-i<=TURN+1) {
                return;
            }
            
            typeij=BP_pair[s[i]][s[j]];
            if (typeij>2) G=P->TerminalAU; else G=0;
            if (typeij && f5[ij]==c[ij]+G) {
                pushblock(i,j,10);
                return;
            }
            if (f5[ij]==f5[ij+1]) {
                pushblock(i+1,j,11);
                return;
            }
            for (k1=i+TURN+1;k1<j;k1++) {
                if (f5[ij]==f5[indx2[k1]+i]+f5[indx2[j]+k1+1]) {
                    pushblock(i,k1,11);
                    pushblock(k1+1,j,11);
                    return;
                }
            }
            break;
        case 12: //decompose fML
            if (j-i<=TURN+1) return;
            
            typeij=BP_pair[s[i]][s[j]];
            
            if (typeij && fML[ij]==P->MLintern[typeij]+c[ij]) {
                pushblock(i,j,10);
                return;
            }
            
            if (fML[ij]==fML[ij+1]+P->MLbase) {
                pushblock(i+1,j,12);
                return;
            }
            
            if (fML[ij]==fML[indx2[j-1]+i]+P->MLbase) {
                pushblock(i,j-1,12);
                return;
            }
            for (k1=T->i+TURN+1;k1<T->j-TURN-1;k1++) {
                if (fML[ij]==fML[indx2[k1]+i]+fML[indx2[j]+k1+1]) {
                    pushblock(i,k1,12);
                    pushblock(k1+1,j,12);
                    return;
                }
            }
            break;
        default : break;
    }
    
    printf( "Error! Do not match the energy, please check your program carefully!\n");
    printf( "%d %d \n", T->type, c[ij]);
}


int fold_sec (char *seq, int * ptable)
{
    int n, i, ene;
    int *s;
    block *T;
    
    n = strlen(seq);
    
    initial_sec(n);
    
    update_fold_params();
    P = scale_parameters();
    s = encode_seq(seq);
    
    fill_array_sec(seq, s);
    
    
    List = NULL;
    pushblock(1, n, 11);
    
    for (i=1;i<=n;i++) ptable[i] = 0;
    ptable[0] = n;
    
    while (List !=NULL) {
        T = List;
        List = List->Next;
        
        de_I(T, seq, s, ptable);
    }
    
    free(s);
    ene = f5[indx2[n]+1];
    free_sec();
    return ene;
}

int check_available (int p, int *pat, int paired)
{
    if (pat[p] == -2) return 1;
    if (pat[p] == -1) {
        if (paired > 0) return 1;
        else return 0;
    }
    if (pat[p] == 0) {
        if (paired == 1) return 0;
        else if (paired == 0) return 1;
    }
    if (pat[p]>0) {
        if (pat[p] == paired) return 1;
        else return 0;
    }
    return 0;
}


void partition(char *seq, int *pat)
{
    int i,j, p, q, k1, r;
    int *s;
    int length,type1,type2,G,type, available;
    
    length=strlen(seq);
    for (i=1;i<=length;i++) {
        for (j=i;j<i+TURN;j++) {
            Qc[indx2[j]+i]=0;
            Qf5[indx2[j]+i]=1;
            Qf51[indx2[j]+i]=0;
            QfML[indx2[j]+i]=0;
            QfML1[indx2[j]+i]=0;
        }
    }
    
    s = encode_seq(seq);
    
    for (i=length-TURN;i>=1;i--) {
        for (j=i+TURN;j<=length;j++) {
            
            type = BP_pair[s[i]][s[j]];
            if (type && check_available(i, pat, j) && check_available(j, pat, i)) {
                G=HairpinE(j-i-1, type, s[i+1], s[j-1], seq+i-1);
                Qc[indx2[j]+i]=EXP(G);
                for (p = i+1; p <= j-2-TURN; p++) {
                    for (q = p+TURN; q < j; q++) {
                        available = 1;
                        
                        for (r=i+1; r<p; r++) if (!check_available(r, pat, 0)) available = 0;
                        for (r=q+1; r<j; r++) if (!check_available(r, pat, 0)) available = 0;
                        if (!check_available(p, pat, q)) available = 0;
                        if (!check_available(q, pat, p)) available = 0;
                        
                        type2 = BP_pair[s[p]][s[q]];
                        if (type2==0) available = 0;
                        if (!available) continue;
                        
                        type2 = rtype[type2];
                        
                        G = LoopEnergy(p-i-1, j-q-1, type, type2, s[i+1], s[j-1], s[p-1], s[q+1]);
                        Qc[indx2[j]+i]+=EXP(G)*Qc[indx2[q]+p];
                        //contribution from int, bulge, helix
                    }
                }
                for (k1=i+1;k1<j-1;k1++) {
                    G=P->MLintern[type]+P->MLclosing;
                    Qc[indx2[j]+i]+=QfML1[indx2[k1]+i+1]*QfML[indx2[j-1]+k1+1]*EXP(G);
                    //contribution from mul
                }
            } else {
                Qc[indx2[j]+i]=0;
            }
            
            //compute Qc[i,j]
            
            Qf51[indx2[j]+i]=0;
            for (k1=i;k1<j;k1++) {
                available = 1;
                for (r=i; r<k1; r++) if (!check_available(r, pat, 0)) available = 0;
                type = BP_pair[s[k1]][s[j]];
                if (!type) available = 0;
                if (!check_available(k1, pat, j) || !check_available(j, pat, k1)) available = 0;
                if (type>2) G=P->TerminalAU; else G=0;
                if (available) Qf51[indx2[j]+i]+=Qc[indx2[j]+k1]*EXP(G);
            }
            //compute Qf51[i,j]
            
            
            Qf5[indx2[j]+i]=Qf51[indx2[j]+i]+1;
            for (k1=i;k1<j;k1++) {
                Qf5[indx2[j]+i]+=Qf51[indx2[k1]+i]*Qf5[indx2[j]+k1+1];
            }
            //compute Qf5[i,j]
            
            QfML1[indx2[j]+i]=0;
            for (k1=i;k1<j;k1++) {
                available = 1;
                type=BP_pair[s[k1]][s[j]];
                if (!type) available = 0;
                if (!check_available(k1, pat, j) || !check_available(j, pat, k1)) available = 0;
                for (r=i; r<k1; r++) if (!check_available(r, pat, 0)) available = 0;
                if (available) {
                    G=P->MLintern[type]+(k1-i)*P->MLbase;
                    QfML1[indx2[j]+i]+=Qc[indx2[j]+k1]*EXP(G);
                }
            }
            //compute QfML1[i,j]
            
            QfML[indx2[j]+i]=QfML1[indx2[j]+i];
            for (k1=i;k1<j;k1++) {
                available = 1;
                for (r=k1+1; r<=j; r++) if (!check_available(r, pat, 0)) available = 0;
                if (available) {
                    QfML[indx2[j]+i]+=QfML1[indx2[k1]+i]*(QfML[indx2[j]+k1+1]);
                    G=(j-k1)*P->MLbase;
                    QfML[indx2[j]+i]*=EXP(G);
                }
            }
            //compute QfML[i,j]
        }
    }
    free(s);
}
void pfunc(char *seq, int *pat)
{
    int length;
    length=strlen(seq);
    update_fold_params();
    P = scale_parameters();
    
    pfunc_initial(length);
    
    partition(seq, pat);
    //	printf("PF: %le\n", Qf5[indx2[length]+1]);
}

//Partion function done




void deblock(char *seq, int *ptable) {
    //assume pf is computed
    
    int type, type2, l, *s, G, p, q, k1, done;
    double seed, Qpro, Q;
    block *T;
    
    
    s = encode_seq(seq);
    l = s[0];
    ptable[0] = l;
    
    List = NULL;
    pushblock(1, l, 11);
    
    while (List !=NULL) {
        T = List;
        List = List->Next;
        Q = 0;
        done = 0;
        
        // printf("(%d, %d): %d\n", T->i, T->j, T->type);
        
        
        if (T->type == 10) {
            seed=1-drand48();
            Qpro = seed * Qc[indx2[T->j]+T->i];
            ptable[T->i] = T->j;
            ptable[T->j] = T->i;
            
            type = BP_pair[s[T->i]][s[T->j]];
            G=HairpinE(T->j-T->i-1, type, s[T->i+1], s[T->j-1], seq+T->i-1);
            Q = EXP(G);
            if (Q > Qpro) {
                done = 1;
            }
            if (done) continue;
            
            for (p = T->i+1; p <= T->j-2-TURN; p++) {
                if (done) break;
                for (q= p+TURN; q < T->j; q++) {
                    type2 = BP_pair[s[p]][s[q]];
                    if (type2==0) continue;
                    type2 = rtype[type2];
                    
                    G = LoopEnergy(p-T->i-1, T->j-q-1, type, type2, s[T->i+1], s[T->j-1],
                                   s[p-1], s[q+1]);
                    Q += EXP(G) * Qc[indx2[q]+p];
                    
                    if (Q > Qpro) {
                        pushblock(p, q, 10);
                        done = 1;
                        break;
                    }
                }
            }
            if (done) continue;
            
            
            for (k1=T->i+1;k1<T->j-1;k1++) {
                G=P->MLintern[type]+P->MLclosing;
                Q += QfML1[indx2[k1]+T->i+1]*QfML[indx2[T->j-1]+k1+1]*EXP(G);
                if (Q > Qpro) {
                    pushblock(T->i+1, k1, 15);
                    pushblock(k1+1, T->j-1, 12);
                    done = 1;
                    break;
                }
            }
            // Qc
        } else if (T->type == 11) {
            seed=1-drand48();
            Qpro = seed * Qf5[indx2[T->j]+T->i];
            
            Q = 1;
            if (Q > Qpro) {
                done = 1;
            }
            if (done) continue;
            
            Q += Qf51[indx2[T->j]+T->i];
            if (Q > Qpro) {
                pushblock(T->i, T->j, 14);
                done = 1;
            }
            if (done) continue;
            
            for (k1=T->i;k1<T->j;k1++) {
                Q += Qf51[indx2[k1]+T->i] * Qf5[indx2[T->j]+k1+1];
                if (Q > Qpro) {
                    pushblock(T->i, k1, 14);
                    pushblock(k1+1, T->j, 11);
                    done = 1;
                    break;
                }
            }
            //Qf5
        } else if (T->type == 12) {
            seed=1-drand48();
            Qpro = seed * QfML[indx2[T->j]+T->i];
            
            Q =QfML1[indx2[T->j]+T->i];
            if (Q > Qpro) {
                pushblock(T->i, T->j, 15);
                done = 1;
            }
            if (done) continue;
            
            for (k1=T->i;k1<T->j;k1++) {
                G=(T->j-k1)*P->MLbase;
                Q += QfML1[indx2[k1]+T->i]*(QfML[indx2[T->j]+k1+1]);
                if (Q > Qpro) {
                    pushblock(T->i, k1, 15);
                    pushblock(k1+1, T->j, 12);
                    done = 1;
                    break;
                }
                Q += QfML1[indx2[k1]+T->i]*EXP(G);
                if (Q > Qpro) {
                    pushblock(T->i, k1, 15);
                    done = 1;
                    break;
                }
            }
            
            // if (!done) printf("Error! %d\n", T->type);
            
            //QfML
        } else if (T->type == 14) {
            seed=1-drand48();
            Qpro = seed * Qf51[indx2[T->j]+T->i];
            
            for (k1=T->i;k1<T->j;k1++) {
                type = BP_pair[s[k1]][s[T->j]];
                if (type>2) G=P->TerminalAU; else G=0;
                Q += Qc[indx2[T->j]+k1]*EXP(G);
                if (Q > Qpro) {
                    pushblock(k1, T->j, 10);
                    done = 1;
                    break;
                }
            }
            //Qf51
            
        } else if (T->type == 15) {
            seed=1-drand48();
            Qpro = seed * QfML1[indx2[T->j]+T->i];
            
            for (k1=T->i;k1<T->j;k1++) {
                type=BP_pair[s[k1]][s[T->j]];
                if (type) {
                    G=P->MLintern[type]+(k1-T->i)*P->MLbase;
                    Q += Qc[indx2[T->j]+k1]*EXP(G);
                    if (Q > Qpro) {
                        pushblock(k1, T->j, 10);
                        done = 1;
                        break;
                    }
                }
            }
            //QfML1
        }
        
        if (!done) {
            printf("Error! (%d %d)%d\n", T->i, T->j, T->type);
        }
    }
    
    free(s);
}




int *sample_structure(char *seq) //assume pf is computed
{
    int *ptable, l, i;
    
    l = strlen(seq);
    
    ptable = (int *) calloc (l+2, sizeof(int));
    for (i=0; i<=l; i++) ptable[i] = 0;
    
    deblock(seq, ptable);
    
    return ptable;
}

patten *struc_sampler(char *seq, int mul, int *pat)
{
    int i, n = strlen(seq);
    int *ptable;
    char *struc;
    patten *List = NULL;
    
    pfunc(seq, pat);
    
    for (i=1; i<=mul; i++) {
        ptable = sample_structure(seq);
        
        // order = BPair_order(ptable);
        // G = energy_eval(seq, ptable, order);
        
        struc = pair2structure(ptable);
        List = patten_detector(struc, 0, n-1, 1, List);
        
        //free(order);
        //free(ptable);
    }
    
    pfunc_freevar();
    return List;
}



