//
//  stat.c
//  Inv_seq
//
//  Created by Fenix Huang on 8/28/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "stat.h"
#include "PF.h", 
#include "pfunc.h"


void ratio (char *seq, int *rat)
{
    int i, n;
    
    n = strlen(seq);
    
    for (i=0; i<5; i++) rat[i]=0;
    for (i=0; i<n; i++){
        switch (seq[i]) {
            case 'A':
                rat[1]++;
                break;
            case 'C':
                rat[2]++;
                break;
            case 'G':
                rat[3]++;
                break;
            case 'U':
                rat[4]++;
                break;
            default:
                rat[0]++;
                break;
        }
    }
}

int string_compare (char *a, char *b) {

    int l, i;
    if (strlen(a)!=strlen(b)) return 0;
    l = strlen(a);
    for (i=0; i<l; i++) {
        if (a[i]!=b[i]) return 0;
    }
    return 1;
}

int p_exist(char *p, patten *List, int add)
{
    patten *indx;
    
    indx = List;
    while (indx != NULL) {
        if (string_compare(indx->pat, p)) {
            if (add) { //if add then add the patten to the list and free p
                indx -> cnt += add;
            }
            return 1;
        }
        indx = indx->next;
    }
    
    return 0;
    
}


patten *patten_detector(char *seq, int i, int j, int frequency, patten *List)
{
    patten *indx, *new;
    char *str;
    int r;
    
    str = (char *) space (sizeof(char)*(j-i+2));
    for (r=0; r<j-i+1; r++) {
        str[r] = seq[i+r];
    }
    str[j-i+1] = 0;
    
    indx = List;
    if (List == NULL) {
        new = (patten *) calloc (1, sizeof(patten));
        new -> pat = str;
        new -> cnt = 1;
        new -> next = List;
        return new;
    }
    while (indx!=NULL) {
        if (!p_exist(str, List, frequency)) {
            new = (patten *) calloc (1, sizeof(patten));
            new -> pat = str;
            new -> cnt = frequency;
            new -> next = List;
            return new;
        } else {
            free(str);
            break;
        }
        indx = indx->next;
    }
    return List;
}

int triplet(char *a, char *b, char *c, int start, int end)
{
    int dab=0, dac=0, dbc=0, i;
    
    for (i = start-1; i<end; i++) {
        if (a[i]!=b[i]) dab++;
        if (b[i]!=c[i]) dbc++;
        if (a[i]!=c[i]) dac++;
    }
    if (dab <= 1 && dbc <=1 && dac <=1) {
        return 2;
    } else if ((dab<=1 && dbc <= 1) || (dbc <= 1 && dac <= 1) || (dab <= 1 && dac <= 1)) {
        return 1;
    } else {
        return 0;
    }
}
//return: 0: not a triplet
//retrun: 1: a triplet but not close
//return: 2: a close triplet


float clustering (char **set, int start, int end, int n)
{
    int i, j, k, indicator;
    float value=0, weight=0, weight_close=0;
    for (i=1; i<=n-2; i++) {
        for (j=i+1; j<=n-1; j++) {
            for (k=j+1; k<=n; k++) {
                indicator = triplet(set[i],set[j],set[k],start, end);
                if (indicator == 2) {
                    weight_close+=1;
                    weight+=1;
                } else if (indicator == 1) {
                    weight+=1;
                }
            }
        }
    }
    if (weight) {
        return weight_close / weight;
    } else
        return 0;
}


double *initial_secondary_structure (int l) // l is the length
{
    int i, j;
    double *Secondary;
    
    Secondary = (double *) space(sizeof(double)*(l+2));
    
    for (i=0;i<=l;i++) Secondary[i]=0;
    Secondary[0]=1;
    Secondary[1]=1;
    Secondary[2]=1;
    Secondary[3]=1;
    Secondary[4]=1;
    
    for (i=5;i<=l;i++) {
        Secondary[i]+=Secondary[i-1]; // allow unpaired vertices
        for (j=3; j<=i-2; j++) {
            Secondary[i]+=Secondary[j]*Secondary[i-2-j];
        }
    }
    //for (i=1;i<=l;i++) printf("%d: %lf\n", i, Secondary[i]);
    return Secondary;
}

void Sample_secodnary_structure(int i, int j, double *Secondary, int *ptable)
{
    int l=j-i+1, r;
    double Qtemp=0, Qv;
    
    if (Secondary[l]==1) return;
    
    Qv = (double) ((1-drand48()) * Secondary[l]);
    Qtemp = Secondary[l-1];
    if (Qv < Qtemp) {
        Sample_secodnary_structure(i+1, j, Secondary, ptable);
        return;
    }
    for (r=3; r<=l-2; r++) {
        Qtemp += Secondary[r] * Secondary[l-2-r];
        if (Qv < Qtemp) {
            ptable[i] = i+r+1;
            ptable[i+r+1] = i;
            Sample_secodnary_structure(i+1, i+r, Secondary, ptable);
            Sample_secodnary_structure(i+r+2, j, Secondary, ptable);
            return;
        }
    }
}

int *SSampler(int l)
{
    int *ptable, i;
    double *Secondary;
    
    ptable = (int *) space(sizeof(int)*(l+2));
    ptable[0]=l;
    for (i=1;i<=l;i++) ptable[i]=0;
    Secondary = initial_secondary_structure (l);
    Sample_secodnary_structure(1, l, Secondary, ptable);
    free(Secondary);
    
    return ptable;
}

void sortList (patten *List)
{
    patten *index, *loop;
    char *str;
    int num;
    
    if (List == NULL) return;
    index = List;
    while (index != NULL) {
        loop = index->next;
        while (loop != NULL) {
            if (index->cnt < loop->cnt) {
                str = index->pat;
                index->pat = loop->pat;
                loop->pat = str;
                
                num = index->cnt;
                index->cnt = loop->cnt;
                loop->cnt=num;
            }
            loop = loop->next;
        }
        index = index->next;
    }
}
