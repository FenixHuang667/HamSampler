//
//  simc.c
//  Inv_seq
//
//  Created by Fenix Huang on 4/28/17.
//  Copyright (c) 2017 Fenix Huang. All rights reserved.
//

#include "simc.h"
#include "pfunc.h"
#include "PF.h"
#include "stat.h"

void Compare_Boltzmann(char *seqA, char *seqB, int *struc_pat, FILE *out);



void Overlap(char *seqA, char *seqB, char *sequence_patten, int *struc_pat, FILE *out)
{
    patten *struc_list, *seq_list, *temp_list, *index, *index2;
    int threshold, energy;
    int *order, *ptable, n = strlen(seqA);
    
    struc_list = struc_sampler(seqA, 1000, struc_pat);
    
    sortList(struc_list);
    ptable = structure2pair(struc_list->pat);
    order = BPair_order(ptable);
    energy = energy_eval(seqA, ptable, order);
    
    free(ptable);
    free(order);
    
    threshold = energy*0.95;
    
    seq_list = NULL;
    index = struc_list;
    while (index!=NULL) {
        ptable = structure2pair(index->pat);
        order = BPair_order(ptable);
        energy = energy_eval(seqA, ptable, order);
        
        if (energy < threshold) {
            fprintf(out, "%s %ld %d\n\n", index->pat, index->cnt, energy);
            
            temp_list = seq_sampler_ham(index->pat, seqA, sequence_patten, 5, 1000, out);
            
            index2 = temp_list;
            while (index2!=NULL) {
                seq_list = patten_detector(index2->pat, 0, n-1, index2->cnt, seq_list);
                index2 = index2->next;
            }
        }
        free(ptable);
        free(order);
        index = index->next;
    }
    
    index = seq_list;
    while (index!=NULL) {
        fprintf(out, "%s %ld\n", index->pat, index->cnt);
        index = index->next;
    }
}

void Compare_Boltzmann(char *seqA, char *seqB, int *struc_pat, FILE *out)
{
    patten *struc_listA, *struc_listB, *index;
    int *order, *ptable, n = strlen(seqA), energy, threshold;
    
    
    struc_listA = struc_sampler(seqA, 1000, struc_pat);
    struc_listB = struc_sampler(seqB, 1000, struc_pat);
    
    sortList(struc_listA);
    ptable = structure2pair(struc_listA->pat);
    order = BPair_order(ptable);
    energy = energy_eval(seqA, ptable, order);
    
    threshold = energy*0.95;
    
    index = struc_listA;
    while (index!=NULL) {
        ptable = structure2pair(index->pat);
        order = BPair_order(ptable);
        energy = energy_eval(seqA, ptable, order);
        
        if (energy < threshold) {
            fprintf(out, "%s %ld %d\n\n", index->pat, index->cnt, energy);
        }
        free(ptable);
        free(order);
        index = index->next;
    }
    
    fprintf(out, "----------------------------\n");
    
    index = struc_listB;
    while (index!=NULL) {
        ptable = structure2pair(index->pat);
        order = BPair_order(ptable);
        energy = energy_eval(seqB, ptable, order);
        
        if (energy < threshold) {
            fprintf(out, "%s %ld %d\n\n", index->pat, index->cnt, energy);
        }
        free(ptable);
        free(order);
        index = index->next;
    }
}


void Inverse_fold_rate(char *struc, char *ref_seq, char *sequence_patten, int distance, int mul, FILE *out)
{
    patten * seqList = NULL, *index, *shift_index, *shift_seq, *rand_seq, *rand_seq_index, *rand_seq_ham, *rand_seq_ham_index;
    int energy, ptable[1000], ref_ptable[1000], ref_energy;
    int i, cnt, num_seq, over_cnt;
    char *new_struc, *ref_struc;
    float ref_rate, rate;
    int neighbor = 0;
    
    ref_energy = fold_sec(ref_seq, ref_ptable);
    
    seqList = seq_sampler_ham(struc, ref_seq, sequence_patten, distance, mul, NULL);
    // ref_struc = pair2structure(ref_ptable);
    //fprintf(out, "%s\n", ref_struc);
    
    
    // seq_list = seq_sampler_ham(structure, reference_sequence, sequence_patten, 10, 1000, Output);
    cnt = 0;
    index = seqList;
    while (index!=NULL) {
        energy = fold_sec(index->pat, ptable);
        if (idetical_struc(ptable, ref_ptable)) {
            // fprintf(out, "%s\n", index->pat);
            cnt+=index->cnt;
        }
        
        //new_struc = pair2structure(ptable);
        //fprintf(out, "%s\n", new_struc);
        index = index->next;
    }
    
    ref_rate = ((float) cnt / mul);
    fprintf(out, "EMR of the input sequence at distance %d: %d/%d = %f\n", distance, cnt, mul, ref_rate);
    
    cnt = 0;
    rand_seq = sample_compatible_seq(struc, 100);
    rand_seq_index = rand_seq;
    while (rand_seq_index!=NULL) {
        energy = fold_sec(rand_seq_index->pat, ptable);
        if (idetical_struc(ptable, ref_ptable)) {
            fprintf(out, "%s\n", rand_seq_index->pat);
            cnt+=rand_seq_index->cnt;
        }

        rand_seq_index = rand_seq_index->next;
    }
    
    fprintf(out, "EMR of the random sampled sequences at distance %d: %d/%d = %f\n", distance, cnt, mul, ((float) cnt / mul));

    
    if (neighbor) {
        shift_seq = seq_sampler_ham(struc, ref_seq, sequence_patten, 5, 50, NULL);
    
        num_seq = 0;
        over_cnt = 0;
        shift_index = shift_seq;
        while (shift_index!=NULL) {
            num_seq++;
            fprintf(out, "Seq: %d, fre: %d\n", num_seq, shift_index->cnt);
            seqList = seq_sampler_ham(struc, shift_index->pat, sequence_patten, distance, mul, NULL);
        
            cnt = 0;
            index = seqList;
            while (index!=NULL) {
                energy = fold_sec(index->pat, ptable);
                if (idetical_struc(ptable, ref_ptable)) cnt+=index->cnt;
                index = index->next;
            }
        
            rate = ((float) cnt / mul);
            fprintf(out, "%d/%d = %f\n", cnt, mul, rate);
            if (rate > ref_rate) over_cnt+=shift_index->cnt;
            shift_index = shift_index->next;
        }
        fprintf(out, "Over: %d\n", over_cnt);
    }
    
}



