//
//  main.c
//  Inv_seq
//
//  Created by Fenix Huang on 8/19/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include "head.h"


#include "utils.h"
#include "stat.h"
#include "PF.h"
#include "pfunc.h"
#include "simc.h"



char Inputfilename[500];
char **Infile;
int sample, heat, mul;
FILE *Input;
FILE *Output;


void usage()
{
    printf("usage:\n"
           "Smapler [-OPTIONS}\n"
           "[-i inputfilename] Specify the input file (default input.in)\n"
           "[-o outputfilename] Specify the output file (default output.out)\n"
           "[-m number] Specify the sample size (defaout:1000)\n"
           "[-d] Specify the hamming distance\n"
           "[-r] enable to compute EMR-spectrum for <=h (default disable, default h=5)\n"
           //"[-showS] show the sampled sequences"
           //"[-showF] show the folded structures"
           //"[-showE] show the energy"
    );
    exit(0);
}

int main(int argc, const char * argv[]) {
    int i,r;
    char buff[1000];
    int *ptable, *seq_pat, *str_pat, *ref_seq, distance;
    patten *seq_list, *index;
    char structure[1000], sequence_patten[1000], structure_patten[1000], reference_sequence[1000];
    
    srand48(time(NULL)); //random seed
    
    strcpy(Inputfilename, "./input.in");
    strcpy(Outputfilename, "./output.out");
    //defaut input and output
    
    
    sample = 1;
    heat = 0;
    mul = 1000;
    Interval = 7;
    compute_ratio = 1;
    refold = 0;
    show_sequence = 1;
    show_fold_structure =0;
    show_energy=0;
    eval = 1;
    distance = 5;
    
    for (i=1; i<argc; i++) {
        if (argv[i][0]=='-') {
            switch (argv[i][1]) {
                case 'i':
                    if (i==argc-1) usage();
                    Infile = argv[++i];
                    strcpy(Inputfilename, Infile);
                    break;
                case 'o':
                    if (i==argc-1) usage();
                    Outfile = argv[++i];
                    strcpy(Outputfilename, Outfile);
                    break;
                case 'm': if (i==argc-1) usage();
                    sample =1;
                    heat = 0;
                    r=sscanf(argv[++i],"%d", &mul);
                    if (!r) usage();
                    break;
                case 'r':
                    refold = 1;
                    break;
                case 'd':
                    if (i==argc-1) usage();
                    r = sscanf(argv[++i],"%d", &distance);
                    if (!r) usage();
                    break;
                case 's':
                    refold = 1;
                    if (strcmp(argv[i], "-showE")==0) show_energy = 1;
                    if (strcmp(argv[i], "-showF")==0) show_fold_structure = 1;
                    if (strcmp(argv[i], "-showS")==0) show_sequence = 1;
                    break;
                default: usage();
            }
        }
    }
    
    Input=fopen(Inputfilename,"r");
    if (Input==NULL) {
        printf("Input file error!\n");
        exit(0);
    }
    
    Output=fopen(Outputfilename,"w");
    if (Output==NULL) {
        printf("Output file error!\n");
        exit(0);
    }

    while(!feof(Input)) {
        fscanf(Input, "%s", structure);

        fscanf(Input, "%s", reference_sequence);
    }
    
    //ptable = structure2pair(structure);
    //seq_pat = seq2code(sequence_patten);
    //ref_seq = seq2code(reference_sequence);
    
    //printf("%s\n", structure);
    //printf("%s\n", reference_sequence);
    // printf("%d\n", distance);
    
    for (i=0; i<1000; i++) sequence_patten[i] = '_';
    
    if (!refold) {
        seq_list = seq_sampler_ham(structure, reference_sequence, sequence_patten, distance, mul, NULL);
    
        index = seq_list;
        while (index!=NULL) {
            if (index->cnt>0) fprintf(Output, "%s %ld\n", index->pat, index->cnt);
            index = index->next;
        }
    } else if (refold) {
        for (i=1; i<=distance; i++) {
            Inverse_fold_rate(structure, reference_sequence,  sequence_patten, i, mul, Output);
        }
    }
    
    
    fclose(Input);
    fclose(Output);
    
    printf("done!\n");
    
    return 0;
}
