/* Header file for utils.c */
#ifdef HAVE_CONFIG_H
#endif
#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include <stdio.h>
#include <stdlib.h>
#include "dmalloc.h"
#define space(S) calloc(1,(S))
#else
extern /*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
                            /* allocate space safely */ 
extern /*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

extern /*@exits@*/ void nrerror(const char message[]);  /* die with error message */
extern void   init_rand(void);                /* make random number seeds */
extern unsigned short xsubi[3];               /* current 48bit random number */
extern double urn(void);                      /* random number from [0..1] */
extern int    int_urn(int from, int to);      /* random integer */
extern /*@observer@*/ char  *time_stamp(void);               /* current date in a string */
extern /*@only@*/ /*@notnull@*/ char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
extern int    hamming(const char *s1, const char *s2);
/* calculate hamming distance */

extern char *pack_structure(const char *struc);
/* pack secondary secondary structure, 5:1 compression using base 3 encoding */
extern char *unpack_structure(const char *packed);
/* unpack sec structure packed with pack_structure() */
extern short *make_pair_table(const char *structure);
/* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
   0 if i is unpaired, table[0] contains the length of the structure. */

extern int bp_distance(const char *str1, const char *str2);
/* dist = {number of base pairs in one structure but not in the other} 
   same as edit distance with open-pair close-pair as move-set */

extern void randomseq(int n, char *res); 
extern char *pair2structure(int *pair);
extern int *structure2pair(char *s);
extern int compare(int *A, int *B);
extern int number(int *A);
extern int *read_patten(char *pat);
extern int * seq2code(char *seq);
extern char code2nucleotide (int a);
extern int nucleotide2code (char a);

extern int idetical_struc(int *pa, int *pb);
