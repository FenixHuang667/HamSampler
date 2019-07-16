/* Header file for gfold.c */



extern int *symbFold(char *seq, int *FoldEnergy);
extern int LoopEnergy(int n1, int n2, int type, int type_2,
                      int si1, int sj1, int sp1, int sq1);

extern int HairpinE(int size, int type, int si1, int sj1, const char *string);
extern void update_fold_params(void);
extern long double pfunc(char *seq); 