#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "fold_vars.h"
#include "pair_mat.h"
#include "energy_par.h"
#include "params.h"
#include "utils.h"
#include "stat.h"
#include "energy_const.h"

#include "gfold.h"


#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif

#define R GASCONST

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAXENG 100000

#define PUBLIC

int *ptable;

typedef struct block block;

struct block
{
    int i;
    int j;
    int r;
    int s;
    int type;
    // 0:Ibc
    // 1:I5 2:IML //3: PKnot
    
    // 10:c
    // 11:f5 12:fML
    
    // 20:GBond 21:GBond_ML
    // 22:GR_T 23:GR_ML 24:GR_PL
    // 25:GL_T 26:GL_ML 27:GL_PL
    int value;
    block *Next;
};


 int *c,*fML, *fPL,*indx2, *f5, *sym;
 int *DMLi, *DMLi1;
 int *I5, *IML, *IPL, *Ibc, *PKnot;
 unsigned long MEM=0;
 int *cand, *cand_index;

 block *initialblock(void);
 void pushblock(int i, int r, int s, int j, int type, int value);
 block* popblock();

long double *Qc, *Qf5, *Qf51, *QfML, *QfML1, *QI51, *QI5;

PUBLIC block *Bstack;

PUBLIC int *symbFold(char *seq, int *FoldEnergy);

PUBLIC void encode_seq(const char *seq);

INLINE int HairpinE(int size, int type, int si1, int sj1, const char *string);
INLINE int LoopEnergy(int n1, int n2, int type, int type_2,
		      int si1, int sj1, int sp1, int sq1);
void update_fold_params(void);


 void initial (int n)
{
  ptable=(int *) space (sizeof(int)*(n+2));
}

 void initial_sec (int n)
{
	int i,k;

	indx2 = (int *) space(sizeof(int)*(n+5));
	for (i = 1; i <= n; i++)
	  indx2[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */

	c     = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	f5    = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	DMLi    = (int *) space(sizeof(int)*(n+3));
	DMLi1    = (int *) space(sizeof(int)*(n+3));
	sym    = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	fML   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	//fPL   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	
	MEM+=4*(sizeof(int)*((n*(n+1))/2+2))+2*(sizeof(int)*(n+3))+(sizeof(int)*(n+5));
}


 void freevar()
{
  free(ptable);
}

 void freevar_sec()
{
  free (f5);
  free (c);
  free (DMLi);
  free (DMLi1);
  free (fML);
  free (sym);
}

 void fill_array_sec(const char *seq)
{
	int i,j,k1,k2,p,q,ij;
	int length,type2,G,type;
	int new_c, new_ML, new_P, new_f5;
	int *FF;

	length=strlen(seq);
	
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

	    type = BP_pair[S[i]][S[j]];
	    ij=indx2[j]+i;

	    if (type) {
	      
	      new_c = HairpinE(j-i-1, type, S[i+1], S[j-1], seq+i-1); //Hairpin loop

	      for (p = i+1; p <=  MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
	        int minq = j-i+p-MAXLOOP-2;
	        if (minq<p+1+TURN) minq = p+1+TURN;
	        for (q = minq; q < j; q++) {
	          type2 = BP_pair[S[p]][S[q]];

	          if (type2==0) continue;
	          type2 = rtype[type2];

	          G = LoopEnergy(p-i-1, j-q-1, type, type2, S[i+1], S[j-1], S[p-1], S[q+1]);
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
	      type2 = BP_pair[S[i]][S[k1]];
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
	      //if (!sym[indx2[j]+k1+1]) continue;
	      type2=BP_pair[S[k1+1]][S[j]];
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



 void de_I(block *T, const char *seq)
{
	int p,q,k1,k2,minq,maxp,G,V,i,j,r,s,ij;
	int typeij, type_l, type_r;

	if (T->value==0) return;
	
	i=T->i;
	j=T->j;
	r=T->r;
	s=T->s;
	ij=indx2[j]+i;

	switch (T->type) {
	  
	case 0:   //decompose Ibc
	  typeij = BP_pair[S[i]][S[j]];
	  if (typeij) {
	    ptable[i]=j;
	    ptable[j]=i;
	    
	    if (T->value==HairpinE(j-i-1, typeij, S[i+1], S[j-1], seq+i-1)) {
	      return;
	    }
	    for (p = i+1; p <=  MIN2(j-2-TURN,i+MAXLOOP+1); p++) { //Interior loop
	      minq = j-i+p-MAXLOOP-2;
	      if (minq<p+1+TURN) minq = p+1+TURN;
	      for (q = minq; q < j; q++) {
	        type_l = BP_pair[S[q]][S[p]];
	        if (type_l==0) continue;
	        G = LoopEnergy(p-i-1, j-q-1, typeij, type_l, S[i+1], S[j-1], S[p-1], S[q+1]);
	        if (T->value==G+Ibc[indx2[q]+p]) {
		  pushblock(p,0,0,q,0,Ibc[indx2[q]+p]);
		  return;
		}
	      }
	    }
	    
	    if (T->value==P->MLintern[typeij]+P->MLclosing+IML[indx2[j-1]+i+1]) {
	      pushblock(i+1,0,0,j-1,2,IML[indx2[j-2]+i+1]);
	      return;
	    }
	  } else {
	    printf("Error! Do not match the energy, please check your program carefully!\n");
	    return;
	  }
	  break;
	case 1:  //decompose I5
          if (j-i<=TURN+1) {
	    return;
	  }
	  
	  if (T->value==I5[ij+1]) {
	    pushblock(i+1,0,0,j,1,I5[ij+1]);
	    return;
	  }
	  typeij = BP_pair[S[i]][S[j]];
	  if (typeij>2) G=P->TerminalAU; else G=0;
	  if (T->value==G+Ibc[ij]) {
	    pushblock(i,0,0,j,0,Ibc[ij]);
	    return;
	  }
	  if (T->value==PKnot[ij]) {
	    pushblock(i,0,0,j,3,PKnot[ij]);
	    return;
	  }
	  for (k1=i+TURN+1;k1<j;k1++) {
	    type_l = BP_pair[S[i]][S[k1]];
	    if (type_l>2) G=P->TerminalAU; else G=0;
	    if (type_l && T->value==Ibc[indx2[k1]+i]+I5[indx2[j]+k1+1]+G) {
	      pushblock(i,0,0,k1,0,Ibc[indx2[k1]+i]);
	      pushblock(k1+1,0,0,j,1,I5[indx2[j]+k1+1]);
	      return;
	    }
	  }
	  for (k1=i+TURN+1;k1<j;k1++) {
	    if (T->value==PKnot[indx2[k1]+i]+I5[indx2[j]+k1+1]) {
	      pushblock(i,0,0,k1,3,PKnot[indx2[k1]+i]);
	      pushblock(k1+1,0,0,j,1,I5[indx2[j]+k1+1]);
	      return;
	    }
	  }
	  break;
	case 2:  //decompose IML
          if (j-i<=TURN+1) {
	    return;
	  }
	  
	  if (T->value==IML[ij+1]+P->MLbase) {
	    pushblock(i+1,0,0,j,2,IML[ij+1]);
	    return;
	  }
	  if (T->value==IML[indx2[j-1]+i]+P->MLbase) {
	    pushblock(i,0,0,j-1,2,IML[indx2[j-1]+i]);
	    return;
	  }
	  typeij = BP_pair[S[i]][S[j]];
	  if (T->value==Ibc[ij]+P->MLintern[typeij]){
	    pushblock(i,0,0,j,0,Ibc[ij]);
	    return;
	  }
	  for (k1=i+TURN+1;k1<j-TURN-1;k1++) {
	    if (T->value==IML[indx2[k1]+i]+IML[indx2[j]+k1+1]) {
	      pushblock(i,0,0,k1,2,IML[indx2[k1]+T->i]);
	      pushblock(k1+1,0,0,j,2,IML[indx2[j]+k1+1]);
	      return;
	    }
	  }
	  break;
	case 10: //decompose c
	  typeij=BP_pair[S[i]][S[j]];
	  if (typeij) {
	    ptable[i]=j;
	    ptable[j]=i;
	    
	    if (T->value==HairpinE(j-i-1, typeij, S[i+1], S[j-1], seq+i-1)) {
	      return;
	    }
	    
	    for (p=i+1;p<=MIN2(j-2-TURN,i+MAXLOOP+1);p++) {
	      minq = j-i+p-MAXLOOP-2;
	      if (minq<p+1+TURN) minq = p+1+TURN;
	      for (q=j-1; q>=minq; q--) {
	        type_l=BP_pair[S[q]][S[p]];
	        if (type_l==0) continue;
	        G=LoopEnergy(p-i-1, j-q-1, typeij, type_l, S[i+1], S[j-1], S[p-1], S[q+1]);
	        if (T->value==G+c[indx2[q]+p]) {
	          pushblock(p,0,0,q,10,c[indx2[q]+p]);
	          return;
	        }
	      }
	    }
	    
          for (k1=j-1; k1>i+TURN; k1--) {
              for (k2=k1-TURN; k2>i; k2--) {
                  if (!BP_pair[S[k2]][S[k1]]) continue;
                  type_l = BP_pair[S[k2]][S[k1]];
                  if (T->value==P->MLintern[typeij]+P->MLclosing+fML[indx2[k2-1]+i+1]+
                      P->MLintern[type_l]+c[indx2[k1]+k2]+(j-k1-1)*P->MLbase) {
                  pushblock(i+1,0,0,k2-1,12,fML[indx2[k2-1]+i+1]);
                  pushblock(k2,0,0,k1,10,c[indx2[k1]+k2]);
                  return;
                }
              }
	    }
	  } else {
	    printf( "Error! Do not match the energy, please check your program carefully!\n");
	    return;
	  }
	  break;
	  	  
	case 11: //decompose f5
          if (j-i<=TURN+1) {
	    return;
	  }
	  
          typeij=BP_pair[S[i]][S[j]];
	  if (typeij>2) G=P->TerminalAU; else G=0;
	  if (typeij && T->value==c[ij]+G) {
	    pushblock(i,0,0,j,10,c[ij]);
	    return;
	  }
	  if (T->value==f5[ij+1]) {
	    pushblock(i+1,0,0,j,11,f5[ij+1]);
	    return;
	  }
	  for (k1=i+TURN+1;k1<j;k1++) {
	    if (T->value==f5[indx2[k1]+i]+f5[indx2[j]+k1+1]) {
	      pushblock(i,0,0,k1,11,f5[indx2[k1]+i]);
	      pushblock(k1+1,0,0,j,11,f5[indx2[j]+k1+1]);
	      return;
	    }
	  }
	  break;
	case 12: //decompose fML
          if (j-i<=TURN+1) return;
	  
          typeij=BP_pair[S[i]][S[j]];

	  if (typeij && T->value==P->MLintern[typeij]+c[ij]) {
	    pushblock(i,0,0,j,10,c[ij]);
	    return;
	  }

	  if (T->value==fML[ij+1]+P->MLbase) {
	    pushblock(i+1,0,0,j,12,fML[ij+1]);
	    return;
	  }

	  if (T->value==fML[indx2[j-1]+i]+P->MLbase) {
	    pushblock(i,0,0,j-1,12,fML[indx2[j-1]+i]);
	    return;
	  }
	  for (k1=T->i+TURN+1;k1<T->j-TURN-1;k1++) {
	    if (T->value==fML[indx2[k1]+i]+fML[indx2[j]+k1+1]) {
	      pushblock(i,0,0,k1,12,fML[indx2[k1]+i]);
	      pushblock(k1+1,0,0,j,12,fML[indx2[j]+k1+1]);
	      return;
	    }
	  }
	  break;  
	default : break;
	}

	printf("Error! Do not match the energy, please check your program carefully!\n");
	printf("%d\n", T->type);
	printf("%d\n", T->value);
}


int * symbFold(char *seq, int *FoldEnergy)
{	
  int length,i,j,r,s,si,sj,sr,ss,correct, total,total_nat, cand=0;
  int score=0,G,type1,type2,new_score;
  long index;
  char *struc;
  block *T, *IN;

  length=strlen(seq);
  initial(length);

  update_fold_params();
  P = scale_parameters();
  encode_seq(seq);
  


    MEM=0;
    initial_sec(length);
    fill_array_sec(seq);

    pushblock(1,0,0,length,11,f5[indx2[length]+1]);
    while (Bstack!=NULL) {
      T=popblock();
      de_I(T,seq);
/*
 	  printf("-----------------------------------\n");
 	  IN=Bstack;
 	  while (IN!=NULL) {
 	    printf("(%d %d %d %d) %d %d\n", IN->i, IN->j, IN->r, IN->s, IN->type,IN->value);
 	    IN=IN->Next;
 	  }
*/
    }
    ptable[0]=length;
    //printf("%d\n", f5[indx2[length]+1]);
    *FoldEnergy = f5[indx2[length]+1]; 
        
    freevar_sec();
    return ptable;
}

INLINE int HairpinE(int size, int type, int si1, int sj1, const char *string) {
  int energy;
  energy = (size <= 30) ? P->hairpin[size] :
    P->hairpin[30]+(int)(P->lxc*log((size)/30.));
  if (tetra_loop)
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
	energy += P->TETRA_ENERGY[(ts - P->Tetraloops)/7];
    }
  if (size == 3) {
    char tl[6]={0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(P->Triloops, tl)))
      energy += P->Triloop_E[(ts - P->Triloops)/6];

    if (type>2)  /* neither CG nor GC */
      energy += P->TerminalAU; /* penalty for closing AU GU pair */
  }
  else  /* no mismatches for tri-loops */
    energy += P->mismatchH[type][si1][sj1];

  return energy;
}


INLINE int LoopEnergy(int n1, int n2, int type, int type_2,
		      int si1, int sj1, int sp1, int sq1) {
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];    /* stack */

  if (ns==0) {                       /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                             /* interior loop */
    if (ns==1) {
      if (nl==1)                     /* 1x1 loop */
	return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                   /* 2x1 loop */
	if (n1==1)
	  energy = P->int21[type][type_2][si1][sq1][sj1];
	else
	  energy = P->int21[type_2][type][sq1][si1][sp1];
	return energy;
      }
    }
    else if (n1==2 && n2==2)         /* 2x2 loop */
      return P->int22[type][type_2][si1][sp1][sq1][sj1];
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->F_ninio[2]);

      energy += P->mismatchI[type][si1][sj1]+
	P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

void update_fold_params(void)
{
  P = scale_parameters();
  make_pair_matrix();
  //if (init_length < 0) init_length=0;
}
block *initialblock(void)
{
	block *new;

	new=(block *) space (sizeof(block));
	new->i=0; new->j=0;
	new->r=0; new->s=0;
	new->type=-1;
	new->value=MAXENG;
	new->Next=NULL;
	return new;
}
void pushblock(int i, int r, int s, int j, int type, int value)
{
	block *T;

	T=initialblock();
	T->i=i; T->j=j;
	T->r=r; T->s=s;
	T->type=type;
	T->value=value;
	T->Next=Bstack;
	Bstack=T;
	return;
}
block* popblock()
{
	block *T;
	T=Bstack;
	Bstack=Bstack->Next;
	T->Next=NULL;
	return T;
}


PUBLIC void encode_seq(const char *seq) {
  unsigned int i,l;

  l = strlen(seq);
  S = (int *) space(sizeof(int)*(l+2));
//  S1= (short *) space(sizeof(short)*(l+2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(seq[i-1]));
  }
    //S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  /* for circular folding add first base at position n+1 */
  //S[l+1] = S[1]; S1[l+1]=S1[1];
}




long double *Qc, *Qf5, *Qf51, *QfML, *QfML1, *QI51, *QI5;
int *indx2;


long double EXPP(int G)
{
    return (long double) exp(-1*(10*(G))/(R*(P->temperature+K0)));
}


void pfunc_initial (int n)
{
    int i;
    
    indx2 = (int *) space(sizeof(int)*(n+5));
    for (i = 1; i <= n; i++)
        indx2[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */
    
    Qc     = (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    Qf5    = (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    Qf51   = (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    QfML   = (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    QfML1  = (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    QI51= (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    QI5= (long double *) space(sizeof(long double)*((n*(n+1))/2+2));
    
}

void pfunc_freevar()
{
    free(Qc); free(QfML); free(QfML1);
    free(Qf5); free(Qf51);
    free(QI51); free(QI5);
    free(indx2);
}


void partition(char *seq)
{
    int i,j,k1,p,q;
    int length,type2,G,type;
    
    length=strlen(seq);
    for (i=1;i<=length;i++)
        for (j=i;j<i+TURN;j++) {
            Qc[indx2[j]+i]=0;
            Qf5[indx2[j]+i]=1;
            Qf51[indx2[j]+i]=0;
            QfML[indx2[j]+i]=0;
            QfML1[indx2[j]+i]=0;
        }
    //#pragma omp parallel for
    for (i=length-TURN;i>=1;i--)
        for (j=i+TURN;j<=length;j++) {
            
            // 	for (gap=TURN;gap<length;gap++)
            //         #pragma omp parallel default (shared) private \
            //  (i,j,r,s,k1,k2,p,q,type1, type2,typek,G,type,left,right,index,indexk)
            //         #pragma omp for schedule (static)
            // 	  for (i=1;i<length;i++) {
            // 	    j=i+gap;
            
            // printf("%d %d\n", i,j);
            type = BP_pair[S[i]][S[j]];
            if (type) {
                G=HairpinE(j-i-1, type, S[i+1], S[j-1], seq+i-1);
                Qc[indx2[j]+i]=EXPP(G);
                for (p = i+1; p <= j-2-TURN; p++)
                    for (q = p+TURN; q < j; q++) {
                        type2 = BP_pair[S[p]][S[q]];
                        
                        if (type2==0) continue;
                        type2 = rtype[type2];
                        
                        G = LoopEnergy(p-i-1, j-q-1, type, type2,
                                       S[i+1], S[j-1], S[p-1], S[q+1]);
                        Qc[indx2[j]+i]+=EXPP(G)*Qc[indx2[q]+p];
                    }
                for (k1=i+1;k1<j-1;k1++) {
                    G=P->MLintern[type]+P->MLclosing;
                    Qc[indx2[j]+i]+=QfML1[indx2[k1]+i+1]*QfML[indx2[j-1]+k1+1]*EXPP(G);
                }
            } else
                Qc[indx2[j]+i]=0;
            
            Qf51[indx2[j]+i]=0;
            for (k1=i;k1<j;k1++) {
                type = BP_pair[S[k1]][S[j]];
                if (type>2) G=P->TerminalAU; else G=0;
                if (type) Qf51[indx2[j]+i]+=Qc[indx2[j]+k1]*EXPP(G);
                
            }
            
            Qf5[indx2[j]+i]=Qf51[indx2[j]+i]+1;
            for (k1=i;k1<j;k1++) {
                Qf5[indx2[j]+i]+=Qf51[indx2[k1]+i]*Qf5[indx2[j]+k1+1];
            }
            
            QfML1[indx2[j]+i]=0;
            for (k1=i;k1<j;k1++) {
                type=BP_pair[S[k1]][S[j]];
                if (type) {
                    G=P->MLintern[type]+(k1-i)*P->MLbase;
                    QfML1[indx2[j]+i]+=Qc[indx2[j]+k1]*EXPP(G);
                }
            }
            
            QfML[indx2[j]+i]=QfML1[indx2[j]+i];
            for (k1=i;k1<j;k1++) {
                G=(j-k1)*P->MLbase;
                QfML[indx2[j]+i]+=QfML1[indx2[k1]+i]*(QfML[indx2[j]+k1+1]+EXPP(G));
            }
            // partition for secondary structure S done!
            
        }
}
long double pfunc(char *seq)
{
    int length;
    long double PF;
    
    length=strlen(seq);
    update_fold_params();
    P = scale_parameters();
    encode_seq(seq);
    
    pfunc_initial(length);
    
    partition(seq);
    PF = Qf5[indx2[length]+1];
    
    pfunc_freevar();
    return PF;
    
}

