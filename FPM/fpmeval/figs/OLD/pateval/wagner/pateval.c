#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "analyzestat.h"

// we assume a maximum size itemset and items sorted
// eventually it should be replaced to read the output file directly

#define MAXISET 10
#define NULLITEM -1

typedef struct names{
int code;
char * name; 
} names_t, *ptr_names_t;

names_t namesvet[]={
{11, "sl1"},
{12, "sl2"},
{13, "sl3"},
{21, "sw1"},
{22, "sw2"},
{23, "sw3"},
{31, "pl1"},
{32, "pl2"},
{33, "pl3"},
{41, "pw1"},
{42, "pw2"},
{43, "pw3"},
{51, "cl1"},
{52, "cl2"},
{53, "cl3"},
{-1, " "}
};

typedef struct itemset{
  int nitems;
  int items[MAXISET];
  int support;
} itemset_t, * ptr_itemset_t;

itemset_t isinput[] ={
#include "iset.db"
{-1, {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},-1}
};

typedef struct htnode{
  struct htnode * next;
  struct htnode * children;
  int count;
  int item;
} htnode_t, * ptr_htnode_t;


// global variables
ptr_htnode_t globalt;
int numtrans;

char * getname(int code){
  int i;
  i = 0;
  while (namesvet[i].code>0){
    if (namesvet[i].code == code)
      return namesvet[i].name;
    i++;
  }
  return namesvet[i].name;
}

int printht(ptr_htnode_t t, int level, ptr_itemset_t piaux){
  ptr_htnode_t paux;
  int i;
  char mask[20];
  mask[0]=0;
  for(i=0;i<=level;i++){
    strcat(mask," ");
  }
  printf("%s",mask);
  for (i=0;i<piaux->nitems;i++){ 
    printf("I%d ",piaux->items[i]);
  }
  printf("I%d (%d)\n",t->item,t->count);
  piaux->items[piaux->nitems] = t->item;
  piaux->nitems++;
  paux = t->children;
  while (paux){
    printht(paux,level+1,piaux);
    paux=paux->next;
  }  
  piaux->nitems--;
  return 0;
}

void printv(int mask[], int n) {
	int i;
	printf("{ ");
	for (i = 0; i < n; ++i)
		if (mask[i])
			printf("%d ", i + 1); /*i+1 is part of the subset*/
	printf("\b }\n");
}

int complsubset(int mask[], int n, int compl[]){
  int i;
  for (i = 0; (i < n); ++i)
    compl[i] = 1 - mask[i];
}  

int nextsubset(int mask[], int n) {
  int i;
  for (i = 0; (i < n) && mask[i]; ++i)
    mask[i] = 0;

  if (i < n) {
    mask[i] = 1;
    return 1;
  }
  return 0;
}

// print the itemsets
int printitemsets(ptr_itemset_t piaux, int mask[MAXISET], int compl[MAXISET], int n ){
  int i;
  if (n>piaux->nitems){
    fprintf(stderr,"Error: printitemsets nitems=%d n=%d\n", piaux->nitems, n);
    exit(1);
  }
  printf("Mask { ");
  for (i = 0; i < n; ++i)
    if (mask[i]) printf(" %d %s",piaux->items[i+1],getname(piaux->items[i+1]));
  printf("} Compl{");
  for (i = 0; i < n; ++i)
    if (compl[i]) printf(" %d %s",piaux->items[i+1],getname(piaux->items[i+1]));
  printf("}\n");
  return 1; 
} 

// create an itemset in piout according to piaux and mask
int mountitemset(ptr_itemset_t piout, ptr_itemset_t piaux, int mask[MAXISET], int n){
  int i,cc=0;
  for (i = 0; i < n; ++i)
    if (mask[i]){
     piout->items[cc] = piaux->items[i+1];
     cc++;
    }
  piout->nitems = cc;
  //printf("mountitemset %d\n",cc);
  return cc;
}


int getsupport(ptr_htnode_t t, ptr_itemset_t pi, int level){
  ptr_htnode_t paux;
  ptr_htnode_t found;
  int i;
#ifdef DEBUGGETSUPP
        printf("getsupport: Recursion %d \n",level);
#endif

  // locate the node, if it exists
  paux = t->children;
  found = (ptr_htnode_t) 0;
  while (paux){
    if (level==0){ // root. it has just one node
      found = paux; 
    } else{ 
#ifdef DEBUGGETSUPP
        printf("Testing Level %d Item %d \n",level,paux->item);
#endif
      if (paux->item == pi->items[level-1]){
        found = paux; 
#ifdef DEBUGGETSUPP
        printf("Found Level %d Item %d \n",level,paux->item);
#endif
      }
    }
    paux = paux->next;
  } 

  // found now has the item. continue recursively
  if (level < pi->nitems){
#ifdef DEBUGGETSUPP
    printf("Recursion Level %d\n",level);
#endif
    return getsupport(found, pi, level+1);
#ifdef DEBUGGETSUPP
    printf("Back from Recursion Level %d\n",level);
#endif
  }

  // we reached a itemset
  if (level == pi->nitems){
    return found->count;
#ifdef DEBUGGETSUPP
    printf("Reached a leaf Level %d Item %d Count %d\n",level,found->item,found->count);
#endif
  }
  return -1; 
}

int bitson(int mask[MAXISET], int n){
  int i,cc=0;
  for (i=0; i<n; i++){
    cc += (mask[i]>0);
  }
  return cc;
}

int traverseht(ptr_htnode_t t, int level, ptr_itemset_t piaux, int support){
  ptr_htnode_t paux;
  itemset_t pimask, picompl;
  int masksupp, complsupp, subsetmask[MAXISET], subsetcompl[MAXISET];
  float probt, probmask, probcompl;
  int i;
  hist_t hist;
  float auxvet[1000];
  int szauxvet;
  float avgvet,varvet;
  char buff[1000];

  if (t->count>0){
    piaux->items[piaux->nitems] = t->item;
    piaux->nitems++;
  }
  sprintf(buff, " %d supp ",t->count);
  // at this point, we do have an itemset in piaux, time to check.
  for (i=0;i<piaux->nitems;i++){ 
    printf("I%d %s ",piaux->items[i],getname(piaux->items[i]));
    sprintf(buff, "%s I%d %s ",buff,piaux->items[i],getname(piaux->items[i]));
  }
  printf("(%d)\n",t->count);
  probt = t->count*1.0/numtrans;
  // first generate the subsets
  if (level >=1) {
    for (i = 0; i < MAXISET; ++i) subsetmask[i] = 0;
    
    szauxvet = 0; 

    // empty set is the first
    //printf("mask "); printv(subsetmask, level-1);
    //printf("comp "); printv(subsetcompl, level-1);
    
    if (bitson(subsetmask,level-1) && bitson(subsetcompl,level-1)) {
      printitemsets(piaux, subsetmask, subsetcompl, level-1);
      mountitemset(&pimask, piaux, subsetmask, level-1);
      mountitemset(&picompl, piaux, subsetcompl, level-1);
      masksupp = getsupport(globalt, &pimask,0);
      complsupp = getsupport(globalt, &picompl,0);
      //printf("supports mask %d compl %d\n",masksupp, complsupp); 
      probmask = masksupp*1.0/numtrans; 
      probcompl = complsupp*1.0/numtrans; 
      printf("lift = %f / (%f X %f)= %f\n", probt, probmask, probcompl, probt/(probmask*probcompl));
      auxvet[szauxvet] = probt/(probmask*probcompl);
      szauxvet++;
    }
  
    while (nextsubset(subsetmask, level-1)){
      complsubset(subsetmask, level-1, subsetcompl);
      //printf("mask "); printv(subsetmask, level-1);
      //printf("comp "); printv(subsetcompl, level-1);
      if (bitson(subsetmask,level-1) && bitson(subsetcompl,level-1)) {
        printitemsets(piaux, subsetmask, subsetcompl, level-1);
        mountitemset(&pimask, piaux, subsetmask, level-1);
        mountitemset(&picompl, piaux, subsetcompl, level-1);
        masksupp = getsupport(globalt, &pimask,0);
        complsupp = getsupport(globalt, &picompl,0);
        //printf("supports mask %d compl %d\n",masksupp, complsupp); 
        probmask = masksupp*1.0/numtrans; 
        probcompl = complsupp*1.0/numtrans; 
        printf("lift = %f / (%f X %f)= %f\n", probt, probmask, probcompl, probt/(probmask*probcompl));
        auxvet[szauxvet] = probt/(probmask*probcompl);
        szauxvet++;
      }
    }
    if (szauxvet>1){
      analyzevector(szauxvet,auxvet, &avgvet, &varvet, &hist);
      printf("%f avg %f var %d sz %s\n", avgvet, varvet, szauxvet, buff); 
    }
  }
  // then, for each subset, calculate the complement
  // calculate the measures for both
  // calculate the statistics for them

  paux = t->children;
  while (paux){
    traverseht(paux,level+1,piaux,support);
    paux=paux->next;
  }  
  piaux->nitems--;
  return 0;
}

int insertht(ptr_htnode_t t, ptr_itemset_t pi, int level){
  ptr_htnode_t paux;
  ptr_htnode_t found;
  int i;
  

#ifdef DEBUGINS
  printf("Enter insertht Level %d n %d (", level, pi->nitems);
  for(i=0; i<pi->nitems; i++){
    printf(" %d ",pi->items[i]);
  }
  printf(") %d Locate %d\n", pi->support, pi->items[level-1]);
#endif

  // locate the node, if it exists
  paux = t->children;
  found = (ptr_htnode_t) 0;
  while (paux){
    if (level==0){ // root. it has just one node
      found = paux; 
    } else{ 
      if (paux->item == pi->items[level-1]){
        found = paux; 
#ifdef DEBUGINS
        printf("Found Level %d Item %d \n",level,paux->item);
#endif
      }
    }
    paux = paux->next;
  } 

  // item not in the list. insert it
  if (!found){
    paux = t->children;
    found = (ptr_htnode_t) malloc(sizeof(htnode_t));
    found->next = paux;
    t->children = found;
    if (level>0)
      found->item = pi->items[level-1];
    else
      found->item = NULLITEM;
    found->count = 0;
    found->children = (ptr_htnode_t) 0;
#ifdef DEBUGINS
    printf("Inserting Level %d Item %d\n",level,found->item);
#endif
  }

  // found now has the item. continue recursively
  if (level < pi->nitems){
#ifdef DEBUGINS
    printf("Recursion Level %d\n",level);
#endif
    insertht(found, pi, level+1);
#ifdef DEBUGINS
    printf("Back from Recursion Level %d\n",level);
#endif
  }

  // we reached a leaf
  if (level == pi->nitems){
    found->count = pi->support;
#ifdef DEBUGINS
    printf("Reached a leaf Level %d Item %d Count %d\n",level,found->item,found->count);
#endif
  }
#ifdef DEBUGINS
  printf("Exit insertht Level %d n %d (", level, pi->nitems);
  for(i=0; i<pi->nitems; i++){
    printf(" %d ",pi->items[i]);
  }
  printf(") %d\n", pi->support);
#endif
  return 0;
}



int main(){

itemset_t iaux;
int i,support=3;
htnode_t htree;

htree.children = (ptr_htnode_t) 0;
htree.next = (ptr_htnode_t) 0;
htree.count = -1;
htree.item = NULLITEM; 
globalt = &htree;


// build hash tree
i = 0;
while (isinput[i].nitems!=-1){
  insertht(&htree,&(isinput[i]),0);
  iaux.nitems = 0;
#ifdef DEBUGINS
  printf("After inserting %d\n",i);
  printht(&htree,0,&iaux);
#endif
  i++;
}

printf("After inserting \n");
printht(&htree,0,&iaux);

printf("Number of transactions: %d\n", htree.children->count);
numtrans=htree.children->count;

printf("Traversal\n");

traverseht(&htree,0,&iaux,support);

// generate the subsets and calculate the metrics
// we should traverse the subtree and for each frequent itemset
// that has more than 2 itemsets, we should calculate the statistics

// traversal is similar to the print.

return 0;

}
