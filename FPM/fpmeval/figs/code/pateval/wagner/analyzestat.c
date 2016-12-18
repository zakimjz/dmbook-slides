#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "analyzestat.h"

static int cmpfloat(const void *p1, const void *p2) {
  float f1=(float)(*(float * const )p1), f2=(float)(*(float * const) p2);
  if (f1<f2) return -1;
  else if (f1>f2) return 1;
  else return 0;
}

int analyzevector(int size, float * vec, float *avg, float *var, hist_t * hist){
  int i;
  float aux=0.0;
  float binwidth;

  qsort(vec, size, sizeof(float), cmpfloat);
  hist->bins = 5;
  hist->inf = vec[0];
  hist->sup = vec[size-1];
  binwidth=(hist->sup-hist->inf)/hist->bins;
  (*avg)=0;
  for (i=0;i<MAXBINS;i++)
    hist->count[i]=0;
  // first calculate the average
  for (i=0;i<size;i++){
    (*avg) += vec[i];
    hist->count[(int)((vec[i]-hist->inf)/binwidth)]++; 
  }
  (*avg) /= size;
  // then calculate the variance
  for (i=0; i<size; i++)
    aux += ((*avg)-vec[i])*((*avg)-vec[i]);
  (*var) = aux/(size-1);
  return size;
}



