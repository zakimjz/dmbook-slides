#ifndef ANALYZESTATH 

#define  ANALYZESTATH

// Constants

#define MAXBINS 100

typedef struct hist{
   float inf, sup;
   int bins;
   int count[MAXBINS];
} hist_t, *ptr_hist_t;

#endif
