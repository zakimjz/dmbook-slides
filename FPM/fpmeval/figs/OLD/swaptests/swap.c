
#include <stdlib.h>
#include "mex.h"

int randint(int N) {
  /* Return random integer between 0 and N-1 */
  return (int)(((double)rand()/(double)(RAND_MAX + 1.0))*(double)N);
}

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
 int m, n, i, j, swap;
 double *IN, *A, *Counts, *realswaps, *IndX, *IndY;
 
 /* Input data matrix */
 IN = mxGetPr(prhs[0]);

 /* Input data matrix dimensions */
 m = mxGetM(prhs[0]);
 n = mxGetN(prhs[0]);

 /* Create an mxArray for the output */
 plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
 /* Create an array for the output's data */
 A = mxGetPr(plhs[0]);

 /* return as second output parameter the number of real swaps */
 plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
 realswaps = mxGetPr(plhs[1]);

 /* Third output parameter: 
    Create a second mxArray output array to record the 
    number of swaps in each position */
 plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
 Counts = mxGetPr(plhs[2]);

 /* Put data in the array */
 for (j = 0; j < m*n; j++) {
   A[j] = IN[j];
   Counts[j] = 0; 
 }
 
 /* Count number of edges and check if matrix is 0-1 */
 int E = 0;
 for (i=0; i<m; i++)
   for (j=0; j<n; j++)
     if (A[j*m + i] == 1)
       E++; 
     else if (A[j*m + i] != 0)
       mexErrMsgTxt("Input matrix should be 0-1"); 
 
 /* IndX and IndY keep a sparse representation of matrix A, 
    therefore they should be updated simultaneously with A 
    they are needed in order to sample edges fast */
  
 IndX = (double *) mxMalloc(m*n * sizeof(double));
 IndY = (double *) mxMalloc(m*n * sizeof(double));

 int k = 0; 
 for (i=0; i<m; i++)
   for (j=0; j<n; j++)
     if (A[j*m + i] == 1) {
       IndX[k] = i; 
       IndY[k] = j; 
       k++; 
     }

 int S; 
 if (nrhs==1) 
   S = E; 
 else 
   S = (int)(mxGetScalar(prhs[1]));
 /* printf("Attempted swaps: %d\n", S); */

 (*realswaps) = 0.0; 
 for (swap = 0; swap<S; swap++) {
   int e1 = randint(E);
   int e2 = randint(E);
   int x1 = IndX[e1]; 
   int y1 = IndY[e1]; 
   int x2 = IndX[e2]; 
   int y2 = IndY[e2]; 
   
   if ((A[y2*m + x1]==0) && (A[y1*m + x2]==0) && (x1!=x2) && (y1!=y2)) {
     (*realswaps)++;
     
     /* swap in the full matrix */
     A[y2*m + x1] = 1; 
     A[y1*m + x2] = 1; 
     A[y1*m + x1] = 0; 
     A[y2*m + x2] = 0; 

     /* swap in the sparse matrix */ 
     IndY[e1] = y2; 
     IndY[e2] = y1; 

     /* update the counts */
     Counts[y2*m + x1]++; 
     Counts[y1*m + x2]++; 
     Counts[y1*m + x1]++; 
     Counts[y2*m + x2]++; 
   }
 }
 /* printf("Swaps happened: %d\n", (int)(*realswaps)); */

 /* Free memory */

 mxFree(IndX); 
 mxFree(IndY); 
}
