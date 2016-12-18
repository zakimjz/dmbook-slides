This directory contains methods for assessing data mining results
using swap randomization as described in

Aristides Gionis, Heikki Mannila, Taneli Mielikäinen, and Panayiotis
Tsaparas. Assessing data mining results via swap randomization. In
Mark Craven and Dimitrios Gunopulos (Eds.): The Twelfth Annual SIGKDD
International Conference on Knowledge Discovery and Data Mining (KDD
2006). ACM, 2006.


Swap randomization for transaction databases and frequent itemsets
------------------------------------------------------------------

Perl scripts swaptdb.pl and swapfreq.pl produce swap-randomized
transaction databases and frequent itemsets.  freqchanges.pl compares
the frequencies of frequent itemsets in the original and
swap-randomized transaction databases.

Usage: perl swaptdb.pl tdbfile tdbprefix stepsize steps

 parameters
  tdbfile   : the name of the file consisting the transaction database
  tdbprefix : the prefix of the filenames of the swap-randomized
              transaction databases produced by the program
  stepsize  : the number of attempted swaps between the stored 
              transaction databases
  steps     : the number of swap-randomized transaction databases
              produced

The format of the transaction database is the same as used in FIMI
repository (http://fimi.cs.helsinki.fi), i.e., the items in the
transaction databases are positive integers and each row in the file
is one transaction represented as a list of items in ascending order.


Usage: perl swapfreq.pl tdbfile freqprefix tmpfile minsupp stepsize steps

 parameters
  tdbfile    : the name of the file consisting the transaction database
  freqprefix : the prefix of the filenames of the swap-randomized
               transaction databases produced by the program
  tmpfile    : the name of the temporary file used to store the
               swap-randomized transaction databases
  minsupp    : the minimum support threshold for frequent itemset mining
  stepsize   : the number of attempted swaps between the stored 
               transaction databases
  steps      : the number of swap-randomized transaction databases
               produced

The format of the transaction database is the same as for swaptdb.pl.
The script need a program 'fim_all' for mining frequent itemsets, as
given in the FIMI repository.


Usage: perl freqchanges.pl tdbfile tmpfile1 tmpfile2 minsupp minsupporig \
                           minsuppswap stepsize

 parameters
  tdbfile     : the name of the file consisting the transaction database
  freqprefix  : the prefix of the filenames of the swap-randomized
               transaction databases produced by the program
  tmpfile1    : the name of the temporary file used to store the frequent
                itemsets in the swap-randomized transaction databases
  tmpfile2    : the name of the temporary file used to store the
                swap-randomized transaction databases
  minsupp     : the minimum support threshold for items to be considered
  minsupporig : the minimum support threshold for frequent itemset mining
                in the original transaction database
  minsuppswap : the minimum support threshold for frequent itemset mining
                in the swap-randomized transaction database
  stepsize    : the number of attempted swaps for producing the 
                swap-randomized transaction database

The format of the transaction database is the same as with swaptdb.pl
and swapfreq.pl, and the some implementation of 'fim_all' e.g. from
the FIMI repository is needed.

The columns of the output are the following:
 freqorig freeqswap relerrorig relerrswap avgcorrorig avgcorrswap \
 freqswap/freqorig freqorig/freqswap ratioorig ratioswap ; itemset

 freqorig    : the support of the itemset in the original data
 freqswap    : the support of the itemset in the swap-randomized data
 relerrorig  : (freqorig-freqswap)/freqorig
 relerrswap  : (freqorig-freqswap)/freqorig
 avgcorrorig : the average correlation between the items of the itemset
               in the original transaction database
 avgcorrswap : the average correlation between the items of the itemset
               in the swap-randomized transaction database
 swaplift    : freqswap/freqorig
 liftswapped : freqorig/freqswap
 ratioorig   : the lift in the original transaction database
 ratioswap   : the lift in the swap-randomized transaction database
 itemset     : the list of the items in the itemset in ascending order


Mex file for swaping a 0-1 matrix
---------------------------------

C code that compiles and runs through Matlab.

File: swap.c

In order to compile follow the following steps inside Matlab:

1. set the compiler option:

>> mex -setup

(in our setup it is option 2 for gcc)

2. Compile the file:

>> mex swap.c

Usage of function:

  Y = swap(X) is the default call
  X is the input binary matrix.
  Y is the output swapped matrix.
  The number of swaps attempted is the number of ones of matrix X.

  Y = swap(X,k) attempts to perform k swaps.

  [Y,t] = swap(X) records t, the actual number of swaps that took
  place.



Matlab function to test significance of clustering
--------------------------------------------------

Compute clustering error on the original data and the swapped data

File: clusteringtest.m

Usage:

  [orig,permuted,sw] = clusteringtest(Dinp,k,nosamples,L);

  Input
    Dinp: input 0-1 matrix
    k: number of clusters to be tested
    nosamples: number of samples (swapped datasets) to be drawn
    L: number of swaps to be attempted in order to generate each 
    sample

  Output
    orig: error on the original data
    permuted: vector of errors on permuted data
    sw: vector of actual swaps that took place in each sample


For example, generate artificial data that contain clustered structure
(use files gendata.m and randgd.m provided in this directory)

>> A = gendata(100,20,5,0.1);  % 100 points, 20 clusters, 5 clusters, 0.1 noise

and run

>> [or,perm,sw] = clusteringtest(A,5,100);




