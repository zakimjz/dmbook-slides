/*  routines used FIM algorithms */
/* 2003/9/1 Takeaki Uno */
/* This program is available for only academic use.
   Neither commercial use, modification, nor re-distribution is allowed */

#ifndef _freqset_c_
#define _freqset_c_

#include<stdio.h>

#define   malloc2(f,a,b,c,d)     if(!(f=(a *)malloc(sizeof(a)*(b)))){printf("%s:memory error %s",c,d);exit(1);}

#define   fopen2r(f,a,c)     if(!(f=fopen(a,"r+"))){printf("%s:file open error %s\n",c,a);exit(1);}
#define   fopen2w(f,a,c)     if(!(f=fopen(a,"w+"))){printf("%s:file open error %s\n",c,a);exit(1);}

/* structure for bit-array */
typedef struct {
  unsigned long *a;
  int end;
} BARRAY;

/* macros for bit-array operations */
#define BARRAY_SET(A,x) ((A).a[(x)/32]|=BITMASK_1[(x)%32])
#define BARRAY_RESET(A,x) ((A).a[(x)/32]&=BITMASK_31[(x)%32])
#define BARRAY_BIT(A,x) ((A).a[(x)/32]&BITMASK_1[(x)%32])
#define BARRAY_01(A,x) (((A).a[(x)/32]&BITMASK_1[(x)%32])/BITMASK_1[(x)%32])

/* structure for integer queue */
typedef struct {
  int s;
  int t;
  int end;
  int *q;
} QUEUE;

/* structure for graph represented adjacency list, list is represented by queue */
typedef struct {
  int node_end;    
  int edge_end;    
  int arc_end;     
  int node_num;    
  int edge_num;    
  int arc_num;     
  int *buf;
  int *buf_w;
  QUEUE *edge;     
  QUEUE *in;       
  QUEUE *out;      
  int node1_num;   
  int *node_w; 
  QUEUE *edge_w; 
  QUEUE *in_w;  
  QUEUE *out_w;  
} SGRAPH;

/* macros for QUEUE operation */
#define QUEUE_LENGTH(Q) ((Q).t-(Q).s)
#define QUEUE_F_LOOP_(Q,i)  for((i)=(Q).s;(i)<(Q).t;(i)++)
#define QUEUE_FE_LOOP_(Q,i,x)  for((i)=(Q).s,x=(Q).q[i];(i)<(Q).t;(i)++,x=(Q).q[i])
#define QUEUE_B_LOOP_(Q,i)  for((i)=(Q).t-1;(i)>=(Q).s;(i)--)
#define QUEUE_BE_LOOP_(Q,i,x)  for((i)=(Q).t-1,x=(Q).q[i];(i)>=(Q).s;(i)--,x=(Q).q[i])
#define QUEUE_RMALL(Q) ((Q).t=(Q).s)


/* constants for bit mask */
int BITMASK_UPPER1[32] = { 0xffffffff, 0xfffffffe, 0xfffffffc, 0xfffffff8,
                           0xfffffff0, 0xffffffe0, 0xffffffc0, 0xffffff80,
                           0xffffff00, 0xfffffe00, 0xfffffc00, 0xfffff800,
                           0xfffff000, 0xffffe000, 0xffffc000, 0xffff8000,
                           0xffff0000, 0xfffe0000, 0xfffc0000, 0xfff80000,
                           0xfff00000, 0xffe00000, 0xffc00000, 0xff800000,
                           0xff000000, 0xfe000000, 0xfc000000, 0xf8000000,
                           0xf0000000, 0xe0000000, 0xc0000000, 0x80000000 };
int BITMASK_UPPER1_[32] = { 0xfffffffe, 0xfffffffc, 0xfffffff8, 0xfffffff0,
                            0xffffffe0, 0xffffffc0, 0xffffff80, 0xffffff00,
                            0xfffffe00, 0xfffffc00, 0xfffff800, 0xfffff000,
                            0xffffe000, 0xffffc000, 0xffff8000, 0xffff0000,
                            0xfffe0000, 0xfffc0000, 0xfff80000, 0xfff00000,
                            0xffe00000, 0xffc00000, 0xff800000, 0xff000000,
                            0xfe000000, 0xfc000000, 0xf8000000, 0xf0000000,
                            0xe0000000, 0xc0000000, 0x80000000, 0x00000000 };

int BITMASK_LOWER1[32] = { 0x00000000, 0x00000001, 0x00000003, 0x00000007,
                           0x0000000f, 0x0000001f, 0x0000003f, 0x0000007f,
                           0x000000ff, 0x000001ff, 0x000003ff, 0x000007ff,
                           0x00000fff, 0x00001fff, 0x00003fff, 0x00007fff,
                           0x0000ffff, 0x0001ffff, 0x0003ffff, 0x0007ffff,
                           0x000fffff, 0x001fffff, 0x003fffff, 0x007fffff,
                           0x00ffffff, 0x01ffffff, 0x03ffffff, 0x07ffffff,
                           0x0fffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff };
int BITMASK_LOWER1_[32] = { 0x00000001, 0x00000003, 0x00000007, 0x0000000f,
                            0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff,
                            0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff,
                            0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff,
                            0x0001ffff, 0x0003ffff, 0x0007ffff, 0x000fffff,
                            0x001fffff, 0x003fffff, 0x007fffff, 0x00ffffff,
                            0x01ffffff, 0x03ffffff, 0x07ffffff, 0x0fffffff,
                            0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff };

int BITMASK_1[32] = { 0x00000001, 0x00000002, 0x00000004, 0x00000008,
                      0x00000010, 0x00000020, 0x00000040, 0x00000080,
                      0x00000100, 0x00000200, 0x00000400, 0x00000800,
                      0x00001000, 0x00002000, 0x00004000, 0x00008000,
                      0x00010000, 0x00020000, 0x00040000, 0x00080000,
                      0x00100000, 0x00200000, 0x00400000, 0x00800000,
                      0x01000000, 0x02000000, 0x04000000, 0x08000000,
                      0x10000000, 0x20000000, 0x40000000, 0x80000000 };
int BITMASK_31[32] = { 0xfffffffe, 0xfffffffd, 0xfffffffb, 0xfffffff7,
                       0xffffffef, 0xffffffdf, 0xffffffbf, 0xffffff7f,
                       0xfffffeff, 0xfffffdff, 0xfffffbff, 0xfffff7ff,
                       0xffffefff, 0xffffdfff, 0xffffbfff, 0xffff7fff,
                       0xfffeffff, 0xfffdffff, 0xfffbffff, 0xfff7ffff,
                       0xffefffff, 0xffdfffff, 0xffbfffff, 0xff7fffff,
                       0xfeffffff, 0xfdffffff, 0xfbffffff, 0xf7ffffff,
                       0xefffffff, 0xdfffffff, 0xbfffffff, 0x7fffffff };

/* Global variables */
FILE *fp;
int  pflag=2, maxflag = 0;
int  v1bound = 1, v2bound=1000000000;
unsigned int cc=0, count = 0, c1=0, c2=0, c3=0;  /* counters */
int *el;
SGRAPH G;     /* graph */
int *perm=NULL;    /* permutation */
QUEUE jump, jump2;   /* jumplist */
BARRAY *BA=NULL;   /* adjacent matrix */
int vend, BAratio = 1000; /* deg/#node ratio */
int maxd1, maxd2, maxd2_;  /* maxdegree of V1, V2 */
QUEUE *jQ=NULL, incQ;   /* jQ:store X(cliq \cup {e})  */
QUEUE cliq, cliq_; /* cliq */
int *sc, scmax = -1; /* score */
QUEUE *ed, cliqtmp, stmp;

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* allocate integer array and set value */
int *intarray_malloc_const ( int end, int c ){
  int i, *a;
  malloc2 ( a, int, end, "intarray_malloc_const", "a");
  for ( i=0 ; i<end ; i++ ) a[i]=c;
  return (a);
}
/* quick sort of integers and its relationship function */
int qsort_int_cmp ( const void *x, const void *y ){
  if ( *((int *)x) < *((int *)y) ) return (-1);
  else return ( *((int *)x) > *((int *)y) );
}
int qsort_int_cmp_ ( const void *x, const void *y ){
  if ( *((int *)x) > *((int *)y) ) return (-1);
  else return ( *((int *)x) < *((int *)y) );
}
void qsort_int ( int *a, int siz ){
    qsort ( a, siz, sizeof(int), qsort_int_cmp );
}
void qsort_int_ ( int *a, int siz ){
    qsort ( a, siz, sizeof(int), qsort_int_cmp_ );
}
/* give an identical permutation */
int *init_perm ( int end ){
  int i, *p;
  malloc2 ( p, int, end, "init_perm", "p" );
  for ( i=0 ; i<end ; i++ ) p[i] = i;
  return ( p );
}
/* give an inverse of a permutation */
int *inverse_perm ( int *perm, int end ){
  int i, *p;
  malloc2 ( p, int, end, "inverse_perm", "p" );
  for ( i=0 ; i<end ; i++ ) p[i] = -1;
  for ( i=0 ; i<end ; i++ )
      if ( perm[i]>=0 && perm[i]<end ) p[perm[i]] = i; 
  return ( p );
}
/* permute array of struct according to permutation */
void perm_struct ( void *a, int unit, int *perm, int siz ){
  int i;
  char *s, *ss = a;
  malloc2 ( s, char, unit*siz, "perm_struct", "s" );
  memcpy ( s, ss, unit*siz );
  for ( i=0 ; i<siz ; i++ )
    memcpy ( ss + unit*perm[i], s + unit*i, unit );
  free ( s );
}
/* radix sort */
int *radix_sort ( void *a, int siz, int unit, int mm, int m, int *perm, int flag ){ 
  int *ll, *l, k, i, t;
  l = intarray_malloc_const ( m-mm, -1 );
  malloc2 ( ll, int, siz, "radix_sort", "ll");
  for ( i=0 ; i<siz ; i++ ){
    k = (*((int *)(((char *)a) + unit*i ))) - mm;
    ll[i] = l[k];
    l[k] = i;
  }
  if ( perm ){
    for ( k=0,i=0 ; k<m-mm ; k++ ){
      while ( l[k] >= 0 ){
        t = l[k];
        l[k] = ll[t];
        ll[t] = perm[i];
        i++;
      }
    }
    memcpy ( perm, ll, sizeof(int)*siz );
    free ( ll );
    free ( l );
    return ( perm );
  } else {
    for ( k=0,i=0 ; k<m-mm ; k++ ){
      while ( l[k] >= 0 ){
        t = l[k];
        l[k] = ll[t];
        ll[t] = i;
        i++;
      }
    }
    free ( l );
    return ( ll );
  }
}




/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* initialization of bit-array structure */
void BARRAY_init ( BARRAY *A, int siz ){
  A->end = siz;
  malloc2 ( A->a, unsigned long, (siz-1)/32+1, "BARRAY_init", "A->a" );
  A->a[(siz-1)/32] = 0;
}
/* free bit-array structure */
void BARRAY_end ( BARRAY *A ){
  if ( A->a  ){
    free ( A->a );
    A->a = NULL;
  }
}
/* set 0 to interval of bit-array */
void BARRAY_reset_interval ( BARRAY *A, int x1, int x2 ){
  int x=x1/32, i, xx=x2/32;
  if ( x == xx ){
    for ( i=x1%32 ; i<=x2%32 ; i++ )
      BARRAY_RESET (*A, i );
  } else {
    A->a[x]&=BITMASK_LOWER1[x1%32];
    x++;
    while ( x<xx ){
      A->a[x] = 0;
      x++;
    }
    A->a[xx]&=BITMASK_UPPER1_[x2%32];
  }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* print a QUEUE */
void QUEUE_print ( QUEUE *Q ){
  int i;
  for ( i=Q->s ; i!=Q->t ; i++ ){
    printf ("%d,",Q->q[i]);
  }
  printf ("\n");
}
/* initialize QUEUE structure */
void QUEUE_init ( QUEUE *Q, int siz ){
  Q->s = 0;
  Q->t = 0;
  Q->end = siz+1;
  malloc2 ( Q->q, int, siz+1, "QUEUE_init", "Q->q");
}
/* free QUEUE */
void QUEUE_end ( QUEUE *Q ){
  if ( Q->q ){
    free ( Q->q );
    Q->q = 0;
  }
}
/* insert an element to the end of QUEUE */
void QUEUE_ins_ ( QUEUE *Q, int e ){
  Q->q[Q->t] = e;
  Q->t++;
}
/* remove jth element and move the last element to jth position */
void QUEUE_rm_ ( QUEUE *Q, int j ){
  Q->t--;
  Q->q[j] = Q->q[Q->t];
}
/* concatinate Q2 to Q1 */
void QUEUE_concat_ ( QUEUE *Q1, QUEUE *Q2 ){
  memcpy ( &(Q1->q[Q1->t]), &(Q2->q[Q2->s]), QUEUE_LENGTH(*Q2)*sizeof(int));
  Q1->t += Q2->t-Q2->s;
}
/* copy interval starting from s2 of Q2 to Q1 of position s1 with length l */
void QUEUE_subcpy_ ( QUEUE *Q1, int s1, QUEUE *Q2, int s2, int l ){
  memcpy ( &(Q1->q[s1]), &(Q2->q[s2]), (l-s2)*sizeof(int));
}
/* copy Q2 to Q1 */
void QUEUE_cpy_ ( QUEUE *Q1, QUEUE *Q2 ){
  QUEUE_RMALL (*Q1);
  QUEUE_concat_ ( Q1, Q2 );
}
/* set Q1 to Q1 and Q2. both Q1 and Q2 have to be sorted */
/* resulted Q1 is sorted */
void QUEUE_and_ ( QUEUE *Q1, QUEUE *Q2 ){
  int i=Q1->s, i2 = Q2->s, ii=Q1->s;
  while ( i != Q1->t && i2 != Q2->t){
    if ( Q1->q[i] > Q2->q[i2] ) i2++;
    else {
      if ( Q1->q[i] == Q2->q[i2] ){
        Q1->q[ii] = Q1->q[i];
        ii++;
      }
      i++;
    }
  }
  Q1->t = ii;
}
/* set Q1 to Q1 union Q2. both Q1 and Q2 have to be sorted */
/* resulted Q1 is sorted */
void QUEUE_merge_ ( QUEUE *Q1, QUEUE *Q2 ){
  int i=Q1->t-1, j=Q2->t-1, t=i+j-Q2->s+1;
  int ei, ej;
  if ( i+1 == Q1->s || j+1 == Q2->s ){
    QUEUE_concat_ ( Q1, Q2 );
    return;
  }
  Q1->t = t+1;
  ei = Q1->q[i];
  ej = Q2->q[j];
  while (1){
    if ( ei > ej ){
      Q1->q[t] = ei;
      if ( i == Q1->s ){
        QUEUE_subcpy_ ( Q1, Q1->s, Q2, Q2->s, j+1 );
        return;
      }
      i--;
      ei = Q1->q[i];
    } else {
      Q1->q[t] = ej;
      if ( j == Q2->s ) return;
      j--;
      ej = Q2->q[j];
    }
    t--;
  }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* max degree of graph */
/* flag &3 == 0:all 1:V1 2:V2, flag&12 == 0:normal 4:directed, 8:in, 12:out */
int SGRAPH_max_degree ( SGRAPH *G, int flag ){
  int i, m=0, k;
  for ( i=((flag&3)==2?G->node1_num:0) ; i<((flag&3)==1?G->node1_num:G->node_end) ; i++ ){
    k = QUEUE_LENGTH (G->edge[i]);
    if ( k > m ) m = k;
  }
  return ( m );
}

/* replace node i by perm[i] */
void SGRAPH_replace_index ( SGRAPH *G, int *perm ){
  int i, v, j;
  for ( i=0 ; i<G->node_end ; i++ ){
    if ( G->edge[i].q ){
      QUEUE_FE_LOOP_ ( G->edge[i], j, v ){
        G->edge[i].q[j] = perm[v];
      }
    }
  }
  if ( G->edge ) perm_struct ( G->edge, sizeof(QUEUE), perm, G->node_end );
  if ( G->in ){
    perm_struct ( G->in, sizeof(QUEUE), perm, G->node_end );
    perm_struct ( G->out, sizeof(QUEUE), perm, G->node_end );
  }
  if ( G->node_w ){
    perm_struct ( G->node_w, sizeof(int), perm, G->node_end );
  }
}
/* sort edges incident to each node  */
void SGRAPH_sort_incident_edges ( SGRAPH *G, int flag ){
  int i, j, v, *x, z;
  malloc2 ( x, int, G->node_end, "SGRAPH_sort_incident_edge_", "x");
  for ( i=0 ; i<G->node_end ; i++ ){
    x[i] = G->edge[i].t;
    QUEUE_FE_LOOP_ ( G->edge[i], j, v )
        if ( v <= i ) QUEUE_ins_ ( &(G->edge[v]), i );
    QUEUE_RMALL ( G->edge[i]);
  }
  for ( i=0 ; i<G->node_end ; i++ ){
    z = x[i]-QUEUE_LENGTH (G->edge[i]);
    QUEUE_BE_LOOP_ ( G->edge[i], j, v )
        G->edge[i].q[j+z] = v;
    G->edge[i].t = G->edge[i].s;
  }
  for ( i=0 ; i<G->node_end ; i++ ){
    for ( j=G->edge[i].t ; j<x[i] ; j++ ){
      v = G->edge[i].q[j];
      if ( v > i ) QUEUE_ins_ ( &(G->edge[v]), i );
    }
    G->edge[i].t = x[i];
  }
  free (x);
}

/* sort nodes by degrees */
int *SGRAPH_sort_node_degree ( SGRAPH *G, int flag ){
  int *perm, *q;
  perm = init_perm ( G->node_end );
  radix_sort ( &(G->edge[G->node1_num].t), G->node_end-G->node1_num,
              sizeof(QUEUE), 0, G->node_end, &(perm[G->node1_num]), flag);
  if ( G->node1_num > 0 )
      radix_sort ( &(G->edge[0].t), G->node1_num,
                 sizeof(QUEUE), 0, G->node_end, perm, flag);
  SGRAPH_replace_index ( G, perm );
  q = inverse_perm (perm, G->node_end );
  free ( perm );
  return ( q );
}

/* convert to simple graph ( after sort indiced edges ) */
void SGRAPH_simple ( SGRAPH *G ){
  int i, j, u, jj;
  for ( i=0 ; i<G->node_end ; i++ ){
    if ( QUEUE_LENGTH ( G->edge[i] ) == 0 ) continue;
    jj = G->edge[i].s;
    for ( j=jj=G->edge[i].s+1 ; j<G->edge[i].t ; j++ ){
      if ( G->edge[i].q[j] != G->edge[i].q[j-1] && G->edge[i].q[j] != i ){
        G->edge[i].q[jj] = G->edge[i].q[j];
        jj++;
      }
    }
    G->edge[i].t = jj;
  }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* FAST file output */
#define FPRINTP_BUFSIZ 100000
char sss[FPRINTP_BUFSIZ];
void QUEUE_perm_fprintp ( QUEUE *S, FILE *fp, int *perm ){
  int i, e, j=FPRINTP_BUFSIZ-1;
  QUEUE_BE_LOOP_ ( *S, i, e ){
    e = perm[e];
    sss[j--] = ' ';
    while ( e>0 ){
      sss[j--] = e%10+'0';
      e /= 10;
    }
  }
  fwrite ( &(sss[j+1]), FPRINTP_BUFSIZ-j-1, 1, fp );
}

#define intarray2_BUFSIZ 1048576
int  **buff=NULL, bc=-1, bcc=intarray2_BUFSIZ, prebcc, pbcc, ppbcc;  /* original buff, its (sub) counter, stored subcounter */ 


/* init 2D-buffer */
void intarray2_init_buf (){
  prebcc = bc = -1;
  bcc = intarray2_BUFSIZ;
  malloc2 ( buff, int *, 65536, "intarray2", "buff" );  /* for allocating buffers */
  pbcc = prebcc;
  ppbcc = prebcc = bcc;
}
/* free 2D-buffer */
void intarray2_free_buf (){
  for ( ; bc>=0 ; bc-- ) free ( buff[bc] );
  free ( buff );
}


/* write an integer to 2D-buffer */
void intarray2_write_buf ( int a ){
  if ( bcc >= intarray2_BUFSIZ ){
    bc++;
    malloc2 ( buff[bc], int, intarray2_BUFSIZ, "intarray2_write_buf", "buff[bc]" );  /* malloc large memory */
    if ( prebcc < bcc ){
      memcpy ( buff[bc], &(buff[bc-1][prebcc]), sizeof(int)*(bcc-prebcc) );
    } 
    bcc = intarray2_BUFSIZ - prebcc;
    prebcc = 0;
  }
  buff[bc][bcc++] = a;
}

/* push current positions of block */
void intarray2_push (){
  intarray2_write_buf ( pbcc );  /* store the start/end position in buff */
  intarray2_write_buf ( ppbcc );
  pbcc = prebcc;
  ppbcc = prebcc = bcc;
}

/* pop previous positions of block */
int intarray2_pop (){
  prebcc = buff[bc][bcc-2];
  bcc = buff[bc][bcc-1];
  if ( pbcc == 0 ) bc--;
  pbcc = prebcc;
  return ( pbcc );
}



/* read transaction datas from file */
/* allocate large memory, and devide it into small peices each of which corresponds to the queue of items of a transaction */
/* *q1&1:alloc #line+#ele &2:add #line to each element, &4:sort each line */
/* &8:transpose  */
QUEUE *QUEUE_load2 ( char *fname, int *q1, int *q2, int flag ){
  int i=0, item, bctmp;
  unsigned char c;
  QUEUE *D;
  FILE *fp;

  intarray2_init_buf ();
  *q1 = -1;
  fopen2r ( fp, fname, "intarray_load");
  
  while (1){
    do {
      for ( bctmp=0,item=0 ; 1 ; bctmp++,item=item*10 +(int)(c-'0') ){
        if ( feof(fp) ) goto READ_END;
        c = getc(fp);
        if ( (c < '0') || (c > '9')) break;
      }
      if ( bctmp>0 ){
        if ( *q1 < item ) *q1 = item;  /* *q1 := max number */
        intarray2_write_buf ( item );
      }
    } while (c != '\n');
    if ( bcc != prebcc ){
      intarray2_push ();
      i++;
    }
  }

READ_END:;
  (*q1)++;  /* #items */
  *q2 = i;
  malloc2 ( D, QUEUE, i+((flag&1)?*q1:0), "QUEUE_load2", "D" );
             /* queues storing items of transactions */
  bctmp=bc;
  prebcc = pbcc;
  if ( flag&8 ) i += *q1;
  do {
    i--;
    D[i].q = &(buff[bc][pbcc]);
    D[i].s = 0;
    D[i].t = bcc-pbcc-2;
    D[i].end = bcc-pbcc-1;
    if ( flag&2 ) QUEUE_F_LOOP_ ( D[i], item ) D[i].q[item] += *q2;
    if ( flag&4 ) qsort_int ( D[i].q, D[i].t );
    intarray2_pop ();
  } while ( i != 0 );
  fclose ( fp );
  bc = bctmp;

  return ( D );
}

/* set edges of V2 vertices */
/* flag&1: count only, &2:mk v1edge */
void SGRAPH_mk_v2edge ( SGRAPH *G, int flag ){
  int i, j, e, i1=G->node1_num, i2=G->node_end, j1=0, j2=G->node1_num;
  if ( flag&2) { e=i1, i1=j1; j1=e; e=i2, i2=j2; j2=e; } 
  for ( i=i1 ; i<i2 ; i++ ) G->edge[i].t =0;
  for ( i=j1 ; i<j2 ; i++ )
      QUEUE_FE_LOOP_ ( G->edge[i], j, e ) G->edge[e].t++;
  if ( flag ) return;
  for ( i=i1 ; i<i2 ; i++ )
      QUEUE_init ( &(G->edge[i]), G->edge[i].t );
  for ( i=j1 ; i<j2 ; i++ )
      QUEUE_FE_LOOP_ ( G->edge[i], j, e ) QUEUE_ins_ ( &(G->edge[e]), i );
}


/* output current set */
void FIM_output ( QUEUE *S, int siz ){
  int i, e, s = QUEUE_LENGTH(*S);
  count++;
  if ( scmax < s ) scmax = s;
  sc[s]++;
  if ( !fp ) return;
  QUEUE_FE_LOOP_ ( *S, i, e ) fprintf (fp, "%d ", e-G.node1_num );
  fprintf (fp, "(%d)\n", siz );
}
/* output permutated current set */
void FIM_perm_output ( QUEUE *S, int siz ){
  int i, e, s = QUEUE_LENGTH(*S), j, f;
  count++;
  if ( scmax < s ) scmax = s;
  sc[s]++;
  if ( !fp ) return;
#if FREQSET_OUTPUT==1
  QUEUE_FE_LOOP_ ( *S, i, e ) fprintf (fp, "%d ", perm[e] );
  fprintf (fp, "(%d)\n", siz );
#endif
#if FREQSET_OUTPUT==2
#ifdef BMACE_DEBUG
  QUEUE_cpy_ ( &cliqtmp, S );
  QUEUE_FE_LOOP_ ( cliqtmp, i, e ){
    for ( j=i-1; j>=0 ; j-- ){
      f = cliqtmp.q[j];
      if ( perm[e] < perm[f] ){
        cliqtmp.q[j+1] = f;
        cliqtmp.q[j] = e;
      }
    }
  }
  
  QUEUE_perm_fprintp ( &cliqtmp, fp, perm );
#else
  QUEUE_perm_fprintp ( S, fp, perm );
#endif
  fprintf (fp, "(%d)\n", siz );
#endif
}

/* output #outputs of each size */
void FIM_output2 (){
  int i;
  //  if ( !fp ) return;
#if FREQSET_PROBLEM==1
  if(fp) fprintf (fp,"(%d)\n", G.node1_num );
  printf ("1\n");
#endif
#if FREQSET_PROBLEM==2
  if ( maxd2 == G.node1_num ) printf ("0\n", G.node1_num );
  else { 
    if(fp) fprintf (fp,"(%d)\n", G.node1_num );
    printf ("1\n");
  }
#endif
#if FREQSET_PROBLEM==3
  if ( maxd2 >= v1bound ) printf ("0\n", G.node1_num );
  else { 
    if(fp) fprintf (fp, "(%d)\n", G.node1_num );
    printf ("1\n");
  }
#endif
  for ( i=1 ; i<=scmax ; i++ ){
    printf ("%d\n", sc[i] );
  }
}


/* Compute X(K) of K */
int FREQSET_compute_XK ( QUEUE *K, QUEUE *S ){
  int i;
  if ( QUEUE_LENGTH(*K) <= 0 ){
    for ( i=0 ; i<G.node1_num ; i++ )
      QUEUE_ins_ ( S, i );
  } else {
    QUEUE_cpy_ (S, &(G.edge[K->q[0]]));
    for ( i=1 ; i<K->t ; i++ )
        QUEUE_and_ ( S, &(G.edge[K->q[i]]));
  }
  return ( QUEUE_LENGTH(*S) );
}

/* Set K to greatest maximal frequent set including K by using S=X(K) */
void FREQSET_mk_maxfreq ( QUEUE *K, QUEUE *S ){
  int u, v, vv, uu, j, i, c;

  QUEUE_RMALL ( *K );
  for ( v=G.node1_num ; v<G.node_end ; v++ ){
    el[v] = -1; G.edge[v].s = G.edge[v].end = 0;
  }
  QUEUE_FE_LOOP_ ( *S, j, u ){
    vv = G.edge[u].q[0];
    el[u] = el[vv];
    el[vv] = u;
  }

  for ( j=0, v=G.node1_num ; v<G.node_end ; v++ ){
    for ( c=0,u=el[v] ; (u!=-1)&&(G.edge[u].s==j) ; u=el[u] ) c++;
    if ( c >= v1bound ){
      for ( u=el[v] ; (u!=-1)&&(G.edge[u].s==j) ; c++,u=uu ){
        uu = el[u];
        G.edge[u].s = j+1;
        G.edge[u].end++;
        vv = G.edge[u].q[G.edge[u].end];
        el[u] = el[vv];
        el[vv] = u;
      }
      j++;
      QUEUE_ins_ ( K, v );
    } else {
      for ( u=el[v] ; (u!=-1)&&(G.edge[u].s==j) ; c++,u=uu ){
        uu = el[u];
        G.edge[u].end++;
        vv = G.edge[u].q[G.edge[u].end];
        el[u] = el[vv];
        el[vv] = u;
      }
    }
  }

  for ( v=G.node1_num ; v<G.node_end ; v++ ){
    G.edge[v].s = 0;
    G.edge[v].end = G.edge[v].t+1;
  }
  QUEUE_FE_LOOP_ ( *S, j, u ){
    G.edge[u].s = 0;
    G.edge[u].end = 0;
  }
}



/* initialize */
/* flag&1: read both direction */
/* flag&2: minus G.node1_num from perm */
/* flag&4: allocate memory for jQ */
/* flag&8: remove unnecessary items from transactions */
void FREQSET_init ( int flag ){
  int i, j, u, w;

  sc = intarray_malloc_const ( G.node_end-G.node1_num, 0 );
  malloc2 ( jQ, QUEUE, G.node_end, "FREQSET_init", "jQ");

  SGRAPH_mk_v2edge ( &G, (flag&1)?0:1 );

#ifdef BMACE_DEBUG
  SGRAPH_mk_v2edge ( &G, 0 );
  QUEUE_init ( &stmp, G.node_end );
  QUEUE_init ( &cliqtmp, G.node_end-G.node1_num );
  QUEUE_ins_ ( &cliqtmp, 3+G.node1_num );
  QUEUE_ins_ ( &cliqtmp, 6+G.node1_num );
  QUEUE_ins_ ( &cliqtmp, 11+G.node1_num );
  QUEUE_ins_ ( &cliqtmp, 27+G.node1_num );
  QUEUE_ins_ ( &cliqtmp, 86+G.node1_num );
  FREQSET_compute_XK ( &cliqtmp, &stmp );
  printf ("stmp = %d\n", QUEUE_LENGTH ( stmp ));
#endif

  perm = SGRAPH_sort_node_degree ( &G, 0 ); 
  maxd1 = G.edge[G.node1_num-1].t;
  maxd2 = G.edge[G.node_end-1].t;
  maxd2_ = G.edge[G.node_end-2].t;
  if ( flag&2 ) for (i=G.node1_num ; i<G.node_end ; i++) perm[i]-=G.node1_num;
  if ( flag&8 ){
    for ( i=0 ; i<G.node1_num ; i++ ){
      QUEUE_FE_LOOP_ ( G.edge[i], j, u ){
        if ( G.edge[u].t < v1bound ) G.edge[i].s++;
      } 
    }
  }
  if ( flag&1 ){
    SGRAPH_sort_incident_edges ( &G, 0 );
  } else {
    for ( i=0 ; i<G.node1_num ; i++ ){
      qsort_int ( G.edge[i].q, G.edge[i].t );
      G.edge[i].q[G.edge[i].t] = G.node_end;
    }
    for ( i=G.node1_num ; i<G.node_end ; i++ ){
      if ( flag&4 ) QUEUE_init ( &jQ[i], G.edge[i].t );
      G.edge[i].t = 0;
      G.edge[i].q = NULL;
    }
  }
  SGRAPH_simple ( &G );
  for ( i=0 ; i<G.node_end ; i++ ){
    G.edge[i].end = 0;
    jQ[i].end = 0;
    jQ[i].s = 0;
    jQ[i].t = 0;
    jQ[i].q = NULL;
  }
  for ( i=0 ; i<G.node1_num ; i++ ){
    G.edge[i].q[G.edge[i].t] = G.node_end; /* loop stopper */
  }
  if ( flag&1 ) el = NULL;
  else { 
    el = intarray_malloc_const  ( G.node_end, -1 );
    for ( i=0 ; i<G.node1_num ; i++ ){
      j = G.edge[i].q[0];
      el[i] = el[j];
      el[j] = i;
    }
  }
  QUEUE_init ( &incQ, maxd2 );
  QUEUE_init ( &cliq, maxd1+2 );
  QUEUE_init ( &cliq_, maxd1+2 );
  QUEUE_init ( &jump, (G.node_end-G.node1_num+2)*2 );
  QUEUE_init ( &jump2, (G.node_end-G.node1_num+2)*2 );
  ed = G.edge;
/*  QUEUE_init ( &cliqtmp, maxd1+2 );
  QUEUE_init ( &stmp, maxd2+2 ); */
}

/* ending operation */
void FREQSET_end ( int flag ){
  int i;
  QUEUE_end ( &incQ );
  QUEUE_end ( &cliq );
  QUEUE_end ( &cliq_ );
  QUEUE_end ( &jump );
  QUEUE_end ( &jump2 );
  for ( i=G.node1_num ; i<G.node_end ; i++ ){
    if ( flag&1 ) QUEUE_end ( &G.edge[i] );
    if ( flag&4 ) QUEUE_end ( &jQ[i] );
  }
  free ( G.edge );
  free ( jQ );
  if ( perm ) free ( perm );
  if ( sc ) free ( sc );
/*  QUEUE_end ( &cliqtmp );
  QUEUE_end ( &stmp ); */
}

#endif
