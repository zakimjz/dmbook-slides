/*  BMACEa: Bipartite MAximal Clique Enumerator +alpha */
/* 2003/9/1 Takeaki Uno */
/* This program is available for only academic use.
   Neither commercial use, modification, nor re-distribution is allowed */

/*  BMACE enumerates all 1. frequent sets, 2. closed item sets, 3. maximal frequent sets of input transactions */
/* usege: BMACE filename [fcm] threshold [flag]
   f: frequent set, c: closed item set, m: maximal frequent set
   threshold: output (frequent, closed item, maximal frequent) sets included
    in at least "threshold" transactions.
   q: no output but count the number of output
*/

#define errr(a,b) if(a){printf("okasii %s\n",b);exit(1);}
/* #define BMACE_DEBUG */

/* FREQSET_PROBLEM:specify the problem, 
    1:Frequent set, 2:Closed set 3:Maximal Frequent Set  */
/* #define FREQSET_PROBLEM 1 */

/* FREQSET_OUTPUT: specify output way
    1:usual, 2:fast version, 3:no output */
/* #define FREQSET_OUTPUT 2 */

int BMACE_mode = FREQSET_PROBLEM;

#include"freqset_.c"


  /* return 0 if no cliq\cup j, j<i is freqeunt, and 1 otherwise */
int BMACEmax_maximal ( QUEUE *S, int i ){
  int u, *v, j, flag = 0, js = jump.s, tt;
  jump.s = jump.t;
  QUEUE_FE_LOOP_ ( *S, j, u ){
    for ( v=ed[u].q ; (*v)<i ; v++ )
      if ( jQ[*v].s++ == 0 ) QUEUE_ins_ ( &jump, *v );
  }
/*  QUEUE_BE_LOOP_ ( jump, j, u ){
    if ( jQ[u].s >= v1bound && G.edge[u].end == 0 ) { flag = 1; break; }
    jQ[u].s = 0;
  } 
  for ( ; j<jump.t ; j++ ) jQ[jump.q[j]].s = 0; */
  QUEUE_BE_LOOP_ ( jump, j, u ){
    if ( jQ[u].s >= v1bound && G.edge[u].end == 0 ) flag = 1;
    jQ[u].s = 0;
  } 

  jump.t = jump.s;
  jump.s = js; 
  return ( flag );
}


  /* return least item <i whose addition does not change X(cliq) */
  /* return -1 if not exist */
int BMACEclosed_maximal ( QUEUE *S, int i ){
  int j, u, jj=0, t, *q=S->q, m=-1, uu=-1, v, *vv, vvv;  
  S->q[S->t] = G.node_end;
  for ( j=0 ; (u=G.edge[q[0]].q[j])<i ; j++ ){
    if ( G.edge[u].end == 1 ) continue;
    for ( jj=0 ; (v=q[jj])<vend ; jj++ ){
      while ( (vvv=G.edge[v].q[jQ[v].s]) < u ) { jQ[v].s++; }
      if ( vvv != u ) goto END1;
    }
    for ( vv=&q[jj]; *vv<G.node_end ; vv++ ){
      if ( !BARRAY_BIT(BA[*vv], u-G.node1_num )) goto END1;
    }
    for ( j=0 ; (u=q[j])<vend ; j++ ) jQ[u].s = 0;
    return (u);    /* non-maximal if all are added */
    END1:;
    if ( m<jj ) m=jj;
  }
  END2:;
  for ( jj=0 ; jj<=m ; jj++ ) jQ[q[jj]].s = 0;
  return ( -1 );
}


/* make jumplist and jQ for S and indices >i */
int BMACE_mk_jump ( QUEUE *S, int i, QUEUE *jump, QUEUE *jump_, QUEUE *cl ){
  int u, *v, j, c=0, cc=0, ccc=0;
  QUEUE_FE_LOOP_ ( *jump, j, u ){
    if ( u <= i ) break;
    jQ[u].t = 0;
    jQ[u].end = 0;
  }
  QUEUE_FE_LOOP_ ( *S, j, u ){
    for ( v=&(ed[u].q[ed[u].t-1]) ; *v>i ; v-- ){
      if ( jQ[*v].end == 0 ) jQ[*v].q[jQ[*v].t++] = u;
    }
  }
  QUEUE_FE_LOOP_ ( *jump, j, u ){
    cc += jQ[u].t;
    if ( u <= i ) break;
    if ( jQ[u].t >= v1bound ){
      if ( jQ[u].t == QUEUE_LENGTH(*S) ){
        QUEUE_ins_ ( cl, u );
        jQ[u].end = 1;
      } else {
        QUEUE_ins_ ( jump_, u );
        c += QUEUE_LENGTH(*S) - jQ[u].t;
        ccc++;
      }
    } else jQ[u].end = 1;
  }
  c = (c+cc)*10 / (c+1);
  if ( c>= 100 ) c = 99;
  return ( c + ccc*100);
}

/* malloc jQ for S and indices >i */
void BMACE_jQ_malloc ( QUEUE *S, int i, QUEUE *cl ){
  int u, v, j, jj;
  jump.t = jump.s = 0;
  QUEUE_FE_LOOP_ ( *S, j, u ){
    for ( jj=G.edge[u].t-1 ; (jj>=0)&&(v=G.edge[u].q[jj])>i ; jj-- ){
      jQ[v].t++;
      jQ[v].end = 1;
      if ( jQ[v].t == v1bound ){
        QUEUE_ins_ ( &jump, v );
      }
    }
  }
  jj=0; QUEUE_FE_LOOP_ ( jump, j, u ){
    jump.q[jj++] = u;
    QUEUE_init ( &jQ[u], jQ[u].t+2 );
  }
  jump.t = jj;
  qsort_int_ ( jump.q, jump.t );
}

/* malloc jQ_ for S and indices >i */
QUEUE *BMACE_malloc_mk_jQ_ ( QUEUE *S, int i, QUEUE *jump, int uu ){
  int u, v, j, jj, t;
  QUEUE *IjQ_;
  jump2.s = jump2.t = jump2.end = 0;
  malloc2 ( IjQ_, QUEUE, jump->t*2, "BMACE_malloc_mk_jQ", "IjQ_" );

    /* IjQ_[t] := X(i\cup v) - X(i) , jump2 := {v| i\cup v is frequent} */
  QUEUE_FE_LOOP_ ( *jump, jj, v ){
    if ( v <= uu ) break;
    QUEUE_init ( &IjQ_[jump2.t], QUEUE_LENGTH(*S)-jQ[v].t+2 );
    t=0; 
    QUEUE_FE_LOOP_ ( *S, j, u ){
      if ( jQ[v].q[t] == u ) t++;
      else QUEUE_ins_ ( &IjQ_[jump2.t], u );
    }
    QUEUE_ins_ ( &jump2, v );
  }
  return ( IjQ_ );
}

/* malloc jQ_ for S and indices >i */
QUEUE *BMACEclosed_malloc_mk_jQ_ ( QUEUE *S, int i, int flag ){
  int u, v, j, jj, t, f;
  QUEUE *IjQ_;
  jump2.s = jump2.t = 0;
  
    /* compute |X(i\cup v)| for each v */
  QUEUE_FE_LOOP_ ( *S, j, u ){
    for ( jj=0 ; jj<G.edge[u].t ; jj++ ){
      v = G.edge[u].q[jj];
      jQ[v].s++;
      if ( jQ[v].s == 1 ) QUEUE_ins_ ( &jump2, v );
    }
  }
  malloc2 ( IjQ_, QUEUE, jump2.t*2, "malloc_mk_jQ", "IjQ_" );
  qsort_int_ ( jump2.q, jump2.t );

    /* malloc IjQ_[t], jump2 := {v| i\cup v is frequent} */
  t=0; QUEUE_FE_LOOP_ ( jump2, j, u ){
    if ( jQ[u].s >= v1bound ){
      if ( jQ[u].s < S->t ){
        jump2.q[t] = u;
        QUEUE_init ( &IjQ_[t], S->t - jQ[u].s+2 );
        IjQ_[t].end = 0;
        t++;
      }
    }
    jQ[u].s = 0;
  }
  jump2.t = t;

    /* IjQ_[t] := X(i\cup v) - X(i)  */
  QUEUE_FE_LOOP_ ( *S, j, u ){
    t = G.edge[u].t-1;
    f = 0;
    QUEUE_FE_LOOP_ ( jump2, jj, v ){
      while ( G.edge[u].q[t] > v ){
        if ( t == 0 ) break;
        t--;
      }
      if ( G.edge[u].q[t] != v ){
        if ( (f==0) && (v<i) ){
          if ( flag == 0 ){
            IjQ_[jj].end = -1;
          } else {
            IjQ_[jj].end++;
          }
        } else {
          QUEUE_ins_ ( &IjQ_[jj], u );
          f = 1;
        }
      }
    }
  }
  if ( flag == 0 ){
    t = 0; QUEUE_FE_LOOP_ ( jump2, jj, v ){
      if ( IjQ_[jj].end != -1 ){
        jump2.q[t] = v;
        IjQ_[t] = IjQ_[jj];
        t++;
      }
    }
    jump2.t = t;
  }
  return ( IjQ_ );
}



/* update inverse jQ_ */
void BMACE_update_jump_ ( QUEUE *jump, int i, QUEUE *IjQ, QUEUE *jump_, QUEUE *IjQ_, int fq_, QUEUE *cl ){
  int j, jj, v, u, tt, *x, *y;

  QUEUE_FE_LOOP_ ( IjQ[i], j, u ){
    G.edge[u].s = -1;
  }
  jump_->s = jump_->t = 0;
  for ( jj=jump->s ; jj<i ; jj++ ){
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    y = IjQ[jj].q + IjQ[jj].t;
    for ( x=IjQ[jj].q ; x<y ; x++ )
        if ( G.edge[*x].s == 0 ) QUEUE_ins_ ( &IjQ[tt], *x );
    if ( IjQ[tt].t <= fq_ ){
      if ( IjQ[tt].t == 0 ){
        QUEUE_ins_ ( cl, jump->q[jj] );
      } else {
        IjQ_[jump_->t] = IjQ[tt];
        IjQ_[jump_->t].end = fq_ +v1bound - IjQ[tt].t;
        QUEUE_ins_ ( jump_, jump->q[jj] );
      }
    }
  }
  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = 0;
}

/* update inverse jQ_ */
int BMACEclosed_update_jump_ ( QUEUE *jump, int i, QUEUE *IjQ, QUEUE *jump_, QUEUE *IjQ_, int fq_ ){
  int j, jj, v, u, tt, *x, *y, f;

  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = -G.node_end;
  jump_->s = jump_->t = 0;

  for ( jj=jump->s ; jj<i ; jj++ ){
    f = IjQ[jj].end;
    if ( f > fq_ ) goto END0;
    for ( y=IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ) if ( f++ >= fq_ ) goto END0;
    }
    if ( f <= fq_ ){
      if ( f == 0 ){
        QUEUE_ins_ ( &cliq, jump->q[jj] );
      } else {
        QUEUE_ins_ ( jump_, jump->q[jj] );
        for ( f=0, y=IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ )
            G.edge[*x].s++;
        goto END0;
      }
    }
    END0:;
  }

  for ( jj=i+1 ; jj<jump->t ; jj++ ){
    if ( IjQ[jj].end > 0 ) goto END1;
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    for ( y = IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ){
        if ( G.edge[*x].s == 0 ) goto END1;
        QUEUE_ins_ ( &IjQ[tt], *x );
      }
    }
    if ( IjQ[tt].t <= fq_ ){
      if ( IjQ[tt].t == 0 ){
        for ( jj=0 ; jj<=i ; jj++ )
          for ( y=IjQ[jj].q+IjQ[jj].t,x=IjQ[jj].q ; x<y ; x++ ) G.edge[*x].s=0;
        return (1);
      } else {
        IjQ[tt].end = 0;
        IjQ_[jump_->t] = IjQ[tt];
        for ( y=IjQ[tt].q+IjQ[tt].t,x=IjQ[tt].q ; x<y ; x++ )
          if ( G.edge[*x].s >0 ) G.edge[*x].s=2;
        QUEUE_ins_ ( jump_, jump->q[jj] );
      }
    }
    END1:;
  }

  for ( j=0, jj=jump->s ; jj<i ; jj++ ){
    if ( jump->q[jj] != jump_->q[j] ) continue;
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    f = IjQ[jj].end;
    for ( y = IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ){
        if ( G.edge[*x].s == 1 ){
          f++;
          G.edge[*x].s = 0;
        } else {
          QUEUE_ins_ ( &IjQ[tt], *x );
        }
      }
    }
    IjQ_[j] = IjQ[tt];
    IjQ_[j].end = f;
    j++;
  }

  for ( jj=0 ; jj<j ; jj++ ){
    for ( y=IjQ_[jj].q+IjQ_[jj].t, x=IjQ_[jj].q ; x<y ; x++ ) G.edge[*x].s=0;
  }
  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = 0;
  return ( 0 );
}

/* update inverse jQ_ */
int BMACEmax_update_jump_ ( QUEUE *jump, int i, QUEUE *IjQ, QUEUE *jump_, QUEUE *IjQ_, int fq_ ){
  int j, jj, v, u, tt, *x, *y, flag=0, f;

  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = -G.node_end;
  jump_->s = jump_->t = 0;

  for ( jj=jump->s ; jj<i ; jj++ ){
    f = IjQ[jj].end;
    if ( f > fq_ ) goto END0;
    for ( y=IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ) if ( f++ >= fq_ ) goto END0;
    }
    if ( f <= fq_ ){
      if ( f == 0 ){
        QUEUE_ins_ ( &cliq, jump->q[jj] );
      } else {
        QUEUE_ins_ ( jump_, jump->q[jj] );
        for ( f=0, y=IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ )
            G.edge[*x].s++;
        goto END0;
      }
    }
    END0:;
  }

  for ( jj=i+1 ; jj<jump->t ; jj++ ){
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    f = IjQ[jj].end;
    for ( y = IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ){
        if ( G.edge[*x].s == 0 ) f++; 
        else QUEUE_ins_ ( &IjQ[tt], *x );
      }
    }
    if ( f+IjQ[tt].t <= fq_ ){
      if ( f+IjQ[tt].t == 0 ){
        flag = 1;
        for ( jj=0 ; jj<=i ; jj++ )
          for ( y=IjQ[jj].q+IjQ[jj].t,x=IjQ[jj].q ; x<y ; x++ ) G.edge[*x].s=0;
        return (1);
      } else {
        IjQ[tt].end = f;
        IjQ_[jump_->t] = IjQ[tt];
        for ( y=IjQ[tt].q+IjQ[tt].t,x=IjQ[tt].q ; x<y ; x++ )
          if ( G.edge[*x].s >0 ) G.edge[*x].s=2;
        QUEUE_ins_ ( jump_, jump->q[jj] );
      }
    }
    END1:;
  }

  for ( j=0, jj=jump->s ; jj<i ; jj++ ){
    if ( jump->q[jj] != jump_->q[j] ) continue;
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    f = IjQ[jj].end;
    for ( y = IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ){
        if ( G.edge[*x].s == 1 ){
          f++;
          G.edge[*x].s = 0;
        } else {
          QUEUE_ins_ ( &IjQ[tt], *x );
        }
      }
    }
    IjQ_[j] = IjQ[tt];
    IjQ_[j].end = f;
    j++;
  }

  for ( jj=0 ; jj<j ; jj++ ){
    for ( y=IjQ_[jj].q+IjQ_[jj].t, x=IjQ_[jj].q ; x<y ; x++ ) G.edge[*x].s=0;
  }
  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = 0;
  return ( 0 );
}

/* update inverse jQ_ */
int BMACEmax_update_jump_o ( QUEUE *jump, int i, QUEUE *IjQ, QUEUE *jump_, QUEUE *IjQ_, int fq_ ){
  int j, jj, v, u, tt, *x, *y, flag=0, f;

  QUEUE_FE_LOOP_ ( IjQ[i], j, u ) G.edge[u].s = -G.node_end;
  jump_->s = jump_->t = 0;

  for ( jj=jump->s ; jj<i ; jj++ ){
    tt = jj+jump->t;
    IjQ[tt].t = 0;
    for ( y=IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s++ >= 0 ) QUEUE_ins_ ( &IjQ[tt], *x );
    }
    if ( IjQ[tt].t <= fq_ ){
      if ( IjQ[tt].t == 0 ){
        QUEUE_ins_ ( &cliq, jump->q[jj] );
      } else {
        IjQ_[j] = IjQ[tt];
        IjQ_[j].end = 0;
        QUEUE_ins_ ( jump_, jump->q[jj] );
      }
    }
  }

  for ( jj=i+1 ; jj<jump->t ; jj++ ){
    if ( (f=IjQ[jj].end) > fq_ ) continue;
    IjQ[tt=jj+jump->t].t = 0;
    for ( y = IjQ[jj].q+IjQ[jj].t, x=IjQ[jj].q ; x<y ; x++ ){
      if ( G.edge[*x].s >= 0 ){
        if ( G.edge[*x].s == 0 ){
          f++;
        } else {
          QUEUE_ins_ ( &IjQ[tt], *x );
        }
      }
    }
    if ( f+IjQ[tt].t <= fq_ ){
      if ( f+IjQ[tt].t == 0 ){
        flag = 1;
        goto END;
      } else {
        IjQ[tt].end = f;
        IjQ_[jump_->t] = IjQ[tt];
        QUEUE_ins_ ( jump_, jump->q[jj] );
      }
    }
  }


  END:;
  QUEUE_FE_LOOP_ ( *jump_, jj, v ){
    printf ("IJQ %d,%d (%d)\n", IjQ_[jj].end, IjQ_[jj].t, jj );
  }
  for ( jj=0 ; jj<=i ; jj++ )
     for ( y=IjQ[jj].q+IjQ[jj].t,x=IjQ[jj].q ; x<y ; x++ ) G.edge[*x].s=0;
  return ( 0 );
}

/* check closed itemset or not */
int BMACEclosed_maximal_ ( QUEUE *jump, int i, QUEUE *IjQ, QUEUE *jump_, QUEUE *IjQ_, int fq ){
  int j, jj, v, u, tt, flag = 0;

  for ( j=i+1 ; j<jump->t ; j++ ) jQ[jump->q[j]].t = IjQ[j].end;
  QUEUE_FE_LOOP_ ( IjQ[i], j, u ){
    for ( jj=0 ; jj<G.edge[u].t ; jj++ ){
      if ( (v=G.edge[u].q[jj]) >= jump->q[i] ) break;
      jQ[v].t--;
    }
  }
  if ( i==0 ){
    for ( j=1 ; j<jump->t ; j++ ){
      v = jump->q[j];
      if ( jQ[v].t >= v1bound ){
        if ( jQ[v].t == fq ) return (1);
        flag = 2;
      }
    }
  } else {
    for ( j=i+1 ; j<jump->t ; j++ ){
      v = jump->q[j];
      if ( jQ[v].t >= v1bound ){
        if ( jQ[v].t == fq ) return (1);
        flag = 2;
        IjQ_[jump_->t].end = jQ[v].t;
        QUEUE_ins_ ( jump_, v );
      }
    }
  }
  return ( 0 );
}

/* free jQ for S and indices >i */
void BMACE_jQ_free ( QUEUE *S, int i ){
  int u, v, j, jj;
  QUEUE_FE_LOOP_ ( *S, j, u ){
    for ( jj=G.edge[u].t-1 ; (v=G.edge[u].q[jj])>i ; jj-- )
        QUEUE_end ( &jQ[v] );
  }
}

/* add items >i whose addition not changes X(cliq) to cl */
void BMACE_addcliq ( QUEUE *S, QUEUE *cl, int i ){ 
  int ii, u;
  for ( ii=jump.t-1 ; ((u=jump.q[ii])>i)&&ii>=0 ; ii-- );
  for ( ii++ ; ii<jump.t ; ii++ ){
    u=jump.q[ii];
    if ( QUEUE_LENGTH(jQ[u]) >= QUEUE_LENGTH(*S) )
        QUEUE_ins_ ( cl, u );
  }
}

/* make list of transactions including v by using el */
void BMACE_mk_inclist ( int v ){
  int u, w, uu;
  QUEUE_RMALL ( incQ );
  for ( u=el[v] ; u!=-1 ; u=uu ){
    QUEUE_ins_ ( &incQ, u );
    uu = el[u];
    if ( G.edge[u].end < G.edge[u].t-1 ){
      G.edge[u].end++;
      w = G.edge[u].q[G.edge[u].end];
      el[u] = el[w];
      el[w] = u;
    }
  }
}

/* make column of adjacency matrix for high-degree vertices */
BARRAY *SGRAPH_sparse_badjmat ( SGRAPH *G, int BAratio ){
  BARRAY *BA;
  int i, j, u;
  malloc2 ( BA, BARRAY, G->node_end, "BMACEclosed", "BA" ); /* alloc memory for adjacency matrix */
  for ( vend=G->node_end,i=0 ; i<G->node1_num ; i++ ){
    if ( QUEUE_LENGTH(G->edge[i]) >= (G->node_end-G->node1_num)/BAratio ){
      BARRAY_init ( &BA[i], G->node_end-G->node1_num );
      BARRAY_reset_interval ( &BA[i], 0, G->node_end-G->node1_num-1 );
      QUEUE_FE_LOOP_ ( G->edge[i], j, u )
          BARRAY_SET ( BA[i], u-G->node1_num );
      if ( vend == G->node_end ) vend = i;
    } else BA[i].a = NULL;
  }
  return ( BA );
}

/*********************************************************************/
/******************      all frequent item sets      *****************/
/*********************************************************************/
/*********************************************************************/

  /* output current frequent set */
void BMACEfreq_output ( int siz, int s ){
  int i, e;
  if ( s == cliq_.t ){ FIM_perm_output ( &cliq, siz ); return; }
  BMACEfreq_output ( siz, s+1 );
  QUEUE_ins_ ( &cliq, cliq_.q[s] );
  BMACEfreq_output ( siz, s+1 );
  cliq.t--;
}

int iters=0;
/* Enumeration of Frequent Item Set (inverse version): Iteration */
void BMACEfreq_iter_ ( QUEUE *S, QUEUE *jump, QUEUE *IjQ, int fq ){
  int j, u, t, v, jj;
  QUEUE jump_, *IjQ_=0;
  iters++;

  BMACEfreq_output ( fq, 0 );
  if ( QUEUE_LENGTH(*jump) == 0 ) return;       /* If no cliq\cup j is frequent, output and return */
    /* Output current solution, and current solution \cup j where j is the largest s.t. cliq\cup j is frequent */
  QUEUE_ins_ ( &cliq, jump->q[0] );
  BMACEfreq_output ( fq-IjQ[jump->s].t, 0 );
  cliq.t --;
  if ( QUEUE_LENGTH(*jump) == 1 ) return; /* return no more j s.t. cliq \cup j is frequent */
     /* allocate memory for array */
  QUEUE_init ( &jump_, jump->t );
  malloc2 ( IjQ_, QUEUE, jump->t*2, "freq_iter", "IjQ_" );

  for ( j=jump->s+1 ; j<jump->t ; j++ ){
    QUEUE_init ( &IjQ[j-1+jump->t], QUEUE_LENGTH(IjQ[j-1]) );
      /* update variables for cliq \cup jump->q[j] */
    t = cliq_.t;  
    BMACE_update_jump_ ( jump, j, IjQ, &jump_, IjQ_, fq-IjQ[j].t-v1bound, &cliq_);
    QUEUE_ins_ ( &cliq, jump->q[j] );  /* add i to cliq(=current solution) */
      BMACEfreq_iter_ ( &IjQ[j], &jump_, IjQ_, fq-IjQ[j].t ); /* recursive call for cliq+{jump->q[j]}*/
    cliq.t--;   /* recover current solution */
    cliq_.t = t;
  }
     /* clear used variables */
  for ( j=jump->s ; j<jump->t-1 ; j++ ) QUEUE_end ( &IjQ[j+jump->t] );
  free ( IjQ_ );
  QUEUE_end ( &jump_ );
}

  /* Enumeration of Frequent Sets: iteration */
void BMACEfreq_iter ( QUEUE *S, int i, QUEUE *jump ){
  int j, u, th, t=cliq_.t;
  QUEUE *IjQ_, jump_;
  QUEUE_ins_ ( &cliq, i );  /* add i to cliq(=current solution) */
  QUEUE_init ( &jump_, jump->t );
  th = BMACE_mk_jump ( S, i, jump, &jump_, &cliq_ );   /* set jQ[j],j>i to X(cliq + {i}) */
  if ( (th%100 > 50) && ((th/100) > 2 ) ){
    IjQ_ = BMACE_malloc_mk_jQ_ ( S, i, &jump_, i );  /* malloc necessary jQ_ */
    BMACEfreq_iter_ ( S, &jump2, IjQ_, QUEUE_LENGTH(*S)); /* recursive call for {i} */
    QUEUE_F_LOOP_ ( jump2, u ) QUEUE_end ( &IjQ_[u] );
    free ( IjQ_ );
  } else {
    BMACEfreq_output ( QUEUE_LENGTH(*S), 0 ); /* output all cliq\subseteq K \subseteq cliq' */
    QUEUE_FE_LOOP_ ( jump_, j, u ){
      BMACEfreq_iter (&jQ[u], u, &jump_); /* recursive call for cliq + u */
    }
  }
  QUEUE_end ( &jump_ );
  cliq_.t = t;
  cliq.t--;   /* remove the item added at the begining of this iteration */
}

/* Enumeration of Frequent Item Set (inverse version): main */
void BMACEfreq (){
  int u, i, t;
  QUEUE_init ( &cliqtmp , G.node_end-G.node1_num);
  for ( i=G.node1_num ; i<G.node_end ; i++ ){ 
    BMACE_mk_inclist ( i );                    /* incQ := X({i}) */
    if ( QUEUE_LENGTH(incQ) >= v1bound ){      /* if {i} is frequent */
      cliq.t = cliq_.t = 0;
      jump.s = jump.t = 0;
      for ( u=i ; u<G.node_end ; u++ ) jQ[u].t = 0;
      BMACE_jQ_malloc ( &incQ, i, &cliq_ );  /* malloc necessary jQ */
      BMACEfreq_iter ( &incQ, i, &jump );   /* recursive call for {i} */
      BMACE_jQ_free ( &incQ, i );    /* free jQ */
    }
  }
}


/*********************************************************************/
/******************     frequent closed item sets    *****************/
/*********************************************************************/
/*********************************************************************/

/* Enumeration of Frequent Item Set (inverse version): Iteration */
void BMACEclosed_iter_ ( QUEUE *S, QUEUE *jump, QUEUE *IjQ, int fq, int vv ){
  int j, u, t, v, jj;
  QUEUE jump_, *IjQ_=0;
  iters++;
  FIM_perm_output ( &cliq, fq );
  if ( QUEUE_LENGTH(*jump)==0 || jump->q[0]<=vv ) return;
  QUEUE_init ( &jump_, jump->t );
  malloc2 ( IjQ_, QUEUE, jump->t*2, "closed_iter_", "IjQ_" );
  for ( j=0 ; j<jump->t ; j++ )
      QUEUE_init ( &IjQ[j+jump->t], QUEUE_LENGTH(IjQ[j]) );
  for ( j=0 ; j<jump->t ; j++ ){
    if ( (v=jump->q[j]) <= vv ) break;
    t = cliq.t;
       /* update variables for cliq \cup jump->q[j] */
    if ( BMACEclosed_update_jump_ ( jump, j, IjQ, &jump_, IjQ_, fq-IjQ[j].t-IjQ[j].end-v1bound ) != 1 ){
      QUEUE_ins_ ( &cliq, jump->q[j] );  /* add i to cliq(=current solution) */
      BMACEclosed_iter_ ( &IjQ[j], &jump_, IjQ_, fq-IjQ[j].t-IjQ[j].end, v ); /* recursive call for cliq+{v}*/
    }
    cliq.t = t;
  }
     /* clear used variables */
  for ( j=0 ; j<jump->t ; j++ ) QUEUE_end ( &IjQ[j+jump->t] );
  free ( IjQ_ );
  QUEUE_end ( &jump_ );
}

  /* Enumeration of Frequent Sets: iteration */
void BMACEclosed_iter ( QUEUE *S, int i, QUEUE *jump ){
  int j, u, th, t=cliq.t;
  QUEUE *IjQ_, jump_;

  QUEUE_ins_ ( &cliq, i );  /* add i to cliq(=current solution) */
  QUEUE_init ( &jump_, jump->t );
  th = BMACE_mk_jump ( S, i, jump, &jump_, &cliq );   /* set jQ[j],j>i to X(cliq + {i}) */
  for ( j=t ; j<cliq.t ; j++ ) G.edge[cliq.q[j]].end = 1;  /* mark items of cliq */
  if ( (th%100 > 50) && ((th/100) > 5 ) && t == 0 ){
    IjQ_ = BMACEclosed_malloc_mk_jQ_ ( S, i, 0 );  /* malloc necessary jQ_ */
    BMACEclosed_iter_ ( S, &jump2, IjQ_, QUEUE_LENGTH(*S), i); /* recursive call for {i} */
    QUEUE_F_LOOP_ ( jump2, u ) QUEUE_end ( &IjQ_[u] );
    free ( IjQ_ );
  } else {
    iters++;
    FIM_perm_output ( &cliq, QUEUE_LENGTH(*S) ); /* output all cliq\subseteq K \subseteq cliq' */
    QUEUE_FE_LOOP_ ( jump_, j, u ){
      if ( u<i ) break;
      if ( BMACEclosed_maximal(&jQ[u], u) == -1 ){ /* closed item set? */
        BMACEclosed_iter (&jQ[u], u, &jump_); /* recursive call for cliq + u */
      }
    }
  }
  QUEUE_end ( &jump_ );
  for ( j=t ; j<cliq.t ; j++ ) G.edge[cliq.q[j]].end = 0; /* delete marks */
  cliq.t = t;   /* remove the item added at the begining of this iteration */
}

/* Enumeration of Frequent Item Set (inverse version): main */
void BMACEclosed (){
  int u, i, t, j;
  QUEUE *jQ_;
  BA = SGRAPH_sparse_badjmat ( &G, BAratio );
  for ( i=G.node1_num ; i<G.node_end ; i++ ){ 
    BMACE_mk_inclist ( i );                    /* incQ := X({i}) */
    if ( QUEUE_LENGTH(incQ) >= v1bound ){      /* if {i} is frequent */
      cliq.t=0;
      jump.s = jump.t = 0;
      for ( u=i ; u<G.node_end ; u++ ) jQ[u].t = 0;
      qsort_int ( incQ.q, incQ.t );
      if ( BMACEclosed_maximal(&incQ, i) == -1 ){ /* closed item set? */
        BMACE_jQ_malloc ( &incQ, i, &cliq );  /* malloc necessary jQ */
        BMACEclosed_iter ( &incQ, i, &jump );   /* recursive call for {i} */
        BMACE_jQ_free ( &incQ, i );    /* free jQ */
      }
    }
//    if ( iters>0 ) printf ("iter %d\n", iters);
  }
     /* clear used variables */
  for ( i=vend ; i<G.node1_num ; i++ ) BARRAY_end ( &BA[i] );
  free ( BA );
}




/***************************************************************************/
/********************      maximal frequent item sets      *****************/
/***************************************************************************/
/***************************************************************************/


/* Enumeration of Frequent Item Set (inverse version): Iteration */
void BMACEmax_iter_ ( QUEUE *S, QUEUE *jump, QUEUE *IjQ, int fq, int vv ){
  int j, u, t, v, jj;
  QUEUE jump_, *IjQ_=0;
  iters++;
  if ( QUEUE_LENGTH(*jump)==0 ){
    FIM_perm_output ( &cliq, QUEUE_LENGTH(*S) ); /* output cliq */
  }
  if ( jump->q[0]<=vv ) return;
  QUEUE_init ( &jump_, jump->t );
  malloc2 ( IjQ_, QUEUE, jump->t*2, "closed_iter_", "IjQ_" );
  for ( j=0 ; j<jump->t ; j++ ){
    QUEUE_init ( &IjQ[j+jump->t], QUEUE_LENGTH(IjQ[j]) );
  }
  for ( j=0 ; j<jump->t ; j++ ){
    if ( (v=jump->q[j]) <= vv ) break;
    t = cliq.t;
       /* update variables for cliq \cup jump->q[j] */
    if ( BMACEmax_update_jump_ ( jump, j, IjQ, &jump_, IjQ_, fq-IjQ[j].t-IjQ[j].end-v1bound ) != 1 ){
      QUEUE_ins_ ( &cliq, jump->q[j] );  /* add i to cliq(=current solution) */
      BMACEmax_iter_ ( &IjQ[j], &jump_, IjQ_, fq-IjQ[j].t-IjQ[j].end, v ); /* recursive call for cliq+{v}*/
    }
    cliq.t = t;
  }
     /* clear used variables */
  for ( j=0 ; j<jump->t ; j++ ) QUEUE_end ( &IjQ[j+jump->t] );
  free ( IjQ_ );
  QUEUE_end ( &jump_ );
}




  /* Enumeration of Frequent Sets: iteration */
void BMACEmax_iter ( QUEUE *S, int i, QUEUE *jump ){
  int j, u, th, t=cliq.t;
  QUEUE *IjQ_, jump_;

  QUEUE_ins_ ( &cliq, i );  /* add i to cliq(=current solution) */
  QUEUE_init ( &jump_, jump->t );
  th = BMACE_mk_jump ( S, i, jump, &jump_, &cliq );  /* set jQ[j],j>i to X(cliq + {i}) */
  if ( (th%100 > 50) && ((th/100) > 5 ) && t == 0){
    IjQ_ = BMACEclosed_malloc_mk_jQ_(S, i, 1 ); /* malloc necessary jQ_ */
    BMACEmax_iter_ ( S, &jump2, IjQ_, QUEUE_LENGTH(*S), i); /* recursive call for {i} */
    QUEUE_F_LOOP_ ( jump2, u ) QUEUE_end ( &IjQ_[u] );
    free ( IjQ_ );
  } else {
    iters++;
    if ( ( jump_.t==0 ) && BMACEmax_maximal(S, i) == 0){
      FIM_perm_output ( &cliq, QUEUE_LENGTH(*S) ); /* output cliq */
      goto END;
    }
    for (j=t ;j<cliq.t ; j++) G.edge[cliq.q[j]].end=1; /* mark items of cliq */
    QUEUE_FE_LOOP_ ( jump_, j, u ){
      if ( u < i ) break;
      if ( BMACEclosed_maximal(&jQ[u], u) == -1 ){ /* closed item set? */
        BMACEmax_iter (&jQ[u], u, &jump_ ); /* recursive call for cliq + u */
      }
    }
    for ( j=t ; j<cliq.t ; j++ ) G.edge[cliq.q[j]].end = 0; /* delete marks */
  }
  END:;
  QUEUE_end ( &jump_ );
  cliq.t = t;   /* remove the item added at the begining of this iteration */
}


/* Enumeration of Frequent Item Set (inverse version): main */
void BMACEmax (){
  int u, i, t, j, flag;
  QUEUE *jQ_;

  BA = SGRAPH_sparse_badjmat ( &G, BAratio );
  for ( i=G.node1_num ; i<G.node_end ; i++ ){ 
    BMACE_mk_inclist ( i );                    /* incQ := X({i}) */
    if ( QUEUE_LENGTH(incQ) >= v1bound ){      /* if {i} is frequent */
      cliq.t = 0;
      jump.s = jump.t = 0;
      for ( u=i ; u<G.node_end ; u++ ) jQ[u].t = 0;
      qsort_int ( incQ.q, incQ.t );
      if ( BMACEclosed_maximal (&incQ, i) == -1 ){ /* closed item set? */
        BMACE_jQ_malloc ( &incQ, i, &cliq );  /* malloc necessary jQ */
        BMACEmax_iter ( &incQ, i, &jump );   /* recursive call for {i} */
        BMACE_jQ_free ( &incQ, i );    /* free jQ */
      }
    }
//    if ( iters>0 ) printf ("iter %d(%d)\n", iters, count);
  }
     /* clear used variables */
  for ( i=vend ; i<G.node1_num ; i++ ) BARRAY_end ( &BA[i] );
  free ( BA );
}




/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/



main ( int argc, char *argv[] ){
  int i, x=0, c=3;
  if ( argc <3 ) goto ERROR;
  v1bound = atoi(argv[2]);
  G.edge = QUEUE_load2 ( argv[1], &(G.node_end), &(G.node1_num), 3 );
  G.node_end += G.node1_num;
  FREQSET_init ( 6); 
  if ( argc == 3 ) fp = NULL;
  else fopen2w ( fp, argv[3], "BMACE_outputfile");
#if FREQSET_PROBLEM == 1
  BMACEfreq ();
#endif
#if FREQSET_PROBLEM == 2
  BMACEclosed ();
#endif
#if FREQSET_PROBLEM == 3
  BMACEmax ();
#endif
#if FREQSET_OUTPUT == 3
  printf ("%d\n", count);
#else 
  FIM_output2 ();
#endif
  if ( fp ) fclose ( fp );
  FREQSET_end ( 6 );
  intarray2_free_buf ();
  exit (0);

  ERROR:;
  printf ("usege: BMACE input-filename threshold output-filename\n");
  exit (1);
}

