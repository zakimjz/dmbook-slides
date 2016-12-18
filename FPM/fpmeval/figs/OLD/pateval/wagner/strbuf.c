#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "strbuf.h"
#include "error.h"

char strbuf[STRBUFSZ];
int laststr = 0;
char * strbufentry[STRBUFENTRY];
int laststrentry=0;

void printstrbuf(){
  int i;
  char * aux;
  for (i=0; i<laststr; i+=strlen(aux)+1){
    aux = &(strbuf[i]);
    printf("%d %s\n",i,aux);
  }
}

void printstrbuffile(FILE * out){
  int i;
  char * aux;
  for (i=0; i<laststrentry; i++){
    fprintf(out,"%d %s\n",i,strbufentry[i]);
  }
}

char * checkstr(str)
char * str;
{
  int i;
  char * aux;
  for (i=0; i<laststr; i+=strlen(aux)+1){
    aux = &(strbuf[i]);
    if (!strcmp(aux,str)) return  &(strbuf[i]);
  }
  return (char*) 0;
}

char * checkstr_entry(str,pos)
char * str;
int * pos;
{
  int i;
  char * aux;
  for (i=0; i<laststrentry; i++){
    aux = strbufentry[i];
    if (!strcmp(aux,str)) {
      (*pos) = i; 
      return aux;
    }
  }
  (*pos) = -1;
  return (char*) 0;
}

char * getstr_entry(pos)
int pos;
{
  if ((pos < 0) || (pos >= laststrentry))
    fatal ("getstr_entry:entry out of range (%d - %d..%d)",pos, 0,laststrentry);
  return strbufentry[pos];
}

char * addstr(str)
char * str;
{
  int start,i;
  start = laststr;
  if (laststr+strlen(str)+1>STRBUFSZ){
    printstrbuf();
    fatal ("String buffer overflow");
  }
  laststr += strlen(str)+1;
  for (i=0; i<strlen(str); i++)
    strbuf[i+start] = str[i];
  strbuf[laststr-1] = 0;
  return  &(strbuf[start]);
}


char * addstr_entry(str,pos)
char * str;
int * pos;
{
  int start,i;
  start = laststr;
  if (laststr+strlen(str)+1>STRBUFSZ){
    printstrbuf();
    fatal ("String buffer overflow");
  }
  laststr += strlen(str)+1;
  for (i=0; i<strlen(str); i++)
    strbuf[i+start] = str[i];
  strbuf[laststr-1] = 0;
  strbufentry[laststrentry] = &(strbuf[start]);
  (*pos) = laststrentry;
  laststrentry++;
  if (laststrentry > STRBUFENTRY) fatal("String buffer entry overflow");
  return  &(strbuf[start]);
}

