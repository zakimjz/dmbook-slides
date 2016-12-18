#ifndef STRBUFH

#define STRBUFH

#define STRBUFSZ 1000000
#define STRBUFENTRY 50000 /* average of 20 char per entry */
#define VARNAME(i) (strbuf+i)

extern char strbuf[STRBUFSZ];
extern char * strbufentry[STRBUFENTRY];
extern int laststr;
extern int laststrentry;

void printstrbuf(void);
void printstrbuffile(FILE *);
char * checkstr(char*);
char * addstr(char*);
char * checkstr_entry(char*,int*);
char * addstr_entry(char*,int*);
char * getstr_entry(int);

#endif
