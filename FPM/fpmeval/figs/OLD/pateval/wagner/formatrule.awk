{
  printf("%s ",$1)
  str=$2
  last = $1
  i=3;
  while ($i != "conf:"){
    str = str " " $i;
    last = $i
    #print i, $i, str
    i++;
  }
  printf(" supp: %.3f ",last/150.0);
  for (;i<=NF;i++){
    printf("%s ",$i)
  } 
  print str
}
