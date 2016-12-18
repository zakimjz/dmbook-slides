BEGIN{
  numitem = 100
  maxitem = 10
}
{
  for (i=1; i<=numitem; i++){
    present[i] = 0
  }
  for (i=1; i<NF; i++){
    present[$i] = 1
  }
  printf("\{%d,\{",NF-1);
  for (i=1; i<=numitem;i++){
    if (present[i]==1){
      printf("%d,",i);
    }
  } 
  for (i=NF; i<maxitem;i++){
      printf("-1,");
  } 
  printf("\},%d\},\n",$NF);
}


