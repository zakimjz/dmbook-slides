BEGIN{
  numitem = 5
}
{
  for (i=1; i<=numitem; i++){
    present[i] = 0
  }
  for (i=1; i<NF; i++){
    present[$i] = 1
  }
  for (i=1; i<=numitem;i++){
    if (present[i]==1){
      printf("%d ",i);
    }
  } 
  printf("%s\n",$NF);
}
