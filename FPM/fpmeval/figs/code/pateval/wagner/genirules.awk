{
   i=NF-9;
   j=NF-7;
   printf ("{%f, %f,\"",$i/150.0,$j);
   for (k=2; k<=i; k++){
     printf ("%s",$k);
   }
   printf("\"},\n");
}

