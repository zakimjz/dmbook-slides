BEGIN{
  maxsz = -1;
}
{
  sz = NF-1
  num[sz]++
  count[sz,num[sz]] =$NF  
  if (sz>maxsz){
    maxsz = sz
  }
  print sz, num[sz], count[sz,num[sz]]
}
END{
   for (i= 1; i<=maxsz; i++){
     avg[i] = 0
     sd[i] = 0;
     for (j=1;j<=num[i];j++){
       avg[i] += count[i,j]
     }
     avg[i] = avg[i]/num[i]
     for (j=1;j<=num[i];j++){
       sd[i] = sd[i] + (count[i,j]-avg[i])*(count[i,j]-avg[i])
     } 
     sd[i] = sd[i]/num[i]
   }
    
   for (i=1; i<=maxsz;i++){
     print i, num[i], avg[i], sd[i]
   }
}
