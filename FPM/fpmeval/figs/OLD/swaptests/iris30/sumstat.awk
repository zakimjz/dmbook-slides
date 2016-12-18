BEGIN{
  maxszitemset = -1
} 
{
  szitemset = $1
  numitemset = $2
  avgsupp = $3
  avgstd = $4
  
  nitem[szitemset] += numitemset
  asupp[szitemset] += avgsupp
  astd[szitemset] += avgstd 
  count[szitemset]++

  if (szitemset>maxszitemset){
    maxszitemset = szitemset
  }

}
END{
  for(i=1; i<=maxszitemset; i++){
    print i, count[i], nitem[i]/count[i], asupp[i]/count[i], astd[i]/count[i]
  }
}
