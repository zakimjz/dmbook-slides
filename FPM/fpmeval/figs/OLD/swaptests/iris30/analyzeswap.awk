{
  freqorig[NR] = $1
  freqswap[NR] = $2
  relerrorig[NR] = $3
  relerrswap[NR] = $4
  avgcorrorig[NR] = $5
  avgcorrswap[NR] = $6
  swaplift[NR] = $7
  liftswapped[NR] = $8
  ratioorig[NR] = $9
  ratioswap[NR] = $10
  itemsetsz[NR] = NF - 11
  itemset[NR] = ""
  for (i=12; i<=NF; i++){
    itemset[NR] = itemset[NR] " " $i 
  }
  print freqorig[NR],freqswap[NR],relerrorig[NR],relerrswap[NR],avgcorrorig[NR],avgcorrswap[NR],swaplift[NR],liftswapped[NR],ratioorig[NR],ratioswap[NR],itemsetsz[NR],itemset[NR]
}
