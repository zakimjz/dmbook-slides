BEGIN{ 
  prev = -1
  count = 0
} 
{
  if ((prev != -1) && (prev != $1)){
   print prev " " count
   count = 0
  }
  prev = $1
  count = count + 1
}
