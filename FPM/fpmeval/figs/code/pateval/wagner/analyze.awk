BEGIN{
bm["sl1"] = 1
bm["sl2"] = 2
bm["sl3"] = 4
bm["sw1"] = 8
bm["sw2"] = 16
bm["sw3"] = 32
bm["pl1"] = 64
bm["pl2"] = 128
bm["pl3"] = 256
bm["pw1"] = 512
bm["pw2"] = 1024
bm["pw3"] = 2048
bm["cl1"] = 4096
bm["cl2"] = 8192
bm["cl3"] = 16384
numis=0
}
{
  bmval = 0
  for (i=2; i<=NF; i++){
    if ($i in bm){
      bmval += bm[$i]
    }
    if ($i == "conf:"){
     j = i-1;
     supp = $j/150.0
     j = i+1;
     conf = $j
    }
    if ($i == "lift:"){
     j = i+1;
     lift = $j
    }
    if ($i == "lev:"){
     j = i+1;
     lev = $j
    }
    if ($i == "conv:"){
     j = i+1;
     conv = $j
    }
  }  
  print bmval, supp, conf, lift, lev, conv, $0
  for (i=1; i<=numis; i++){
    if (id[i] == bmval){
      pos = i
    } else {
      numis++;
      id[numis] = bmval
      pos = numis;
    }
  }
  count[pos] ++;
  vsupp[pos,count[pos]] = supp 
  vconf[pos,count[pos]] = conf 
}
END{
   for (i=1; i<=numis; i++){
     is = id[i]
     for (j=1; j<=count[i];i++){
       
     }
   }
}
