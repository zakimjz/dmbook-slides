#!/usr/bin/env perl

use strict;

print "#cluster 1\n";
print "0 3 0\n";
print "0 3 2\n";
print "0 5 2\n";
print "0 5 0\n";
print "0 3 0\n";

print "\n\n";
print "4 3 0\n";
print "4 3 2\n";
print "4 5 2\n";
print "4 5 0\n";
print "4 3 0\n";
print "\n\n";

print "0 3 0\n";
print "4 3 0\n\n\n";
print "0 3 2\n";
print "4 3 2\n\n\n";
print "0 5 2\n";
print "4 5 2\n\n\n";
print "0 5 0\n";
print "4 5 0\n\n\n";

print "#cluster 2\n";
print "0 0 2\n";
print "1 0 2\n";
print "1 1 2\n";
print "0 1 2\n";
print "0 0 2\n";
print "\n\n";
print "0 0 4\n";
print "1 0 4\n";
print "1 1 4\n";
print "0 1 4\n";
print "0 0 4\n";
print "\n\n";

print "0 0 2\n";
print "0 0 4\n\n\n";
print "1 0 2\n";
print "1 0 4\n\n\n";
print "1 1 2\n";
print "1 1 4\n\n\n";
print "0 1 2\n";
print "0 1 4\n\n\n";


print "#cluster1 (0,3,0) to (4,5,2)\n";
for (my $i=0; $i < 40; $i++){
  my ($x, $y, $z) = (rand, rand, rand);
  $x = $x*4;
  $y = $y*2+3;
  $z = $z*2;
  print "$x $y $z\n";
}

print "\n\n#cluster2 (0,0,2) to (1,1,4)\n";
for (my $i=0; $i < 20; $i++){
  my ($x, $y, $z) = (rand, rand, rand);
  $z = $z*2+2;
  print "$x $y $z\n";
}
