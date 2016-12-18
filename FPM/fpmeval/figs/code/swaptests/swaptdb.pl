# A program for computing swap-randomized transaction databases

use strict;
use warnings;

die "parameters: tdbfile tdbprefix iterlen steps\n" unless scalar(@ARGV)==4;

my ($tdbfile,$tdbprefix,$iterlen,$steps)=@ARGV;

print STDERR "tdbfile: $tdbfile\n";
print STDERR `du -hs $tdbfile`;
print STDERR "tdbprefix: $tdbprefix\n";
print STDERR "iterlen: $iterlen steps: $steps\n";

#a subroutine for swapping
sub swapedges(\%\@\@$) {
    my ($eref,$jref,$iref,$n)=@_;

    my $a=int(rand($n));
    my $b=int(rand($n));

    my ($aj,$ai)=($$jref[$a],$$iref[$a]);
    my ($bj,$bi)=($$jref[$b],$$iref[$b]);

    #edges ($aj,$ai) and ($bj,$bi) are swappable iff there are no edges ($aj,$bi) and ($bj,$ai)
    if(not defined $$eref{"$aj,$bi"} and 
       not defined $$eref{"$bj,$ai"}) {
	#delete edges ($aj,$ai) and ($bj,$bi)
	delete $$eref{"$aj,$ai"};
	delete $$eref{"$bj,$bi"};
	#add edges ($aj,$bi) and ($bj,$ai)
	$$eref{"$aj,$bi"}=$a;
	$$eref{"$bj,$ai"}=$b;
	#replace $ai with $bi and $bi with $ai
	$$iref[$a]=$bi;
	$$iref[$b]=$ai;
	#inform about successful swap
	return 1;
    }
    #inform about a failed swap
    return 0;
}

#
# THE MAIN PROGRAM
#
print STDERR "reading data...\n";
my %e=();
my $n=0;
my $rows=0;
my @j=();
my @i=();
#my %items=();

open FIN,"$tdbfile";
#read the data
while(<FIN>) {
    chomp;
    my @row=split " ";

    for my $i (@row) {
	my $j=$rows; 
	$e{"$j,$i"}=$n; 
	$j[$n]=$j;
	$i[$n]=$i;
	$n++;
#	$items{$i}=0 unless defined $items{$i};
#	$items{$i}++;

    }
    $rows++;
}
close FIN;

print STDERR "data read.\n";

my $swaps=0;
srand;

#start swapping
print STDERR "starting to swap...\n";
for(my $i=0; $i<=$iterlen*$steps; $i++) {
    #if the time is right...
    if(($i%$iterlen)==0) {
	my $k=0; my @row=();
	
# writing the swapped data
	
	$tdbfile=sprintf("%s.%d.dat",$tdbprefix,$i/$iterlen);
	
	open FOUT, ">$tdbfile";
	for(my $l=0; $l<=$n; $l++) {
	    if($l<$n and $k==$j[$l]) {
		push @row,$i[$l];
	    } else {
		@row=sort {$a <=> $b} @row;
		print FOUT join(" ",@row),"\n";
		
		if($l<$n) {
		    @row=($i[$l]); $k=$j[$l];
		}
	    }
	}
	close FOUT;
	
# print statistics

	if($i>0) {
	    print "$i $swaps ",$swaps/$i, " ",$swaps/$n,"\n";
	} else {
	    print "0 0 0 0\n";
	}
    }
    
    #try to swap regardless of the iteration
    $swaps+=swapedges(%e,@j,@i,$n);
}
