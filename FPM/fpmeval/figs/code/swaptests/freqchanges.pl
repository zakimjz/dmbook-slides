use strict;
use warnings;

die "tdb tmpfile1 tmpfile2 minsupp minsupporig minsuppswapped stepsize\n" unless scalar(@ARGV)==7;

my ($tdbfile,$tmpfile1,$tmpfile2,$minsupp,$minsupporig,$minsuppswapped,$stepsize)=@ARGV;

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

sub count($\%) {
    my ($set,$rref)=@_;
    my @set=split " ",$set;

    my $supp=0;
    for my $row (keys %$rref) {
	my @row=split " ",$row;

	next if scalar(@row)<scalar(@set);

	my $i=0; my $j=0;

	while($i<scalar(@set) && $j<scalar(@row)) {
	    if($set[$i]==$row[$j]) {
		$i++; $j++;
	    } elsif($set[$i]>$row[$j]) {
		$j++;
	    } else {
		$j=scalar(@row);
	    }
	}

	$supp+=$$rref{$row} if $i==scalar(@set);
    }
    return $supp;
}

sub loadfreq($\%) {
    my ($fin,$fref)=@_;

    open FIN, "$fin";
    while(<FIN>) {
	chomp;
	my @row=split " ";
	my $supp=pop @row;
	$supp=substr($supp,1,length($supp)-2);
	$$fref{join(" ",sort {$a <=> $b} @row)}=$supp;
    }
    close FIN;
}

srand;

#
# THE MAIN PROGRAM
#
print STDERR "reading data...";
my %e=();
my $n=0;
my $rows=0;
my @j=();
my @i=();
my %items=();
my @rows=();
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
	$items{$i}=0 unless defined $items{$i};
	$items{$i}++;
    }

    push @rows,join(" ",sort {$a <=> $b} @row);

    $rows++;
}
close FIN;

for my $item (keys %items) {
    delete $items{$item} if $items{$item}<$minsupp;
}

my %rowsorig=();
for(my $i=0; $i<scalar(@rows); $i++) {
    my @row=split " ",$rows[$i];
    my @newrow=();
    for my $item (@row) {
	push @newrow,$item if defined $items{$item};
    }
    my $row=join(" ",@newrow);
    if(not defined $rowsorig{$row}) { $rowsorig{$row}=1; } else { $rowsorig{$row}++; }
}

print STDERR "$rows rows, $n ones, ",scalar(keys %rowsorig)," different rows w.r.t. frequent items\n";

print STDERR "computing frequent sets in the original data...";
`./fim_all $tdbfile $minsupporig $tmpfile1`;
my %freqorig=(); loadfreq($tmpfile1,%freqorig);
print STDERR scalar(keys %freqorig), " sets\n";

print STDERR "swapping...";
my $swaps=0;
for(my $i=0; $i<$stepsize; $i++) {
    $swaps+=swapedges(%e,@j,@i,$n);
}

# writing the swapped data
@rows=();
my @row=(); my $k=0;
open FOUT, ">$tmpfile2";
for(my $l=0; $l<=$n; $l++) {
    if($l<$n and $k==$j[$l]) {
	push @row,$i[$l];
    } else {
	@row=sort {$a <=> $b} @row;
	push @rows, join(" ",@row);
	print FOUT "$rows[$#rows]\n";
	if($l<$n) {
	    @row=($i[$l]); $k=$j[$l];
	}
    }
}
close FOUT;

my %rowsswapped=();
for(my $i=0; $i<scalar(@rows); $i++) {
    my @row=split " ",$rows[$i];
    my @newrow=();
    for my $item (@row) {
	push @newrow,$item if defined $items{$item};
    }
    my $row=join(" ",@newrow);
    if(not defined $rowsswapped{$row}) { $rowsswapped{$row}=1; } else { $rowsswapped{$row}++; }
}

print STDERR "$stepsize attempted, $swaps successful swaps, ", scalar(keys %rowsswapped)," different rows w.r.t. frequent items\n";


print STDERR "computing frequent sets in the swapped data...";
`./fim_all $tmpfile2 $minsuppswapped $tmpfile1`;

my %freqswapped=(); loadfreq($tmpfile1,%freqswapped);
print STDERR scalar(keys %freqswapped), " sets\n";

print STDERR "adding missing sets to swapped...\n";
my @missing=();
for my $set (keys %freqorig) {
    next if $freqorig{$set}<$minsupp;
    push @missing, $set unless defined $freqswapped{$set};
}
for(my $i=0; $i<scalar(@missing); $i++) {
    print STDERR "$i / ",scalar(@missing),"\r" if ($i%10)==0;
    $freqswapped{$missing[$i]}=count($missing[$i],%rowsswapped);
}

print STDERR "adding missing sets to orig...\n";
@missing=();
for my $set (keys %freqswapped) {
    next if $freqswapped{$set}<$minsupp;
    push @missing, $set unless defined $freqorig{$set};
}
for(my $i=0; $i<scalar(@missing); $i++) {
    print STDERR "$i / ",scalar(@missing),"\r" if ($i%10)==0;
    $freqorig{$missing[$i]}=count($missing[$i],%rowsorig);
}

my $total=$freqorig{''};
for my $set (keys %freqorig) {
    next if $freqorig{$set}<$minsupp;
    my @set=split " ",$set;
    next if scalar(@set)<2;
    
    my @pairs=();
    for(my $i=0; $i<scalar(@set); $i++) {
	for(my $j=$i+1; $j<scalar(@set); $j++) {
	    push @pairs,"$set[$i] $set[$j]";
	}
    }
    
    my $ratioorig=$freqorig{$set}/$total;
    my $ratioswapped=$freqswapped{$set}/$total;
    
    for my $item (@set) {
	$ratioorig/=$freqorig{$item}/$total;
	$ratioswapped/=$freqorig{$item}/$total;
    }
    
    my $avgcorrorig=0;
    my $avgcorrswapped=0;
    
    for my $pair (@pairs) {
	my ($i,$j)=split " ", $pair;
	$avgcorrorig+=($total*$freqorig{$pair}-$freqorig{$i}*$freqorig{$j})/sqrt(($total*$freqorig{$i}-$freqorig{$i}*$freqorig{$i})*($total*$freqorig{$j}-$freqorig{$j}*$freqorig{$j}));
	
	$avgcorrswapped+=($total*$freqswapped{$pair}-$freqswapped{$i}*$freqswapped{$j})/sqrt(($total*$freqswapped{$i}-$freqswapped{$i}*$freqswapped{$i})*($total*$freqswapped{$j}-$freqswapped{$j}*$freqswapped{$j}));
    }
    if(scalar(@pairs)>0) {
	$avgcorrorig/=scalar(@pairs);
	$avgcorrswapped/=scalar(@pairs);
    }

    print $freqorig{$set}, " ", $freqswapped{$set}, " ",
    ($freqorig{$set}>0 ? ($freqorig{$set}-$freqswapped{$set})/$freqorig{$set} : 'Inf'), " ",
    ($freqswapped{$set}>0 ? ($freqswapped{$set}-$freqorig{$set})/$freqswapped{$set} : 'Inf'), " ",
    $avgcorrorig, " ", $avgcorrswapped, " ",
    ($freqorig{$set}>0 ? $freqswapped{$set}/$freqorig{$set} : 'Inf'), " ",
    ($freqswapped{$set}>0 ? $freqorig{$set}/$freqswapped{$set} : 'Inf'), " ",
    $ratioorig, " ", $ratioswapped, " ", 
    
    "; $set\n";
}
#print STDERR "The columns were: freqorig freqswapped relerrorig relerrswapped avgcorrorig avgcorrswapped freqswapped/freqorig freqorig/freqswapped ratioorig ratioswapped ; itemset\n";
