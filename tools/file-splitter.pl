#!/usr/bin/perl
open(OUT, "> /dev/null");
my $count = 0;
my $write = 0;
print "For every iteration, enter 'e'. For every hundreth iteration, enter 'h'. For every thousandth iteration, enter 't'.";
my $input = <STDIN>;
print "Enter number of processors used:";
my $np = <STDIN>;
my $itcount = 0;
while (<>) {
    my $number;
    if ($input =~ m/e/) {
	$number = 1;
    } elsif ($input =~ m/h/) {
	$number = 100;
    } elsif ($input =~ m/t/) {
	$number = 1000;
    }
    if ( m/\# iteration/ and ($itcount%$np==0)) {
#    if ( m/\# iteration/) {
	close OUT;
	$count++;
	($it) = m/iteration (\d+)/;
	my ($time)= m/time (\d+\.\d*)/;
	my $filename;
	if ($ARGV =~ m/rho/) {
	    $filename = sprintf("density_%09u.dat", $it);
	} elsif ($ARGV =~ m/hamiltonian/) {
	    $filename = sprintf("hamiltonian_%09u.dat", $it);
	} elsif ($ARGV =~ m/metric/) {
            $filename = sprintf("metric_%09u.dat", $it);
	} elsif($ARGV =~ m/curv/) {
            $filename = sprintf("curv_%09u.dat", $it);
	} elsif($ARGV =~ m/momentum/) {
	    $filename = sprintf("momentum_%09u.dat",$it);
	} elsif($ARGV =~ m/vel/) {
	    $filename = sprintf("velocity_%09u.dat",$it);
        } else {
	    $filename = sprintf("iteration_%09u.dat", $it);
	}
	if (($count%$number)==0 or $count==1) {
	    print "ITERATION IS $it, writing to $filename\n";
	    open (OUT, ">", "$filename") or die "cannot open $filename $!";
	    print OUT "$time\n";
	    $write = 1;
	} else {
	    $write = 0;
	}
    }
    if ( m/\# iteration/ ) {
	    $itcount++;
	}
    
    if ($write==1) {
	$line=$_;
	if ($line =~ m/^\s*\d+/) {
	    my @column = split ( /\s+/, $line );
	    chomp($line);
	    print OUT "$line\n";
	} else {
	    print OUT "$line";
	}
    }
}
