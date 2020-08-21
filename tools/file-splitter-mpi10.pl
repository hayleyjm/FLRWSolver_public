#!/usr/bin/perl
my $it = -1;
print "IT is $it\n";
my $count = 0;
my $blocks = 0;
my @block = ();

print "For every iteration, enter 'e'. For every hundreth iteration, enter 'h'. For every thousandth iteration, enter 't'. For every ten thousandth iteration, enter '10'.";
my $input = <STDIN>;
my $time;
## splitter for 10 nodes

while (<>) {
#    my @block = ();
    my $number = 0;
    my $ittemp;
    my $interval;
    if ($input =~ m/e/) {
	$interval = 1;
    } elsif ($input =~ m/h/) {
	$interval = 100;
    } elsif ($input =~ m/t/) {
	$interval = 1000;
    }elsif ($input =~ m/10/) {
	$interval = 10000;
    }

    if (m/\# iteration/) {
	($ittemp) = m/iteration (\d+)/;
	($time)= m/time (\d+\.\d*)/;
#	print "TIME is $time\n";
	if ($ittemp==$it) {
	    # do nothing
	} else {
	    $it = $ittemp;
        }
    }
    
    while ( m/^\d.*/) {
	$number = 1;
	push @block,$_;
	$_ = <>;
    }

    if ($number==1) {
	$blocks++;
#	print "BLOCKS is $blocks\n";
    }

    if ($blocks==2) {
	$count++;
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


	if (($count%$interval)==0 or $count==1) {
	    print "ITERATION is $it, TIME is $time, writing to $filename\n";
	    open (OUT, ">", "$filename") or die "cannot open $filename $!";
	    print OUT "$time\n";
	    print OUT @block;
	    close (OUT)
	}	    
	$number = 0;
	$blocks = 0;
	@block = ();
    }

}
