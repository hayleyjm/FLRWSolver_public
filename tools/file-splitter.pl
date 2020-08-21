#!/usr/bin/perl
open(OUT, "> /dev/null");
my $count = 0;
my $write = 0;
print "For every iteration, enter 'e'. For every hundreth iteration, enter 'h'. For every thousandth iteration, enter 't'.";
my $input = <STDIN>;
while (<>) {
    if ($input =~ m/e/) {
	if ( m/\# iteration/ ) {
	    close OUT;
	    $count++;
	    ($it) = m/iteration (\d+)/;
	    my ($time)= m/time (\d+\.\d*)/;
	    my $densitystring = "rho";
	    my $hamstring = "hamiltonian";
	    my $filename;
	    if (index($ARGV, $densitystring) != -1) {
		$filename = sprintf("density_%09u.dat", $it);
	    } elsif (index($ARGV, $hamstring) != -1) {
		$filename = sprintf("hamiltonian_%09u.dat", $it);
	    }
	    print "ITERATION IS $it, writing to $filename\n";
	    open (OUT, ">", "$filename") or die "cannot open $filename $!";
	    print OUT "$time\n";
	    $write = 1;
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

    elsif ($input =~ m/h/) {
	if ( m/\# iteration/ ) {
            close OUT;
            $count++;
            ($it) = m/iteration (\d+)/;
            my ($time)= m/time (\d+\.\d*)/;
            my $densitystring = "rho";
            my $hamstring = "hamiltonian";
            my $filename;
            if (index($ARGV, $densitystring) != -1) {
                $filename = sprintf("density_%09u.dat", $it);
            } elsif (index($ARGV, $hamstring) != -1) {
                $filename = sprintf("hamiltonian_%09u.dat", $it);
            }
	    if (($count%100)==0) {
		print "ITERATION IS $it, writing to $filename\n";
		open (OUT, ">", "$filename") or die "cannot open $filename $!";
		print OUT "$time\n";
		$write = 1;
	    } else {
		$write = 0;
	    }
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

    elsif ($input =~ m/t/) {
        if ( m/\# iteration/ ) {
            close OUT;
            $count++;
            ($it) = m/iteration (\d+)/;
            my ($time)= m/time (\d+\.\d*)/;
            my $densitystring = "rho";
            my $hamstring = "hamiltonian";
            my $filename;
            if (index($ARGV, $densitystring) != -1) {
                $filename = sprintf("density_%09u.dat", $it);
            } elsif (index($ARGV, $hamstring) != -1) {
                $filename = sprintf("hamiltonian_%09u.dat", $it);
            }
	    if (($count%1000)==0) {
                print "ITERATION IS $it, writing to $filename\n";
		open (OUT, ">", "$filename") or die "cannot open $filename $!";
                print OUT "$time\n";
                $write = 1;
            } else {
                $write = 0;
            }
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
}
