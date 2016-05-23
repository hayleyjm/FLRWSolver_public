#!/usr/bin/perl
my $count = 0;
while (<>) {
    my @block = ();
    my $number = 0;
    while ( m/^\d.*/) {
	$number = 1;
	push @block,$_;
        $_ = <>;
    }

    if ($number==1) {
	$count++;
#	print "Hello\n";
#	print @block;

	$filename = sprintf("iteration_%09u.dat", $count);
	print "ITERATION IS $count, writing to $filename\n";
	open (OUT, ">", "$filename") or die "cannot open $filename $!";
	print OUT @block;
	$number = 0;
	close (OUT)
    }
    
}
