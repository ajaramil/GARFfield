#!/usr/bin/perl


$fname = $ARGV[0];
$out = $ARGV[1];

open (IN,"<$fname");
open (OUT,">$out");

while (<IN>) {

	if ( /^\s{1,}\d/ ) {
		my @tks = split();
		print OUT "$tks[12]\n" if defined($tks[12]);
		}

}

close IN;
close OUT;
