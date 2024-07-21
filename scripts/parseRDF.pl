#! /usr/bin/perl

print "Parsing RDFs ...\n";

$infile = $ARGV[0];
$outfile = $ARGV[1];

open(IN,"<$infile");
open(OUT,">$outfile");

print "-- parseRDF opened $outfile\n";

while(<IN>){
    if (/{(.*?)} {(.*?)}.*/){
    @r = split(' ',$1);
    @g = split(' ',$2);
    for $i (1 .. $#r){
       print OUT "$r[$i]\t$g[$i]\n";
    }
  } 
}

#system "rm temp.txt";
