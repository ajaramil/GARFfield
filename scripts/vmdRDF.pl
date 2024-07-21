#!/usr/bin/perl

@name = qw/TOT A MB1 EB MA Me OH MB2/;

open(TCL, ">rdf.tcl");

$trj = $ARGV[0];

print TCL "package require pbctools
mol new $trj
pbc set {69.8428 69.8428 69.8428} -all\n";

print TCL "set sel1 [atomselect top \"all\"]
set foo [measure gofr \$sel1 \$sel1 delta rmax 10 usepbc TRUE]
set bar [open temp00.txt w]
puts \$bar \$foo
exec [./parseRDF.pl temp00.txt ../test/cg/TOT]\n";

for ($i = 1; $i <= 7; $i++){
   for ($j = $i; $j <= 7; $j++){
      $outfile = "$name[$i]$name[$j]";
      $tempfile = "temp$i$j.txt";
      print TCL "set sel1 [atomselect top \"type $i\"]
set sel2 [atomselect top \"type $j\"]
set foo [measure gofr \$sel1 \$sel2 rmax 10 usepbc TRUE] 
set bar [open $tempfile w]
puts \$bar \$foo
exec [./parseRDF.pl $tempfile $outfile]\n";
   }
}

close TCL;

#system("/project/DOW-latex/bin/vmd -dispdev text < rdf.tcl > temp.txt");
system("/Applications/VMD_1.9.1.app/Contents/vmd/vmd_MACOSXX86 -dispdev text < rdf.tcl > temp.txt");
#system ("./vmd  -dispdev text < rdf.tcl > temp.txt");
#system("/project/DOW-latex/optimize/vmd/bin/vmd -dispdev text < rdf.tcl > temp.txt");

