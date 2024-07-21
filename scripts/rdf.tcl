package require pbctools
mol new CG.nvt.lammpstrj
pbc set {69.8428 69.8428 69.8428} -all
set sel1 [atomselect top "all"]
set foo [measure gofr $sel1 $sel1 delta rmax 10 usepbc TRUE]
set bar [open temp00.txt w]
puts $bar $foo
exec [./parseRDF.pl temp00.txt TOT]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 1"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp11.txt w]
puts $bar $foo
exec [./parseRDF.pl temp11.txt AA]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 2"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp12.txt w]
puts $bar $foo
exec [./parseRDF.pl temp12.txt AMB1]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 3"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp13.txt w]
puts $bar $foo
exec [./parseRDF.pl temp13.txt AEB]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 4"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp14.txt w]
puts $bar $foo
exec [./parseRDF.pl temp14.txt AMA]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 5"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp15.txt w]
puts $bar $foo
exec [./parseRDF.pl temp15.txt AMe]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp16.txt w]
puts $bar $foo
exec [./parseRDF.pl temp16.txt AOH]
set sel1 [atomselect top "type 1"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp17.txt w]
puts $bar $foo
exec [./parseRDF.pl temp17.txt AMB2]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 2"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp22.txt w]
puts $bar $foo
exec [./parseRDF.pl temp22.txt MB1MB1]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 3"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp23.txt w]
puts $bar $foo
exec [./parseRDF.pl temp23.txt MB1EB]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 4"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp24.txt w]
puts $bar $foo
exec [./parseRDF.pl temp24.txt MB1MA]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 5"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp25.txt w]
puts $bar $foo
exec [./parseRDF.pl temp25.txt MB1Me]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp26.txt w]
puts $bar $foo
exec [./parseRDF.pl temp26.txt MB1OH]
set sel1 [atomselect top "type 2"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp27.txt w]
puts $bar $foo
exec [./parseRDF.pl temp27.txt MB1MB2]
set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 3"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp33.txt w]
puts $bar $foo
exec [./parseRDF.pl temp33.txt EBEB]
set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 4"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp34.txt w]
puts $bar $foo
exec [./parseRDF.pl temp34.txt EBMA]
set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 5"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp35.txt w]
puts $bar $foo
exec [./parseRDF.pl temp35.txt EBMe]
set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp36.txt w]
puts $bar $foo
exec [./parseRDF.pl temp36.txt EBOH]
set sel1 [atomselect top "type 3"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp37.txt w]
puts $bar $foo
exec [./parseRDF.pl temp37.txt EBMB2]
set sel1 [atomselect top "type 4"]
set sel2 [atomselect top "type 4"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp44.txt w]
puts $bar $foo
exec [./parseRDF.pl temp44.txt MAMA]
set sel1 [atomselect top "type 4"]
set sel2 [atomselect top "type 5"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp45.txt w]
puts $bar $foo
exec [./parseRDF.pl temp45.txt MAMe]
set sel1 [atomselect top "type 4"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp46.txt w]
puts $bar $foo
exec [./parseRDF.pl temp46.txt MAOH]
set sel1 [atomselect top "type 4"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp47.txt w]
puts $bar $foo
exec [./parseRDF.pl temp47.txt MAMB2]
set sel1 [atomselect top "type 5"]
set sel2 [atomselect top "type 5"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp55.txt w]
puts $bar $foo
exec [./parseRDF.pl temp55.txt MeMe]
set sel1 [atomselect top "type 5"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp56.txt w]
puts $bar $foo
exec [./parseRDF.pl temp56.txt MeOH]
set sel1 [atomselect top "type 5"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp57.txt w]
puts $bar $foo
exec [./parseRDF.pl temp57.txt MeMB2]
set sel1 [atomselect top "type 6"]
set sel2 [atomselect top "type 6"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp66.txt w]
puts $bar $foo
exec [./parseRDF.pl temp66.txt OHOH]
set sel1 [atomselect top "type 6"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp67.txt w]
puts $bar $foo
exec [./parseRDF.pl temp67.txt OHMB2]
set sel1 [atomselect top "type 7"]
set sel2 [atomselect top "type 7"]
set foo [measure gofr $sel1 $sel2 rmax 10 usepbc TRUE] 
set bar [open temp77.txt w]
puts $bar $foo
exec [./parseRDF.pl temp77.txt MB2MB2]
