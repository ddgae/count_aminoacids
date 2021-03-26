#writing a function using TCL
#usage: vmd -pdb file.pdb -psf file.psf -dcd file.dcd
#usage: open tcl prompt and source momentofI.tcl
# david gae

#using vmd's measure rgyr
puts "source to get the answer: source momentofI.tcl"
puts "radius of gyration per frame:  gyr_radius top"
proc gyr_radius top {
set sel [atomselect top all]
set num_steps [ molinfo top get numframes]
for {set i 0} {$i < $num_steps} {incr i} {
    $sel frame $i
    set ans [ measure rgyr $sel]
    puts $ans
      }
}


# dalke et al. Using Tcl for molecular visualization and analysis, Pac symp Biocomput. 1997: 85-96
set I 0
set ans 0
set sel1 [ atomselect top all]
set com_ref [ measure center $sel1 weight mass ]
foreach m [$sel1 get mass ] coord [ $sel1 get { x y z } ] {
	# radius of gyration sqrt{ M* ri - center(r1)^2) /N }
	set I [ expr $I + $m * [veclength2 [vecsub $coord $com_ref ]]]
	set ans [ expr $ans + $m ]
}
set r1 [expr sqrt ($I/$ans)]
puts "moment of inertia per x,y,z position: $r1"

# angular inertia
set I 0
set ans 0
set v 0
set outfile [ open "data.csv" w]

set sel2 [ atomselect top all]
set com_ref [ measure center $sel2 weight mass ]
foreach m [$sel2 get mass ] coord [ $sel2 get { x y z } ] {
	# radius of gyration [M * (r - com_ref)^2]
	#http://hyperphysics.phy-astr.gsu.edu/hbase/mi.html
	set I [ expr $I + $m * [veclength2 [vecsub $coord $com_ref ]]]
	#puts "$I"
	# L = Iw (  L = mrv, v = wr, v = 2*pi * r )
    set r1 [veclength2 [vecsub $coord $com_ref]]
    #puts "$r1"
    set r2 [ veclength [ vecsub $coord $com_ref]]
    #puts "$r2"
	set v [ expr  2 * 3.14 * $r1 ]
	#puts $v
	set L [ expr  $I + $m * ($v / $r1) ]
	 #puts $L
	set ans [ expr $ans + $m]
	#symmetry axis
	set r3 [ expr sqrt($I/$ans)]
	#solid cylinder
	set r4 [ expr 0.5 * sqrt($I/$ans)]
	#solid sphere
	set r5 [ expr 0.4 * sqrt($I/$ans)]
	#rod about center
	set r6 [ expr 0.08 * sqrt($I/$ans)]
	#solid cylinder, center diameter
	set r7 [ expr 0.4 * sqrt($I/$ans)+ 0.08 * sqrt($L/$ans)]
	set r8 [ expr 0.67 * sqrt($I/$ans) ]
	set r9 [ expr 0.3 * sqrt($L/$ans) ]
}
puts  "moment of inertia, $r3"
puts  "solid cylinder, $r4"
puts  "solid sphere, $r5"
puts  "rod about center, $r6"
puts  "solid cylinder, $r7"
puts  "thin shell, $r8"
puts  "Rod about end, $r9"

puts $outfile "moment of inertia, $r3"
puts $outfile "solid cylinder, $r4"
puts $outfile "solid sphere, $r5"
puts $outfile "rod about center, $r6"
puts $outfile "solid cylinder, $r7"
puts $outfile "thin shell, $r8"
puts $outfile "Rod about end, $r9"
close $outfile

