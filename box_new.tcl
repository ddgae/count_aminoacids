#based on Analysis scripts tutorial
#https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node204.html
#select only the oxygen from a WAT file, since hydrogen can built with reduce, psfgen, etc. 
set sel1 [atomselect top "name O"]
#determine the center of the box. 
#set sel2 [measure center $sel1]

#function collect x,y,z coordinates to create a waterbox.  
proc box { sel1 } {
	#save the output filename into x,y,z coordinates 
	set output1 [ open box1.xyz "w"]
	#get x,y,z from a file
	set coords [ $sel1 get { x y z } ]
	#look at the set sel2
	set coords2 [ list -0.14582006633281708 0.1844363510608673 0.018788645043969154]
	#lvarpop identifies an element in array and assigns it. Can remove by giving a zero. 	
	set coord [ lvarpop coords ]
	lassign $coord minx miny minz
	lassign $coord maxx maxy maxz
	lassign $coords2 x1 y1 z1
	#loop through $coords of x,y,z coordinates.	
	foreach coord $coords {
		# assign a new array with current x,y,z		
		lassign $coord   x y z 
		lassign $coords2 x1 y1 z1
		#select only the sorted and unique x,y,z
		set X1 [ lsort -unique "$x"]
		set Y1 [ lsort -unique "$y"]
		set Z1 [ lsort -unique "$z"]
		#center at origin
		set X2 [ lsort -unique "$x1"]
		set Y2 [ lsort -unique "$y1"]
		set Z2 [ lsort -unique "$z1"]
		# [x, y, z ----|center.
		#corner 1
		if { $X1 > $X2 }  { set d1 [format %.4f [expr $X1]]}
		if { $Y1 > $Y2 }  { set d2 [format %.4f [expr $Y1]]}
		if { $Z1 > $Z2 }  { set d3 [format %.4f [expr $Z1]]} 
		#corner 2
		if { -$X1 < $X2 } { set d4 [format %.4f [expr -$X1]]}
		if { -$Y1 < $Y2 } { set d5 [format %.4f [expr -$Y1]]}
		if { -$Z1 < $Z2 } { set d6 [format %.4f [expr -$Z1]]}
		#corner 3
		if { -$X1 < $X2 } { set d7 [format %.4f [expr -$X1]]}
		if { -$Y1 < $Y2 } { set d8 [format %.4f [expr -$Y1]]}
		if {  $Z1 > $Z2 } { set d9 [format %.4f [expr  $Z1]]}
		#corner 4
		if {  $X1 > $X2 } { set d10 [format %.4f [expr  $X1]]}
		if { -$Y1 < $Y2 } { set d11 [format %.4f [expr -$Y1]]}
		if { -$Z1 < $Z2 } { set d12 [format %.4f [expr -$Z1]]}
		#corner 5
		if {  $X1 > $X2 } { set d13 [format %.4f [expr  $X1]]}
		if { -$Y1 < $Y2 } { set d14 [format %.4f [expr -$Y1]]}
		if {  $Z1 > $Z2 } { set d15 [format %.4f [expr  $Z1]]}
		#corner 6
		if { -$X1 < $X2 } { set d19  [format %.4f [expr -$X1]]}
		if {  $Y1 > $Y2 } { set d20  [format %.4f [expr  $Y1]]}
		if { -$Z1 < $Z2 } { set d21  [format %.4f [expr -$Z1]]}
		#corner 7
		if { -$X1 < $X2 } { set d22  [format %.4f [expr -$X1]]}   
		if {  $Y1 > $Y2 } { set d23  [format %.4f [expr  $Y1]]}
		if {  $Z1 > $Z2 } { set d24  [format %.4f [expr  $Z1]]}
		#corner 8
		if { $X1 > $X2 }  { set d25   [format %.4f [expr  $X1]]}   
		if { $Y1 > $Y2 }  { set d26   [format %.4f [expr  $Y1]]}
		if {-$Z1 < $Z2 }  { set d27   [format %.4f [expr -$Z1]]}
		#check results
		puts  "$d1      $d2       $d3"
		#minbox
		#set outside1  [lsort -unique [lindex "O         $d1        $d2        $d3"]]
		#set outside2  [lsort -unique [lindex "O         $d4        $d5        $d6"]]
		#set outside4  [lsort -unique [lindex "O         $d10       $d11       $d12"]]
		#set outside5  [lsort -unique [lindex "O         $d13       $d14       $d15"]]
		set outside1  [lindex "O         $d1        $d2        $d3"]
		set outside2  [lindex "O         $d4        $d5        $d6"]
		set outside3  [lindex "O         $d7        $d8        $d9"]
		set outside4  [lindex "O         $d10       $d11       $d12"]
		set outside5  [lindex "O         $d13 	    $d14       $d15"]
		set outside7  [lindex "O 	 $d19	    $d20       $d21"]
		set outside8  [lindex "O 	 $d22	    $d23       $d24"]
		set outside9  [lindex "O	 $d25 	    $d26       $d27"]
                #output to file
		puts  $output1 "$outside1"
		puts  $output1 "$outside2"
		puts  $output1 "$outside3"
		puts  $output1 "$outside4"
		puts  $output1 "$outside5"
		puts  $output1 "$outside7"
		puts  $output1 "$outside8"
		puts  $output1 "$outside9"
		
	}
	close $output1 
}
