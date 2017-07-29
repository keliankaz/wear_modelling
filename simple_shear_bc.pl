#! /usr/bin/perl
use POSIX;
# This script does the following tasks:
#    1. writes a fric2d input file that has displacement boundary conditions that simulate uniform stress field in  
#	rectangle undergoing simple shear and uniform shortening.
#
#	usage: simple_shear_bc.pl height length elesize  
#
#	length = the length of the model in meters 
#       size = size of the element in meters
#   
#   2017-04-20 adapted by Kelian Dascher-Cousineau:
#   added flexibility  for input (in order to better integrate with matlab
#   workflow)
#  
#   NEW USAGE: 
#   simple_shear_pc.pl heght length elesize yougn poisson sigxx sigyy shear


($height, $length, $size, $E, $nu, $sigxx, $sigyy, $shear) = @ARGV;

# here are a bunch of parameters that we are not changing at the moment
# For test we use  0.06 0.15 0.0006/0.0002



# commented out - kelian
# $E = 30000;    #Young's Modulus MPa
# $nu = 0.25;    #Poisson's ratio
# $sigxx = 100;  #normal stress parallel to faults MPa
# $sigyy = 100;  #normal stress perpendicular to faults MPa
# $shear = 4000000;   #shear stress on faults  MPa

$name = 'input';

$elenum=0;  #counts the number of elements

open(BOUND, ">boundaries.txt");

#step through the four sides of the model starting with the top

#TOP OF THE MODEL
$num = floor($length/$size+0.5);
$elenum += $num;
for ($i=0; $i<$num; $i++){
	$xbeg = $i*$size;
	$xend = $xbeg + $size;
	$ybeg = $height/2;
        $yend = $height/2;

	#The shear displacement includes the component due to simple shear and the component
	#due to lateral shortening
	$us = -($xbeg+$size/2)/$E*((1-$nu*$nu)*$sigxx-$nu*(1+$nu)*$sigyy) + $shear*$height*2*(1+$nu)/$E;
	#The normal displacement is due to vertical shortening
	$un = -$height/$E*( (1-$nu*$nu)*$sigyy-$nu*(1+$nu)*$sigxx);	
	print BOUND "1	$xbeg	$ybeg	$xend	$yend	2	0	0	$us	$un\n";	
}#for  


#RIGHT SIDE OF THE MODEL
$num = floor($height/$size+0.5);
$elenum += $num;
for ($i=0; $i<$num; $i++){
        $xbeg = $length;
        $xend = $length;
        $ybeg = $height/2-$i*$size;
        $yend = $ybeg-$size;

        #The shear displacement includes the component due to vertical shortening                 
        $us = ($ybeg-$size/2+$height/2)/$E*( (1-$nu*$nu)*$sigyy-$nu*(1+$nu)*$sigxx);
        #The normal displacement includes the component due to simple shear and the component
        #due to lateral shortening
        $un = -$length/$E*((1-$nu*$nu)*$sigxx-$nu*(1+$nu)*$sigyy) + $shear*($ybeg-$size/2+$height/2)*2*(1+$nu)/$E;
        print BOUND "1	$xbeg	$ybeg	$xend	$yend	2	0	0	$us	$un\n";         
}#for


#BOTTOM OF THE MODEL
$num = floor($length/$size+0.5);
$elenum += $num;
for ($i=0; $i<$num; $i++){
        $xbeg = $length - $i*$size;
        $xend = $xbeg - $size;
        $ybeg = -$height/2;
        $yend = -$height/2;

        #The shear displacement includes the component due to lateral shortening
        $us = ($xbeg-$size/2)/$E*((1-$nu*$nu)*$sigxx-$nu*(1+$nu)*$sigyy);
        #The normal displacement is zero 
        $un = 0;
        print BOUND "1	$xbeg	$ybeg	$xend	$yend	2	0	0	$us	$un\n";         
}#for


#LEFT SIDE OF THE MODEL
$num = floor($height/$size+0.5);
$elenum += $num;
for ($i=0; $i<$num; $i++){
        $xbeg = 0;
        $xend = 0;
        $ybeg = -$height/2+$i*$size;
        $yend = $ybeg+$size;

        #The shear displacement includes the component due to vertical shortening                          
        $us = -($ybeg+$size/2+$height/2)/$E*( (1-$nu*$nu)*$sigyy-$nu*(1+$nu)*$sigxx);
        #The normal displacement includes the component due to simple shear 
	$un = -$shear*($ybeg+$size/2+$height/2)*2*(1+$nu)/$E;
        print BOUND "1	$xbeg	$ybeg	$xend	$yend	2	0	0	$us	$un\n";
}#for

close(BOUND);

#------------ 
#PRINT THE INPUT FILE INFORMATION 
#------------ 

#Consider releasing steps of 5km underlap, zero underlap and 5km overlap
	
		open(FILE, ">$name");
		open(TMP, ">tmp.txt"); 

		print TMP "*------------------------                           \n";                                                                    
		print TMP "*Problem        Variables                           \n";                                                            
		print TMP "*------------------------                           \n\n";                                                                    
		print TMP "*Titles   \n";                                                                                      
		print TMP "title1  =       \"Releasing stepover in sandstone\" \n";
		print TMP "title2  =       \" friction = 0; stepover = $j; overlap = $overlap\"             \n\n";                                                           
		print TMP "*Rock   Properties                                                                           \n";          
		print TMP "pr    =     $nu                  *(possion's     ratio)                                            \n";              
		print TMP "e     =     $E                   *(young's       modulus MPa)                                                    \n";
		print TMP "k1c   =     1.0                  *(mode  I       fracture toughness,     MPa sqrt(m), from rel.in) \n";                                      
		print TMP "tensileStrength = 5e8    *(tensile strength,     MPa, from rel.in) \n\n";                                      
		print TMP "*Symmetry            \n";                                                                                   
		print TMP "ksym    =       1        \n\n";                                                                       
		print TMP "*Gravitational  stresses     \n";                                                                                   
		print TMP "gravity =       0           *(apply gravity:        0=no,   1=yes)  \n";                                        
		print TMP "kratio  =       1           *(ratio of      horizontal:vertical     gravitation     stress) \n";                                
		print TMP "density =       2500    *(density       of      rock    kg/m^3)           \n\n";                              
		print TMP "*Number of      Observation     Lines   and     Boundary        Lines \n";                                          
		print TMP "nolines =       0           *(number        of      observation     lines,  must    mach    #       lines   designated      below)\n";
		print TMP "nblines =       $elenum     *(number        of      boundary        lines,  must    mach    #       lines   designated      below)\n\n";
		print TMP "*Options        Commands                                                                             \n";    
		print TMP "nsteps  =       1       *(number        of      loading steps)                                          \n";
		print TMP "niters  =       1       *(number        of      opening-mode    fracture        propagation     iterations)     \n";                        
		print TMP "pdfiter =       1       *(print data    from    fracture        prop    iteration)                          \n";    
		print TMP "bdata   =       1       *(print data    along   boundaries:     0=no    1=yes)                          \n";
		print TMP "fdata   =       1       *(print data    along   (new    fractures:      0=no    1=yes)                  \n";
		print TMP "faultData       =       1       *(print data    along   faults: 0=no    1=yes)                          \n";
		print TMP "odata   =       0          *(print observation     line    data:   0=no    1=yes)                          \n";
		print TMP "tolerance       =       1.00E-06     *(test  for     frictinal       slip    convergence)                                    \n";
		print TMP "tolCounter      =       100          *(by-pass       tolerance       test    and     exit    convergence     at      this    iteration)\n";      
		print TMP "end     *(terminate     list)      \n\n";
		print TMP "*Observation    Lines              \n";                                                                     
		print TMP "*------------------                \n";                                                                             
		print TMP "*xbeg   ybeg    xend    yend    numpb\n\n";                                                           
		print TMP "*------------------   \n";                                                                                          
		print TMP "*Boundary       Lines  \n";                                                                                 
		print TMP "*------------------    \n";                                                                                         
		print TMP "*(bvs   =       shear;  bvn     =       normal) STATIC  MONOTONIC \n";                                      
		print TMP "*num    xbeg    ybeg    xend    yend    knode   bvs     bvn     bvs     bvn \n";                    
		print TMP "*---    ----    ----    ----    ----    -----   ---     ---     ---     --- \n";

		system "cat tmp.txt boundaries.txt > $name";
		close(FILE);
		open(FILE, ">>$name");	


                print FILE "\n*---------------------- \n";                                                                                        
                print FILE "*Fault  Conditions         \n";                                                                             
                print FILE "*----------------------         \n";                                                                                
                print FILE "*fault  name    grow_tails?     from_end1?      from_end2?  \n";                                                             
                print FILE "*num    xbeg    ybeg    xend    yend    stiffS          stiffN          ten-str init-coh slid-coh stat-fric dy-fric crit-slip \n";
                print FILE "*----   ----    ----    ----    ----    ------          ------          ---     ---     -----     -------   ------  --------- \n";
                print FILE "fault   rough    yes     no      no                                                 \n";             




	
		close(FILE);

system "rm tmp.txt boundaries.txt";

