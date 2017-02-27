# hbonds
The code can count the hbond based on geometric criteria (customizable),
Meantime generating POSCARs for VASP runs with both full water solvation
and partial preserved solvation(non-hbonded removed). 

02/27/17 Update: Added support for TI4P water model support, also support CxHyOz system on Pt, xyz>0.

06/27 Update:   Added the function to calculate the dipole moment for each hydrogen bond found.

04/05 Update:   Now the code works fast enough so that you don't have to wait minutes for a single frame hbond searching.
                Default criteria, AD<3.5 ADH<30 AHD>120.

03/10/16 Update: Output improment, now tells the donar oxygen(s) and the acceptor oxgen(s).

                 
12/12 Update:	Improvement in calculation speed by cleaning out some legacy float-string conversions.
				Now can recenter the system to align the bottom Pt to the bottom of box.
				Formatting improvement.

12/02 Update:	Reconfigured PBC, now works with all 8 neighboring cells.
				Known issue: computational time doubled.

11/12 Update:	Now compatible with C3HxOx adsorbates.

10/22 Update: 	Reconfigured outputting algorithm, now can export with POSCARS
				partial fixed flags on for NEB-dimer/geometry optimization
				or,
				simply all fixed coordinates for single point energy calc.

*****************************IMPORTANT***********************************

NOTE: for **************   COUNTING   ************* lammps frames
 Turned off manual atom type input for speed
   (Can turn it back on when needed).

   Default settings:

rC = [2,3,4,5]         
rHS = [10,11]
rOS = [6,7]  
rHW = [9]
rOW = [8]  

fixedposcar = 1

Future improvement,
1. Update neighbor list to fasten hbonding counting.