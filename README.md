# hbonds
Based on the honds currently have,
want to develop the code also capable of
generating POSCAR with both explicit water solvation
and implicit water solvation. 

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

NOTE: for **************   EXPORTING   ************* POSCAR files
*****************************IMPORTANT***********************************

For now the code can EXPORT C3HxOx related POSCARS.
   i.e. No correct number of C or O atoms will show up
        in the POSCAR files IF the system is not for C3HxO3.

NOTE: for **************   COUNTING   ************* lammps frames
 Turned off manual atom type input for speed
   (Can turn it back on when needed).

   Default settings:
rHS = ['10','11']
rOS = ['6']
rHW = ['9']
rOW = ['8']

Future improvement,
1. Update neighbor list to fasten hbonding counting.
2. Reduce input by reading the data file for automated identification of atom types.