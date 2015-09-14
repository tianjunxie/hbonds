# hbonds
Based on the honds currently have,
want to develop the code also capable of
generating POSCAR with both explicit water solvation
and implicit water solvation. 

NOTE: for **************   EXPORTING   ************* POSCAR files
*****************************IMPORTANT***********************************

For now the code can EXPORT C3HxO3 related POSCARS.
   i.e. No correct number of C or O atoms will show up
        in the POSCAR files if the system is not for C3HxO3.

NOTE: for **************   COUNTING   ************* lammps frames
 Turned off manual atom type input for speed
   (Can turn it back on when needed).

   Default settings:
**********************************rHS = ['8','9']************************   
**********************************rOS = ['4','5']************************   
**********************************rHW = ['7']    ************************
**********************************rOW = ['6']    ************************

Future improvement,
1. add a switch for PBC
2. add C O variation capabablity(outputting POSCAR)
3. for (1, 2), needed to improve the neighbor list
