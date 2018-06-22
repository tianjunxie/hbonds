#!/usr/bin/env python
import string,sys
import numpy as np
from math import *
import os
import timeit
# from decimal import Decimal
# from array import *
temp = lwat = lads = output = xp = xn = yp = yn = coords = []
h1 = h2 = o1 = o2 =  link = [None] * 10000
A  = 0
flag = nlines = 0

def dist(X,Y):
    x1 = X[1]
    x2 = Y[1]
    y1 = X[2]
    y2 = Y[2]
    z1 = X[3]
    z2 = Y[3]
    return 2*(x1-x2)*(y1-y2)*a*b*cosgamma + 2*(y1-y2)*(z1-z2)*b*c*cosalpha + 2*(x1-x2)*(z1-z2)*a*c*cosbeta + (x1-x2)*(x1-x2)*a*a + (y1-y2)*(y1-y2)*b*b + (z1-z2)*(z1-z2)*c*c


def init():
    x = [None] * 50000
    file = open(sys.argv[1],'r')
    for i in range(9):
        file.readline() # skip first 9 lines
    for i in range(N): # read in x, y, z coordinates 
        line = file.readline()
        words = [float(n) for n in line.split()]
        x.append(words)
    file.close()
    x = [a for a in x if a != None]
    dz = x[0][4]
    for i in range(len(x)):
        x[i][0] = int(x[i][0])
        x[i][1] = int(x[i][1])      
        # if x[i][4] < dz:
        #     x[i][4] = 1-(x[i][4]+dz)
        # else:
        #     x[i][4] = x[i][4] - dz
        x[i][4] = x[i][4] - dz    
    return (x)

def ang(d1,d2,d3):# Note: square of the distances
    return degrees( acos((d1 + d2 - d3)/(2.0 * sqrt(d1) * sqrt(d2))) )

def find(y,c):
    x = grab= grabxp = grabxn = grabyp = grabyn = grabxpyp = grabxnyn = grabxnyp = grabxpyn = [None] * 50000
    for i in range(len(c)):
        grab = ((c[i][1],c[i][2],c[i][3],c[i][4],c[i][0]))
        grabxp = ((c[i][1],c[i][2]+1,c[i][3],c[i][4],c[i][0]))
        grabxn = ((c[i][1],c[i][2]-1,c[i][3],c[i][4],c[i][0]))
        grabyp = ((c[i][1],c[i][2],c[i][3]+1,c[i][4],c[i][0]))
        grabyn = ((c[i][1],c[i][2],c[i][3]-1,c[i][4],c[i][0]))
        grabxpyp = ((c[i][1],c[i][2]+1,c[i][3]+1,c[i][4],c[i][0]))
        grabxpyn = ((c[i][1],c[i][2]+1,c[i][3]-1,c[i][4],c[i][0]))
        grabxnyp = ((c[i][1],c[i][2]-1,c[i][3]+1,c[i][4],c[i][0]))
        grabxnyn = ((c[i][1],c[i][2]-1,c[i][3]-1,c[i][4],c[i][0]))
        for j in range(len(y)):
            if grab[0] == y[j]:
                x[9*i+0]=grab # original
                x[9*i+1]=grabxp # +x
                x[9*i+2]=grabxn # -x
                x[9*i+3]=grabyp # +y
                x[9*i+4]=grabyn # -y
                x[9*i+5]=grabxpyp # +x+y
                x[9*i+6]=grabxpyn # +x-y
                x[9*i+7]=grabxnyp # -x+y
                x[9*i+8]=grabxnyn # -x-y
    x = filter(None, x)
    return(x)

def ss(e):
    if len(e)==2 and e[0]<=e[1]:return e[1]
    return ss(e[:-1]) if e[0]<=e[-1]>=e[1] else ss([e[-1]]+e[:-1])


# if element is found it returns index of element else returns None
def find_element_in_list(element, list_element):
    count = 0
    out = [None]*2
    for i in range(len(list_element)):
        if list_element[i][0] == element:
            out[count] = list_element[i][1]
            count = count +1
    if count ==2:
        return out
###############################################   Manual input starts here               ################################################
# rHS = string.split((raw_input('\natom ID of H_adsorbate : ')))                                                         ################                                                                                                         
# rOS = string.split((raw_input('atom ID of O_adsorbate : ')))                                                           ################                                                                                                       
# rHW = string.split((raw_input('atom ID of H_water     : ')))                                                           ################                                                                                                       
# rOW = string.split((raw_input('atom ID of O_water     : ')))                                                           ################                                                                                                       
                                                                                                                         ################                                         
                                                                                                                    ################                                                                      
rC = [2,3,4,5]                                                                                                           ################                                                       
rHS = [10,11,13]                                                                                                            ################                                                      
rOS = [6,7,12]                                                                                                              ################                                                    
rHW = [9]                                                                                                                ################                                                  
rOW = [8]                                                                                                                ################                                                  
#rC =   [6,7]                                                                                                          ################                                                       
#rHS = [10,11]                                                                                                            ################                                                      
#rOS = [2,3,4,5]                                                                                                             ################                                                    
#rHW = [9]                                                                                                                ################                                                  
#rOW = [8]    
#                                                                                                                          ################                                         
################################              SET your value HERE                                       ################                                                        
fixedposcar = 0   # 1 all fixed, 0 not all fixed

if fixedposcar == 1:
    tail = 'F     F     F'
else:
    tail = 'T     T     T'                                                                                                                      ################
###############################################   Manual input ends here               ##################################################
start = timeit.timeit()
########################################################################################################################
file = open(sys.argv[1],'r')
for i in range(3):
    file.readline() # skip first 3 lines
line = file.readline()
N = int(line)
file.readline() # skip the next line
line = file.readline()
words = string.split(line)
xlo_bound = float(words[0])
xhi_bound = float(words[1])
if len(words)==3:
    xy =  float(words[2])
else:
    xy = 0    # get xy value for the box, now supports both triclinic and orthorgonal
line = file.readline()
words = string.split(line)
ylo_bound =  float(words[0])
yhi_bound =  float(words[1])
if len(words)==3:
    xz =  float(words[2])
else:
    xz = 0    # get xz value for the box, now supports both triclinic and orthorgonal
line = file.readline()
words = string.split(line)
zlo =   float(words[0])
zhi =   float(words[1])
if len(words)==3:
    yz =  float(words[2])
else:
    yz = 0    # get yz value for the box, now supports both triclinic and orthorgonal
xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
ylo = ylo_bound - min(0.0,yz)
yhi = yhi_bound - max(0.0,yz)
lx=xhi-xlo
ly=yhi-ylo
lz=zhi-zlo
a=lx
b=sqrt(ly*ly+xy*xy)
c=sqrt(lz*lz+xz*xz+yz*yz)
cosalpha = (xy*xz+ly*yz)/b/c
cosbeta = xz/c
cosgamma = xy/b
V = sqrt((1-cosalpha*cosalpha-cosbeta*cosbeta-cosgamma*cosgamma+2*cosalpha*cosbeta*cosgamma))*a*b*c
# print xlo, xhi, ylo, yhi,zlo, zhi, xy, xz, yz
print '--------------------------------------------------------------------------------------------------------'
print 'lattice info: ', round(a,6), round(b,6), round(c,6), round(degrees(acos(cosalpha)),1), round(degrees(acos(cosbeta)),1), round(degrees(acos(cosgamma)),1), round(V,4)
bx = b * cosgamma
cx = c * cosbeta
by = b * sqrt(1-cosgamma * cosgamma)
cy = c * (cosalpha-cosbeta * cosgamma)/sqrt(1-cosgamma * cosgamma)
file.close
# ########################################################################################################################################################
  
coords = init() #reads in atoms coords

h1 = find (rHS,coords)
h2 = find (rHW,coords)
o1 = find (rOS,coords)
o2 = find (rOW,coords)
C= find (rC,coords) # find atoms with designated type numbers

NH = (len(h1)+len(h2))/9
NHW = len(h2)/9
NHS = NH - NHW
NO = (len(o1)+len(o2))/9
NOW = len(o2)/9
NOS = NO - NOW
NC = len(C)/9
NPt=N-NC-NO-NH

for i in range(len(h1)):
    for j in range(len(o1)):
        d = dist(h1[i],o1[j])
        if d <= 1.2:
            temp.append( (o1[j][4], o1[j][1] ,o1[j][2] ,o1[j][3] , h1[i][4], h1[i][1] ,h1[i][2] ,h1[i][3] , d)) # pbc link for OH in adsorbate: o1(xyz),h1(xyz),do1h1
lads = temp # find all the OH pairs in adsorbates

temp = []
bond = []
for i in range(len(h2)):
    for j in range(len(o2)):
        d = dist(h2[i],o2[j])
        if d <= 1.2:
            temp.append( (o2[j][4], o2[j][1] ,o2[j][2] ,o2[j][3] , h2[i][4], h2[i][1] ,h2[i][2] ,h2[i][3] , d))  #link for OH in water: o2(xyz),h2(xyz),do2h2
            bond.append( (o2[j][4], h2[i][4]))
lwat = temp # find all the OH pairs in water
bond = set(bond)
bond = list(bond)
bond.sort()
myarray = np.asarray(bond) ## this array contains the all the water OH links
# find_element_in_list(55, myarray)
#########################################################################################################################################################
print '--------------------------------------------------------------------------------------------------------'
if fixedposcar == 0:
    print('\x1b[6;30;42m'+'export ALL FIXED POSCAR coordinates'+'\x1b[0m')
else:
    print('\x1b[6;30;42m'+'export PARTIAL RELAXED POSCAR coordinates'+'\x1b[0m')
#########################################################################################################################################################
####################   To change, fixedposcar = 1    will give    all fixed flags '          ############################################################
# ##################              fixedposcar = 0    will give    partailly fixed flags '    ############################################################
#########################################################################################################################################################
print("Found Total Number of atoms:"+"\033[95m {}\033[00m".format((' %d' % N  )))
print 'Pt C O H'
print NPt, NC, NO, NH
print '_________________________________________________________________________________________________________________________________________'
print('O_acceptor         O_donar        O-O_Dist            A-H_Dist           AHD_Angle            ADH_Angle            Ow_coord ')
print '_________________________________________________________________________________________________________________________________________'
fo1 = fo2 =temp = []

## new working loop
for i in range(len(lads)):
    for j in range(len(lwat)):
        do1o2 = dist( lads[i] , lwat[j] )
        h2_t = lwat[j][4:8]
        h1_t = lads[i][4:8]
        do2h2 = lwat[j][8]
        do1h1 = lads[i][8]
        do1h2 = dist( lads[i][:4], h2_t )
        do2h1 = dist( lwat[j][:4], h1_t )
        if do1o2 <= float(12.25) and (lads[i][0],lwat[j][0]) not in link and do1h2 <= 6.25:
            A1 = ang( do1o2 , do2h2, do1h2)# angle of o1h2
            A2 = ang( do2h2 , do1h2, do1o2)# angle of o1o2    
            if A1 <= 30 and A2 >= 120:
                link.append((lads[i][0],lwat[j][0])) 
                temp.append((lads[i][0] , lwat[j][0] , round(sqrt(do1o2),8), round(sqrt(do1h2),8) , round(A2,8), round(A1,8), lwat[j][1],lwat[j][2],lwat[j][3]  ))
        if do1o2 <= float(12.25) and (lwat[j][0], lads[i][0]) not in link and do2h1<=6.25:
            A1 = ang( do1o2 , do1h1, do2h1 )# angle of o2h1
            A2 = ang( do1h1 ,do2h1, do1o2 ) # angle of o1o2
            if  A1 <= 30 and A2 >= 120:
                link.append((lwat[j][0], lads[i][0]))
                temp.append((lwat[j][0], lads[i][0] ,  round(sqrt(do1o2),8) , round(sqrt(do2h1),8), round(A2,8), round(A1,8), lwat[j][1],lwat[j][2],lwat[j][3]  ))

## legacy loop
# for i in range(len(o1)):
#   for j in range(len(o2)):
#        do1o2 = dist(o1[i] , o2[j])  
#        if do1o2 <= float(12.25):
#           #print i,j
#           for l in range(len(h2)):
#               do2h2 = dist(o2[j],h2[l])
#               if do2h2<=1.3 and (o1[i][4],o2[j][4]) not in link:
#                    do1h2 = dist(o1[i] , h2[l])
#                    A1 = ang( do1o2 , do2h2, do1h2)# angle of o1h2
#                    A2 = ang( do2h2 , do1h2, do1o2) # angle of o1o2
#                    # if A2>=120:  
#                    if A1<=30 and A2>=120 and do1h2<=6.25:                     
#                        link.append((o1[i][4],o2[j][4]))
#                        temp.append((o1[i][4] , o2[j][4] , round(sqrt(do1o2),8), round(sqrt(do1h2),8) , round(A2,8), round(A1,8), o2[i][1],o2[i][2],o2[i][3]   ))
#           for l in range(len(h1)):
#               do1h1 = dist(o1[i],h1[l])
#               if do1h1<=1.3 and (o2[j][4], o1[i][4]) not in link:
#                    do2h1 = dist(o2[j] , h1[l])
#                    # print do1h1,do1o2,do2h1
#                    A1 = ang( do1o2 , do1h1, do2h1 )# angle of o2h1
#                    A2 = ang( do1h1 ,do2h1, do1o2 ) # angle of o1o2
#                    # if  A2>=120:
#                    if  A1<=30 and A2>=120 and do2h1<=6.25:
#                        link.append((o2[j][4], o1[i][4]))
#                        temp.append(( o2[j][4], o1[i][4] ,  round(sqrt(do1o2),8) , round(sqrt(do2h1),8), round(A2,8), round(A1,8), o2[i][1],o2[i][2],o2[i][3]  ))

temp=list(set(temp))
temp=sorted(temp)
temp = [x for x in temp if x != None]
link=list((set(link)))
link=sorted(link)
link = [x for x in link if x != None]
fh = []   ## list of hbond related hydrogens
for i in range(len(link)):
    fh.append(find_element_in_list(link[i][0], myarray) )
    fh.append(find_element_in_list(link[i][1], myarray) )
fh = [x for x in fh if x != None]
fh = np.asarray(fh)
#print fh
#print link
if link!=[]:
    fo1, fo2  = zip(*link)  ## list of hbond oxygens
f1=open('hbonds.log', 'w+')
f1.write(' , '.join('%s-%s' % x for x in link))
for i in range(len(fo1)):
    if fo1.count(fo1[i])==1 and fo2.count(fo2[i])==1:
        print '   ', "%.3d"%temp[i][0], "    I      ", "%.3d"%temp[i][1],  "         %.8f    " %temp[i][2], "     %.8f      " %temp[i][3]," %.8f" %temp[i][4],'     ', "    %.8f     " %temp[i][5] , "%.4f," %temp[i][6], "%.4f," %temp[i][7], "%.4f" %temp[i][8]
    if fo1.count(fo1[i])==2 or fo2.count(fo2[i])==2:
        print '   ', "%.3d"%temp[i][0], "    II     ", "%.3d"%temp[i][1],  "         %.8f    " %temp[i][2], "     %.8f      " %temp[i][3]," %.8f" %temp[i][4],'     ', "    %.8f     " %temp[i][5] , "%.4f," %temp[i][6], "%.4f," %temp[i][7], "%.4f" %temp[i][8]
    if fo1.count(fo1[i])==3 or fo2.count(fo2[i])==3:
        print '   ', "%.3d"%temp[i][0], "    III    ", "%.3d"%temp[i][1],  "         %.8f    " %temp[i][2], "     %.8f      " %temp[i][3]," %.8f" %temp[i][4],'     ', "    %.8f     " %temp[i][5] , "%.4f," %temp[i][6], "%.4f," %temp[i][7], "%.4f" %temp[i][8]      

print '_________________________________________________________________________________________________________________________________________'
print 'Total hydrogen bond(s) btw water and the adsorbate:  %d\n' % (len(link))


## POSCAR outputting
header = """                             
   1.00000000000000000     
     ax                   0.0                      0.0
     bx            by                      0.0
     cx            cy               cz
   Pt   C    O    H 
    NPt     NC     NO     NH
Selective dynamics
Direct"""
header = header.replace("ax" , '%10s'%str(a)).replace("bx" , '%10s'%str(bx)).replace("cx" , '%10s'%str(cx)).replace("by" , '%10s'%str(by)).replace("cy" , '%10s'%str(cy)).replace("cz" , '%10s'%str(c))
# lattice vectors seeded into header and use the template thereafter
header = header.replace("NPt" ,str(NPt))

output1 = open('POSCAR', 'w')
output1.writelines((header.replace("NC" ,str(NC)).replace("NO" ,str(NO)).replace("NH" , str(NH)),'\n'))
output2 = open('POSCAR_expbg', 'w')
output2.writelines((header.replace("NC" ,"").replace("C" ,"").replace("NO" ,str(NOW)).replace("NH" , str(NHW)),'\n'))
output3 = open('POSCAR_imp', 'w')
output3.writelines((header.replace("NC" ,str(NC)).replace("NO" ,str(len(fh)+NOS)).replace("NH" , str(2*len(fh)+NHS)),'\n'))
output4 = open('POSCAR_impbg', 'w')
output4.writelines((header.replace("NC" ,"").replace("C" ,"").replace("NO" ,str(len(fh))).replace("NH" , str(2*len(fh))),'\n'))

for i in range(NPt):
    print >> output1, "%.8f    " % coords[i][2], "%.8f    " % coords[i][3] , "%.8f    " % coords[i][4], '        ' , 'F     F     F'
    print >> output2, "%.8f    " % coords[i][2], "%.8f    " % coords[i][3] , "%.8f    " % coords[i][4], '        ' , 'F     F     F'
    print >> output3, "%.8f    " % coords[i][2], "%.8f    " % coords[i][3] , "%.8f    " % coords[i][4], '        ' , 'F     F     F'
    print >> output4, "%.8f    " % coords[i][2], "%.8f    " % coords[i][3] , "%.8f    " % coords[i][4], '        ' , 'F     F     F'
    
for i in range(NC):
    print >> output1, "%.8f    " % C[i*9][1], "%.8f    " % C[i*9][2] , "%.8f    " % C[i*9][3], '        ' , tail
    print >> output3, "%.8f    " % C[i*9][1], "%.8f    " % C[i*9][2] , "%.8f    " % C[i*9][3], '        ' , tail
    
for i in range(len(o1)/9):
        print >> output1, "%.8f    " % o1[i*9][1], "%.8f    " % o1[i*9][2] , "%.8f    " % o1[i*9][3], '        ' , tail
        print >> output3, "%.8f    " % o1[i*9][1], "%.8f    " % o1[i*9][2] , "%.8f    " % o1[i*9][3], '        ' , tail

for i in range(len(o2)/9):
    if o2[i*9][4] in fo1 or o2[i*9][4] in fo2:
        print >> output1, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , tail
        print >> output2, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , tail
        print >> output3, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , tail
        print >> output4, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , tail    
    else:
        print >> output1, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , 'F     F     F'
        print >> output2, "%.8f    " % o2[i*9][1], "%.8f    " % o2[i*9][2] , "%.8f    " % o2[i*9][3], '        ' , 'F     F     F'
for i in range(len(h1)/9):
        print >> output1, "%.8f    " % h1[i*9][1], "%.8f    " % h1[i*9][2] , "%.8f    " % h1[i*9][3], '        ' , tail
        print >> output3, "%.8f    " % h1[i*9][1], "%.8f    " % h1[i*9][2] , "%.8f    " % h1[i*9][3], '        ' , tail

for i in range(len(h2)/9):
    if h2[i*9][4] in fh:
        print >> output1, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , tail
        print >> output2, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , tail
        print >> output3, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , tail
        print >> output4, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , tail
    else:
        print >> output1, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , 'F     F     F'
        print >> output2, "%.8f    " % h2[i*9][1], "%.8f    " % h2[i*9][2] , "%.8f    " % h2[i*9][3], '        ' , 'F     F     F'

print fh

### legacay output
# for i in range(len(coords)):
#     grab = [coords[i][0],coords[i][2],coords[i][3],coords[i][4]]
#     tailF = 'F     F     F\n'
#     tailT = 'T     T     T\n'
#     if fixedposcar == 0:
#         if i+1 <= NPt:
#             print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output4, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#         elif ( i+1  <= NPt+NC+NOS or  i+1  > N-(N-NPt-NC-NO-NHW)):
#             print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#             print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#         else:
#             if (grab[0] in temp):
#                 print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#                 print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#                 print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#                 print >> output4, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailT,
#             else:
#                 print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#                 print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,

#     elif fixedposcar == 1:
#         if  i+1  <= NPt:
#             print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output4, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,

#         elif ( i+1  <= NPt+NC+NOS or  i+1  > N-(N-NPt-NC-NO-NHW)):
#             print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#             print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,

#         else:
#             if (grab[0] in temp):
#                 print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#                 print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#                 print >> output3, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#                 print >> output4, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,

#             else:
#                 print >> output1, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,
#                 print >> output2, "%.8f    " % grab[1], "%.8f    " % grab[2] , "%.8f    " % grab[3], '        ' , tailF,

output1.close
output2.close
output3.close
output4.close
end = timeit.timeit()
print end - start
#references
# https://ac.els-cdn.com/S0167732201003427/1-s2.0-S0167732201003427-main.pdf?_tid=e84e133b-8606-465d-96be-2d0432d25fc9&acdnat=1529636744_214c9e5968b46c24bb2139c81277ab9f
# https://aip.scitation.org/doi/pdf/10.1063/1.458652
# https://aip.scitation.org/doi/pdf/10.1063/1.464521