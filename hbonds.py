#!/usr/bin/env python
import string,sys
import math
# from decimal import Decimal
# from array import *


temp = output = xp = xn = yp = yn = []
h1 = h2 = o1 = o2  = [None] * 10000

A  = 0
flag = 0
count = 0

def dist(X,Y):
    x1=float(X[1])
    x2=float(Y[1])
    y1=float(X[2])
    y2=float(Y[2])
    z1=float(X[3])
    z2=float(Y[3])
    return 2*(x1-x2)*(y1-y2)*a*b*cosgamma + 2*(y1-y2)*(z1-z2)*b*c*cosalpha + 2*(x1-x2)*(z1-z2)*a*c*cosbeta + (x1-x2)*(x1-x2)*a*a + (y1-y2)*(y1-y2)*b*b + (z1-z2)*(z1-z2)*c*c

def nsth(x,y):
    d=dist(x , y)
    # d=min(do,dxp1,dxp2,dxn1,dxn2,dyp1,dyp2,dyn1,dyn2)
    if d<=float(1.2):
        return(d)
    else:
        return(5)

def uniq(input):
  for x in input:
    if x not in output:
      output.append(x)
  return output

def ang(d1,d2,d3):# Note: square of the distances
    return math.degrees(math.acos((d1 + d2 - d3)/(2.0 * math.sqrt(d1) * math.sqrt(d2))))

def find(y):
    x=grab=grabxp=grabxn=grabyp=grabyn = [None] * 10000
    file = open(sys.argv[1],'r')
    for i in range(9):
        file.readline() # skip first 9 lines
    for i in range(N): # read in x, y, z coordinates 
        line = file.readline()
        words = string.split(line)
        grab = ((words[1],words[2],words[3],words[4],words[0]))
        grabxp = ((words[1],str(float(words[2])+1),words[3],words[4],words[0]))
        grabxn = ((words[1],str(float(words[2])-1),words[3],words[4],words[0]))
        grabyp = ((words[1],words[2],str(float(words[3])+1),words[4],words[0]))
        grabyn = ((words[1],words[2],str(float(words[3])-1),words[4],words[0]))
        for j in range(len(y)):
            if grab[0] == y[j]:
                x[5*i+0]=grab
                x[5*i+1]=grabxp
                x[5*i+2]=grabxn
                x[5*i+3]=grabyp
                x[5*i+4]=grabyn
    file.close()
    return(x)

# rHS = string.split((raw_input('\natom ID of H_adsorbate : ')))
# rOS = string.split((raw_input('atom ID of O_adsorbate : ')))
# rHW = string.split((raw_input('atom ID of H_water     : ')))
# rOW = string.split((raw_input('atom ID of O_water     : ')))
# rHS = []
# rOS = ['3']
# rHW = ['5']
# rOW = ['4']

rHS = ['8','9']
rOS = ['4','5']
rHW = ['7']
rOW = ['6']


file = open(sys.argv[1],'r')
for i in range(5):
    file.readline() # skip first 5 lines
line = file.readline()
words = string.split(line)
xlo_bound = float(words[0])
xhi_bound = float(words[1])
xy =  float(words[2])
line = file.readline()
words = string.split(line)
ylo_bound =  float(words[0])
yhi_bound =  float(words[1])
xz =   float(words[2])
line = file.readline()
words = string.split(line)
zlo =   float(words[0])
zhi =   float(words[1])
yz =    float(words[2])
xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
ylo = ylo_bound - min(0.0,yz)
yhi = yhi_bound - max(0.0,yz)
lx=xhi-xlo
ly=yhi-ylo
lz=zhi-zlo
a=lx
b=math.sqrt(ly*ly+xy*xy)
c=math.sqrt(lz*lz+xz*xz+yz*yz)
cosalpha = (xy*xz+ly*yz)/b/c
cosbeta = xz/c
cosgamma = xy/b
# print xlo, xhi, ylo, yhi,zlo, zhi, xy, xz, yz
# print a, b, c, cosalpha, cosbeta, cosgamma
file.close
#identifying the total number of atoms as N and write corresponding poscar
file = open(sys.argv[1],'r')
output = open('POSCAR', 'w')
header = """C  H  O  Pt                             
   8.41590000000000     
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     0.5000000095058164    0.8660253995413444    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.5775489799070801
   Pt   C    O    H 
    27     3     27     55
Selective dynamics
Direct"""
output.writelines((header,'\n'))
for line in file:
    if not line: continue
    if flag == 1:
        N = int(line)
        print('\n----------------------------------------------------')
        print('Found Total Number of atoms: %d' % N)
        print('     ID            ID             Dist            Angle')
        flag = 0
    elif flag == -1:
        words = string.split(line)
        # print words
        grab = ((words[2],'    ' , words[3] ,'    ' , words[4], '     F    F    F', '\n'))
        output.writelines(grab)
        # output.close()
    if line.startswith('ITEM: NUMBER'):
       flag = 1
    elif line.startswith('ITEM: ATOMS'):
       flag = -1
file.close
output.close
# ########################################################################################################################################################
#reading in coordinates plus atom types
h1 = find (rHS)
h2 = find (rHW)
o1 = find (rOS)
o2 = find (rOW)


h1 = [x for x in h1 if x != None]
h2 = [x for x in h2 if x != None]
o1 = [x for x in o1 if x != None]
o2 = [x for x in o2 if x != None]


# ########################################################################################################################################################
for i in range(len(o1)):
    for j in range(len(o2)):
        for l in range(len(h2)):
            do2h2 = nsth(o2[j],h2[l])
            if do2h2<=1.2:
                # print o2[j][4],h2[l][4],dist(o2[j],h2[l])
                do1o2 = dist(o1[i] , o2[j])
                do1h2 = dist(o1[i] , h2[l])

                A = ang( do2h2 , do1o2 , do1h2)

                if  do1o2 <= float(12.25) and A <= 30:
                    temp.append((o1[i][4] , o2[j][4] , round(math.sqrt(do1o2),8) , round(A,8)))
    for i in range(len(o1)):
        for j in range(len(o2)):
            for l in range(len(h1)):
                do1h1 = nsth(o1[i],h1[l])
                if do1h1<=1.2:

                    do1o2 = dist(o1[i] , o2[j])
                    do2h1 = dist(o2[j] , h1[l])
                    # print do1h1,do1o2,do2h1
                    A = ang( do1h1 , do1o2 , do2h1)

                    if  do1o2 <= float(12.25) and A <= 30:
                        temp.append((o1[i][4] , o2[j][4] , round(math.sqrt(do1o2),8) , round(A,8)))
temp=list(set(temp))
temp=sorted(temp)

for i in range(len(temp)):
    if i>0 and temp[i][1]==temp[i-1][1] and temp[i][3]<temp[i-1][3]:
        print '    ',temp[i][0],'         ',temp[i][1],'     ',temp[i][2], '     ',temp[i][3],'       '
        count=count+1
    elif i==0 or temp[i][1]!=temp[i-1][1]:
        print '    ',temp[i][0],'         ',temp[i][1],'     ',temp[i][2], '     ',temp[i][3],'       '
        count=count+1

print '----------------------------------------------------'
print 'Total hydrogen bond(s) btw water and the adsorbate:   %d\n' % count
