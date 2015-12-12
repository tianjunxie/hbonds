#!/usr/bin/env python
import string,sys
import math
# from decimal import Decimal
# from array import *
temp = lw = la = output = xp = xn = yp = yn = coords = []
h1 = h2 = o1 = o2 = link = [None] * 10000
A  = 0
flag = nlines = 0
count = 0

def dist(X,Y):
    x1 = X[1]
    x2 = Y[1]
    y1 = X[2]
    y2 = Y[2]
    z1 = X[3]
    z2 = Y[3]
    return 2*(x1-x2)*(y1-y2)*a*b*cosgamma + 2*(y1-y2)*(z1-z2)*b*c*cosalpha + 2*(x1-x2)*(z1-z2)*a*c*cosbeta + (x1-x2)*(x1-x2)*a*a + (y1-y2)*(y1-y2)*b*b + (z1-z2)*(z1-z2)*c*c

def nsth(x,y):
    d=dist(x , y)
    # d=min(do,dxp1,dxp2,dxn1,dxn2,dyp1,dyp2,dyn1,dyn2)
    if d<=float(1.2):
        return(d)
    else:
        return(5)

def init():
    x = [None] * 5000
    file = open(sys.argv[1],'r')
    for i in range(9):
        file.readline() # skip first 9 lines
    for i in range(N): # read in x, y, z coordinates 
        line = file.readline()
        words = [float(n) for n in line.split()]
        x.append(words)
    file.close()
    x = [a for a in x if a != None]
    d = x[0][4]
    for i in range(len(x)):
        x[i][0] = int(x[i][0])
        x[i][1] = int(x[i][1])
        x[i][4] = x[i][4]-d
    return (x)

def ang(d1,d2,d3):# Note: square of the distances
    return math.degrees(math.acos((d1 + d2 - d3)/(2.0 * math.sqrt(d1) * math.sqrt(d2))))

def find(y,c):
    x=grab= grabxp = grabxn = grabyp = grabyn = grabxpyp = grabxnyn = grabxnyp = grabxpyn = [None] * 50000
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
                x[9*i+0]=grab
                x[9*i+1]=grabxp
                x[9*i+2]=grabxn
                x[9*i+3]=grabyp
                x[9*i+4]=grabyn
                x[9*i+5]=grabxpyp
                x[9*i+6]=grabxpyn
                x[9*i+7]=grabxnyp
                x[9*i+8]=grabxnyn
    file.close()
    return(x)
###############################################   Manual input starts here               #########################################
# rHS = string.split((raw_input('\natom ID of H_adsorbate : ')))
# rOS = string.split((raw_input('atom ID of O_adsorbate : ')))
# rHW = string.split((raw_input('atom ID of H_water     : ')))
# rOW = string.split((raw_input('atom ID of O_water     : ')))

# rHS = []
# rOS = ['3']
# rHW = ['5']
# rOW = ['4']

rHS = [8,9]
rOS = [4,5]
rHW = [7]
rOW = [6]

# rHS = [10,11]
# rOS = [6]
# rHW = [9]
# rOW = [8]

fixedposcar = 1

###############################################   Manual input ends here               ########################################

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
# ########################################################################################################################################################

coords = init()
h1 = find (rHS,coords)
h2 = find (rHW,coords)
o1 = find (rOS,coords)
o2 = find (rOW,coords)
h1 = [x for x in h1 if x != None]
h2 = [x for x in h2 if x != None]
o1 = [x for x in o1 if x != None]
o2 = [x for x in o2 if x != None]

NH = (len(h1)+len(h2))/9
NHW = len(h2)/9
NHS = NH - NHW
NO = (len(o1)+len(o2))/9
NOW = len(o2)/9
NOS = NO - NOW

for i in range(len(h1)):
    for j in range(len(o1)):
        d = nsth(h1[i],o1[j])
        if d <= 1.2:
            temp.append((o1[j][4],h1[i][4]))
temp=list(set(temp))
la=sorted(temp)
temp = []
for i in range(len(h2)):
    for j in range(len(o2)):
        d = nsth(h2[i],o2[j])
        if d <= 1.2:
            temp.append((o2[j][4],h2[i][4]))
temp=list(set(temp))
lw=sorted(temp)
# for i in range(len(la)):
    # print la[i]
# print len(la)

#########################################################################################################################################################
print('\nFound Total Number of atoms: %d' % N)
print '----------------------------------------------------------------'
if fixedposcar == 1:
    print('Will export ALL FIXED POSCAR coordinates')
else:
    print('Will export PARTIAL RELAXED POSCAR coordinates')
print 'To change, fixedposcar = 1    will give    all fixed flags '
print '           fixedposcar = 0    will give    partailly fixed flags '
print '--------------------------------------------------------------'
print('O_adsorbate     O_water             Dist            Angle')
print '--------------------------------------------------------------'

temp = []
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
                    link.append(o2[j][4])
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
                        link.append(o2[j][4])
                        temp.append((o1[i][4] , o2[j][4] , round(math.sqrt(do1o2),8) , round(A,8)))
temp=list(set(temp))
temp=sorted(temp)
temp = [x for x in temp if x != None]
link=list((set(link)))
link=sorted(link)
link = [x for x in link if x != None]
for i in range(len(temp)):
    if i>0 and temp[i][1]==temp[i-1][1] and temp[i][3]<temp[i-1][3]:
        print '    ', temp[i][0],'         ',temp[i][1],'     ',temp[i][2], '     ',temp[i][3],'       '
        count=count+1
    elif i==0 or temp[i][1]!=temp[i-1][1]:
        print '    ',temp[i][0],'         ',temp[i][1],'     ',temp[i][2], '     ',temp[i][3],'       '
        count=count+1
print 'Total hydrogen bond(s) btw water and the adsorbate:   %d\n' % count


temp = []
for i in range(len(lw)):
    for j in range(len(link)):
        if lw[i][0] == link[j]:
            temp.append(lw[i][0])
            temp.append(lw[i][1])
temp=list(set(temp))
# identifying the total number of atoms as N and write corresponding poscar

output1 = open('POSCAR', 'w')
header1 = """Pt   C   O   H                             
   8.41590000000000     
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     0.5000000095058164    0.8660253995413444    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.5775489799070801
   Pt   C    O    H 
    27     3     NO     NH
Selective dynamics
Direct"""
output1.writelines((header1.replace("NO" ,str(NO)).replace("NH" , str(NH)),'\n'))

output2 = open('POSCAR_expbg', 'w')
header2 = """Pt   C   O   H                             
   8.41590000000000     
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     0.5000000095058164    0.8660253995413444    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.5775489799070801
   Pt     O     H 
    27     NO     NH
Selective dynamics
Direct"""
output2.writelines((header2.replace("NO" ,str(NOW)).replace("NH" , str(NHW)),'\n'))

output3 = open('POSCAR_imp', 'w')
header3 = """Pt   C   O   H                             
   8.41590000000000     
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     0.5000000095058164    0.8660253995413444    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.5775489799070801
   Pt   C    O    H 
    27   3    NO    NH
Selective dynamics
Direct"""
output3.writelines((header3.replace("NO" ,str(len(temp)/3+NOS)).replace("NH" , str(2*len(temp)/3+NHS)),'\n'))

output4 = open('POSCAR_impbg', 'w')
header4 = """Pt   C   O   H                             
   8.41590000000000     
     1.0000000000000000    0.0000000000000000    0.0000000000000000
     0.5000000095058164    0.8660253995413444    0.0000000000000000
     0.0000000000000000    0.0000000000000000    4.5775489799070801
   Pt     O     H 
    27    NO    NH
Selective dynamics
Direct"""
output4.writelines((header4.replace("NO" ,str(len(temp)/3)).replace("NH" , str(2*len(temp)/3)),'\n'))

for i in range(len(coords)):
    grab = [coords[i][1],coords[i][2],coords[i][3],coords[i][4]]
    if fixedposcar != 1:
        tail = 'F     F     F\n'
        if i+1 <= 27:
            print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output4, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
        elif ( i+1  < 34 or  i+1  > N-(N-27-3-NO-NHW)):
            tail = 'T     T     T\n'
            print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
        else:
            if ((i+1) in temp):
                tail = 'T     T     T\n'
                print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output4, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            else:
                tail = 'F     F     F\n'
                print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,

    elif fixedposcar == 1:
        tail = 'F     F     F\n'
        if  i+1  <= 27:
            print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output4, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,

        elif ( i+1  < 34 or  i+1  > N-(N-27-3-NO-NHW)):
            print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
            print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,

        else:
            if ((i+1) in temp):
                print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output3, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output4, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,

            else:
                print >> output1, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,
                print >> output2, "%.8f" % grab[1], "%.8f" % grab[2] , "%.8f" % grab[3], '        ' , tail,

output1.close
output2.close
output3.close
output4.close


