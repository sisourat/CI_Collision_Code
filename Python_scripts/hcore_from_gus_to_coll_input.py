import sys
import numpy as np
from printout import *


# from xyz.txt to input.xml
fxyz = open("xyz.txt","r")
fxml = open("input.xml","w")

print >> fxml, '<GetStaInput>'
ncent = int(fxyz.readline().split()[0])
for i in range(ncent):
 print >> fxml, " "
 print >> fxml, "  <Center>"
 w = fxyz.readline().split()
 print >> fxml, "     <Position  x='"+str(w[2])+"' y='"+str(w[3])+"' z='"+str(w[4])+"' />" 
 print >> fxml, "     <Basisset file='/mnt/nfs/.../aobasis"+str(i+1)+".bas' />"
 print >> fxml, "     <Potential charge='"+str(w[1])+"' exponent='0.0' />" 
 print >> fxml, "  </Center>"
print >> fxml, " "
print >> fxml, '</GetStaInput>'

fxyz.close()
fxml.close()


# from hcore_energy.txt to i1e eigenvalues
fe = open("hcore_energy.txt","r")
nmo = int(fe.readline().split()[0])
feig = open("eigenvalues","w")
e = []
print >> feig, nmo
for i in range(nmo):
 dat = fe.readline().split()
 print >> feig, dat[0], dat[1]
 e.append(float(dat[1]))
feig.close()


# from hcore_eigenvectors.txt to eigenvectors
fmo = open("hcore_eigenvectors.txt","r")
feigv = open("eigenvectors","w")
print >> feigv, " # eigenvectors coming from GUS"

w = fmo.readline().split()
nAO =  int(w[0])
nMO =  int(w[1])
print >> feigv, nMO, nAO
for i in range(nMO):
 print >> feigv, i+1, e[i]
 for j in range(nAO):
   w = fmo.readline().split()
   print >> feigv, '      ',w[2]

fmo.close()
feigv.close()

# from AObasis.txt to aobasis.bas
fbas = open("AObasis.txt","r")
gtos = []
gtop = []
gtod = []
gtof = []
gtog = []
gtoh = []

ts = 'f'
tp = 'f'
td = 'f'
tf = 'f'
tg = 'f'
th = 'f'

ic = 1

for d in fbas:
 l = d.split()

 if(len(l)==1 and l[0]=='END'):
  break

 if(len(l)==2 and l[0]=='E'):
  print printgto(ic,gtos,gtop,gtod,gtof,gtog,gtoh)
  ic+=1
  gtos = []
  gtop = []
  gtod = []
  gtof = []
  gtog = []
  gtoh = []

  ts = 'f'
  tp = 'f'
  td = 'f'
  tf = 'f'
  tg = 'f'
  th = 'f'


 if(len(l)==2 and not(l[1]=='1' or l[1]=='0')):
  print 'only uncontracted basis set'
  sys.exit()

 if(len(l)>1 and ts=='t'):
  gtos.append(l[1])
  ts = 'f'
 if(len(l)>1 and l[0]=='S' and l[1]=='1'):
  ts = 't'

 if(len(l)>1 and tp=='t'):
  gtop.append(l[1])
  tp = 'f'
 if(len(l)>1 and l[0]=='P' and l[1]=='1'):
  tp = 't'

 if(len(l)>1 and td=='t'):
  gtod.append(l[1])
  td = 'f'
 if(len(l)>1 and l[0]=='D' and l[1]=='1'):
  td = 't'

 if(len(l)>1 and tf=='t'):
  gtof.append(l[1])
  tf = 'f'
 if(len(l)>1 and l[0]=='F' and l[1]=='1'):
  tf = 't'

 if(len(l)>1 and tg=='t'):
  gtog.append(l[1])
  tg = 'f'
 if(len(l)>1 and l[0]=='G' and l[1]=='1'):
  tg = 't'

 if(len(l)>1 and th=='t'):
  gtoh.append(l[1])
  th = 'f'
 if(len(l)>1 and l[0]=='H' and l[1]=='1'):
  th = 't'

fbas.close()
  
