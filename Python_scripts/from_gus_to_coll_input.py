import sys
import numpy as np
from printout import *

nmo_occ = int(sys.argv[1])

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
 print >> fxml, "     <Basisset file='/home/.../aobasis"+str(i+1)+".bas' />"
 print >> fxml, "     <Potential charge='"+str(w[1])+"' exponent='0.0' />" 
 print >> fxml, "  </Center>"
print >> fxml, " "
print >> fxml, '</GetStaInput>'

fxyz.close()
fxml.close()

# from HForbenergy.txt to HF eigenvalues
ehf = np.loadtxt("HForbenergy.txt",dtype={'names': ('ind', 'sym', 'energy'), 'formats': ('i4', 'S1', 'f')})
feig = open("HFeigenvalues","w")
print >> feig, len(ehf)-nmo_occ
for i in range(len(ehf)):
 if(i>nmo_occ-1):
  print >> feig, ehf[i][2]
feig.close()

# from h1eMO.txt to i1e eigenvalues
e = np.loadtxt("h1eMO.txt")
e1e = []
nmo = int(e[0])
feig = open("eigenvalues","w")
k=0
print >> feig, nmo
for i in range(nmo):
 for j in range(i+1):
  k+=1
# print k,e[k]
 if(not(k==0)):
  print >> feig, e[k]
  e1e.append(e[k])
feig.close()



# from mocoef.txt to eigenvectors
fmo = open("mocoef.txt","r")
feigv = open("eigenvectors","w")
print >> feigv, " # eigenvectors coming from GUS"

w = fmo.readline().split()
nAO =  int(w[0])
nMO =  int(w[1])
print >> feigv, nMO, nAO
for i in range(nMO):
 print >> feigv, i+1-nmo_occ, ehf[i][2], e1e[i]
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
  
