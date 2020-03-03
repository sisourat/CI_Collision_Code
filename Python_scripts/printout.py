import sys
import numpy as np

def printgto(i, gtos, gtop, gtod, gtof, gtog, gtoh):

 fcol=open('aobasis'+str(i)+'.bas','w')

 if(len(gtos)>0):
  print >> fcol, '! s functions'
  print >> fcol, 'H', len(gtos), len(gtos)
  for i in range(len(gtos)):
   c = np.zeros(len(gtos))
   c[i]=1.0
   print >> fcol, gtos[i],' '.join(map(str, c))
 
 if(len(gtop)>0):
  print >> fcol, ' '
  print >> fcol, '! p functions'
  print >> fcol, 'H', len(gtop), len(gtop)
  for i in range(len(gtop)):
   c = np.zeros(len(gtop))
   c[i]=1.0
   print >> fcol, gtop[i],' '.join(map(str, c))
 
 if(len(gtod)>0):
  print >> fcol, ' '
  print >> fcol, '! d functions'
  print >> fcol, 'H', len(gtod), len(gtod)
  for i in range(len(gtod)):
   c = np.zeros(len(gtod))
   c[i]=1.0
   print >> fcol, gtod[i],' '.join(map(str, c))

 if(len(gtof)>0):
  print >> fcol, ' '
  print >> fcol, '! f functions'
  print >> fcol, 'H', len(gtof), len(gtof)
  for i in range(len(gtof)):
   c = np.zeros(len(gtof))
   c[i]=1.0
   print >> fcol, gtof[i],' '.join(map(str, c))

 if(len(gtog)>0):
  print >> fcol, ' '
  print >> fcol, '! g functions'
  print >> fcol, 'H', len(gtog), len(gtog)
  for i in range(len(gtog)):
   c = np.zeros(len(gtog))
   c[i]=1.0
   print >> fcol, gtog[i],' '.join(map(str, c))

 if(len(gtoh)>0):
  print >> fcol, ' '
  print >> fcol, '! h functions'
  print >> fcol, 'H', len(gtoh), len(gtoh)
  for i in range(len(gtoh)):
   c = np.zeros(len(gtoh))
   c[i]=1.0
   print >> fcol, gtoh[i],' '.join(map(str, c))
 

 fcol.close()
 return 0
