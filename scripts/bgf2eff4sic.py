import os, sys
from string import atof, atoi
from math import sqrt

bgf=open(sys.argv[1],'r')
bgflines=bgf.readlines()

geos={}
types={}
atoms=0
numAtoms={}
for line in bgflines:
  if line.find('DESCRP')==0:
    name=line.split()[1]
    geos[name]={}
    types[name]=[]
  elif line.find('END')==0:  
    numAtoms[name]=atoms
    atoms=0
  elif line.find('HETATM')==0:
    atoms+=1
    tokens=line.split()
    id=atoi(tokens[1])
    #print name,id,type,x,y,z
    geos[name][id]={}
    geos[name][id]['type']=tokens[6]
    geos[name][id]['x']=atof(tokens[3])
    geos[name][id]['y']=atof(tokens[4])
    geos[name][id]['z']=atof(tokens[5])
    if tokens[6]=="Si": geos[name][id]['q']=4.0
    elif tokens[6]=="C": geos[name][id]['q']=4.0
    elif tokens[6]=="H": geos[name][id]['q']=1.0
    if tokens[6] not in types[name]:
      types[name].append(tokens[6])
  else:
    continue

a2b=1.889725989
for key,value in geos.iteritems():
  out=open("data."+key,'w')
  file=""
#  file+="AJB\n\n%d atom types\n%d atoms\n\nMasses\n\n"%(len(types[key])+1,numAtoms[key]*3)
  molecule=value
  nAtoms=numAtoms[key]
  typecount=1
  type={}
  spin={}
  radius={}
  for t in types[key]:
    if t=="Si": 
      file+="%d %4.2f\n"%(typecount,28.0855)
      type["Si"]=typecount
      spin["Si"]=3
      radius["Si"]=1.691398
    elif t=="C": 
      file+="%d %4.2f\n"%(typecount,12.0107)
      type["C"]=typecount
      spin["C"]=3
      radius["C"]=0.621427
    elif t=="H": 
      file+="%d %4.2f\n"%(typecount,1.00794)
      type["H"]=typecount
      spin["H"]=0
      radius["H"]=0.0
    typecount+=1
  file+="%d %4.2f\n%d %4.2f"%(typecount,1.0000,typecount+1,1.0000)
  file+="\n\nAtoms\n\n"
  for id,atom in molecule.iteritems():
    file+="%d %d %4.2f %d %4.2f %4.2f %4.2f %4.2f\n"%(id,type[atom['type']],atom['q'],spin[atom['type']],radius[atom['type']],a2b*atom['x'],a2b*atom['y'],a2b*atom['z'])
  for id1,prop1 in molecule.iteritems():
    if prop1['type'] == "Si":
      for id2,prop2 in molecule.iteritems():
        dx=prop1['x']-prop2['x']
        dy=prop1['y']-prop2['y']
        dz=prop1['z']-prop2['z']
        r=sqrt(dx*dx+dy*dy+dz*dz)
        if r==0: continue
        elif prop2['type'] == "C":
          if r < 2.2:	#2.02
            #print "%s: Si-C "%(key), id1, id2, r
            nAtoms+=1
            file+="%d %d 0.0 %d 1.5 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+1,1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
            nAtoms+=1
            file+="%d %d 0.0 %d 1.5 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+2,-1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
        elif prop2['type'] == "H":
          if r < 1.6:	#1.48
            #print "%s: Si-H "%(key), id1, id2, r
            nAtoms+=1
            file+="%d %d 0.0 %d 1.984 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+1,1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
            nAtoms+=1
            file+="%d %d 0.0 %d 1.984 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+2,-1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
        elif prop2['type'] == "Si":
          if r < 2.5:	#2.33
            #print "%s: Si-Si "%(key), id1, id2, r
            nAtoms+=1
            file+="%d %d 0.0 %d 1.5 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+1,1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
            nAtoms+=1
            file+="%d %d 0.0 %d 1.5 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+2,-1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
    if prop1['type'] == "C":
      for id2, prop2 in molecule.iteritems():
        dx=prop1['x']-prop2['x']
        dy=prop1['y']-prop2['y']
        dz=prop1['z']-prop2['z']
        r=sqrt(dx*dx+dy*dy+dz*dz)
        if r==0: continue
        elif prop2['type'] == "C":
          if r < 1.65:   #1.54
            #print "%s: C-C "%(key), id1, id2, r
            nAtoms+=1
            file+="%d %d 0.0 %d 1.291042 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+1,1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
            nAtoms+=1
            file+="%d %d 0.0 %d 1.291042 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+2,-1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
        elif prop2['type'] == "H":
          if r < 1.2:   #1.09
            #print "%s: C-H "%(key), id1, id2, r
            nAtoms+=1
            file+="%d %d 0.0 %d 1.54937 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+1,1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
            nAtoms+=1
            file+="%d %d 0.0 %d 1.54937 %4.2f %4.2f %4.2f\n"%(nAtoms,len(types[key])+2,-1,a2b*(prop1['x']-dx/2),a2b*(prop1['y']-dy/2),a2b*(prop1['z']-dz/2))
  typestring=""
  for atom,value in type.iteritems():
    if atom=="Si" or atom=="C": typestring+="%d %s "%(value, atom)
  file="eFF-ECP: %s\n\n%d atom types\n%d atoms\n\n-1000.0 1000.0 xlo xhi\n-1000.0 1000.0 ylo yhi\n-1000.0 1000.0 zlo zhi\n\nMasses\n\n"%(typestring,len(types[key])+2,nAtoms)+file
  out.write(file)
  out.close()
