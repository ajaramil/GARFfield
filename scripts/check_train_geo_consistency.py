## Andres Jaramillo-Botero
## Script to eliminate inconsistent/redundant entries in geo and trainset files

import sys,os,re

if len(sys.argv) != 3:
  print "Usage: checkconsistency.py geo_file trainset_file\n"

t=open(sys.argv[2],'r')
g=open(sys.argv[1],'r')
train=t.readlines()
geo=g.readlines()

keys=["CHARGE","ENDCHARGE","GEOMETRY","ENDGEOMETRY","ENERGY","ENDENERGY","CELL", "PARAMETERS","ENDCELL"]
other_keys=["a","b","c","xx","yy","zz","xz","yz","xz"]

def tsplit(string, delimiters):
    """Behaves str.split but supports multiple delimiters."""
    
    delimiters = tuple(delimiters)
    stack = [string,]
    
    for delimiter in delimiters:
        for i, substring in enumerate(stack):
            substack = substring.split(delimiter)
            stack.pop(i)
            for j, _substring in enumerate(substack):
                stack.insert(i+j, _substring)
            
    return stack

def grep(string,list):
    expr = re.compile(string)
    return [elem for elem in list if expr.match(elem)]

def plist(list,message):
  for i in list:
    print "Structure %s not in %s"%(i,message)

# Extract files names from geo
gnames=[]
descrp=grep("DESCRP",geo)
for i in descrp:
  gnames.append(i.split()[1])

# Extract file names called in trainset.in
tnames=[]
entry=0
linesdic={}
for l in train:
  linesdic[entry]=[]
  noc=l.split('!')[0]
  noc2=noc.split('#')[0]	# remove comments starting with # or !
  if l[0] != "#" and len(noc2)>0:
    nops=tsplit(noc2,'+-/,\n')
    for il in nops:
      i=il.split()
      for j in i:
        try: 
          float(j)
        except:
          if j not in keys:
            linesdic[entry].append(j)
          if j not in tnames and j not in keys:
            tnames.append(j)
  entry+=1

intersection_geo_trainset=list(set(tnames).intersection( set(gnames) ))
extra_geo=list(set(intersection_geo_trainset).symmetric_difference(set(gnames)))
extra_train=list(set(gnames).symmetric_difference(set(intersection_geo_trainset)))

#plist(extra_geo,"trainset.in file (remove from geo)")
plist(extra_train,"trainset.in file")

# in trainset.in: comment out lines with file names not in geo
print "Patching trainset.in file into trainset.in.new!"
outfile=open("trainset.in.new",'w')

not_in_geo=[]
flag=False
for k,v in linesdic.iteritems():
  for entry in v:
    if entry not in intersection_geo_trainset and entry not in other_keys:
      flag=True
      break
  if flag:
# Don't write anything to new trainset
#    outfile.write("#"+train[k])
    flag=False
  else:
    t=train[k].split("#")
#    outfile.write(t[0]+"\n")
    outfile.write(t[0])
outfile.close()

outfile=open("geo.new",'w')

print "Cleaning geo file (to match trainset.in) !"

# build clean geo file from trainset.in
geoacc = ""
read_flag = False
for i in geo:
  tmp=i.split("\n")[0]
  tok=tmp.split()
  if len(tok) == 0:
    geoacc += i
    continue
  if tok[0] == 'BIOGRF' or tok[0] == 'XTLGRF':
    geoacc += i
    continue
  if tok[0] == 'END' and read_flag:
    geoacc += i
    outfile.write(geoacc)
    geoacc = ""
    read_flag = False
    continue
  if tok[0] == 'DESCRP':
    if tok[1] in tnames:
      geoacc += i
      read_flag = True
      continue
    else:
      print "WARNING: Removing %s from geo file"%(tok[1])
      geoacc = ""
      read_flag = False
      continue
  if read_flag and len(tok)>0:
#    if tok[0] != '#':
    geoacc += i

outfile.close()
