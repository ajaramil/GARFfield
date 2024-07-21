import sys, os

fin=open(sys.argv[1],'r')
finlines=fin.readlines()

atom={}
for line in finlines:
  if (line.find("CONECT")==0):
    tokens=line.split()
    atom[atoi(token[0])]=[]
    for t in range(1:len(tokens)):
      atom[atoi(token[0])].append(atoi(tokens[t]))
      print ""
