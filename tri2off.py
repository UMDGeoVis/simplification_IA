import sys
import os

theNode = -1;

target = open("complex.off", 'w')
target.write("OFF\n")


val=0

simplexes = []
tuples = []

nodes_string = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        words = line.split()

        if len(words) == 1:
            simplexes.append(words[0]);

        else:
            tuples.append(words);



target.write(""+simplexes[0]+" "+simplexes[1]+" "+repr(0)+"\n")
for ind in range(0,int(simplexes[0])):
    for item in tuples[ind]:
        target.write("%s " % item)
    target.write("\n")


for ind in range(int(simplexes[0]),int(simplexes[1])+int(simplexes[0])):
    target.write("3 ")
    for item in tuples[ind]:
        target.write("%s " % item)
    target.write("\n")

target.close()





target = open("complex.vtk", 'w')

target.write("# vtk DataFile Version 2.0\n\n")
target.write("ASCII\n")
target.write("DATASET UNSTRUCTURED_GRID\n\n")
target.write("POINTS " + simplexes[0] + " float\n")

for ind in range(0,int(simplexes[0])):
    for item in tuples[ind]:
        target.write("%s " % item)
    target.write("\n")

target.write("CELLS " + simplexes[1] + " " + repr(int(simplexes[1])*4) +"\n")

for ind in range(int(simplexes[0]),int(simplexes[1])+int(simplexes[0])):
    target.write("3 ")
    for item in tuples[ind]:
        target.write("%s " % item)
    target.write("\n")


target.write("CELL_TYPES " + simplexes[1] + "\n")

for i in range(int(simplexes[1])):
    target.write(repr(5) + " "),
target.write("\n")
