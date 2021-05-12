import os
import sys
aa = ["ALA", "ASN", "CYS", "GLU", "HIS", "LEU", "MET", "PRO", "THR", "TYR", "ARG", "ASP", "GLN", "GLY", "ILE", "LYS", "PHE", "SER", "TRP", "VAL", "CPR"]
dataset = sys.argv[1]

middle = "CPR"
for i in aa:
    for j in aa:
        os.system("./calculate_triplet_cum %s %s %s %s > out.out" % (i,middle,j, dataset))

i = "CPR"

for middle in aa:
    for j in aa:
        os.system("./calculate_triplet_cum %s %s %s %s > out.out" % (i,middle,j, dataset))

j = "CPR"

for middle in aa:
    for i in aa:
        os.system("./calculate_triplet_cum %s %s %s %s > out.out" % (i,middle,j, dataset))
