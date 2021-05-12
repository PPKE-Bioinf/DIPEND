import os
import sys
aa = ["ALA", "ASN", "CYS", "GLU", "HIS", "LEU", "MET", "PRO", "THR", "TYR", "ARG", "ASP", "GLN", "GLY", "ILE", "LYS", "PHE", "SER", "TRP", "VAL"]
dataset = sys.argv[1]

for middle in aa:
    for i in aa:
        for j in aa:
            os.system("./calculate_triplet_cum %s %s %s %s > out.out" % (i,middle,j, dataset))

