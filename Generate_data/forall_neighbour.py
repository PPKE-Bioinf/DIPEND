import os
import sys
aa = ["ALA", "ASN", "CYS", "GLU", "HIS", "LEU", "MET", "PRO", "THR", "TYR", "ARG", "ASP", "GLN", "GLY", "ILE", "LYS", "PHE", "SER", "TRP", "VAL"]

dataset = sys.argv[1] # conly or TCBIG
nt = sys.argv[2] # l or r

for i in aa:
    for j in aa:
        os.system("./neighbour-template %s %s %s %s > o" % (dataset, nt, i, j))
