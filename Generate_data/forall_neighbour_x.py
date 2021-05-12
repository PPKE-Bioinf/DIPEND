import os
import sys
aa = ["ALA", "ASN", "CYS", "GLU", "HIS", "LEU", "MET", "PRO", "THR", "TYR", "ARG", "ASP", "GLN", "GLY", "ILE", "LYS", "PHE", "SER", "TRP", "VAL", "CPR"]

dataset = sys.argv[1] # conly or TCBIG
nt = sys.argv[2] # l or r

j = "CPR"
for i in aa:
    os.system("./neighbour-template %s %s %s %s > o" % (dataset, nt, i, j))
