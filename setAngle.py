import sys
import chimera
from chimera import runCommand

infn = sys.argv[1]
phi = sys.argv[2]
psi = sys.argv[3]
outfn = sys.argv[4]
aanum = sys.argv[5]

runCommand("open %s\n" % (infn))
try:
    runCommand("setattr r phi %s :%s\n" % (phi, aanum))
except ValueError:
    print("pychimera ValueError: bond is part of a cycle for residue %s, not changing its phi angle" % (aanum))
try:
    runCommand("setattr r psi %s :%s\n" % (psi, aanum))
except ValueError:
    print("pychimera ValueError: bond is part of a cycle for residue %s, not changing its psi angle" % (aanum))
runCommand("write format pdb #0 %s" % (outfn))
