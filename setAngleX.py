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
    runCommand("setattr :%s res phi %s\n" % (aanum, phi))
except BondRotationError:
    print("BondRotationError: bond is part of a cycle for residue %s, not changing its phi angle" % (aanum))
runCommand("setattr :%s res psi %s\n" % (aanum, psi))
runCommand("save %s" % (outfn))
