import os
import sys
import re

num = int(sys.argv[1])
base = sys.argv[2]

os.system("ls %s*_result.pdb > res.del" % (base))

resultfnames = []
missnums = []

with open("res.del", 'r') as res:
    resultfnames=(res.readlines())

for n in range(1,num+1):
    found = 0
    for resultn in resultfnames:
        if re.search(r"%s%s_result.pdb" % (base, n), resultn):
            found = 1
    if found == 0:
        if n not in missnums:
            missnums.append(n)

with open("missing.out", "w") as missing:
    for m in missnums:
        missing.write("%s\n" % (m))

os.system("rm res.del")
