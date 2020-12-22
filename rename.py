import os
import sys

num = int(sys.argv[1])
base = sys.argv[2]
inf = sys.argv[3]

os.system("ls %s?.pdb > b1.del" % (base))
os.system("ls %s?_result.pdb > r1.del" % (base))
if num>9:
    os.system("ls %s??.pdb > b2.del" % (base))
    os.system("ls %s??_result.pdb > r2.del" % (base))
if num>99:
    os.system("ls %s???.pdb > b3.del" % (base))
    os.system("ls %s???_result.pdb > r3.del" % (base))
if num>999:
    os.system("ls %s????.pdb > b4.del" % (base))
    os.system("ls %s????_result.pdb > r4.del" % (base))
os.system("cat b*.del > all_b.del")
os.system("cat r*.del > all_r.del")

list_of_numbers = []

with open(inf, 'r') as infh:
    for f in infh.readlines():
        list_of_numbers.append(int(f.strip()))

with open("all_r.del", 'r') as allr:
    resultnames = allr.readlines()
with open("all_b.del", 'r') as allb:
    basenames = allb.readlines()

for r in resultnames:
    r = r.strip()
    os.system("mv %s bak_%s" % (r, r))
for b in basenames:
    b = b.strip()
    os.system("mv %s bak_%s" % (b, b))

if len(list_of_numbers)>len(resultnames):
    print("Too many numbers!")
else:
    for counter in list_of_numbers:
        n = list_of_numbers.index(counter)
        os.system("mv bak_%s %s%s.pdb" % (basenames[n].strip(), base, str(counter)))
        os.system("mv bak_%s %s%s_result.pdb" % (resultnames[n].strip(), base, str(counter)))

os.system("rm *.del")
