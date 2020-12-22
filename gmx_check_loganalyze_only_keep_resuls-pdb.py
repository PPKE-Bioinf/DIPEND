import os
import re
import sys
import getopt
import random
import math
import subprocess
import time

start_time = time.time()

#######################################################################
# writes to the logfile # attention: it appends! make sure to delete the file before each run
def WriteLog(line, base):
    logFileName = "log_gmx_%s.out" % (base)
    with open(logFileName, "a") as logHandle:
        logHandle.write(line)

#######################################################################
def GmxCheck(i, GromacsPath, GromacsSuffix, base, DataPath):
    # this function runs a short energy minimisation in GROMACS on the resulting Peptides to see if everything is OK
    fname = "%s%s.pdb" % (base, i)
    os.environ['GMX_MAXBACKUP'] = '-1'
    c1 = "%sgmx%s pdb2gmx -f %s -o check.gro -p check.top -ff amber99sb-ildn -water none -ignh" % (GromacsPath, GromacsSuffix, fname)
    os.system(c1)
    c2 = "%sgmx%s editconf -f check.gro -o check_box.gro -bt cubic -d 2" % (GromacsPath, GromacsSuffix)
    if os.path.isfile("check.top"):
        os.system(c2)
    else:
        WriteLog("Trouble with Gromacs check step 1!\n", base)
        sys.exit()
    c3 = "%sgmx%s grompp -f %sem.mdp -c check_box.gro -p check.top -o em2.tpr" % (GromacsPath, GromacsSuffix, DataPath)
    if os.path.isfile("check_box.gro"):
        os.system(c3)
    else:
        WriteLog("Trouble with Gromacs check step 2!\n", base)
        sys.exit()
    c4 = "%sgmx%s mdrun -s em2.tpr -o em2.trr -c em2.gro" % (GromacsPath, GromacsSuffix)
    if os.path.isfile("em2.tpr"):
        os.system(c4)
    else:
        WriteLog("Trouble with Gromacs check step 3!\n", base)
        sys.exit()
    if os.path.isfile("em2.trr"):
        pass
    else:
        WriteLog("Trouble with Gromacs check step 4!\n", peptides)
        sys.exit()

#######################################################################
def GmxLoganalyse(i, GromacsPath, GromacsSuffix, base, DataPath):
    # this function extracts the results of the energy minimisation, whether it was successful
    converged = 0
    outfname = "gmx_summary_%s.out" % (base)
    out = open(outfname, 'a')
    infname = "md.log"
    if os.path.isfile(infname):
        for line in open(infname, 'r').readlines(): # GROMACS log
            if "converged to Fmax" in line:
                converged = 1
                c5 = "echo 0 | %sgmx%s trjconv -f em2.trr -s em2.tpr -o result.pdb -pbc nojump -sep" % (GromacsPath, GromacsSuffix)
                if os.path.isfile("em2.trr"):
                    os.system(c5)
                else:
                    WriteLog("Big trouble with Gromacs optim step 4!\n", peptides)
                    sys.exit()
                os.system("ls result*.pdb | wc > check_num.dat")
                with open("check_num.dat", "r") as numf:
                    nu = numf.readlines()
                    nu_ = nu[0].split()
                    number = nu_[0]
                    opt_num = int(number)-1 # to see how many optim?.pdb files were created and opening the last one of them
                    os.system("mv result%s.pdb %s%s_result.pdb" % (opt_num, base, i))  
                    os.system("rm result*.pdb")
            if "Segmentation fault" in line:
                converged = 0
    if converged == 1:
        out.write("Successful energy minimisation on structure %s\n" % i)
    else:
        out.write("Could not successfully run energy minimisation on structure %s\n" % i)
    out.close()
        
########################################################################
####################### M A I N

def Main():
    base = sys.argv[1]
    os.system("ls *.pdb > pdbs")
    GromacsPath = "/home/gromdev/gromacs-2020/build/bin/" # TODO
    GromacsSuffix = "" # TODO
    DataPath = "/home/zita/Scripts-Research/PEPROB/Data/" # TODO

    with open("pdbs", 'r') as pdbs:
        lines = pdbs.readlines()
        for line in lines:
            line = line.strip()
            i = int(re.match(r'%s(\d*)\.pdb' % (base), line).group(1))
            print(i)
            GmxCheck(i, GromacsPath, GromacsSuffix, base, DataPath)
            GmxLoganalyse(i, GromacsPath, GromacsSuffix, base, DataPath)

Main()
