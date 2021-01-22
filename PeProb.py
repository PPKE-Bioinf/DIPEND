import os
import re
import sys
import getopt
import random
import math
import subprocess
import time

start_time = time.time() # for measuring purposes

####################### T H E    P R O G R A M

def Help():
    # prints helping message
    print("""this script takes the sequence of a disordered Peptides and tries to build it by starting from a long helix, using scwrl4 to build sidechains and optimize it using GROMACS. After that, the program is setting the phi psi angles one by one, checking for CA-CA clashes.
Random mode - RANDOM: phi and psi angles are chosen randomly, regardless of sequence
Left Probability mode - LEFT: phi and psi angles are chosen by the probabilities of the Dunbrack Lab based on the left neighbour alone
Right Probability mode - RIGHT: phi and psi angles are chosen by the probabilities of the Dunbrack Lab based on the right neighbour alone
Derived Triplet Probability mode - DERIVED_TRIPLET: triplet probabilities were calculated from the data of Dunbrack Ting et al. with their instructions (summing the log probabilities, calculating back the probabilities and summing them to create cumulative probabilities), using their bins to choose with a random number

The input parameters:
-b --base: the beginning of the output filenames (optional, empty string by default)
-c --cycle: how many times it goes through trying to build the sequences (optional, 10 by deault)
-d --dataset: which dataset from Ting et al. 2010 to use: Conly or TCBIG (optional, TCBIG by default)
-g --gmxcheck: perform Gromacs check at the end for each successfully generated structure (1 yes by default, set it to any other value to skip Gromacs check)
-h --help: prints this help message and exits
-m --mode: choose the building mode (see text above, optional: DERIVED_TRIPLET by default)
-n --numofstructures: the number of structures to be generated (optional but highly recommended, 1 by default)
-p --proline: whether we want to rotate the phi angels of prolines as well, in order to do so, one has to tweak the code of chimera, 1 by default, 0 if you do not want to tweak your chimera code and you want to leave the prolines alone
-r --remain: does not delete temporary files (0 (deletes temporary files) by default, set it to 1 if you want the files to remain). Beware! From the Gromacs files, only the very last run is retained even with value 1! (I set the GMX_MAXBACKUP to -1 in order to avoid creating a lot of files, modify the Groamcs running sh files if you need all the Gromacs files!)
-s --sequence: the sequence of the Peptides to be generated (Required!)
""")
    sys.exit()

##################### CLASS ########################################
# the Class of the Peptides to handle input parameters
#class Peptides(sequence, base, mode, numberOfStructures, dataSet, length):
class Peptides():

    aminoAcidNames = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}

    # a list of lists to store all the Ramachandran angle combinations calculated
    phiPsiAngles = []
    angles = [int(x) for x in range(-180, 180, 5)]
    for phi in angles:
        for psi in angles:
            a = [phi, psi]
            phiPsiAngles.append(a)


    def __init__(self, sequence, base, mode, numberOfStructures, dataSet, cycle, remain, gmxCheck, phi0, psi0, sx, sy, pro):
        self.sequence = sequence
        self.base = base
        self.mode = mode
        self.numberOfStructures = numberOfStructures
        self.dataSet = dataSet
        self.cycle = cycle
        self.remain = remain
        self.gmxCheck = gmxCheck
        self.optnum = 0 # the frames Gromacs makes with the optimization of the initial helix
        self.angles = [] # phi psi angles for all peptides
        self.proline = pro # whether we want to change the dihedral angles of prolines as well

    def setOptNum(self, num):
        self.optnum = num
    
    def addAngles(self, angles):
        self.angles = angles

##################### CLASS ########################################
# the Class of Install to handle program paths
class Install():
    def __init__(self, GromacsPath, GromacsSuffix, DataPath, WorkingDirectory, ChimeraXPath, Scwrl4Path):
        self.GromacsPath = GromacsPath
        self.GromacsSuffix = GromacsSuffix
        self.DataPath = DataPath
        self.WorkingDirectory = WorkingDirectory
        self.ChimeraXPath = ChimeraXPath
        self.Scwrl4Path = Scwrl4Path

#######################################################################
# writes to the logfile # attention: it appends! make sure to delete the file before each run
def WriteLog(line, peptides):
    base = peptides.base
    logFileName = "log%s.out" % (base)
    with open(logFileName, "a") as logHandle:
        logHandle.write(line)

#######################################################################
# writes to the statusfile # attention: it appends! make sure to delete the file before each run
def WriteStatus(line, peptides):
    statusFileName = "status%s.out" % (peptides.base)
    with open(statusFileName, "a") as statusHandle:
        statusHandle.write(line)

#######################################################################
# chooses a piece based on the left neighbour, coil only Dunbrack probabilities
def ChooseAnglesLeftProb(aminoAcidNumber, peptides, install):
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    draw = random.random() # random value generation between 0 and 1
    phi = ""
    psi = ""
    command = "%sfetch_angles %s l %s %.3f %s %s %s" % (install.DataPath[0:-5], peptides.aminoAcidNames[leftaa], peptides.aminoAcidNames[middleaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
    try:
        rc = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "fetched_one_neighbour.dat" # this is the output file of the fetch_angles C program
    with open(filename, "r") as handle:
        lines = handle.readlines()
        l = lines[0].split(";")
        phi = float(l[1])
        psi = float(l[2])
        ang = [int(math.ceil(phi)), int(math.ceil(psi))]
    return ang

#######################################################################
# chooses a piece based on the right neighbour, coil only Dunbrack probabilities
def ChooseAnglesRightProb(aminoAcidNumber, peptides, install):
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid after the last one is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]
    draw = random.random() # random value generation between 0 and 1
    phi = ""
    psi = ""
    command = "%sfetch_angles %s r %s %.3f %s %s %s" % (install.DataPath[0:-5], peptides.aminoAcidNames[middleaa], peptides.aminoAcidNames[rightaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
    try:
        rc = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "fetched_one_neighbour.dat" # this is the output file of the fetch_angles C program
    with open(filename, "r") as handle:
        lines = handle.readlines()
        l = lines[0].split(";")
        phi = float(l[1])
        psi = float(l[2])
        ang = [int(math.ceil(phi)), int(math.ceil(psi))]
    return ang

######################################################################
# chooses a piece with the derived Dunbrack-Ting probabilities (counting both neighbours with summing up the logarithms of the probabilities)
def ChooseAnglesDerivedProb(aminoAcidNumber, peptides, install):
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid after the last one is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]

    draw = random.random() # random value generation between 0 and 1
    phi = ""
    psi = ""
    command = "%sfetch_both_angles %s %s %s %.3f %s %s %s" % (install.DataPath[0:-5], peptides.aminoAcidNames[leftaa], peptides.aminoAcidNames[middleaa], peptides.aminoAcidNames[rightaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
    try:
        rc = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "fetched_both.dat" # this is the output file of the fetch_angles_for_both C program
    with open(filename, "r") as handle:
        lines = handle.readlines()
        l = lines[-1].split(";") # we want to read the last line (By default the first line is the last data in case the search crashes, the second line is what we seek. If the program crashes we get phi=175 psi=175 and cum=1 by default)
        phi = l[0] # we only write the phi psi and cum with ; separated into the searched_both.dat file
        psi = l[1]
        ang = [int(phi), int(psi)]
    return ang

######################################################################
# Collects the angles derived for the given structure
def CollectAllAngles(peptides, install):
    angles = []
    for n in range(1,len(peptides.sequence)+1):
        if peptides.mode == "LEFT":
            a = ChooseAnglesLeftProb(n, peptides, install)
        elif peptides.mode == "RIGHT":
            a = ChooseAnglesRightProb(n, peptides, install)
        elif peptides.mode == "DERIVED_TRIPLET" or peptides.mode == "TRIPLET":
            a = ChooseAnglesDerivedProb(n, peptides, install)
        else:
            WriteLog("No valid mode given, running with DERIVED_TRIPLET mode\n", peptides)
            peptides.mode = "DERIVED_TRIPLET"
            a = ChooseAnglesDerivedProb(n, peptides, install)
        if a == []:
            WriteLog("Trouble with amino acid %s in mode %s for amino acid number %s\n" % (peptides.sequence[n], peptides.mode, n), peptides)
        angles.append(a)
    if len(angles[0])<2:
        angles.pop(0)
    peptides.addAngles(angles)

#######################################################################
# Builds up the initial structure with sidechains using ChimeraX
def Build(peptides, install):
    with open("build.cxc", "w") as build_script:
        build_script.write('build start peptide "built" %s ' % (peptides.sequence))
        for build in range(1,len(peptides.sequence)+1):
            build_script.write("-65,135 ") # the default initial peptide has a psi of an antiparallel beta sheet, and the phi of the most probable proline conformation
        build_script.write("rotlib Dunbrack\nsave %s/%sinitial.pdb\n" % (install.WorkingDirectory, peptides.base, ))
        build_script.write("q\n")
    os.system("%schimerax --nogui --offscreen %s/build.cxc" % (install.ChimeraXPath, install.WorkingDirectory)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script

#######################################################################
# This optimizes the initial structure with Gromacs
def GmxOptim(peptides, install):
    # this function runs a short energy minimisation in GROMACS on the resulting Peptides and gives the result in a pdb file
    #bash_command = 'bash gmx_optim.sh -g %s -s "%s" -f %sinitial.pdb -d %s' % (install.GromacsPath, install.GromacsSuffix, peptides.base, install.DataPath)
    #os.system(bash_command)
    fname = "%sinitial.pdb" % peptides.base
    os.environ['GMX_MAXBACKUP'] = '-1'
    c1 = "%sgmx%s pdb2gmx -f %s -o optim.gro -p optim.top -ff amber99sb-ildn -water none -ignh" % (install.GromacsPath, install.GromacsSuffix, fname)
    os.system(c1)
    c2 = "%sgmx%s editconf -f optim.gro -o optim_box.gro -bt cubic -d 2" % (install.GromacsPath, install.GromacsSuffix)
    if os.path.isfile("optim.top"):
        os.system(c2)
    else:
        WriteLog("Trouble with Gromacs optim step 1!\n", peptides)
        sys.exit()
    c3 = "%sgmx%s grompp -f %sem.mdp -c optim_box.gro -p optim.top -o em1.tpr" % (install.GromacsPath, install.GromacsSuffix, install.DataPath)
    if os.path.isfile("optim_box.gro"):
        os.system(c3)
    else:
        WriteLog("Trouble with Gromacs optim step 2!\n", peptides)
        sys.exit()
    c4 = "%sgmx%s mdrun -s em1.tpr -o em1.trr -c em1.gro" % (install.GromacsPath, install.GromacsSuffix)
    if os.path.isfile("em1.tpr"):
        os.system(c4)
    else:
        WriteLog("Trouble with Gromacs optim step 3!\n", peptides)
        sys.exit()
    c5 = "echo 0 | %sgmx%s trjconv -f em1.trr -s em1.tpr -o optim.pdb -pbc nojump -sep" % (install.GromacsPath, install.GromacsSuffix)
    if os.path.isfile("em1.trr"):
        os.system(c5)
    else:
        WriteLog("Trouble with Gromacs optim step 4!\n", peptides)
        sys.exit()
    os.system("ls optim*.pdb | wc > optim_num.dat")
    with open("optim_num.dat", "r") as numf:
        nu = numf.readlines()
        nu_ = nu[0].split()
        number = nu_[0]
        opt_num = int(number)-1 # to see how many optim?.pdb files were created and opening the last one of them
    return opt_num

#######################################################################
def CheckContacts(i, peptides, install):
# Checks for steric clashes in the built structure
    rc = -3
    try:
        rc = subprocess.check_output("%scheckcontacts %s %s %s" % (install.DataPath[0:-5], i, peptides.base, install.WorkingDirectory), shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "contacts_%s.dat" % (peptides.base)
    rc = int(open(filename, "r").readlines()[0])
    return rc

#######################################################################
def GmxCheck(i, peptides, install):
    # this function runs a short energy minimisation in GROMACS on the resulting Peptides to see if everything is OK
    fname = "%s%s_sc.pdb" % (peptides.base, i)
    os.environ['GMX_MAXBACKUP'] = '-1'
    c1 = "%sgmx%s pdb2gmx -f %s -o check.gro -p check.top -ff amber99sb-ildn -water none -ignh" % (install.GromacsPath, install.GromacsSuffix, fname)
    os.system(c1)
    c2 = "%sgmx%s editconf -f check.gro -o check_box.gro -bt cubic -d 2" % (install.GromacsPath, install.GromacsSuffix)
    if os.path.isfile("check.top"):
        os.system(c2)
    else:
        WriteLog("Trouble with Gromacs check step 1!\n", peptides)
        sys.exit()
    c3 = "%sgmx%s grompp -f %sem.mdp -c check_box.gro -p check.top -o em2.tpr" % (install.GromacsPath, install.GromacsSuffix, install.DataPath)
    if os.path.isfile("check_box.gro"):
        os.system(c3)
    else:
        WriteLog("Trouble with Gromacs check step 2!\n", peptides)
        sys.exit()
    c4 = "%sgmx%s mdrun -s em2.tpr -o em2.trr -c em2.gro" % (install.GromacsPath, install.GromacsSuffix)
    if os.path.isfile("em2.tpr"):
        os.system(c4)
    else:
        WriteLog("Trouble with Gromacs check step 3!\n", peptides)
        sys.exit()
    if os.path.isfile("em2.trr"):
        pass
    else:
        WriteLog("Trouble with Gromacs check step 4!\n", peptides)
        sys.exit()

#######################################################################
def GmxLoganalyse(i, peptides, install):
    # this function extracts the results of the energy minimisation, whether it was successful and if it is successful, saves the result in a separate pdb file with the ending _result.pdb
    converged = 0
    infname = "md.log"
    if os.path.isfile(infname):
        with open(infname, 'r') as inhandle:
            for line in inhandle.readlines(): # GROMACS log
                if "converged to Fmax" in line:
                    converged = 1
                    c5 = "echo 0 | %sgmx%s trjconv -f em2.trr -s em2.tpr -o result.pdb -pbc nojump -sep" % (install.GromacsPath, install.GromacsSuffix)
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
                        os.system("mv result%s.pdb %s%s_result.pdb" % (opt_num, peptides.base, i))  
                        os.system("rm result*.pdb")
                if "Segmentation fault" in line:
                    converged = 0
    if converged == 1:
        WriteStatus("Successful energy minimisation on structure %s\n" % i, peptides)
    else:
        WriteStatus("Could not successfully run energy minimisation on structure %s\n" % i, peptides)

#######################################################################
def SetAngleForall(i, peptides, install):
    # loops through the structure in question, setting every dihedral angle to the already predefined values
    #copy = "cp optim%s.pdb %s%s.pdb" % (str(peptides.optnum), peptides.base, str(i))
    copy = "cp %s/%sinitial.pdb %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, install.WorkingDirectory, peptides.base, str(i))
    os.system(copy)
    with open("rotate.cxc", "w") as rotating_script:
        with open("rotateProPhi.cxc", "w") as prorot_script:
            rotating_script.write("open %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
            prorot_script.write("open %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i))) # another session due to possible errors
            for build in range(1,len(peptides.sequence)+1):
                phi = peptides.angles[build-1][0]   
                psi = peptides.angles[build-1][1]
                if peptides.sequence[build-1]=="P":                    
                    if peptides.proline==0:
                        WriteLog("Proline res %s phi is not rotated.\n" % (build), peptides)
                        rotating_script.write("setattr :%s res psi %s\n" % (build, psi))
                    else:
                        WriteLog("Trying to rotate Proline res %s phi...\n" % (build), peptides)
                        prorot_script.write("setattr :%s res phi %s\n" % (build, phi)) # another session due to possible errors
                        rotating_script.write("setattr :%s res psi %s\n" % (build, psi))

                else:
                    rotating_script.write("setattr :%s res phi %s\n" % (build, phi))
                    rotating_script.write("setattr :%s res psi %s\n" % (build, psi))

            rotating_script.write("save %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
            prorot_script.write("save %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i))) # another session due to possible errors
            rotating_script.write("q\n")
            prorot_script.write("q\n")

    os.system("%schimerax --nogui --offscreen %s/rotate.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    if peptides.proline==1: # rotating proline phi angles
        os.system("%schimerax --nogui --offscreen %s/rotateProPhi.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script # another session due to possible errors

#######################################################################
# Runs Scwrl4
def RunSCWRL4(i, peptides, install):
    infile = "%s/%s%s.pdb" % (install.WorkingDirectory, peptides.base, str(i))
    outfile = "%s/%s%s_sc.pdb" % (install.WorkingDirectory, peptides.base, str(i))
    scwrl4command = '%s -i %s -o %s\n' % (install.Scwrl4Path, infile, outfile)
    os.system(scwrl4command)

#######################################################################
def ChimeraxClashcheck(i, peptides, install, gmxed): # gmxed means whether it is on the optimized structure
    clashes = -1
    maxOverlap = -1.0
    maxDist = -1.0
    headerNum = -1
    with open("clashcheck.cxc", "w") as clashcheck_script:
        if gmxed == 0:
            clashcheck_script.write("open %s/%s%s_sc.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
        else:
            clashcheck_script.write("open %s/%s%s_result.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
        clashcheck_script.write("clashes #1 saveFile chimerax_clashcheck.dat\n")
        clashcheck_script.write("q\n")
    os.system("%schimerax --nogui --offscreen %s/clashcheck.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    with open("chimerax_clashcheck.dat", "r") as clashcheck_file:
        for lineNum in range(len(clashcheck_file.readlines())):
            line = clashcheck_file.readlines()[lineNum]
            line_list = line.split()
            if len(line_list) == 2:
                try:
                    clashes = int(line_list[0])
                except:
                    WriteLog("Trouble with the ChimeraX clashcheck file read for peptide %s" % (i), peptides)
            if line_list[0]=="atom1" and line_list[1]=="atom2" and line_list[3]=="overlap" and line_list[4]=="distance":
                headerNum = lineNum
        if headerNum == len(clashcheck_file.readlines()): # if there are no clashes at all, arbitrary values
            if clashes!=0:
                WriteLog("Trouble with chimerax clascheck report reading...", peptides)
            maxOverlap = 0.0
            maxDist = 1000.0
        else:
            firstClash = clashcheck_file.readlines[headerNum+1].split()
            maxOverlap = firstClash[6]
            maxDist = firstClash[7]
        clashReport = [clashes, maxOverlap, maxDist]
        WriteLog(clashReport, peptides)
        return clashReport


#######################################################################
def Rotate(i, peptides, install):
    # This is the main function to govern the fate of one structure, first it derives all the phi psi angles, after that it sets them using chimera and afterwards tries to run Gromacs on it (if not told otherwise)
    CollectAllAngles(peptides, install)
    success = -1
    SetAngleForall(i, peptides, install)
    contacts = CheckContacts(i, peptides, install)
    if contacts == 0:
        if peptides.gmxCheck == 1:
            RunSCWRL4(i, peptides, install)
            clash_report= ChimeraxClashcheck(i, peptides, install, 0)
            WriteLog("# Clash report (success?) for structure %s: number of clashes: %s, maximum overlap: %s, distance: %s\n" % (i, clash_report[0], clash_report[1], clash_report[2]), peptides)
            if clash_report[1]<2.2 and clash_report[3]>1: # TODO threshold! 
                success = 1
                WriteLog("# Clash report for structure %s: number of clashes: %s, maximum overlap: %s, distance: %s\n" % (i, clash_report[0], clash_report[1], clash_report[2]), peptides)
                GmxCheck(i, peptides, install)
                GmxLoganalyse(i, peptides, install)
                clash_report= ChimeraxClashcheck(i, peptides, install, 1)
                WriteLog("# Clash report for the gmx optimized structure %s: number of clashes: %s, maximum overlap: %s, distance: %s\n" % (i, clash_report[0], clash_report[1], clash_report[2]), peptides)
    if success == 1:
        WriteLog("# Chosen angles for structure %s:\n" % (i), peptides)
        for k in range(1,len(peptides.sequence)+1):
            WriteLog("%s\t%s\t%s\t%s\n" % (k, peptides.sequence[k-1], peptides.angles[k-1][0], peptides.angles[k-1][1]), peptides)
    return success

#######################################################################
def Rename(peptides):
    # This renames the resulting peptides based on whether the Gromacs optimization was successful (only run if Gromacs optimization has been set)
    WriteLog("Renaming files based on success...\n", peptides)
    num_of_successful = 0
    success = []
    for k in range(1,peptides.numberOfStructures+1):
        fn = "%s%s_result.pdb" % (peptides.base, k)
        if os.path.isfile(fn):
            success.append([k])
            num_of_successful = num_of_successful + 1
            os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
            os.system("mv %s%s_result.pdb bak_%s%s_result.pdb" % (peptides.base, k, peptides.base, k))
        else:
            os.system("mv %s%s.pdb fail_%s%s.pdb" % (peptides.base, k, peptides.base, k))
    WriteLog("Number of successfully optimized structures: %s\n" % (num_of_successful), peptides)
    for l in range(1, num_of_successful+1):
        os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, success[l-1], peptides.base, l))
        os.system("mv bak_%s%s_result.pdb %s%s_result.pdb" % (peptides.base, success[l-1], peptides.base, l))

#######################################################################
def Cleanup(peptides):
    WriteLog("Cleaning up...\n", peptides)
    # deletes temporary files crated by the program
    seql = len(peptides.sequence)-2
    # the assembled pieces by lsqman
    if seql>99:
        os.system("rm p?.pdb")
        os.system("rm p??.pdb")
        os.system("rm p???.pdb")
    elif seql>9:
        os.system("rm p?.pdb")
        os.system("rm p??.pdb")
    else:
        os.system("rm p?.pdb")
    # the temporary outfiles generated by this program
    os.system("rm fetch*.dat")
    os.system("rm optim_num.dat")
    os.system("rm temp.dat")
    os.system("rm temp%s.dat" % peptides.base)
    os.system("rm contacts_%s.dat" % peptides.base)
    # the Gromacs generated files
    os.system("rm check.top")
    os.system("rm check_box.gro")
    os.system("rm check_num.dat")
    os.system("rm check.gro")
    os.system("rm optim.top")
    os.system("rm optim_box.gro")
    os.system("rm optim.gro")
    os.system("rm em*.gro")
    os.system("rm em*.trr")
    os.system("rm em*.tpr")
    os.system("rm mdout.mdp")
    os.system("rm posre.itp")
    os.system("rm md.log")
    os.system("rm ener.edr")
    # the ChimeraX generated files
    os.system ("rm rotate.cxc")
    os.system("rm build.cxc")
    os.system("rm rotateProPhi.cxc")
    os.system("rm chimerax_clashcheck.dat")
    os.system("rm clashcheck.cxc")
    
########################################################################
####################### M A I N

def Main():

    # default values for the inputs
    numberOfStructures = 1
    mode = "DERIVED_TRIPLET"
    base = ""
    sequence = "" # if it is not given, it means trouble!
    dataSet = "TCBIG"
    cycle = 10
    remain = 0
    gmxCheck = 1
    phi0 = -135
    psi0 = 135
    sx = 0.005
    sy = 0.005
    pro = 1

    textForLater = [] # it will be later written to the logfile, but that cannot yet be opened (because it first needs all the input arguments to be set)

    # taking command line lettered arguments 
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]
    unixOptions = "b:c:d:g:m:n:p:r:s:h"  
    gnuOptions = ["base=", "cycle=", "dataset=", "gmxcheck=", "mode=", "numofstructures=", "proline=", "remain=", "sequence=", "help"]  

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        # output error, and return with an error code
        textForLater.append(str(err)+'\n')
        sys.exit(2)

    # evaluate given options for the inputs
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-b", "--base"):
            base = currentValue
        elif currentArgument in ("-c", "--cycle"):
            cycle = int(currentValue)
        elif currentArgument in ("-d", "--dataset"):
            dataSet = currentValue
        elif currentArgument in ("-g", "--gmxcheck"):
            gmxCheck = int(currentValue)
        elif currentArgument in ("-m", "--mode"):
            mode = currentValue
        elif currentArgument in ("-n", "--numofstructures"):
            numberOfStructures = int(currentValue)
        elif currentArgument in ("-p", "--proline"):
            pro = int(currentValue) 
        elif currentArgument in ("-r", "--remain"):
            remain = int(currentValue)
        elif currentArgument in ("-s", "--sequence"):
            sequence = currentValue
        elif currentArgument in ("-h", "--help"):
            Help()
        else: # I think, it is already handled in the try statement above
            textForLater.append("Option %s with value %s is not recognised\n" % (currentArgument, currentValue))

######################

    MyPeptides = Peptides(sequence, base, mode, numberOfStructures, dataSet, cycle, remain, gmxCheck, phi0, psi0, sx, sy, pro)

    textForLater.append("Command given: %s\n" % (" ".join(sys.argv)))

    # writing the things to the logfile
    for textLine in textForLater:
        WriteLog(textLine, MyPeptides)

    GromacsPath = "/home/gromdev/gromacs-2020/build/bin/" # TODO
    GromacsSuffix = "" # TODO
    DataPath = "/home/zita/Scripts-Research/PEPROB/Data/" # TODO
    WorkingDirectory = os.getcwd()
    ChimeraXPath = "/usr/bin/" # TODO

    Scwrl4Path = "/usr/local/bin/scwrl4/Scwrl4" # TODO
    #Scwrl4Path = "/home/harzi/scwrl4/Scwrl4" # TODO

    MyInstall = Install(GromacsPath, GromacsSuffix, DataPath, WorkingDirectory, ChimeraXPath, Scwrl4Path)

    #########
    Build(MyPeptides, MyInstall)
    #opt_num = GmxOptim(MyPeptides, MyInstall)
    #MyPeptides.setOptNum(opt_num)
    status = []
    for i in range(1,MyPeptides.numberOfStructures+1):
        status.append(-1)
    for i in range(1,MyPeptides.numberOfStructures+1):
        status[i-1]=Rotate(i, MyPeptides, MyInstall)
    for j in range(1,MyPeptides.cycle):
        WriteLog("Trying at %s\n" % (j), MyPeptides)
        for i in range(1,MyPeptides.numberOfStructures+1):
            if status[i-1]!=1:
                status[i-1] = Rotate(i, MyPeptides, MyInstall)
    if MyPeptides.gmxCheck == 1:    
        Rename(MyPeptides)
    else:
        for i in range(1,MyPeptides.numberOfStructures+1):
            WriteStatus("%s%s.pdb %s\n" % (base, i, status[i-1]), MyPeptides)
    if MyPeptides.remain == 0:
        Cleanup(MyPeptides)
    elapsed_time = time.time() - start_time
    WriteLog("Time elapsed: %.2f\n" % (elapsed_time), MyPeptides)

Main()
