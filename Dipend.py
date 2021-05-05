import os
import re
import sys
import getopt
import random
import math
import subprocess
import time
from struct import unpack
import numpy

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
-p --proline: whether we want to rotate the phi angels of prolines as well, in order to do so, one has to tweak the code of chimera, 0 by default, 1 if you want to tweak your chimerax code
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


    def __init__(self, sequence, base, mode, numberOfStructures, dataSet, cycle, remain, gmxCheck, pro):
        self.sequence = sequence
        self.base = base
        self.mode = mode
        self.numberOfStructures = numberOfStructures
        self.dataSet = dataSet
        self.cycle = cycle
        self.remain = remain
        self.gmxCheck = gmxCheck
        self.angles = [] # phi psi angles for all peptides
        self.proline = pro # whether we want to change the dihedral angles of prolines as well
        self.success = [] # 3 values for each structure, first is 1 if there are no CA-CA clashes, the second is 1 if successfully optimized and the third is 1 if the all atom clashcheck is successful on the optimized structure
        self.current_cycle = 0 # at which round of trying are we?
        self.distributions = []
        self.weights = []

    def addAngles(self, angles):
        self.angles = angles
    def plusplus_current_cycle(self):
        self.current_cycle = self.current_cycle + 1
    def init_success(self):
        for p in range(1,self.numberOfStructures+1):
            self.success.append([-1, -1, -1])
    def add_success1(self, k, number):
        self.success[k-1][0] = number
    def add_success2(self, k, number):
        self.success[k-1][1] = number
    def add_success3(self, k, number):
        self.success[k-1][2] = number

    def add_distributions(self, distributions):
        self.distributions = distributions

    def add_weights(self, weights):
        self.weights = weights

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
# writes user defined custom distributions 
def WriteDistributions(peptides, install):
    filename = install.DataPath+"distribution.in"
    counter = 0
    distributions = []
    weights = []
    for k in range(len(peptides.sequence)):
        distributions.append(0)
    for k in range(len(peptides.sequence)):
        weights.append(0.0)
    with open(filename, "r") as f:
        for line in f.readlines(): 
            if line[0] == "#": 
                continue # skipping header line
            line = line.strip()
            line_ = line.split()
            counter+= 1
            try:
                phi = int(line_[1])
            except:
                WriteLog("Error with phi in the distribution file! It should be an integer. This line will be ignored.\n", peptides)
                continue
            try:
                psi = int(line_[2])
            except:
                WriteLog("Error with psi in the distribution file! It should be an integer. This line will be ignored.\n", peptides)
                continue
            try:
                stdev = int(line_[3])
            except:
                WriteLog("Error with stdev in the distribution file! It should be an integer. This line will be ignored.\n", peptides)
                continue
            try:
                weight = float(line_[4])
            except:
                WriteLog("Error with weight in the distribution file! It should be a float. This line will be ignored.\n", peptides)
                continue
            if "-" in line_[0]:
                sp = line_[0].split('-')
                try:
                    first = int(sp[0])
                except:
                    WriteLog("Error with the residue ranges in the distribution file! First residue is smaller than 1.\n", peptides)
                    first = 1
                try:
                    last = int(sp[1])
                except:
                    WriteLog("Error with the residue ranges in the distribution file! Last residue is larger than the whole peptide.\n", peptides)
                    last = len(peptides.sequence)
                if first<1:
                    WriteLog("Error with the residue ranges in the distribution file! First residue is smaller than 1.\n", peptides)
                    first = 1
                if last>len(peptides.sequence):
                    WriteLog("Error with the residue ranges in the distribution file! Last residue is larger than the whole peptide.\n", peptides)
                    last = len(peptides.sequence)
                for k in range(first-1, last):
                    distributions[k] = counter
                    weights[k] = weight
            else:
                try:
                    k = int(line_[0])
                except:
                    WriteLog("Error with the residue writing in the distribution file! This line will be ignored.\n", peptides)
                    continue
                distributions[k-1] = counter
                weights[k-1] = weight

            command = "python3 %s/TwodimensionalNormalDistribution.py %s %s %s %s/%sDistribution_%s.out" % (sys.path[0], phi, psi, stdev, install.WorkingDirectory, peptides.base, counter)
            try:
                rc = subprocess.check_output(command, shell=True)
            except subprocess.CalledProcessError:
                WriteLog("subprocess.CalledProcessError\n", peptides)
    peptides.add_distributions(distributions)
    peptides.add_weights(weights)

#######################################################################        
# reads the binary file from the database for a given amino acid pair and extracts the angles and cumulative probabilities for this pair
def ReadDatabaseLR(aminoAcidNumber, peptides, install):
    distr = []
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid after the last one is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]

    if peptides.mode == "WEIGHTED_LEFT":
        filename = install.DataPath+"left-"+peptides.dataSet+'/'+peptides.aminoAcidNames[middleaa]+"l.bin"
    else:
        filename = install.DataPath+"right-"+peptides.dataSet+'/'+peptides.aminoAcidNames[middleaa]+"r.bin"
    small = []
    with open(filename, 'rb') as reader: # reading as binary file
        while True:
            nt_raw = reader.read(5) # reading correspoinding bytes
            a_raw = reader.read(3)
            n_raw = reader.read(3)
            x1_raw = reader.read(1) # padding bytes 1
            phi_raw = reader.read(4)
            psi_raw = reader.read(4)
            x2_raw = reader.read(4) # padding bytes 2
            cum_raw = reader.read(8)

            nt = nt_raw.decode('utf-8') # decode byte data into string
            a = a_raw.decode('utf-8')
            n = n_raw.decode('utf-8')

            try:
                phi =  unpack('<f', phi_raw) # unpacking byte data into float
            except:
                break
            try:
                psi =  unpack('<f', psi_raw)
            except:
                break
            try:
                cum =  unpack('<d', cum_raw) # unpacking byte data into double
            except:
                break
            if not nt_raw or not a_raw or not n_raw or not x1_raw or not x2_raw or not phi_raw or not psi_raw or not cum_raw: # if EOF
                break
            if peptides.mode == "WEIGHTED_LEFT":
                if nt[0]== 'l' and a == peptides.aminoAcidNames[middleaa] and n == peptides.aminoAcidNames[leftaa]: # only retrieving data of the needed amino acid pair
                    struct = [phi[0], psi[0], cum[0]]
                    distr.append(struct)
            else:
                if nt[0]== 'r' and a == peptides.aminoAcidNames[middleaa] and n == peptides.aminoAcidNames[rightaa]: # only retrieving data of the needed amino acid pair
                    struct = [phi[0], psi[0], cum[0]]
                    distr.append(struct)
    return distr

#######################################################################
# reads the binary file from the database for an amino acid triplet and extracts the angles and cumulative probabilities
def ReadDatabaseT(aminoAcidNumber, peptides, install):
    distr = []
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid after the last one is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]

    filename = install.DataPath+peptides.dataSet+'/'+middleaa+'/'+leftaa+middleaa+rightaa+".bin"

    with open(filename, 'rb') as reader: # reading as binary file
        while True:
            phi_raw = reader.read(4) # reading correspoinding bytes
            psi_raw = reader.read(4)
            x = reader.read(8) # padding bytes
            p1_raw = reader.read(16)
            p2_raw = reader.read(16)
            cum_raw = reader.read(16)

            try:
                phi = unpack('<i',phi_raw) # unpacking byte data into integer
                psi = unpack('<i',psi_raw)
                cum = numpy.frombuffer(cum_raw, dtype=numpy.longdouble)
            except:
                break
            if not phi_raw or not psi_raw or not x or not p1_raw or not p2_raw or not cum_raw: # if EOF
                break
            small = [phi[0], psi[0], cum[0]]
            distr.append(small)
    return distr


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
    command = "%s/fetch_angles %s l %s %.3f %s %s %s" % (sys.path[0], peptides.aminoAcidNames[leftaa], peptides.aminoAcidNames[middleaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
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
    command = "%s/fetch_angles %s r %s %.3f %s %s %s" % (sys.path[0], peptides.aminoAcidNames[middleaa], peptides.aminoAcidNames[rightaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
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
    command = "%s/fetch_both_angles %s %s %s %.3f %s %s %s" % (sys.path[0], peptides.aminoAcidNames[leftaa], peptides.aminoAcidNames[middleaa], peptides.aminoAcidNames[rightaa], draw, install.DataPath, peptides.dataSet, install.WorkingDirectory) # searches in the binary file for the first phi psi value which has a cumulative sum greater than the random value (cumulative sums are in an ascending order)
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

#######################################################################
# chooses a piece based on the left neighbour, combining the Dunbrack probabilities and the user defined distribution with the user defined weight
def ChooseAnglesWeightedLeftProb(aminoAcidNumber, peptides, install):
    ang = [175,175]
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    draw = random.random() # random value generation between 0 and 1
    data = ReadDatabaseLR(aminoAcidNumber, peptides, install)
    distrnum = peptides.distributions[aminoAcidNumber-1]
    weight = peptides.weights[aminoAcidNumber-1]
    if distrnum == 0:
        return ChooseAnglesLeftProb(aminoAcidNumber, peptides, install) # no distribution for this amino acid chosen, using the Dunbrack values - same as in LEFT mode
    else:
        filename_u = "%s/%sDistribution_%s.out" % (install.WorkingDirectory, peptides.base, distrnum) # this is the user defined distribution
        with open(filename_u, "r") as handle_u:
            lines_u = handle_u.readlines()
            for d in data:
                phi_d = float(d[0])
                psi_d = float(d[1])
                cum_d = float(d[2])
                for line_u in lines_u:
                    line_u = line_u.strip()
                    line_u_ = line_u.split()
                    phi_u = int(line_u_[0])
                    psi_u = int(line_u_[1])
                    cum_u = float(line_u_[3])
                    if phi_d == phi_u and psi_d == psi_u:
                        weighted_cum = weight*cum_u + (1-weight)*cum_d
                        if weighted_cum>draw:
                            ang = [phi_d, psi_d]
                            return ang
    return ang

######################################################################
# chooses a piece based on the right neighbour, combining the Dunbrack probabilities and the user defined distribution with the user defined weight
def ChooseAnglesWeightedRightProb(aminoAcidNumber, peptides, install):
    ang = [175,175]
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]
    draw = random.random() # random value generation between 0 and 1
    data = ReadDatabaseLR(aminoAcidNumber, peptides, install)
    distrnum = peptides.distributions[aminoAcidNumber-1]
    weight = peptides.weights[aminoAcidNumber-1]
    if distrnum == 0:
        return ChooseAnglesRightProb(aminoAcidNumber, peptides, install) # no distribution for this amino acid chosen, using the Dunbrack values - same as in RIGHT mode
    else:
        filename_u = "%s/%sDistribution_%s.out" % (install.WorkingDirectory, peptides.base, distrnum) # this is the user defined distribution
        with open(filename_u, "r") as handle_u:
            lines_u = handle_u.readlines()
            for d in data:
                phi_d = int(d[0])
                psi_d = int(d[1])
                cum_d = float(d[2])
                for line_u in lines_u:
                    line_u = line_u.strip()
                    line_u_ = line_u.split()
                    phi_u = int(line_u_[0])
                    psi_u = int(line_u_[1])
                    cum_u = float(line_u_[3])
                    if phi_d == phi_u and psi_d == psi_u:
                        weighted_cum = weight*cum_u + (1-weight)*cum_d
                        if weighted_cum>draw:
                            ang = [phi_d, psi_d]
                            return ang
    return ang

######################################################################
# chooses a piece with the derived Dunbrack-Ting probabilities (counting both neighbours with summing up the logarithms of the probabilities) and with the user defined probabilities weighting the two with the user defined weight
def ChooseAnglesWeightedTripletProb(aminoAcidNumber, peptides, install):
    ang = [175,175]
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

    data = ReadDatabaseT(aminoAcidNumber, peptides, install) # data extracted from the binary database
    distrnum = peptides.distributions[aminoAcidNumber-1] # which distribution to use
    weight = peptides.weights[aminoAcidNumber-1] # which weiht to use
    if distrnum == 0:
        return ChooseAnglesDerivedProb(aminoAcidNumber, peptides, install) # no distribution for this amino acid chosen, using the Dunbrack values - same as in TRIPLET mode
    else:
        filename_u = "%s/%sDistribution_%s.out" % (install.WorkingDirectory, peptides.base, distrnum) # this is the user defined distribution
        with open(filename_u, "r") as handle_u: # opening the corresponding user defined distribution file
            lines_u = handle_u.readlines()
            for d in data: # the extracted from binary
                phi_d = int(d[0])
                psi_d = int(d[1])
                cum_d = float(d[2])
                for line_u in lines_u:
                    line_u = line_u.strip()
                    line_u_ = line_u.split()
                    phi_u = int(line_u_[0])
                    psi_u = int(line_u_[1])
                    cum_u = float(line_u_[3])      
                    if phi_d == phi_u and psi_d == psi_u:
                        weighted_cum = weight*cum_u + (1-weight)*cum_d
                        if weighted_cum>draw:
                            ang = [phi_d, psi_d]
                            return ang
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
        elif peptides.mode == "WEIGHTED_LEFT":
            a = ChooseAnglesWeightedLeftProb(n, peptides, install)
        elif peptides.mode == "WEIGHTED_RIGHT":
            a = ChooseAnglesWeightedRightProb(n, peptides, install)
        elif peptides.mode == "WEIGHTED_TRIPLET":
            a = ChooseAnglesWeightedTripletProb(n, peptides, install)
        else:
            WriteLog("No valid mode given, running with DERIVED_TRIPLET mode\n", peptides)
            peptides.mode = "TRIPLET"
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
    os.system("%s --nogui --offscreen %s/build.cxc" % (install.ChimeraXPath, install.WorkingDirectory)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script

#######################################################################
def CheckContacts(i, peptides, install):
# Checks for steric clashes in the built structure
    rc = -3
    command = "%s/checkcontacts %s %s %s" % (sys.path[0], i, peptides.base, install.WorkingDirectory)
    try:
        rc = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "contacts_%s.dat" % (peptides.base)
    rc = int(open(filename, "r").readlines()[0])
    return rc

#######################################################################
def GmxCheck(i, peptides, install):
    # this function runs a short energy minimisation in GROMACS on the resulting Peptides to see if everything is OK
    fname = "%s%s.pdb" % (peptides.base, i)
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
                        os.system("mv result%s.pdb %sgmxed_%s.pdb" % (opt_num, peptides.base, i))  
                        os.system("rm result*.pdb")
                if "Segmentation fault" in line:
                    converged = 0
    return converged

#######################################################################
def SetAngleForall(i, peptides, install):
    # loops through the structure in question, setting every dihedral angle to the already predefined values
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
                        if i == 1 and peptides.current_cycle == 0:
                            WriteLog("Proline res %s phi is not rotated.\n" % (build), peptides)
                        rotating_script.write("setattr :%s res psi %s\n" % (build, psi))
                    else:
                        if i == 1 and peptides.current_cycle == 0:
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

    os.system("%s --nogui --offscreen %s/rotate.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    if peptides.proline==1: # rotating proline phi angles
        os.system("%s --nogui --offscreen %s/rotateProPhi.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script # another session due to possible errors

#######################################################################
# Runs Scwrl4
def RunSCWRL4(i, peptides, install):
    infile = "%s/%s%s.pdb" % (install.WorkingDirectory, peptides.base, str(i))
    outfile = "%s/%s%s.pdb" % (install.WorkingDirectory, peptides.base, str(i))
    scwrl4command = '%s -i %s -o %s\n' % (install.Scwrl4Path, infile, outfile)
    os.system(scwrl4command)

#######################################################################
def ChimeraxClashcheck(i, peptides, install, gmxed): # gmxed means whether it is on the optimized structure
    clashes = -1
    maxOverlap = -1.0
    maxDist = 1000.0
    headerNum = 7
    with open("clashcheck.cxc", "w") as clashcheck_script:
        if gmxed == 0:
            clashcheck_script.write("open %s/%s%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
        else:
            clashcheck_script.write("open %s/%sgmxed_%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
        clashcheck_script.write("clashes #1 saveFile chimerax_clashcheck.dat\n")
        clashcheck_script.write("q\n")
    os.system("%s --nogui --offscreen %s/clashcheck.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    with open("chimerax_clashcheck.dat", "r") as clashcheck_file:
        rl = clashcheck_file.readlines()
        for lineNum in range(len(rl)):
            line = rl[lineNum].strip()
            line_list = line.split()
            if len(line_list) == 2:
                try:
                    clashes = int(line_list[0])
                except:
                    WriteLog("Trouble with the ChimeraX clashcheck file read for peptide %s" % (i), peptides)
            if "distance" in rl[lineNum]:
                if len(line_list) == 4:
                    if line_list[0]=="atom1" and line_list[1]=="atom2" and line_list[2]=="overlap" and line_list[3]=="distance":
                        headerNum = lineNum
                        continue
        if headerNum > len(rl)-2: # if there are no clashes at all, arbitrary default values
            if clashes!=0:
                WriteLog("Trouble with chimerax clascheck report reading...", peptides)
        else:
            firstClash = rl[headerNum+1].split()
            maxOverlap = float(firstClash[-2])
            maxDist = float(firstClash[-1])
        clashReport = [clashes, maxOverlap, maxDist]
        return clashReport


#######################################################################
def Rotate(i, peptides, install):
    # This is the main function to govern the fate of one structure, first it derives all the phi psi angles, after that it sets them using chimera and afterwards tries to run Gromacs on it (if not told otherwise)
    CollectAllAngles(peptides, install)
    SetAngleForall(i, peptides, install)
    contacts = CheckContacts(i, peptides, install)
    if contacts == 0:
        #WriteLog("#Structure %s went through CA-CA 4 A clash check.\n" % (i), peptides)
        peptides.add_success1(i,1)
        if peptides.gmxCheck == 1:
            RunSCWRL4(i, peptides, install)
            #clash_report= ChimeraxClashcheck(i, peptides, install, 0)
            #WriteLog("Clash_report_for_structure,_number_of_clashes,_maximum_overlap,_distance:\t%s\t%s\t%s\t%s\n" % (i, clash_report[0], clash_report[1], clash_report[2]), peptides)
            GmxCheck(i, peptides, install)
            conv = GmxLoganalyse(i, peptides, install)
            if conv == 1:
                peptides.add_success2(i,1)
                clash_report= ChimeraxClashcheck(i, peptides, install, 1)
                WriteLog("Clash_report_for_the_gmx_optimized_structure,_number_of_clashes,_maximum_overlap,_distance:\t%s\t%s\t%s\t%s\n"  % (i, clash_report[0], clash_report[1], clash_report[2]), peptides)
                if clash_report[0] == 0: # if there are no clashes at all
                #if clash_report[1]<2.2 and clash_report[2]>1: # this is a less strict condition
                    peptides.add_success3(i,1)

    if peptides.success[i-1] == [1,1,1]:
        WriteLog("# Chosen angles for structure %s:\n" % (i), peptides)
        for k in range(1,len(peptides.sequence)+1):
            WriteLog("%s\t%s\t%s\t%s\n" % (k, peptides.sequence[k-1], peptides.angles[k-1][0], peptides.angles[k-1][1]), peptides)

#######################################################################
def Rename(peptides):
    # This renames the resulting peptides based on whether the Gromacs optimization was successful (only run if Gromacs optimization has been set)
    WriteLog("Renaming files based on success...\n", peptides)
    num_of_successful = 0
    for k in range(1,peptides.numberOfStructures+1):
        fn = "%sgmxed_%s.pdb" % (peptides.base, k)
        if peptides.gmxCheck == 1 and os.path.isfile(fn) and peptides.success[k-1] == [1,1,1]:
            num_of_successful = num_of_successful + 1
            os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
            os.system("mv %sgmxed_%s.pdb bak_%sgmxed_%s.pdb" % (peptides.base, k, peptides.base, k))
        elif peptides.gmxCheck == 0 and peptides.success[k-1][0] == 1:
            os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
        else:
            os.system("mv %s%s.pdb fail_%s%s.pdb" % (peptides.base, k, peptides.base, k))
            if peptides.gmxCheck == 1:
                os.system("mv %sgmxed_%s.pdb fail_%sgmxed_%s.pdb" % (peptides.base, k, peptides.base, k))
    WriteLog("Number of successfully optimized structures: %s\n" % (num_of_successful), peptides)
    new_num = 1
    for l in range(1,peptides.numberOfStructures+1):
        if peptides.gmxCheck == 1 and peptides.success[l-1] == [1,1,1]:
            os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, l, peptides.base, new_num))
            os.system("mv bak_%sgmxed_%s.pdb %sgmxed_%s.pdb" % (peptides.base, l, peptides.base, new_num))
            WriteLog("Renaming structure %s to structure %s\n" % (l, new_num), peptides)
            new_num = new_num + 1
        if peptides.gmxCheck == 0 and peptides.success[l-1][0] == 1:
            os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, l, peptides.base, new_num))
            WriteLog("Renaming structure %s to structure %s\n" % (l, new_num), peptides)
            new_num = new_num + 1


#######################################################################
def Cleanup(peptides):
    WriteLog("Cleaning up...\n", peptides)
    # deletes temporary files crated by the program
    os.system("rm fetch*.dat")
    os.system("rm contacts_%s.dat" % peptides.base)
    # the Gromacs generated files
    os.system("rm check.top")
    os.system("rm check_box.gro")
    os.system("rm check_num.dat")
    os.system("rm check.gro")
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
    base = "structure_"
    sequence = "" # if it is not given, it means trouble!
    dataSet = "tcbig"
    cycle = 10
    remain = 0
    gmxCheck = 1
    pro = 0

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

    MyPeptides = Peptides(sequence, base, mode, numberOfStructures, dataSet, cycle, remain, gmxCheck, pro)

    textForLater.append("Command given: %s\n" % (" ".join(sys.argv)))

    # writing the things to the logfile
    for textLine in textForLater:
        WriteLog(textLine, MyPeptides)

    WorkingDirectory = os.getcwd()
    DataPath = sys.path[0]+"/Data/"

    pathfname = sys.path[0]+"/paths.txt"

    with open(pathfname, "r") as paths: # please specify the paths in paths.txt
        for line in paths:
            line = line.strip()
            line_ = line.split()
            if len(line_)>1:
                if line_[0] == "GromacsPath":
                    GromacsPath = line_[1]
                elif line_[0] == "GromacsSuffix":
                    GromacsSuffix = line_[1]
                elif line_[0] == "ChimeraXPath":
                    ChimeraXPath = line_[1]
                elif line_[0] == "Scwrl4Path":
                    Scwrl4Path = line_[1]

    MyInstall = Install(GromacsPath, GromacsSuffix, DataPath, WorkingDirectory, ChimeraXPath, Scwrl4Path)

    #########
    Build(MyPeptides, MyInstall)
    if MyPeptides.mode == "WEIGHTED_LEFT" or MyPeptides.mode == "WEIGHTED_RIGHT" or MyPeptides.mode == "WEIGHTED_TRIPLET":
        WriteDistributions(MyPeptides, MyInstall)
    MyPeptides.init_success()
    for i in range(1,MyPeptides.numberOfStructures+1):
        Rotate(i, MyPeptides, MyInstall)
    for j in range(1,MyPeptides.cycle):
        WriteLog("Trying at %s\n" % (MyPeptides.current_cycle), MyPeptides)
        MyPeptides.plusplus_current_cycle()
        for i in range(1,MyPeptides.numberOfStructures+1):
            if MyPeptides.success[i-1][2]!=1:
                Rotate(i, MyPeptides, MyInstall)
    if MyPeptides.gmxCheck == 1:    
        Rename(MyPeptides)
    else:
        for i in range(1,MyPeptides.numberOfStructures+1):
            WriteStatus("%s%s.pdb %s\n" % (base, i, MyPeptides.success[i-1][0]), MyPeptides)
    if MyPeptides.remain == 0:
        Cleanup(MyPeptides)

    #for i in range(1,MyPeptides.numberOfStructures+1):
        #WriteLog("Structure %s, success 1: %s, succcs 2: %s, success 3: %s\n" % (i, MyPeptides.success[i-1][0], MyPeptides.success[i-1][1], MyPeptides.success[i-1][2]), MyPeptides)

    elapsed_time = time.time() - start_time
    WriteLog("Time elapsed: %.2f\n" % (elapsed_time), MyPeptides)

Main()
