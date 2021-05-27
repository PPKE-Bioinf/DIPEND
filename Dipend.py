import os
import re
import sys
import getopt
import random
import math
import subprocess
import time
import struct
import numpy as np

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

    aminoAcidNames = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "X": "CPR", "Y": "TYR", "V": "VAL"} # X and CPR are for CYS Proline

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
        self.filenames = []

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

    def add_filenames(self, filenames):
        self.filenames = filenames

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
# makes a two dimensional normal distribution with the given parameters
def TwoDimensionalNormalDistribution(peptides, install, phi_c, psi_c, stdev_c, fn):
    # splot (1/(3*sqrt(6.28)))*exp(-0.5*((x-0)/3)**2)*(1/(3*sqrt(6.28)))*exp(-0.5*((y-0)/3)**2) w pm3d
    res=5 # resolution in degrees
    grid=int(1+360/res)
    # set the desired Ramachandran angles and their deviation here
    #phi_avg=-60 
    #psi_avg=-45
    #std=15
    # shifting the averages to make sure to use the closest grid point 
    phi_avg = phi_c+(res/2)
    psi_avg = psi_c+(res/2)
    std = int(stdev_c)
    outfname = fn
    # this will be the array with the final values
    values=np.zeros((grid,grid))
    # a function to take care of the circular nature of the Ramachandran map
    def inrama (a):
        if a < -180:
            return a+360
        if a >= 180:
            return a-360
        return a
    # a sum of all values for scaling at the end (to add up to a probability of 1)
    allvalues=0
    offset=0
    # Calculating the distribution as a product of two 1D distributions
    for phi in range(-180,180,res):
        for psi in range (-180,180,res):
            z= (1/(std*math.sqrt(2*math.pi)))*math.exp(-0.5*(phi/std)**2)
            z*=(1/(std*math.sqrt(2*math.pi)))*math.exp(-0.5*(psi/std)**2)
            #print(x,y,z)
            # adjusting the angles to the desired region
            #phi_out=inrama(phi-phi_avg)
            #psi_out=inrama(psi-psi_avg)
            phi_out=inrama(phi+phi_avg)
            psi_out=inrama(psi+psi_avg)
            # calculating the coordinates in the grid
            x=math.floor(phi_out/res)
            y=math.floor(psi_out/res)
            #x=int(phi_out/res)
            #y=int(psi_out/res)
            # for very shallow distributions (large stdev) the values might not be 0 at the map boundaries.
            # we will subtract this from all values in order to use only the differences between map regions
            if phi == -180 and psi == -180:
                #print ("z at edge:",z)
                offset=z
            values[x][y]=z-offset
            allvalues += z-offset
    cumulative_probability = 0.0
    # writing the distribution to a file
    # note that we write the map in the conventional format
    # from -180 to 180 but use the adjusted values
    #fgz=open("proba2.dat","w")

    with open(outfname, "wb") as f:
        for phi in range(-180,180,res):
            for psi in range (-180,180,res):
                # this is unnecessary for the current values but makes the code
                # more robust if the numbers are modified.
                x=math.floor(phi/res)
                y=math.floor(psi/res)              
                #x=int(phi/res)
                #y=int(psi/res)
                # using the scaled values ensures that we have a cumulative robability of 1
                probability = values[x][y]/allvalues
                cumulative_probability +=  probability
                packed_values = struct.pack('<ddii',cumulative_probability, probability, phi,psi)
                f.write(packed_values)
                #fgz.write('{0:4d}  {1:4d}  {2:.2e}  {3:.2e} | {4:4d} {5:4d}\n'.format(phi,psi,probability,cumulative_probability,x,y))
            #fgz.write("\n")
    #fgz.close()

#######################################################################
# writes user defined custom distributions 
def WriteDistributions(peptides, install):
    filename = install.DataPath+"distribution.in"
    counter = 0
    distributions = []
    weights = []
    filenames = []
    for k in range(len(peptides.sequence)):
        distributions.append(0)
    for k in range(len(peptides.sequence)):
        weights.append(0.0)
    for k in range(len(peptides.sequence)):
        filenames.append("")
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
            try:
                current_filename = line_[5]
            except:
                current_filename = ""
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
                    filenames[k] = current_filename
            else:
                try:
                    k = int(line_[0])
                except:
                    WriteLog("Error with the residue writing in the distribution file! This line will be ignored.\n", peptides)
                    continue
                distributions[k-1] = counter
                weights[k-1] = weight
                filenames[k-1] = current_filename

            if current_filename == "":
                fn = "%s/%sDistribution_%s.bin" % (install.WorkingDirectory, peptides.base, counter)
            else:
                fn = current_filename
            TwoDimensionalNormalDistribution(peptides, install, phi, psi, stdev, fn)
    peptides.add_distributions(distributions)
    peptides.add_weights(weights)
    peptides.add_filenames(filenames)

######################################################################
# chooses angles according to the correspoinding weighted sum of the Dunbrack and the user given distribution
def ChooseAngles(aminoAcidNumber, peptides, install):
    middleaa = peptides.sequence[aminoAcidNumber-1]
    if aminoAcidNumber == 1:
        leftaa = 'A' # avoiding a boundary problem, the amino acid before the first is set to alanine
    else:
        leftaa = peptides.sequence[aminoAcidNumber-2]
    if aminoAcidNumber == len(peptides.sequence):
        rightaa = 'A' # avoiding a boundary problem, the amino acid after the last one is set to alanine
    else:
        rightaa = peptides.sequence[aminoAcidNumber]

    if leftaa == 'X': # cis proline is not defined as neighbour, only as central amino acid... using trans proline instead
        leftaa = 'P'
    if rightaa == 'X': # cis proline is not defined as neighbour, only as central amino acid... using trans proline instead
        rightaa = 'P'

    ang = [175,175]
    draw = random.random() # random value generation between 0 and 1

    # default weight for non.weighted modes (Dunbrack P will be calculated with (1-weight))
    weight=0
    # distribution for weighted modes
    if (peptides.mode.startswith("WEIGHTED")):   
        distrnum = peptides.distributions[aminoAcidNumber-1]
        weight = peptides.weights[aminoAcidNumber-1]
        filename_u = peptides.filenames[aminoAcidNumber-1]
        if filename_u == "": # if no custom filename given
            filename_u = "%sDistribution_%s.bin" % (peptides.base, peptides.distributions[aminoAcidNumber-1])

        # opening the file here
        if distrnum>0:
            reader_u=open(filename_u,'rb')
        
    # choosing the proper distribution files according to left, right or both neigbors as set
    if peptides.mode == "LEFT" or peptides.mode == "WEIGHTED_LEFT":
        filename_d = install.DataPath+peptides.dataSet+"-l/"+middleaa+"/l"+leftaa+".bin"
    if peptides.mode == "RIGHT" or peptides.mode == "WEIGHTED_RIGHT":
        filename_d = install.DataPath+peptides.dataSet+"-r/"+middleaa+"/r"+rightaa+".bin"
    if peptides.mode == "TRIPLET" or peptides.mode == "WEIGHTED_TRIPLET":
        filename_d = install.DataPath+peptides.dataSet+"/"+middleaa+"/"+leftaa+middleaa+rightaa+".bin"

    reader_d=open(filename_d, 'rb')

    found = False

    while (not found):
        if peptides.mode == "TRIPLET" or peptides.mode == "DERIVED_TRIPLET" or peptides.mode == "WEIGHTED_TRIPLET":
            prob_raw_d = reader_d.read(8) # extra data only for triplet
            lnprob_raw_d = reader_d.read(8) # extra data only for triplet
            if not prob_raw_d or not lnprob_raw_d:
                 break

        cum_raw_d = reader_d.read(8)
        phi_raw_d = reader_d.read(4)
        psi_raw_d = reader_d.read(4)
        try:
            cum_d = struct.unpack('<d',cum_raw_d)[0]
            phi_d = struct.unpack('<i',phi_raw_d)[0] # unpacking byte data into integer
            psi_d = struct.unpack('<i',psi_raw_d)[0] # unpacking byte data into integer
        except:
            break
        
        # default for external dist (=no external distribution)
        cum_u=0
        if (peptides.mode.startswith("WEIGHTED")) and distrnum>0: # if there is an external distribution given for this residue
            cum_raw_u = reader_u.read(8)
            prob_raw_u = reader_u.read(8)
            phi_raw_u = reader_u.read(4)
            psi_raw_u = reader_u.read(4)
            try:
                cum_u = struct.unpack('<d',cum_raw_u)[0]
                phi_u = struct.unpack('<i',phi_raw_u)[0] # unpacking byte data into integer
                psi_u = struct.unpack('<i',psi_raw_u)[0]
            except:
                break
            # give a warning if the two angle pairs differ 8should not happen):    
            if phi_d != phi_u or psi_d != psi_u:
                WriteLog("Ooops: %5.3f <=> %5.3f | %5.3f <=> %5.3f" % (phi_d,phi_u,psi_d,psi_u), peptides)

       # calculating final propbability (defaults to cum_d for unweighted modes)     
        weighted_cum = weight*cum_u + (1-weight)*cum_d
        if weighted_cum>draw:
           ang = [phi_d, psi_d] 
           found=True
    return ang

######################################################################
# Collects the angles derived for the given structure
def CollectAllAngles(peptides, install):
    angles = []
    for n in range(1,len(peptides.sequence)+1):
        a = ChooseAngles(n, peptides, install)
        if a == [] or a == None:
            WriteLog("Trouble with amino acid %s in mode %s for amino acid number %s\n" % (peptides.sequence[n], peptides.mode, n), peptides)
        angles.append(a)
    if len(angles[0])<2:
        angles.pop(0)
    peptides.addAngles(angles)

#######################################################################
# Builds up the initial structure with sidechains using ChimeraX
def Build(peptides, install):
    newseq = ""
    for l in peptides.sequence:
        if l == "X": # Chimerax obviously will not recognize X
            newseq = newseq+'P'
        else:
            newseq = newseq+l
    with open("build.cxc", "w") as build_script:
        build_script.write('build start peptide "built" %s ' % (newseq))
        for build in range(1,len(peptides.sequence)+1):
            build_script.write("-65,135 ") # the default initial peptide has a psi of an antiparallel beta sheet, and the phi of the most probable proline conformation
        build_script.write("rotlib Dunbrack\nsave %s/%sinitial.pdb\n" % (install.WorkingDirectory, peptides.base, ))
        build_script.write("q\n")
    os.system("%s --nogui --silent %s/build.cxc" % (install.ChimeraXPath, install.WorkingDirectory)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script

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
                        os.system("mv result%s.pdb %smin_%s.pdb" % (opt_num, peptides.base, i))  
                        #os.system("rm result*.pdb")
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
                if peptides.sequence[build-1]=="P" or peptides.sequence[build-1]=="X":               
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

    os.system("%s --nogui --silent %s/rotate.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    if peptides.proline==1: # rotating proline phi angles
        os.system("%s --nogui --silent %s/rotateProPhi.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script # another session due to possible errors

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
            clashcheck_script.write("open %s/%smin_%s.pdb\n" % (install.WorkingDirectory, peptides.base, str(i)))
        clashcheck_script.write("clashes #1 saveFile chimerax_clashcheck.dat\n")
        clashcheck_script.write("q\n")
    os.system("%s --nogui --silent %s/clashcheck.cxc" % (install.ChimeraXPath, install.WorkingDirectory, )) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
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
        fn = "%smin_%s.pdb" % (peptides.base, k)
        if peptides.gmxCheck == 1 and os.path.isfile(fn) and peptides.success[k-1] == [1,1,1]:
            num_of_successful = num_of_successful + 1
            os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
            os.system("mv %smin_%s.pdb bak_%smin_%s.pdb" % (peptides.base, k, peptides.base, k))
        elif peptides.gmxCheck == 0 and peptides.success[k-1][0] == 1:
            os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
        else:
            os.system("mv %s%s.pdb fail_%s%s.pdb" % (peptides.base, k, peptides.base, k))
            if peptides.gmxCheck == 1:
                os.system("mv %smin_%s.pdb fail_%smin_%s.pdb" % (peptides.base, k, peptides.base, k))
    WriteLog("Number of successfully optimized structures: %s\n" % (num_of_successful), peptides)
    new_num = 1
    for l in range(1,peptides.numberOfStructures+1):
        if peptides.gmxCheck == 1 and peptides.success[l-1] == [1,1,1]:
            os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, l, peptides.base, new_num))
            os.system("mv bak_%smin_%s.pdb %smin_%s.pdb" % (peptides.base, l, peptides.base, new_num))
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
    #os.system("rm fetch*.dat")
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

    # declaring parameters
    numberOfStructures = -99
    mode = ""
    base = ""
    sequence = ""
    dataSet = ""
    cycle = -99
    remain = -99
    gmxCheck = -99
    pro = -99

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

    # precedence is given for the command line parameters, if it is given for a parameter, the value for that parameter in the parameter file is ignored
    parameterfilename = sys.path[0]+"/Data/parameters.in"
    with open(parameterfilename, "r") as par:
        parline = par.readlines()
        for p_ in parline:
            p_ = p_.strip()
            p = p_.split()
            if p[0]=="b" and base=="":
                base = p[1]
            elif p[0]=="c" and cycle==-99:
                try:
                    cycle = int(p[1])
                except:
                    cycle = 10
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="d" and dataSet=="":
                dataSet = p[1]
            elif p[0]=="g" and gmxCheck==-99:
                try:
                    gmxCheck = int(p[1])
                except:
                    gmxCheck = 1
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="m" and mode=="":
                mode = p[1]
            elif  p[0]=="n" and numberOfStructures==-99:
                try:
                    numberOfStructures = int(p[1])
                except:
                    numberOfStructures = 1
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="p" and pro==-99:
                try:
                    pro = int(p[1])
                except:
                    pro = 0
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="r" and remain==-99:
                try:
                    remain = int(p[1])
                except:
                    remain = 0
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="s" and sequence=="":
                sequence = p[1]
    
    MyPeptides = Peptides(sequence, base, mode, numberOfStructures, dataSet.upper(), cycle, remain, gmxCheck, pro)

    #textForLater.append("Command given: %s\n" % (" ".join(sys.argv)))

    # writing the things to the logfile
    for textLine in textForLater:
        WriteLog(textLine, MyPeptides)

    WriteLog("Given parameters:\n \tbase: %s\n \tcycle: %s\n \tdataset: %s\n \tgromacs optimization: %s\n \tmode: %s\n \tnumber of structures: %s\n \tproline phi rotation: %s\n \tremaining temporary files: %s\n \tsequence: %s\n" % (MyPeptides.base, MyPeptides.cycle, MyPeptides.dataSet, MyPeptides.gmxCheck, MyPeptides.mode, MyPeptides.numberOfStructures, MyPeptides.proline, MyPeptides.remain, MyPeptides.sequence), MyPeptides)

    WorkingDirectory = os.getcwd()
    DataPath = sys.path[0]+"/Data/"

    pathfname = sys.path[0]+"/Data/paths.in"

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
