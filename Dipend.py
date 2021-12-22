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
    print("""Dipend.py 

DIPEND is a Python-based pipeline to build random structures for intrinsically disordered protein segments, making use of neighbor-dependent backbone conformation preferences described in Ting et al. PLoS Comput Biol 6(4): e1000763. 

It uses ChimeraX, Scwrl4 and GROMACS as external tools.

Input parameteres can be supplied as options or listed in the Data/parameters.in file.

-a --angletoadd: the angle to perturb the dihedral angles in case of a clash (30 degrees by default)
-b --base: the beginning of the output filenames (optional, structure_ by default)
-c --cycle: how many times it goes through trying to build the sequences (optional, 10 by deault)
-d --dataset: which dataset from Ting et al. 2010 to use: Conly or TCBIG (optional, TCBIG by default)
-g --gmxcheck: perform Gromacs check at the end for each successfully generated structure (1 yes by default, set it to any other value to skip Gromacs check)
-h --help: prints this help message and exits
-k --keep: does not delete temporary files (0 (deletes temporary files) by default, set it to 1 if you want the files to remain). Beware! From the Gromacs files, only the very last run is retained even with value 1! (I set the GMX_MAXBACKUP to -1 in order to avoid creating a lot of files, modify the Gromacs sh files if you need all the Gromacs files!)
-m --mode: choose the building mode (see text above, optional: TRIPLET by default)
     LEFT/RIGHT/TRIPLET modes use the left- or right-neighbor based or the combined probabnilites described in Ting et al. 
     WEIGHTED_LEFT/WEIGHTED_RIGHT/WEIGHTED_TRIPLET combines the above probabilities with user-specified residue-specific distributions 
     defined in the file Data/distributions.in
-n --numofstructures: the number of structures to be generated (optional but highly recommended, 1 by default)
-p --proline: whether we want to rotate the phi angles of prolines as well, in order to do so, one has to tweak the code of chimerax, 0 by default, 1 if you want to tweak your chimerax code
-r --random: an angle to randomly perturb the dihedral angles (5 degrees by default)
-s --sequence: the sequence of the protein to be generated (required!)
-t --threadnum: the number of the threads (parallel running processes), 1 by default
-u --unknotmax: the maximum number angle pairs to combinatorically perturb (6 by default) 
-x --range from (for parallel)
-y --range to (for parallel)
""")
    sys.exit()

##################### CLASS ########################################
# the Class of the Peptides to handle input parameters
#class Peptides(sequence, base, mode, numberofstructures, dataset, length):
class Peptides():

    aminoAcidNames = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "X": "CPR", "Y": "TYR", "V": "VAL"} # X and CPR are for CYS Proline

    # a list of lists to store all the Ramachandran angle combinations calculated
    phiPsiAngles = []
    angles = [int(x) for x in range(-180, 180, 5)]
    for phi in angles:
        for psi in angles:
            a = [phi, psi]
            phiPsiAngles.append(a)


    def __init__(self, sequence, base, mode, numberofstructures, dataset, cycle, keep, gmxcheck, pro, random, angleadd, unknot, threadnum, range_from, range_to):
        self.sequence = sequence
        self.base = base
        self.mode = mode
        self.numberofstructures = numberofstructures
        self.dataset = dataset
        self.cycle = cycle
        self.keep = keep
        self.gmxcheck = gmxcheck
        self.angles = [] # phi psi angles for all peptides
        self.proline = pro # whether we want to change the dihedral angles of prolines as well
        self.success = [] # 3 values for each structure, first is 1 if there are no CA-CA clashes, the second is 1 if successfully optimized and the third is 1 if the all atom clashcheck is successful on the optimized structure
        self.currentcycle = 0 # at which round of trying are we?
        self.distributions = []
        self.weights = []
        self.filenames = []
        self.addrandom = random
        self.angletoadd = angleadd
        self.unknotmax = unknot
        self.threadnum = threadnum
        self.range_from = range_from
        self.range_to = range_to

    def addAngles(self, angles):
        self.angles = angles
    def plusplus_currentcycle(self):
        self.currentcycle = self.currentcycle + 1
    def init_success(self):
        for p in range(1,self.numberofstructures+1):
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
    def __init__(self, gromacspath, gromacssuffix, datapath, workingdirectory, chimeraxpath, scwrl4path):
        self.gromacspath = gromacspath
        self.gromacssuffix = gromacssuffix
        self.datapath = datapath
        self.workingdirectory = workingdirectory
        self.chimeraxpath = chimeraxpath
        self.scwrl4path = scwrl4path

#######################################################################
# writes to the logfile # attention: it appends! make sure to delete the file before each run
def WriteLog(line, peptides):
    base = peptides.base
    logFileName = "log%s%s.out" % (base, peptides.range_from)
    with open(logFileName, "a") as logHandle:
        logHandle.write(line)

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
    filename = install.datapath+"distribution.in"
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
                    WriteLog("Error with the residue ranges in the distribution file! Last residue is larger than the lenght of the sequence.\n", peptides)
                    last = len(peptides.sequence)
                if first<1:
                    WriteLog("Error with the residue ranges in the distribution file! First residue is smaller than 1.\n", peptides)
                    first = 1
                if last>len(peptides.sequence):
                    WriteLog("Error with the residue ranges in the distribution file! Last residue is larger than the length of the sequence.\n", peptides)
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
                fn = "%s/%sDistribution_%s.bin" % (install.workingdirectory, peptides.base, counter)
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
        filename_d = install.datapath+peptides.dataset+"-l/"+middleaa+"/l"+leftaa+".bin"
    if peptides.mode == "RIGHT" or peptides.mode == "WEIGHTED_RIGHT":
        filename_d = install.datapath+peptides.dataset+"-r/"+middleaa+"/r"+rightaa+".bin"
    if peptides.mode == "TRIPLET" or peptides.mode == "WEIGHTED_TRIPLET":
        filename_d = install.datapath+peptides.dataset+"/"+middleaa+"/"+leftaa+middleaa+rightaa+".bin"

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
        build_script.write("rotlib Dunbrack\nsave %s/%sinitial.pdb\n" % (install.workingdirectory, peptides.base, ))
        build_script.write("q\n")
    os.system("%s --nogui --silent %s/build.cxc" % (install.chimeraxpath, install.workingdirectory)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script

#######################################################################
def CheckContacts(i, peptides, install):
# Checks for steric clashes in the built structure
    rc = -3
    command = "%s/checkcontacts %s %s %s" % (sys.path[0], i, peptides.base, install.workingdirectory)
    try:
        rc = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        WriteLog("subprocess.CalledProcessError\n", peptides)
    filename = "contacts_%s%s.dat" % (peptides.base, i)
    rc = int(open(filename, "r").readlines()[0])
    return rc

#######################################################################
def AnalyzeContacts(i, peptides, install):
    filename = "contacts_%s%s.lst" % (peptides.base, i)
    # for faster processing, we will make uese of:
    # all lines will contain exactly 3 columns
    # residue numbers are ordered (increasing)
    contact_min_num=0
    resid_a=[-1 for k in range(100)]
    resid_b=[-1 for k in range(100)]
    mindist=[5 for k in range(100)]
    resid_m=[-1 for k in range(100)]
    #resid_m=[]
    linenum=0
    for line in open(filename, "r").readlines():
      resa,resb,distance=line.split()
      resa=int(resa)
      resb=int(resb)
      distance=float(distance)
      if linenum == 0:
          contact_min_num=0
          resid_a[0] = resa
          resid_b[0] = resb
          mindist[0] = distance
          
      elif resa-resid_a[contact_min_num] <= 1 or resb-resid_b[contact_min_num] <= 1: 
          #print("Ezmi?")               
          if distance < mindist[contact_min_num]:
              mindist[contact_min_num] = distance

      else:
         contact_min_num+=1
         if contact_min_num >= 100:
             contact_min_num=-1
             break
         resid_a[contact_min_num] = resa
         resid_b[contact_min_num] = resb
         mindist[contact_min_num] = distance   

      #print ("cm: ",contact_min_num," a:",resa," am: ",resid_a[contact_min_num]," b: ",resb, " bm: ",resid_b[contact_min_num])
      linenum += 1

      
    for cm in range (0,contact_min_num+1):
         #resid_m.append(int((resid_a+resid_b)/2))
  
         resid_m[cm]=int((resid_a[cm]+resid_b[cm])/2)
         #print ("CONTACT: ",cm, " residues: ",resid_a[cm],resid_b[cm]," resid_m ",resid_m[cm]," dist: ",mindist[cm])

    return resid_m

#######################################################################
def GmxCheck(i, peptides, install):
    # this function runs a short energy minimisation in GROMACS on the resulting Peptides to see if everything is OK
    fname = "%s%s.pdb" % (peptides.base, i)
    os.environ['GMX_MAXBACKUP'] = '-1'
    c1 = "%sgmx%s pdb2gmx -f %s -o check_%s.gro -p check_%s.top -ff amber99sb-ildn -water none -ignh" % (install.gromacspath, install.gromacssuffix, fname, i, i)
    os.system(c1)
    c2 = "%sgmx%s editconf -f check_%s.gro -o check_box_%s.gro -bt cubic -d 2" % (install.gromacspath, install.gromacssuffix, i, i)
    if os.path.isfile("check_%s.top" % (i)):
        os.system(c2)
    else:
        WriteLog("Trouble with Gromacs check step 1!\n", peptides)
        sys.exit()
    c3 = "%sgmx%s grompp -f %sem.mdp -c check_box_%s.gro -p check_%s.top -o em2_%s.tpr" % (install.gromacspath, install.gromacssuffix, install.datapath, i, i, i)
    if os.path.isfile("check_box_%s.gro" % (i)):
        os.system(c3)
    else:
        WriteLog("Trouble with Gromacs check step 2!\n", peptides)
        sys.exit()
    c4 = "%sgmx%s mdrun -s em2_%s.tpr -o em2_%s.trr -c em2_%s.gro -e ener_%s.edr -g md_%s.log" % (install.gromacspath, install.gromacssuffix, i, i, i, i, i)
    if os.path.isfile("em2_%s.tpr" % (i)):
        os.system(c4)
    else:
        WriteLog("Trouble with Gromacs check step 3!\n", peptides)
        sys.exit()
    if os.path.isfile("em2_%s.trr" % (i)):
        pass
    else:
        WriteLog("Trouble with Gromacs check step 4!\n", peptides)
        sys.exit()

#######################################################################
def GmxLoganalyse(i, peptides, install):
    # this function extracts the results of the energy minimisation, whether it was successful and if it is successful, saves the result in a separate pdb file with the ending _result.pdb
    converged = 0
    infname = "md_%s.log" % (i)
    if os.path.isfile(infname):
        with open(infname, 'r') as inhandle:
            for line in inhandle.readlines(): # GROMACS log
                if "converged to Fmax" in line:
                    converged = 1
                    c5 = "echo 0 | %sgmx%s trjconv -f em2_%s.trr -s em2_%s.tpr -o result.pdb -pbc nojump -sep" % (install.gromacspath, install.gromacssuffix, i, i)
                    if os.path.isfile("em2_%s.trr" % (i)):
                        os.system(c5)
                    else:
                        WriteLog("Big trouble with Gromacs optim step 4!\n", peptides)
                        sys.exit()
                    os.system("ls result*.pdb | wc > check_num_%s.dat" % (i))
                    with open("check_num_%s.dat" % (i), "r") as numf:
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
    random_addition = peptides.addrandom
    copy = "cp %s/%sinitial.pdb %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, install.workingdirectory, peptides.base, str(i))
    os.system(copy)
    with open("rotate_%s.cxc" % (i), "w") as rotating_script:
        with open("rotateProPhi_%s.cxc" % (i), "w") as prorot_script:
            rotating_script.write("open %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))
            prorot_script.write("open %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i))) # another session due to possible errors
            for build in range(1,len(peptides.sequence)+1):
                phi = peptides.angles[build-1][0]   #+random_add
                psi = peptides.angles[build-1][1]   #+random_add
                if peptides.sequence[build-1]=="P" or peptides.sequence[build-1]=="X":               
                    if peptides.proline==0:
                        if i == 1 and peptides.currentcycle == 0:
                            WriteLog("Proline res %s phi is not rotated.\n" % (build), peptides)
                        rotating_script.write("setattr :%s res psi %s\n" % (build, psi))
                    else:
                        if i == 1 and peptides.currentcycle == 0:
                            WriteLog("Trying to rotate Proline res %s phi...\n" % (build), peptides)
                        prorot_script.write("setattr :%s res phi %s\n" % (build, phi)) # another session due to possible errors
                        rotating_script.write("setattr :%s res psi %s\n" % (build, psi))

                else:
                    rotating_script.write("setattr :%s res phi %s\n" % (build, phi))
                    rotating_script.write("setattr :%s res psi %s\n" % (build, psi))
                # taking care of cis Pro
                if peptides.sequence[build-1]=="X":
                    rotating_script.write("setattr :%s res omega 0\n" % (build))

            rotating_script.write("save %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))
            prorot_script.write("save %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i))) # another session due to possible errors
            rotating_script.write("q\n")
            prorot_script.write("q\n")

    os.system("%s --nogui --silent %s/rotate_%s.cxc" % (install.chimeraxpath, install.workingdirectory, i)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    if peptides.proline==1: # rotating proline phi angles
        os.system("%s --nogui --silent %s/rotateProPhi_%s.cxc" % (install.chimeraxpath, install.workingdirectory, i)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script # another session due to possible errors

#######################################################################
# Runs Scwrl4
def RunSCWRL4(i, peptides, install):
    infile = "%s/%s%s.pdb" % (install.workingdirectory, peptides.base, str(i))
    outfile = "%s/%s%s.pdb" % (install.workingdirectory, peptides.base, str(i))
    scwrl4command = '%s -i %s -o %s\n' % (install.scwrl4path, infile, outfile)
    os.system(scwrl4command)

#######################################################################
def ChimeraxClashcheck(i, peptides, install, gmxed): # gmxed means whether it is on the optimized structure
    clashes = -1
    maxOverlap = -1.0
    maxDist = 1000.0
    headerNum = 7
    with open("clashcheck_%s.cxc" % (i), "w") as clashcheck_script:
        if gmxed == 0:
            clashcheck_script.write("open %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))
        else:
            clashcheck_script.write("open %s/%smin_%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))
        clashcheck_script.write("clashes #1 saveFile chimerax_clashcheck_%s.dat\n" % (i))
        clashcheck_script.write("q\n")
    os.system("%s --nogui --silent %s/clashcheck_%s.cxc" % (install.chimeraxpath, install.workingdirectory, i)) # unfortunately, it seems to me, that ChimeraX comes with its own Python 53.7, and I did not find a way to integrate its functions into the python used by my kernel, so I decided to invoke it from bash with a chimerax command script
    with open("chimerax_clashcheck_%s.dat" % (i), "r") as clashcheck_file:
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
    # trying to generate structures with defined angles and added noise until
    # there are no contacts or the max number of tries is reached
    try_unknot=1
    random_tries=0
    angle_to_add=peptides.angletoadd # should be put to the options properly
    # note that the actual number of tries will be 2^(max_unknot_angles)
    max_unknot_angles = peptides.unknotmax # -> 64 tries allowed

    SetAngleForall(i, peptides, install)
    contacts = CheckContacts(i, peptides, install)

    WriteLog("i: %s, contacts: %s\n" % (i, contacts), peptides)
    if contacts > 0:
        angle_to_unknot_num=0
        max_tries=1
        anglenam_to_unknot = ["" for r in range(0,max_unknot_angles+3)]
        angleres_to_unknot = [0 for r in range(0,max_unknot_angles+3)]        
        angleval_to_unknot = [0 for r in range(0,max_unknot_angles+3)]
        for residm in AnalyzeContacts(i, peptides, install):
            if residm >= 0 and angle_to_unknot_num <= max_unknot_angles:
                if peptides.sequence[residm]=="P" or peptides.sequence[residm]=="X":
                    anglenam_to_unknot[angle_to_unknot_num]="psi"
                    angleres_to_unknot[angle_to_unknot_num]=residm+1
                    angleval_to_unknot[angle_to_unknot_num]=peptides.angles[residm][1]
                    angle_to_unknot_num += 1
                    max_tries *= 2
                else:
                    anglenam_to_unknot[angle_to_unknot_num]="phi"
                    angleres_to_unknot[angle_to_unknot_num]=residm+1
                    angleval_to_unknot[angle_to_unknot_num]=peptides.angles[residm][0]
                    anglenam_to_unknot[angle_to_unknot_num+1]="psi"
                    angleres_to_unknot[angle_to_unknot_num+1]=residm+1
                    angleval_to_unknot[angle_to_unknot_num+1]=peptides.angles[residm][1]
                    angle_to_unknot_num += 2
                    max_tries *= 4
                    
        if angle_to_unknot_num <= max_unknot_angles:
            # systematically generate a binary code
            WriteLog("angles_to_unknot: %s, max_tries: %s\n" % (angle_to_unknot_num, max_tries),peptides)
            code=[0 for r in range(0,angle_to_unknot_num)]
            try_unknot=0
            contacts=1
            while try_unknot < max_tries and contacts > 0:

                with open("unknot_%s.cxc" % (i), "w") as unknot_script:
                    unknot_script.write("open %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))
                    for p in range(0,angle_to_unknot_num):
                        if code[p] == 1:
                            angle_val=angleval_to_unknot[p]-angle_to_add
                        else:
                            angle_val=angleval_to_unknot[p]+angle_to_add
                        unknot_script.write("setattr :%s res %s %s\n" % (angleres_to_unknot[p], anglenam_to_unknot[p],angle_val))
                    unknot_script.write("save %s/%s%s_%s.pdb\n" % (install.workingdirectory, peptides.base, str(i), str(try_unknot)))
                    unknot_script.write("save %s/%s%s.pdb\n" % (install.workingdirectory, peptides.base, str(i)))               
                    unknot_script.write("q\n")
                #print("invoking unknot.cxc")
                os.system("%s --nogui --silent %s/unknot_%s.cxc" % (install.chimeraxpath, install.workingdirectory, i)) 
                contacts = CheckContacts(i, peptides, install)
                WriteLog("code: %s contacts: %s\n" % (code, contacts),peptides) 
                try_unknot += 1
                # increase code
                pos=angle_to_unknot_num-1
                if code[pos]==0:
                    code[pos]+=1
                else:
                    while code[pos] == 1:
                        code[pos]=0
                        pos -= 1
                    code[pos]=1
 
       
    # continuing only when no contacts   
    if contacts == 0:
        #WriteLog("#Structure %s went through CA-CA 4 A clash check.\n" % (i), peptides)
        peptides.add_success1(i,1)
        if peptides.gmxcheck == 1:
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
    # This renames the resulting peptides based on whether there are clashes or whether the Gromacs optimization was successful

    WriteLog("Renaming files based on success...\n", peptides)
    # renaming every pdb file - inserting bak_ before it in order not to overwrite something later
    for k in range(peptides.range_from, peptides.range_to):
        fn = "%smin_%s.pdb" % (peptides.base, k)
        os.system("mv %s%s.pdb bak_%s%s.pdb" % (peptides.base, k, peptides.base, k))
        if os.path.isfile(fn):
            os.system("mv %smin_%s.pdb bak_%smin_%s.pdb" % (peptides.base, k, peptides.base, k))

    # renaming the successful ones so that their numbering will be consecutive from 1
    num_of_successful = peptides.range_from
    for m in range(peptides.range_from, peptides.range_to):
        fn_ = "bak_%smin_%s.pdb" % (peptides.base, m)
        if peptides.gmxcheck == 1:
            if peptides.success[m-1] == [1,1,1]:
                WriteLog("Trying to rename structure %s to structure %s\n" % (m, num_of_successful), peptides)
                if os.path.isfile(fn_):
                    if not os.path.isfile("%smin_%s.pdb" % (peptides.base, num_of_successful)):
                        os.system("mv bak_%smin_%s.pdb %smin_%s.pdb" % (peptides.base, m, peptides.base, num_of_successful))
                    else:
                        WriteLog("Trouble renaming structure %s to structure %s\n" % (m, num_of_successful), peptides)
                else:
                    WriteLog("Trouble executing command: mv bak_%smin_%s.pdb %smin_%s.pdb, file not found!\n" % (peptides.base, m, peptides.base, num_of_successful), peptides)
                if not os.path.isfile("%s%s.pdb" % (peptides.base, num_of_successful)):
                    os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, m, peptides.base, num_of_successful))
                else:
                     WriteLog("Trouble renaming structure %s to structure %s\n" % (m, num_of_successful), peptides)
                num_of_successful += 1
                
        else: # if not gmx
            if peptides.success[m-1][0] == 1:
                WriteLog("Trying to rename structure %s to structure %s\n" % (m, num_of_successful), peptides)
                if not os.path.isfile("%s%s.pdb" % (peptides.base, num_of_successful)):
                    os.system("mv bak_%s%s.pdb %s%s.pdb" % (peptides.base, m, peptides.base, num_of_successful))
                else:
                    WriteLog("Trouble renaming structure %s to structure %s\n" % (m, num_of_successful), peptides)
                num_of_successful += 1
               
    
    WriteLog("Number of successfully optimized structures: %s\n" % (num_of_successful-1), peptides) # we incremented it once more at the end


    # renameing the failing ones too, continuing the numbering where we left off
    rename_num = num_of_successful
    for o in range(peptides.range_from, peptides.range_to):
        fn_ = "bak_%smin_%s.pdb" % (peptides.base, o)
        if peptides.gmxcheck == 1:
            if peptides.success[o-1] != [1,1,1]:
                WriteLog("Trying to rename structure %s to structure %s\n" % (o, rename_num), peptides)
                if os.path.isfile(fn_):
                    if not os.path.isfile("fail_%smin_%s.pdb" % (peptides.base, rename_num)):
                        os.system("mv bak_%smin_%s.pdb fail_%smin_%s.pdb" % (peptides.base, o, peptides.base, rename_num))
                if not os.path.isfile("fail_%s%s.pdb" % (peptides.base, rename_num)):
                    os.system("mv bak_%s%s.pdb fail_%s%s.pdb" % (peptides.base, o, peptides.base, rename_num))
                else:
                    WriteLog("Trouble renaming structure %s to structure %s\n" % (o, rename_num, peptides))
                rename_num += 1
        else: # if not gmx
            if peptides.success[o-1][0] != 1:
                WriteLog("Trying to rename structure %s to structure %s\n" % (o, rename_num), peptides)
                if not os.path.isfile("fail_%s%s.pdb" % (peptides.base, rename_num)):
                    os.system("mv bak_%s%s.pdb fail_%s%s.pdb" % (peptides.base, o, peptides.base, rename_num))
                rename_num += 1
    if rename_num!=peptides.range_to:
        WriteLog("Error! Not all structures were renamed.", peptides)


#######################################################################
def Cleanup(peptides):
    WriteLog("Cleaning up...\n", peptides)
    # deletes temporary files crated by the program
    #os.system("rm fetch*.dat")
    for n in range(peptides.range_from,peptides.range_to):
        os.system("rm contacts_%s%s.dat" % (peptides.base, n))
        os.system("rm contacts_%s%s.lst" % (peptides.base, n))
        if peptides.gmxcheck==1:
            # the Gromacs generated files
            os.system("rm check_%s.top" % (n))
            os.system("rm check_box_%s.gro" % (n))
            os.system("rm check_num_%s.dat" % (n))
            os.system("rm check_%s.gro" % (n))
            os.system("rm em2_%s.gro" % (n))
            os.system("rm em2_%s.trr" % (n))
            os.system("rm em2_%s.tpr" % (n))
            os.system("rm mdout.mdp")
            os.system("rm posre.itp")
            os.system("rm md_%s.log" % (n))
            os.system("rm ener_%s.edr" % (n))
            os.system("rm chimerax_clashcheck_%s.dat" % (n))
            os.system("rm clashcheck_%s.cxc" % (n))
        os.system ("rm rotate_%s.cxc" % (n))
        os.system("rm rotateProPhi_%s.cxc" % (n))
        os.system("rm unknot_%s.cxc" % (n))
        os.system("rm %s%s_*.pdb" % (peptides.base, n))
    
########################################################################
def ConcatEnsemble(peptides):
    ens_filname = "%sensemble_%s_%s.pdb" % (peptides.base, peptides.range_from, peptides.range_to)
    WriteLog("Concatenating ensemble of successfully built structures to %s\n" % ens_filname, peptides)
    with open(ens_filname, "w") as ens:
        for m in range(peptides.range_from,peptides.range_to):
            if peptides.gmxcheck == 1:
                fn_ = "%smin_%s.pdb" % (peptides.base, m)
            else:
                fn_ = "%s%s.pdb" % (peptides.base, m)
            if os.path.isfile(fn_):
                modelnum = m 
                found = 0
                with open(fn_, "r") as struc:
                    lines = struc.readlines()
                    for l_ in lines:
                        if l_[0:6]=="HEADER " or l_[0:6]=="TITLE " or l_[0:6]=="EXPDTA" or l_[0:6]=="AUTHOR" or l_[0:6]=="REMARK"  or l_[0:6]=="CRYST1":
                            ens.write(l_)
                        elif l_[0:6]=="MODEL ":
                            found = 1
                            if modelnum<10:
                                ens.write("MODEL        %s\n" % modelnum)
                            elif modelnum>10 and modelnum<100:
                                ens.write("MODEL       %s\n" % modelnum)
                            elif modelnum>100 and modelnum<1000:
                                ens.write("MODEL      %s\n" % modelnum)
                            elif modelnum>1000 and modelnum<10000:
                                ens.write("MODEL     %s\n" % modelnum)
                            else:
                                ens.write("MODEL    %s\n" % modelnum)
                        elif l_[0:6]=="ATOM  " or l_[0:6]=="HETATM":  
                            if found == 0: # not found a MODEL line before the first atom line
                                if modelnum<10:
                                    ens.write("MODEL        %s\n" % modelnum)
                                elif modelnum>10 and modelnum<100:
                                    ens.write("MODEL       %s\n" % modelnum)
                                elif modelnum>100 and modelnum<1000:
                                    ens.write("MODEL      %s\n" % modelnum)
                                elif modelnum>1000 and modelnum<10000:
                                    ens.write("MODEL     %s\n" % modelnum)
                                else:
                                    ens.write("MODEL    %s\n" % modelnum)
                                found = 1
                            ens.write(l_)
                        elif l_=="END\n":
                            ens.write("ENDMDL\n")
                        else:
                            ens.write(l_)

########################################################################

####################### M A I N

def Main():

    # declaring parameters
    numberofstructures = -99
    mode = ""
    base = ""
    sequence = ""
    dataset = ""
    cycle = -99
    keep = -99
    gmxcheck = -99
    pro = -99
    addrandom = -99
    angletoadd = -99
    unknotmax = -99
    threadnum = -99
    range_from = -99
    range_to = -99

    textForLater = [] # it will be later written to the logfile, but that cannot yet be opened (because it first needs all the input arguments to be set)

    # taking command line arguments 
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]
    unixOptions = "a:b:c:d:g:k:m:n:p:s:r:t:u:x:y:h"  
    gnuOptions = ["angletoadd=", "base=", "cycle=", "dataset=", "gmxcheck=", "keep=", "mode=", "numofstructures=", "proline=", "sequence=", "random=", "threads", "unknotmax", "xfrom=", "yto=", "help"]  

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        # output error, and return with an error code
        textForLater.append(str(err)+'\n')
        sys.exit(2)

    # evaluate given options for the inputs
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-a", "--angletoadd"):
            angletoadd = int(currentValue)
        elif currentArgument in ("-b", "--base"):
            base = currentValue
        elif currentArgument in ("-c", "--cycle"):
            cycle = int(currentValue)
        elif currentArgument in ("-d", "--dataset"):
            dataset = currentValue
        elif currentArgument in ("-g", "--gmxcheck"):
            gmxcheck = int(currentValue)
        elif currentArgument in ("-k", "--keep"):
            keep = int(currentValue)
        elif currentArgument in ("-m", "--mode"):
            mode = currentValue
        elif currentArgument in ("-n", "--numofstructures"):
            numberofstructures = int(currentValue)
        elif currentArgument in ("-p", "--proline"):
            pro = int(currentValue) 
        elif currentArgument in ("-s", "--sequence"):
            sequence = currentValue
        elif currentArgument in ("-r", "--random"):
            addrandom = int(currentValue)
        elif currentArgument in ("-t", "--threadnum"):
            addrandom = int(currentValue)
        elif currentArgument in ("-u", "--unknotmax"):
            unknotmax = int(currentValue)
        elif currentArgument in ("-x", "--xfrom"):
            range_from = int(currentValue)
        elif currentArgument in ("-y", "--yto"):
            range_to = int(currentValue)
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
            if p[0]=="#":
                continue # skip comment lines
            elif p[0]=="a" and angletoadd==-99:
                try:
                    angletoadd = int(p[1])
                except:
                    angletoadd = 30
                    print("trouble with %s<" % (p[1]))       
            elif p[0]=="b" and base=="":
                base = p[1]
            elif p[0]=="c" and cycle==-99:
                try:
                    cycle = int(p[1])
                except:
                    cycle = 10
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="d" and dataset=="":
                dataset = p[1]
            elif p[0]=="g" and gmxcheck==-99:
                try:
                    gmxcheck = int(p[1])
                except:
                    gmxcheck = 1
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="k" and keep==-99:
                try:
                    keep = int(p[1])
                except:
                    keep = 0
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="m" and mode=="":
                mode = p[1]
            elif  p[0]=="n" and numberofstructures==-99:
                try:
                    numberofstructures = int(p[1])
                except:
                    numberofstructures = 1
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="p" and pro==-99:
                try:
                    pro = int(p[1])
                except:
                    pro = 0
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="r" and addrandom==-99:
                try:
                    addrandom = int(p[1])
                except:
                    addrandom = 5
                    print("trouble with %s<" % (p[1]))
            elif p[0]=="t" and threadnum==-99:
                try:
                    threadnum = int(p[1])
                except:
                    threadnum = 1
                    print("trouble with %s<" % (p[1]))   
            elif p[0]=="u" and unknotmax==-99:
                try:
                    unknotmax = int(p[1])
                except:
                    unknotmax = 6
                    print("trouble with %s<" % (p[1]))       
            elif p[0]=="s" and sequence=="":
                sequence = p[1]
    
    #if there is no parallel
    if range_from == -99:
        range_from = 1
    if range_to == -99:
        range_to = numberofstructures+1

    MyPeptides = Peptides(sequence, base, mode, numberofstructures, dataset.upper(), cycle, keep, gmxcheck, pro, addrandom, angletoadd, unknotmax, threadnum, range_from, range_to)

    #textForLater.append("Command given: %s\n" % (" ".join(sys.argv)))

    # writing the things to the logfile
    for textLine in textForLater:
        WriteLog(textLine, MyPeptides)

    WriteLog("Given parameters:\n \tangle to add: %s\n \tbase: %s\n \tcycle: %s\n \tdataset: %s\n \tgromacs optimization: %s\n \tmode: %s\n \tnumber of structures: %s\n \tproline phi rotation: %s\n \tremaining temporary files: %s\n \trandom addition: %s\n \tunkot max: %s\n \tsequence: %s\n \tthreadnum: %s\n \trange from: %s\n \trange to: %s\n" % (MyPeptides.angletoadd, MyPeptides.base, MyPeptides.cycle, MyPeptides.dataset, MyPeptides.gmxcheck, MyPeptides.mode, MyPeptides.numberofstructures, MyPeptides.proline, MyPeptides.keep, MyPeptides.addrandom, MyPeptides.unknotmax, MyPeptides.sequence, MyPeptides.threadnum, MyPeptides.range_from, MyPeptides.range_to), MyPeptides)

    workingdirectory = os.getcwd()
    datapath = sys.path[0]+"/Data/"

    pathfname = sys.path[0]+"/Data/paths.in"

    with open(pathfname, "r") as paths: # please specify the paths in paths.txt
        for line in paths:
            line = line.strip()
            line_ = line.split()
            if line_[0][0]=="#":
                continue # skip comment lines
            if len(line_)>1:
                if line_[0] == "GromacsPath":
                    gromacspath = line_[1]
                elif line_[0] == "GromacsSuffix":
                    gromacssuffix = line_[1]
                elif line_[0] == "ChimeraXPath":
                    chimeraxpath = line_[1]
                elif line_[0] == "Scwrl4Path":
                    scwrl4path = line_[1]

    MyInstall = Install(gromacspath, gromacssuffix, datapath, workingdirectory, chimeraxpath, scwrl4path)

    #########
    Build(MyPeptides, MyInstall)
    if MyPeptides.mode == "WEIGHTED_LEFT" or MyPeptides.mode == "WEIGHTED_RIGHT" or MyPeptides.mode == "WEIGHTED_TRIPLET":
        WriteDistributions(MyPeptides, MyInstall)
    MyPeptides.init_success()
    #
    for i in range(MyPeptides.range_from, MyPeptides.range_to):
    #for i in range(1,MyPeptides.numberofstructures+1):
        Rotate(i, MyPeptides, MyInstall)
    for j in range(1,MyPeptides.cycle):
        WriteLog("Trying at %s\n" % (MyPeptides.currentcycle), MyPeptides)
        MyPeptides.plusplus_currentcycle()
        for i in range(MyPeptides.range_from, MyPeptides.range_to):
        #for i in range(1,MyPeptides.numberofstructures+1):
            if MyPeptides.success[i-1][2]!=1:
                Rotate(i, MyPeptides, MyInstall)
    Rename(MyPeptides)
    if MyPeptides.keep == 0:
        Cleanup(MyPeptides)

    ConcatEnsemble(MyPeptides)
    elapsed_time = time.time() - start_time
    WriteLog("Time elapsed: %.2f\n" % (elapsed_time), MyPeptides)

Main()
