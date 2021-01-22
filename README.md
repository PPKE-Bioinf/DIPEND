Installation:
Please make sure the paths are set!
Modify them at PeProb.py from line 505, for example:
    GromacsPath = "/path/gromacs/bin/"
    GromacsSuffix = "_mpi"
    LsqmanPath = "/home/user/Lsqman/bin32/lsqman"
    Scwrl4Path = "/home/user/scwrl4/Scwrl4"
    DataPath = "/home/user/PeProb/Data/"
    ChimeraXPath = "/usr/bin/" 

Requirements:
python 3.5
UCSF Chimerax 1.1
gromacs 2020 or 4.5

Example running of the program:
nohup python3 ./PeProb.py -b test_extended_refactored_ -c 100 -d tcbig -g 1 -m LEFT -n 3 -p 1 -r 0 -s QGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHG > out.out 2>&1 &

Call --help for more info.
