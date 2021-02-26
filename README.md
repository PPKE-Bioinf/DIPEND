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
Linux operating system
python 3.5
UCSF Chimerax 1.1
gromacs 2020 or 4.5


You can download Chimerax here:
https://www.rbvi.ucsf.edu/chimerax/download.html

You can donwnload GROMACS here:
https://manual.gromacs.org/documentation/

Example running of the program:
nohup python3 ./PeProb.py -b test_ -c 100 -d tcbig -g 1 -m LEFT -n 3 -p 0 -r 0 -s QGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHG > out.out 2>&1 &

Call --help for more info.

If you would like to compile the c++ executables from the sources instead of using the given executables, please use the following commands:
g++ fetch_angles.cpp -o fetch_angles -std=c++11
g++ fetch_both_angles.cpp -o fetch_both_angles -std=c++11
g++ checkcontacts.cpp -o checkcontacts -std=c++11

Please do not hesitate to contact me in case of any questions at the following e-mail address:
harmat.zita@itk.ppke.hu

