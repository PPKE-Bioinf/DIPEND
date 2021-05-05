Installation:
Please make sure the paths are set!
Modify them at paths.txt as seen in that file!

Requirements:
Linux operating system
python 3.5 or higher
UCSF Chimerax 1.1 or higher
gromacs 4.5 or higher

You can download Chimerax here:
https://www.rbvi.ucsf.edu/chimerax/download.html

You can donwnload GROMACS here:
https://manual.gromacs.org/documentation/

Example running of the program:
nohup python3 /PATH/DIPEND/Dipend.py -b test_ -c 100 -d tcbig -g 1 -m LEFT -n 3 -p 0 -r 0 -s QGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHG > out.out 2>&1 &
where /PATH/ is where you cloned this repository.

Call --help for more info.

If you would like to compile the c++ executables from the sources instead of using the given executables, please use the following commands:
g++ fetch_angles.cpp -o fetch_angles -std=c++11
g++ fetch_both_angles.cpp -o fetch_both_angles -std=c++11
g++ checkcontacts.cpp -o checkcontacts -std=c++11

If you would like to make a non-binary, human readable copy of one of the database files, use retrieve_distribution for the left or right neighbours and retrieve_distribution_triplet for the triplet neighbours. Set the paths accordingly. They are located in the src directory.

Please do not hesitate to contact me in case of any questions at the following e-mail address:
harmat.zita@itk.ppke.hu

