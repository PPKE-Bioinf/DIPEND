Installation:
Please make sure the paths are set in the paths.in file in the Data folder!

Requirements:
Linux operating system, little-endian byte structure
python 3.5 or higher
UCSF Chimerax 1.1 or higher
gromacs 4.5 or higher

You can download Chimerax here:
https://www.rbvi.ucsf.edu/chimerax/download.html

You can download Scwrl4 here:
http://dunbrack.fccc.edu/SCWRL3.php/

You can donwnload GROMACS here:
https://manual.gromacs.org/documentation/

Example running of the program:
nohup python3 /PATH/DIPEND/Dipend.py -b test_ -c 100 -d tcbig -g 1 -m LEFT -n 3 -p 0 -r 0 -s QGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHG > out.out 2>&1 &
where /PATH/ is where you cloned this repository.

Call --help for more info.

Parameters can also be set in the parameters.in file in the Data folder. Parameters given on the command line are always prioritized over parameters in the parameter file.

When you use the WEIGHTED_LEFT, WEIGHTED_RIGHT or WEIGHTED_TRIPLET modes, you can set Gaussian  distributions for each amino acid, setting a phi psi angle as centers and a standard deviation and a weight between 0 and 1 to set how much it is considered relative to the Dunbrack distribution. If you set 0, only the Dunbrack distribution will be consiedered, if you set 1, only your custom distributin will be considered. The distributions are set in the distribution.in file in the Data folder. You can define amino acids line by line or using the '-' character to define amino acid ranges. You can also leave out amino acids. For them, only the Dunbrack distributions will be considered.

In the sequence, P means trans proline (which is most of the cases). If you would like to include cys proline specifically, please use the letter X in the sequence. Note, that cys proline will only be considered as central amino acid, not as neighbour, because the Dunbrack datasets do not regard such a case. Trans proline will be considered as neighbour.

If you would like to compile the c++ executables from the sources instead of using the given executables, please use the following commands:
g++ checkcontacts.cpp -o checkcontacts -std=c++11

The probabilities are stored in binary files. Angles are in a 4 byte signed int format, probabilities are in an 8 byte double precision format. The byte structure is little-endian. If you would like to make a non-binary, human readable copy of one of the Dunbrack database files, use retrieve_distribution for the left or right neighbours and retrieve_distribution_triplet for the triplet neighbours in the Generate_data directory. For the weighted user defined distributions, the retrieve_distribution_user or retrieve_distribution_user_bare file makes the conversion. Give them the full path and filename of the binary file by running any of these scripts.

Please do not hesitate to contact me in case of any questions at the following e-mail address:
harmat.zita@itk.ppke.hu

