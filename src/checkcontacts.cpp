#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

// this function checks whether any non-neighbour amino acid CA atoms are closer than a thershold (presence of a "knot") or whether amino acid CA atoms any non-neighbor but less than 4 amino acid apart CA atoms are further than a threshold (presence of a "tear")

int main(int argc, char* argv[]) {
    string i = argv[1];
    string base = argv[2]; 
    string pwd = argv[3]; // the working directory

    int value = 3;
    int cutoff = 4;
    int upperCutoff = 20;


    ofstream datfile,lstfile;
    string datfilename = pwd+"/contacts_"+base+i+".dat";
    string lstfilename = pwd+"/contacts_"+base+i+".lst"; 
    datfile.open(datfilename);
    lstfile.open(lstfilename); 

    string line;
    ifstream infile; 
    string ifilename = pwd+"/"+base+i+".pdb";
    infile.open(ifilename); 
    int n = 0;
    while (std::getline(infile, line))
{
    if (line.substr(0,4)=="ATOM" && line.substr(13,2)=="CA"){
            n = n+1;
    }
}
    infile.close();

    string line2;
    char residuepairs[50]; // GZ
    ifstream infile2;
    infile2.open(ifilename);

    vector< vector<double> > CA;
    CA.resize(n);
    for(int i = 0 ; i < n ; ++i)
    {
        CA[i].resize(3);
    }
    int tear = 0;
    int knot = 0;
    vector< vector<double> > distMatrix;
    distMatrix.resize(n);
    for(int i = 0; i < n ; ++i)
    {
        distMatrix[i].resize(n);
    }
    int m = 0;
    while (std::getline(infile2, line))
    {
    if (line.substr(0,4)=="ATOM" && line.substr(13,2)=="CA"){
            //CA[m][0]=stod(line.substr(31,6));
            //CA[m][1]=stod(line.substr(39,6));
            //CA[m][2]=stod(line.substr(47,6));
            CA[m][0]=stod(line.substr(30,8));
            CA[m][1]=stod(line.substr(38,8));
            CA[m][2]=stod(line.substr(46,8));
            m = m+1;
        }
    }

    infile2.close();

    for (int a = 0; a < n ; a++ ){
        for (int b = 0; b<a-1; b++){
            double xa = CA[a][0];
            double ya = CA[a][1];
            double za = CA[a][2];
            double xb = CA[b][0];
            double yb = CA[b][1];
            double zb = CA[b][2];
            double dx = abs(xa-xb);
            double dy = abs(ya-yb);
            double dz = abs(za-zb);
            double d = sqrt(dx*dx+dy*dy+dz*dz);
            distMatrix[a][b] = d;
            if (d<=cutoff && d!=0){ // this would mean the peptide chain twisting on itself in a "knot", a problem
                knot = 1;
                value = 1;
            	sprintf(residuepairs,"%d %d %7.3f\n",b,a,d); // GZ
		lstfile << residuepairs;
            }
            if (a-b<4 && d>upperCutoff && CA[a][0]!=0){ // this would mean a discontinuation of the peptide chain, a serious problem
                tear = 1;
                value = 2;
            }
            if (knot == 0 and tear == 0){ // # hopefully everything is fine
                value = 0;
            }
        }
    }

    datfile << value << endl;

    datfile.close();
    lstfile.close();

    return 0;

}
