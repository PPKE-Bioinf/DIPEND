#include <iostream>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {

struct Line{
    char neighbour[5];
    char firstaa[3];
    char secondaa[3];
    float phi;
    float psi;
    double cum;
    };

    ofstream alals ("./conly-left/ALAl.bin", ios::out |   ios::binary);
    if(!alals) {
        cout << "Cannot open file!" << endl;
        return 1;
   }
    
    ofstream argls ("./conly-left/ARGl.bin", ios::out |   ios::binary);
    if(!argls) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream asnls ("./conly-left/ASNl.bin", ios::out |   ios::binary);
    if(!asnls) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream aspls ("./conly-left/ASPl.bin", ios::out |   ios::binary);
    if(!aspls) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream cysls ("./conly-left/CYSl.bin", ios::out |   ios::binary);
    if(!cysls) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    string ifilename_conly = "/home/zita/Research-new_sm_v3/GKAP/de-novo-3D-prediction-GKAP-rat/Dunbrack-Ting-probabilities/NDRD_Conly.txt";
    ifstream infile_conly; 
    string line_conly;
    infile_conly.open(ifilename_conly);
    while (std::getline(infile_conly, line_conly))
    {
        if (line_conly.size()>63 && line_conly.substr(0,1)!="#"){
            Line  l;
            string neigh = line_conly.substr(4,5);
            int i;
            for (i=0; i<neigh.length(); i++){
                l.neighbour[i]=neigh[i];
            }
            string faa = line_conly.substr(0,3);
            int j;
            for (j=0; j<faa.length(); j++){
                l.firstaa[j]=faa[j];
            }
            string saa = line_conly.substr(10,3);
            int k;
            for (k=0; k<saa.length(); k++){
                l.secondaa[k]=saa[k];
            }
            l.phi = stof(line_conly.substr(14,6));
            l.psi = stof(line_conly.substr(20,6));
            l.cum = stod(line_conly.substr(52,12));
            if (l.neighbour[0]=='l' && l.neighbour[1]=='e' && l.neighbour[2]=='f' && l.neighbour[3]=='t'){
                if (l.firstaa[0]=='A' && l.firstaa[1]=='L' && l.firstaa[2]=='A'){alals.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='A' && l.firstaa[1]=='R' && l.firstaa[2]=='G'){argls.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='A' && l.firstaa[1]=='S' && l.firstaa[2]=='N'){asnls.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='A' && l.firstaa[1]=='S' && l.firstaa[2]=='P'){aspls.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='C' && l.firstaa[1]=='Y' && l.firstaa[2]=='S'){cysls.write((char *) &l, sizeof(Line));}
            } // end if left
            } // end if line length
        } // end while getline

    alals.close();
    argls.close();
    asnls.close();
    aspls.close();
    cysls.close();

    infile_conly.close();

    return 0;

}

