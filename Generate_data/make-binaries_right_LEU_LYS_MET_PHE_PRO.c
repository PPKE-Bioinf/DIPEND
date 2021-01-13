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

    ofstream leurs ("./conly-right/LEUr.bin", ios::out |   ios::binary);
    if(!leurs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream lysrs ("./conly-right/LYSr.bin", ios::out |   ios::binary);
    if(!lysrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream metrs ("./conly-right/METr.bin", ios::out |   ios::binary);
    if(!metrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream phers ("./conly-right/PHEr.bin", ios::out |   ios::binary);
    if(!phers) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream prors ("./conly-right/PROr.bin", ios::out |   ios::binary);
    if(!prors) {
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

            if (l.neighbour[0]=='r' && l.neighbour[1]=='i' && l.neighbour[2]=='g' && l.neighbour[3]=='h' && l.neighbour[4]=='t'){
                if (l.firstaa[0]=='L' && l.firstaa[1]=='E' && l.firstaa[2]=='U'){leurs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='L' && l.firstaa[1]=='Y' && l.firstaa[2]=='S'){lysrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='M' && l.firstaa[1]=='E' && l.firstaa[2]=='T'){metrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='P' && l.firstaa[1]=='H' && l.firstaa[2]=='E'){phers.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='P' && l.firstaa[1]=='R' && l.firstaa[2]=='O'){prors.write((char *) &l, sizeof(Line));}
            } // end if right
            
            

            } // end if line length
        } // end while getline

    leurs.close();

    lysrs.close();

    metrs.close();
    
    phers.close();

    prors.close();

    infile_conly.close();

    return 0;

}

