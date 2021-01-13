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

    ofstream serrs ("./conly-right/SERr.bin", ios::out |   ios::binary);
    if(!serrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream thrrs ("./conly-right/THRr.bin", ios::out |   ios::binary);
    if(!thrrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream trprs ("./conly-right/TRPr.bin", ios::out |   ios::binary);
    if(!trprs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream tyrrs ("./conly-right/TYRr.bin", ios::out |   ios::binary);
    if(!tyrrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream valrs ("./conly-right/VALr.bin", ios::out |   ios::binary);
    if(!valrs) {
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
            //cout << "neighb: ";
            //int p;
            //for (p=0; p<sizeof(l.neighbour); p++){cout << l.neighbour[p];}
            //cout << " faa: ";
            //int m;
            //for (m=0; m<sizeof(l.firstaa); m++){cout << l.firstaa[m];}
            //cout << " saa: ";
            //int o;
            //for (o=0; o<sizeof(l.secondaa); o++){cout << l.secondaa[o];}
            //cout << " phi: " << to_string(l.phi) << " psi: " << to_string(l.psi) << " cum: " << to_string(l.cum) << endl;

            if (l.neighbour[0]=='r' && l.neighbour[1]=='i' && l.neighbour[2]=='g' && l.neighbour[3]=='h' && l.neighbour[4]=='t'){
                if (l.firstaa[0]=='S' && l.firstaa[1]=='E' && l.firstaa[2]=='R'){serrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='T' && l.firstaa[1]=='H' && l.firstaa[2]=='R'){thrrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='T' && l.firstaa[1]=='R' && l.firstaa[2]=='P'){trprs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='T' && l.firstaa[1]=='Y' && l.firstaa[2]=='R'){tyrrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='V' && l.firstaa[1]=='A' && l.firstaa[2]=='L'){valrs.write((char *) &l, sizeof(Line));}
            } // end if right
            
            

            } // end if line length
        } // end while getline

    serrs.close();

    thrrs.close();

    trprs.close();

    tyrrs.close();

    valrs.close();

    infile_conly.close();

    return 0;

}

