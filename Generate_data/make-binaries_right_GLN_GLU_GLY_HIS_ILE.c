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

    ofstream glnrs ("./conly-right/GLNr.bin", ios::out |   ios::binary);
    if(!glnrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream glurs ("./conly-right/GLUr.bin", ios::out |   ios::binary);
    if(!glurs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream glyrs ("./conly-right/GLYr.bin", ios::out |   ios::binary);
    if(!glyrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream hisrs ("./conly-right/HISr.bin", ios::out |   ios::binary);
    if(!hisrs) {
        cout << "Cannot open file!" << endl;
        return 1;
   }

    ofstream ilers ("./conly-right/ILEr.bin", ios::out |   ios::binary);
    if(!ilers) {
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
                if (l.firstaa[0]=='G' && l.firstaa[1]=='L' && l.firstaa[2]=='N'){glnrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='G' && l.firstaa[1]=='L' && l.firstaa[2]=='U'){glurs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='G' && l.firstaa[1]=='L' && l.firstaa[2]=='Y'){glyrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='H' && l.firstaa[1]=='I' && l.firstaa[2]=='S'){hisrs.write((char *) &l, sizeof(Line));}
                if (l.firstaa[0]=='I' && l.firstaa[1]=='L' && l.firstaa[2]=='E'){ilers.write((char *) &l, sizeof(Line));}
            } // end if right
            
            

            } // end if line length
        } // end while getline

    glnrs.close();

    glurs.close();

    glyrs.close();

    hisrs.close();

    ilers.close();

    infile_conly.close();

    return 0;

}

