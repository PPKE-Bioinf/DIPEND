#include <iostream>
#include <string>
#include <fstream>
#include <map>

using namespace std;

int main(int argc, char* argv[]) {

    map<string, string> aa_names;

    aa_names["ALA"]="A";
    aa_names["ARG"]="R";
    aa_names["ASN"]="N";
    aa_names["ASP"]="D";
    aa_names["CYS"]="C";
    aa_names["GLN"]="Q";
    aa_names["GLU"]="E";
    aa_names["GLY"]="G";
    aa_names["HIS"]="H";
    aa_names["ILE"]="I";
    aa_names["LEU"]="L";
    aa_names["LYS"]="K";
    aa_names["MET"]="M";
    aa_names["PHE"]="F";
    aa_names["PRO"]="P";
    aa_names["SER"]="S";
    aa_names["THR"]="T";
    aa_names["TRP"]="W";
    aa_names["TYR"]="Y";
    aa_names["VAL"]="V";
    aa_names["CPR"]="X"; // cys proline

struct Line{
    double cum;
    signed int phi;
    signed int psi;
    };

//dataset = "conly";
string dataset = argv[1];
//neighbourtype = "l";
string neighbourtype = argv[2];
//amino_acid = "ALA";
string amino_acid = argv[3];
//neighbour = "ALA";
string neighbour = argv[4];

    ofstream out ("/home/zita/Scripts-Research/DIPEND/Data/"+dataset+"-"+neighbourtype+"/"+aa_names[amino_acid]+"/"+neighbourtype+aa_names[neighbour]+".bin", ios::out |   ios::binary);
    if(!out) {
        cout << "Cannot open file!" << endl;
        return 1;
   }   

    string ifilename = "/home/zita/Research-new_sm_v3/GKAP/de-novo-3D-prediction-GKAP/Dunbrack-Ting-probabilities/NDRD_"+dataset+".txt";
    ifstream infile; 
    string line;
    infile.open(ifilename);
    while (std::getline(infile, line))
    {
        if (line.size()>63 && line.substr(0,1)!="#"){
            Line  l;
            string neigh = line.substr(4,5);
            string faa = line.substr(0,3);
            string saa = line.substr(10,3);
            l.phi = stof(line.substr(14,6));
            l.psi = stof(line.substr(20,6));
            l.cum = stod(line.substr(52,12));
            if (neigh[0]=='l'){
                if (faa[0]==amino_acid[0] && faa[1]==amino_acid[1] && faa[2]==amino_acid[2]){
                    if (saa[0]==neighbour[0] && saa[1]==neighbour[1] && saa[2]==neighbour[2]){
                        out.write((char *) &l, sizeof(Line));
                        } // end if neighbour
                    } // end if amino_acid
            } // end if left
            } // end if line length
        } // end while getline

    out.close();
    infile.close();

    return 0;

}

