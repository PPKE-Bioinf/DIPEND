#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>

using namespace std;

int main(int argc, char* argv[]) {

struct Line{
    double prob;
    double lnprob;
    double cum;
    signed int phi;
    signed int psi;
    string leftOrRight;
    string res;
    string neighbaa;
    };

struct Write{
    double prob;
    double lnprob;
    double cum;
    signed int phi;
    signed int psi;
    };

    string firstAminoAcid = argv[1]; // the left one - three letter code!
    string secondAminoAcid = argv[2]; // the central one - three letter code!
    string thirdAminoAcid = argv[3]; // the right one - three letter code!
    string dataset = argv[4];

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

	Line* relevant_lines = new Line[20800]; // stores the lines about the triplet in the input data
    Write* data = new Write[5184]; // stores for each bin the derived probabilities

    double Sum = 0;

    string ifilename = "/home/zita/Research-new_sm_v3/GKAP/de-novo-3D-prediction-GKAP/Dunbrack-Ting-probabilities/NDRD_"+dataset+".txt";

    // for making a binary file as output
    string outfn = "/home/zita/Scripts-Research/DIPEND/Data/"+dataset+"/"+aa_names[secondAminoAcid]+"/"+aa_names[firstAminoAcid]+aa_names[secondAminoAcid]+aa_names[thirdAminoAcid]+".bin";
    ofstream out (outfn,ios::binary);
    if(!out) {
        cout << "Cannot open Output file!" << endl;
        return 1;
   }

    ifstream infile; 
    string line;
    infile.open(ifilename);

    int counter = 0; // for adding the lines to the array

    while (std::getline(infile, line))
    {
        if (line.size()>63 && line.substr(0,1)!="#"){
            Line  l;
            l.res = line.substr(0,3);
            l.leftOrRight = line.substr(4,1); // it is enough to store only the first letter of "left" or "right"
            l.neighbaa = line.substr(10,3);
            l.phi = stof(line.substr(14,6));
            l.psi = stof(line.substr(20,6));
            l.prob = stof(line.substr(27,13));
            l.lnprob = stof(line.substr(40,11));
            l.cum = stod(line.substr(52,12));

            if (l.res==secondAminoAcid){ // we have found data about the middle amino acid
                if (l.leftOrRight[0]=='l'){
                    if (l.neighbaa=="ALL"){ // we also want to gather data about the sum of all left neighbours
                        relevant_lines[counter] = l;
                        counter++;
                    }
                    if (l.neighbaa==firstAminoAcid){ // if the left neighbour is the first one of the triplet
                        relevant_lines[counter] = l;
                        counter++;
                    }
                 } // end of left neighbours

                if (l.leftOrRight[0]=='r'){
                    if (l.neighbaa=="ALL"){ // we also want to gather data about the sum of all right neighbours
                        relevant_lines[counter] = l;
                        counter++;
                    }
                    if (l.neighbaa==thirdAminoAcid){ // if the right neighbour is the third one of the triplet
                        relevant_lines[counter] = l;
                        counter++;
                    }
                 } // end of right neighbours
            } // end if middle amino acid

            
            } // end if line length
        } // end while getline

    infile.close();

    signed int count_phi;
    signed int count_psi;
    int count=0;
    for (count_phi = -180; count_phi<180; count_phi=count_phi+5){
        for (count_psi = -180; count_psi<180; count_psi = count_psi+5){ // iterating over each bin
            Write Bin;
            Bin.lnprob = 0.0;
            Bin.prob = 0.0;
            Bin.cum = 0.0;
            for(int i=0; i<20800; i++){
                if (relevant_lines[i].phi==count_phi and relevant_lines[i].psi==count_psi){
                    Bin.phi = count_phi;
                    Bin.psi = count_psi;
                    if (relevant_lines[i].neighbaa=="ALL"){
                        Bin.lnprob = Bin.lnprob-relevant_lines[i].lnprob;
                        }
                    else{
                        Bin.lnprob = Bin.lnprob+relevant_lines[i].lnprob;
                        }
                   }
                }
            Bin.prob = exp(-Bin.lnprob); // Reversing the natural logarithm, mindig the minus sign
            Sum=Sum+Bin.prob;
            data[count]=Bin;
            count++;
            } // end psi
        } // end phi

    //cout << "Sum: " << Sum << endl;

    long double Cumulative = 0;
    for (int count2 = 0; count2<5184; count2++){
        Write b = data[count2];
        b.prob=b.prob/Sum;
        b.lnprob=abs(log(b.prob)); // log is the natural logarithm with base e, we are taking the absolute value, because otherwise it is negative!
        Cumulative = Cumulative+b.prob;
        b.cum = Cumulative;
        out.write((char *) &b, sizeof(Write));
        }

    out.close();

    // garbage collection
	delete[] relevant_lines;
    delete[] data;

    return 0;

}

