#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>

using namespace std;

int main(int argc, char* argv[]) {

struct Line{
    string leftOrRight;
    string res;
    string neighbaa;
    signed int phi;
    signed int psi;
    long double prob;
    long double lnprob;
    long double cum;
    };

struct TripletBin{
    string leftaa;
    string middleaa;
    string rightaa;
    signed int phi;
    signed int psi;
    long double prob;
    long double lnprob;
    long double cum;
    };

    string firstAminoAcid = argv[1]; // the left one - three letter code!
    string secondAminoAcid = argv[2]; // the central one - three letter code!
    string thirdAminoAcid = argv[3]; // the right one - three letter code!

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

	Line* relevant_lines = new Line[20800]; // stores the lines about the triplet in the input data
    TripletBin* data = new TripletBin[5184]; // stores for each bin the derived probabilities

    long double Sum = 0;

    string path = "/home/zita/Store/GKAP/Dunbrack-Ting-derived-triplet-data/";
    //string ext_txt = ".txt";
    string ext_bin = ".bin";

    // for making a binary file as output
    string outfn = path+aa_names[secondAminoAcid]+"/"+aa_names[firstAminoAcid]+aa_names[secondAminoAcid]+aa_names[thirdAminoAcid]+ext_bin;
    ofstream out (outfn,ios::binary);
    if(!out) {
        cout << "Cannot open Output file!" << endl;
        return 1;
   }
    


    // This was the old version where I outputted in a text file
    //string outfnstr = path+aa_names[secondAminoAcid]+"/"+aa_names[firstAminoAcid]+aa_names[secondAminoAcid]+aa_names[thirdAminoAcid]+ext_txt;
    //ofstream outstr (outfnstr, ios::out);
    //if(!out){
        //cout << "Cannot open Output file!" << endl;
        //return 1;
    //}


    //ofstream temporaryFile ("temp.txt", ios::out);

    string ifilename_conly = "/home/zita/Research-new_sm_v3/GKAP/de-novo-3D-prediction-GKAP/Dunbrack-Ting-probabilities/NDRD_Conly.txt"; // for coil data
    ifstream infile_conly; 
    string line_conly;
    infile_conly.open(ifilename_conly);

    int counter = 0; // for adding the lines to the array

    while (std::getline(infile_conly, line_conly))
    {
        if (line_conly.size()>63 && line_conly.substr(0,1)!="#"){
            Line  l;
            l.res = line_conly.substr(0,3);
            l.leftOrRight = line_conly.substr(4,1); // it is enough to store only the first letter of "left" or "right"
            l.neighbaa = line_conly.substr(10,3);
            l.phi = stof(line_conly.substr(14,6));
            l.psi = stof(line_conly.substr(20,6));
            l.prob = stof(line_conly.substr(27,13));
            l.lnprob = stof(line_conly.substr(40,11));
            l.cum = stod(line_conly.substr(52,12));

            if (l.res==secondAminoAcid){ // we have found data about the middle amino acid
                if (l.leftOrRight[0]=='l'){
                    if (l.neighbaa=="ALL"){ // we also want to gather data about the sum of all left neighbours
                        //temporaryFile << l.res  << ";"  << l.leftOrRight[0] << ";" << l.neighbaa << ";" << to_string(l.phi) << ";" << to_string(l.psi) << ";" << to_string(l.prob) << ";" << to_string(l.lnprob) << ";" << to_string(l.cum) << endl;
                        relevant_lines[counter] = l;
                        counter++;
                    }
                    if (l.neighbaa==firstAminoAcid){ // if the left neighbour is the first one of the triplet
                        //temporaryFile << l.res  << ";"  << l.leftOrRight[0] << ";" << l.neighbaa << ";" << to_string(l.phi) << ";" << to_string(l.psi) << ";" << to_string(l.prob) << ";" << to_string(l.lnprob) << ";" << to_string(l.cum) << endl;
                        relevant_lines[counter] = l;
                        counter++;
                    }
                 } // end of left neighbours

                if (l.leftOrRight[0]=='r'){
                    if (l.neighbaa=="ALL"){ // we also want to gather data about the sum of all right neighbours
                        //temporaryFile << l.res  << ";"  << l.leftOrRight[0] << ";" << l.neighbaa << ";" << to_string(l.phi) << ";" << to_string(l.psi) << ";" << to_string(l.prob) << ";" << to_string(l.lnprob) << ";" << to_string(l.cum) << endl;
                        relevant_lines[counter] = l;
                        counter++;
                    }
                    if (l.neighbaa==thirdAminoAcid){ // if the right neighbour is the third one of the triplet
                        //temporaryFile << l.res  << ";"  << l.leftOrRight[0] << ";" << l.neighbaa << ";" << to_string(l.phi) << ";" << to_string(l.psi) << ";" << to_string(l.prob) << ";" << to_string(l.lnprob) << ";" << to_string(l.cum) << endl;
                        relevant_lines[counter] = l;
                        counter++;
                    }
                 } // end of right neighbours
            } // end if middle amino acid

            
            } // end if line length
        } // end while getline

    //temporaryFile.close();
    infile_conly.close();

    signed int count_phi;
    signed int count_psi;
    int count=0;
    for (count_phi = -180; count_phi<180; count_phi=count_phi+5){
        for (count_psi = -180; count_psi<180; count_psi = count_psi+5){ // iterating over each bin
            TripletBin Bin;
            Bin.leftaa = firstAminoAcid;
            Bin.middleaa = secondAminoAcid;
            Bin.rightaa = thirdAminoAcid;
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
        TripletBin b = data[count2];
        b.prob=b.prob/Sum;
        b.lnprob=abs(log(b.prob)); // log is the natural logarithm with base e, we are taking the absolute value, because otherwise it is negative!
        Cumulative = Cumulative+b.prob;
        b.cum = Cumulative;
        //outstr << b.leftaa << ";" << b.middleaa << ";" << b.rightaa << ";" << b.phi << ";" << b.psi << ";" << b.prob << ";" << b.lnprob << ";" << b.cum << endl; 
        out.write((char *) &b, sizeof(TripletBin));
        }

    out.close();

    // garbage collection
	delete[] relevant_lines;
    delete[] data;

    return 0;

}

