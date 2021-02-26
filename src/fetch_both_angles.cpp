#include <iostream>
#include <string>
#include <fstream>
#include <map>

using namespace std;

struct TripletBin{ // the size of it is 160
    signed int phi;
    signed int psi;
    long double prob;
    long double lnprob;
    long double cum;
    };

    void write_line(TripletBin * data, ofstream& outfilehandle) // writing out the found data
    {
        outfilehandle << (*data).phi;
        outfilehandle << ";";
        outfilehandle << (*data).psi;
        outfilehandle << ";";
        outfilehandle << (*data).cum;
        outfilehandle << endl;
    }; // end write_line

int main(int argc, char* argv[]) {

    string firstaa = argv[1];
    string secondaa = argv[2];
    string thirdaa = argv[3];
    double rand = stod(argv[4]);
    string path = argv[5];
    string dataset = argv[6];
    string pwd = argv[7]; // the working directory

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

    string filename = "";
    // now we are looking for the corresponding files on my disc
    filename = path+dataset+'/'+aa_names[secondaa]+"/"+aa_names[firstaa]+aa_names[secondaa]+aa_names[thirdaa]+".bin"; // example: data/conly/G/HGI.bin

    ifstream infile(filename, ios::in | ios::binary);
    if(!infile) {
        cout << "Cannot open file!" << endl;
        return 1;
    }

    string ofname = pwd+"/fetched_both.dat";

    ofstream outtext(ofname, ios::out);

    struct TripletBin* data;
    data = new TripletBin;

    infile.seekg(sizeof(TripletBin)*5183);
    infile.read((char *) data, sizeof (TripletBin)); // read out the very last line
    write_line(data, outtext); // write it out, it will be the default output in every case as first line

    for(int n=0; n<5184; n++){
        infile.seekg(sizeof(TripletBin)*n); // going line by line
        infile.read((char *) data, sizeof (TripletBin));
        if((*data).cum>rand){ // if the cumulative sum is greater than the drawn random value
            write_line(data, outtext); // if everything goes well, the second line is the data we seek
            break;
            } // end if random
    } // end for loop

    infile.close();
    outtext.close();
}
