#include <iostream>
#include <string>
#include <fstream>
#include <map>
using namespace std;

struct TripletBin{ // the size of it is 160
    double prob;
    double lnprob;
    double cum;
    signed int phi;
    signed int psi;
    };

    void write_line(TripletBin * data, ofstream& outfilehandle) // writing out the found data
    {
        outfilehandle << (*data).phi;
        outfilehandle << " ";
        outfilehandle << (*data).psi;
        outfilehandle << " ";
        outfilehandle << (*data).prob;
        outfilehandle << " ";
        outfilehandle << (*data).lnprob;
        outfilehandle << " ";
        outfilehandle << (*data).cum;
        outfilehandle << endl;
    }; // end write_line

int main(int argc, char* argv[]) {

    string filename = argv[1];

    ifstream infile(filename, ios::in | ios::binary);
    if(!infile) {
        cout << "Cannot open file!" << endl;
        return 1;
    }

    string ofname = "distribution.dat";

    ofstream outtext(ofname, ios::out);

    struct TripletBin* data;
    data = new TripletBin;

    int n = 0;
    while(!infile.eof()){
        infile.seekg(sizeof(TripletBin)*n); // going line by line
        infile.read((char *) data, sizeof (TripletBin));
        write_line(data, outtext); // if everything goes well, the second line is the data we seek
        n++;
    } // end while

    infile.close();
    outtext.close();
}
