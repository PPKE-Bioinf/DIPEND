#include <iostream>
#include <string>
#include <fstream>
#include <map>
using namespace std;

struct Data{ // the size of it is 160
    double cum;
    double prob;
    signed int phi;
    signed int psi;
    };

    void write_line(Data * data, ofstream& outfilehandle) // writing out the found data
    {
        //outfilehandle << "Phi: ";
        outfilehandle << (*data).phi;
        //outfilehandle << " Psi: ";
        outfilehandle << "\t";
        outfilehandle << (*data).psi;
        //outfilehandle << " Prob: ";
        outfilehandle << "\t";
        outfilehandle << (*data).prob;
        //outfilehandle << " Cum: ";
        outfilehandle << "\t";
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

    struct Data* data;
    data = new Data;

    int n = 0;
    while(!infile.eof()){
        infile.seekg(sizeof(Data)*n); // going line by line
        infile.read((char *) data, sizeof (Data));
        write_line(data, outtext); // if everything goes well, the second line is the data we seek
        n++;
    } // end while

    infile.close();
    outtext.close();
}
