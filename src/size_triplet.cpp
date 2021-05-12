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

    void write_line(TripletBin * data) // writing out the found data
    {
        cout << "Data: ";
        cout << sizeof((*data));
        cout << " prob: ";
        cout << sizeof((*data).prob);
        cout << " lnprob: ";
        cout << sizeof((*data).lnprob);
        cout << " cum: ";
        cout << sizeof((*data).cum);
        cout << " phi: ";
        cout << sizeof((*data).phi);
        cout << " psi: ";
        cout << sizeof((*data).psi);
        cout << endl;
    }; // end write_line

int main(int argc, char* argv[]) {

    string filename = argv[1];

    ifstream infile(filename, ios::in | ios::binary);
    if(!infile) {
        cout << "Cannot open file!" << endl;
        return 1;
    }

    struct TripletBin* data;
    data = new TripletBin;

    infile.seekg(sizeof(TripletBin)); // going line by line
    infile.read((char *) data, sizeof (TripletBin));
    write_line(data); // if everything goes well, the second line is the data we seek

    infile.close();
}
