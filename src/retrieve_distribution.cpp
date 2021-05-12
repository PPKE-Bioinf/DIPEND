#include <iostream>
#include <string>
#include <fstream>
using namespace std;

struct Line{ // a structure to represent one bin
    double cum; // cumulative sum
    signed int phi;
    signed int psi;

    };

    void write_line(Line * data, ofstream& outfilehandle) // writing out the found data
    {
        int c;
        outfilehandle << (*data).cum;
        outfilehandle << " ";
        outfilehandle << (*data).phi;
        outfilehandle << " ";
        outfilehandle << (*data).psi;
        outfilehandle << endl;
    }; // end write_line

int main(int argc, char* argv[]) {

    string filename = argv[1];

    ifstream infile(filename, ios::in | ios::binary);
    if(!infile) {
        cout << "Cannot open file!" << endl;
        return 1;
    }

    ofstream outtext("distribution.dat", ios::out);

    struct Line* data;
    data = new Line;
    int n = 0;
    while(!infile.eof()){
        infile.seekg(sizeof(Line)*n); // going line by line
        infile.read((char *) data, sizeof (Line));  
        write_line(data, outtext);
        n++;
        
    } // end while

    infile.close();
    outtext.close();
}
