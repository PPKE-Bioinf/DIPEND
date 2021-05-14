#include <iostream>
#include <string>
#include <fstream>
using namespace std;

struct Line{ // a structure to represent one bin
    double cum; // cumulative sum
    signed int phi;
    signed int psi;

    };

    void write_line(Line * data) // writing out the found data
    {
        cout << "Data: ";
        cout << sizeof((*data));
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

    struct Line* data;
    data = new Line;
    infile.seekg(sizeof(Line)); // going line by line
    infile.read((char *) data, sizeof (Line));
    write_line(data);

    infile.close();
}
