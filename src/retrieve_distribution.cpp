#include <iostream>
#include <string>
#include <fstream>
using namespace std;

struct Line{ // a structure to represent one bin
    char neighbour[5]; // left or right
    char firstaa[3];
    char secondaa[3];
    float phi;
    float psi;
    double cum; // cumulative sum
    };

    void write_line(Line * data, ofstream& outfilehandle) // writing out the found data
    {
        int c;
        for (c=0; c<sizeof((*data).firstaa); c++){
            outfilehandle << (*data).firstaa[c];}
        outfilehandle << " ";
        int d;
        for (d=0; d<sizeof((*data).neighbour); d++){
            outfilehandle << (*data).neighbour[d];}
        outfilehandle << " ";
        int e;
        for (e=0; e<sizeof((*data).secondaa); e++){
            outfilehandle << (*data).secondaa[e];}
        outfilehandle << " ";
        outfilehandle << (*data).phi;
        outfilehandle << " ";
        outfilehandle << (*data).psi;
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
