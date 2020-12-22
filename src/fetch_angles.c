#include <iostream>
#include <string>
#include <fstream>

using namespace std;

// usage example: ./fetch_angles TRP r HIS 0.76

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
        //cout << "first amino acid: ";
        //int a;
        //for (a=0; a<sizeof((*data).firstaa); a++){cout << (*data).firstaa[a];}
        //cout << " neighbour: ";
        //int b;
        //for (b=0; b<sizeof((*data).neighbour); b++){cout << (*data).neighbour[b];}
        //cout << " second amino acid: ";
        int c;
        for (c=0; c<sizeof((*data).secondaa); c++){
            //cout << (*data).secondaa[c]; 
            outfilehandle << (*data).secondaa[c];}
        outfilehandle << ";";
        //cout << " phi: ";
        //cout << to_string((*data).phi);
        outfilehandle << to_string((*data).phi);
        outfilehandle << ";";
        //cout << " psi: ";
        //cout << to_string((*data).psi);
        outfilehandle << to_string((*data).psi);
        outfilehandle << ";";
        //cout << " cumul: ";
        //cout << to_string((*data).cum);
        outfilehandle << to_string((*data).cum);
        //cout << endl;
        outfilehandle << endl;
    }; // end write_line

int main(int argc, char* argv[]) {

    string aa = argv[1];
    string neighb = argv[2]; 
    string secondaa = argv[3];
    double rand = stod(argv[4]);
    string path = argv[5];
    string dataset = argv[6];
    string pwd = argv[7]; // the working directory

    string filename = "";
    // now we are looking for the corresponding files on my disc
    if (neighb == "r"){
        //cout << "right\n";
        filename = path+"right-"+dataset+'/'+aa+neighb+".bin"; // example: data/right-conly/GLUl.bin
    }
    if (neighb == "l"){
        //cout << "left\n";
        filename = path+"left-"+dataset+'/'+aa+neighb+".bin"; // example: data/left-conly/GLUl.bin
    }
    
    //cout << "filename: " << filename << endl;

    ifstream infile(filename, ios::in | ios::binary);
    if(!infile) {
        cout << "Cannot open file!" << endl;
        return 1;
    }

    ofstream outtext(pwd+"/fetched_one_neighbour.dat", ios::out);

    struct Line* data;
    data = new Line;
    int n = 0;
    int found = 0;
    while(found==0){
        infile.seekg(sizeof(Line)*n); // going line by line
        infile.read((char *) data, sizeof (Line));
        if((*data).secondaa[0]==secondaa[0] && (*data).secondaa[1]==secondaa[1] && (*data).secondaa[2]==secondaa[2] && (*data).cum>rand){ // if I find the right amino acids and a cumulative sum greater than the drawn random value
            found = 1;
            write_line(data, outtext);
            } // end if random
        if(n>1000000) // this is only a precaution in order not to run into an endless cycle or a very long running time
            {
            found = 1;
            write_line(data, outtext);
            cout << "Uh-oh, the search took too long, getting random value for phi psi!\n";
}
        n = n+1;
    } // end while found==0

    infile.close();
    outtext.close();
}
