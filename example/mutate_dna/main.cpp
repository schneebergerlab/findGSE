// written by Hequan Sun, Max Planck Institude for Plant Breeding Research.
#include <iostream>
#include  <fstream>
#include  <sstream>
#include      <map>
#include   <string>
#include <stdlib.h>
using namespace std;
//
char mutate_subseq(char orich);
//
bool verbose = false;
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\nFunction: mutate dna sequences in fasta with specified mutation rate (only point mutations)";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\nUsage: mutate_dna seq.fasta 0.001 outfile_prefix" << endl;
        cout << "\nPls keep order of input as: sequence_file mutation_rate, output_file_prefix." << endl;
        cout << "\nOutput will be outfile_prefix.fasta." << endl << endl;
        return 1;
    }
    cout << "\nInput DNA sequence is from file: " << argv[1] << "." << endl;
    double pmrate = atof(argv[2]);
    if(pmrate > 1.0)
    {
        cout << "Error. Mutation rate must be less than 1. Exited. " << endl;
        return 1; 
    }
    cout << pmrate*100 << "% bases of each original sequence have been asked to be mutated. " << endl;
    // open seq file
    std::ifstream fp(argv[1]);
    if(!fp.is_open())
    {
        cout << "Cannot open dna seq file " << argv[1] << ". Exited." << endl;
        return 1;
    }
    // initialize a file where mutated sequences will be recorded.
    fstream outfp;
    std::stringstream ss("");
    ss << argv[3] << ".rate" << pmrate << ".fasta";
    outfp.open((ss.str()).c_str(), ios::out);
    if(!outfp.is_open())
    {
        cout << "Cannot open file " << ss.str() << " to write mutated dna info. Exited.\n";
        return 1;
    }
    cout << "\nMutated DNA sequence will be collected in file: " << ss.str() << "." << endl;
    // read sequence and mutate
    std::string line("");
    getline(fp, line);
    while(fp.good())
    {
        if(line.size()==0 || line.find("#")!=std::string::npos)
        {
            getline(fp, line);
            continue;
        }
        
        if(line[0]=='>')
        {
            string seq("");
            string chrseq("");   
            while(fp.good())
            {             
                getline(fp, chrseq);
                if(chrseq.find(">")!=std::string::npos) break; // next seq
                seq    += chrseq;
            }
            
            long mutnum = (long)seq.length()*pmrate;
            int window = (int)1.0/pmrate;
            
            for(size_t i=0; i<seq.length(); i++)
            {
                int pos = 0;
                if(window*(i+1)-1 < seq.length()-1)
                {
                    pos     = (int)rand()%window;
                }
                else
                {
                    pos     = (int)rand()%(seq.length()-window*i);
                }
                //cout << " pos set as " << pos << endl;
                int  oripos = pos + window*i - 1;
                char orich  = seq[oripos];
                char ch     = mutate_subseq(orich);
                
                //cout << "i=" << i << ", base "     << seq[oripos] << " will be changed as ";
                seq[oripos] = ch;
                //cout << seq[oripos] << endl;
                if(window*(i+1)-1 >= seq.length()-1) 
                {
                    break;
                    //cout << " break here at i=" << i << endl;
                }
            }
            outfp << line << endl << seq << endl;
            cout << "Length of original sequence " << line.substr(1) << " is " << seq.length() << ", of which "
                 << mutnum << " have been point-mutated."    << endl;
            line.clear();
            line += chrseq;
        }
        else
        {
            getline(fp, line);
        }
    }
    
    fp.close();
    outfp.close();
}

char mutate_subseq(char orich)
{
    char ch  = 'N';
    
    if(orich == 'A')
    {
        int num = (int)rand()%30;
        if(num <= 16)
        {
            ch = 'G';
        }
        else
        if(num >= 24)
        {
            ch = 'C';
        }
        else
        if(num>=17 && num<=23)
        {
            ch = 'T';
        }
        else ;
    }
    else
    if(orich == 'C')
    {
        int num = (int)rand()%30;
        if(num <= 16)
        {
            ch = 'T';
        }
        else
        if(num >= 24)
        {
            ch = 'A';
        }
        else
        if(num>=17 && num<=23)
        {
            ch = 'G';
        }
        else ;
    }
    else
    if(orich == 'G')
    {
        int num = (int)rand()%30;
        if(num <= 16)
        {
            ch = 'A';
        }
        else
        if(num >= 24)
        {
            ch = 'C';
        }
        else
        if(num>=17 && num<=23)
        {
            ch = 'T';
        }
        else ;
    }
    else
    if(orich == 'T')
    {
        int num = (int)rand()%30;
        if(num <= 16)
        {
            ch = 'C';
        }
        else
        if(num >= 24)
        {
            ch = 'G';
        }
        else
        if(num>=17 && num<=23)
        {
            ch = 'A';
        }
        else ;
    }
    else ;
    
    return ch;    
}
