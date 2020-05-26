//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 2, O = 1;
char normalize;

//Description
string descr;
descr += "This makes the transformation matrix, T,\n";
descr += "to convert from STFT power on a linear frequency scale.\n";
descr += "to power at the cfs of a new frequency scale.\n";
descr += "\n";
descr += "This is done with triangle-shaped weighting functions in T.\n";
descr += "\n";
descr += "T has F columns, one for each positive STFT freq.\n";
descr += "T has B rows, one for each output cf.\n";
descr += "\n";
descr += "Inputs X1 (STFT freqs) and X2 (cfs) are in units of Hz.\n";
descr += "Output T has the same data type and file format as X1.\n";
descr += "\n";
descr += "Include -n (--normalize) to normalize the triangle weights to sum to 1.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ get_spectrogram_T_mat X1 X2 -o T \n";
descr += "$ get_spectrogram_T_mat X1 X2 > T \n";
descr += "$ cat X2 | get_spectrogram_T_mat X1 -n > T \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_lit  *a_nrm = arg_litn("n","normalize",0,1,"normalize output to sum to 1");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (T)");

//Get options

//Get normalize
normalize = (a_nrm->count>0);

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i2.N(); o1.C = i1.N();
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X1, *X2, *T;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { T = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (T)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (ov::get_spectrogram_T_mat_s(T,o1.iscolmajor(),X1,int(i1.N()),X2,int(i2.N()),normalize))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(T),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (T)" << endl; return 1; }
    }
    delete[] X1; delete [] X2; delete[] T;
}

//Finish

