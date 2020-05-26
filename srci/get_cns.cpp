//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 2, O = 1;
int *Y;

//Description
string descr;
descr += "Gets center \"numbers\" for a set of center frequencies\n";
descr += "within a set of STFT frequncies. That is, the output cns\n";
descr += "are the indices of the closest-fit STFT freq for each cf.\n";
descr += "\n";
descr += "Inputs X1 (STFT freqs) and X2 (cfs) are in units of Hz.\n";
descr += "\n";
descr += "Y has the same size and file format as X2 (cfs), but has int data type.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ get_cns X1 X2 -o Y \n";
descr += "$ get_cns X1 X2 > Y \n";
descr += "$ cat X2 | get_cns X1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a vector" << endl; return 1; }
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i2.F!=i1.F) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format" << endl; return 1; }

//Set output header info
o1.F = i2.F; o1.T = (sizeof(int)==4) ? 32 : 64;
o1.R = i2.R; o1.C = i2.C; o1.S = i2.S; o1.H = i2.H;

//Other prep
try { Y = new int[o1.N()]; }
catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }

//Process
if (i1.T==1)
{
    float *X1, *X2;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (ov::get_cns_s(Y,X1,int(i1.N()),X2,int(i2.N()))) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    delete[] X1; delete [] X2;
}

//Finish
if (wo1)
{
    try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
}
delete[] Y;

