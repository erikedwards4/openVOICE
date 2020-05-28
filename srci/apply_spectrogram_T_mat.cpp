//Includes
#include "apply_spectrogram_T_mat.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 2, O = 1;

//Description
string descr;
descr += "This applies the transformation matrix, T, to X,\n";
descr += "such that Y = T*X, or Y = X*T' (depending on orientation).\n";
descr += "\n";
descr += "T must be BxF, as from get_spectrogram_T_mat.\n";
descr += "\n";
descr += "This is used to convert from STFT power on a linear\n";
descr += "frequency scale to power on a new frequency scale.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ apply_spectrogram_T_mat X T -o Y \n";
descr += "$ apply_spectrogram_T_mat X T > Y \n";
descr += "$ cat X | apply_spectrogram_T_mat - T > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,T)");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (T) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a matrix" << endl; return 1; }
if (!i2.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (T) must be a matrix" << endl; return 1; }
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.iscolmajor()!=i2.iscolmajor()) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format" << endl; return 1; }
if (i1.R!=i2.C && i1.C!=i2.C) { cerr << progstr+": " << __LINE__ << errstr << "inputs must be size compatible" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (i1.R==i2.C) ? i2.R : i1.R;
o1.C = (i1.R==i2.C) ? i1.C : i2.R;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *T, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { T = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (T)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(T),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (T)" << endl; return 1; }
    if (ov::apply_spectrogram_T_mat_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),T,int(i2.R),int(i2.C)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete [] T; delete[] Y;
}

//Finish

