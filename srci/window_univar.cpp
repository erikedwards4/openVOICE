//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 2, O = 1;
int dim;
double fs, fr, stp;
int c0, T;
char mn0;

//Description
string descr;
descr += "Takes univariate X and produces a series of windowed frames.\n";
descr += "A window W of length L must also be input (e.g. from the output of hamming).\n";
descr += "\n";
descr += "The output Y has size LxT or TxL, where L is the length of each frame, \n";
descr += "and T is the total number of frames (output time points).\n";
descr += "\n";
descr += "Use -r (--fr) to give the frame rate in Hz [default=100.0].\n";
descr += "This actually has to be in the same units as fs (sample rate),\n";
descr += "but this is usually in Hz (if fs is in kHz, then use kHz for fr).\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is the sample number (within X) at the center of the first frame.\n";
descr += "\n";
descr += "Include -m (--mean0) to subtract the mean from each frame [default=no].\n";
descr += "The mean subtraction occurs after windowing.\n";
descr += "\n";
descr += "Use -d (--dim) to give the output dimension of frames [default=0].\n";
descr += "If d=0, then Y is LxT [default], but if d=1, then Y is TxL.\n";
descr += "\n";
descr += "Output (Y) has the same data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ window_univar X W -o Y \n";
descr += "$ window_univar -d1 X W > Y \n";
descr += "$ cat W | window_univar -r200 X > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,W)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_lit  *a_mn0 = arg_litn("m","mean0",0,1,"include to subtract mean from each frame");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"frame dimension of Y [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get fs
fs = (a_fs->count>0) ? a_fs->dval[0] : 10000.0;
if (fs<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fs (sample rate) must be positive" << endl; return 1; }

//Get fr
fr = (a_fr->count>0) ? a_fr->dval[0] : 100.0;
if (fr<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fr (frame rate) must be positive" << endl; return 1; }
if (fr>fs) { cerr << progstr+": " << __LINE__ << errstr << "fr (frame rate) must be <= fs (sample rate)" << endl; return 1; }

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get c0
if (a_c0->count==0) { c0 = 0; }
else if (a_c0->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be nonnegative" << endl; return 1; }
else { c0 = a_c0->ival[0]; }

//Get mn0
mn0 = (a_mn0->count>0);

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) must be a vector" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) found to be empty" << endl; return 1; }
if (i2.N()>=i1.N()) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < length of X" << endl; return 1; }
if (c0>=int(i1.N())) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < length of X" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = 1 + int(floor((int(i1.N())-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? i2.N() : uint32_t(T);
o1.C = (dim==0) ? uint32_t(T) : i2.N();
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *W, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { W = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (W)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(W),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (W)" << endl; return 1; }
    if (ov::window_univar_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),W,dim,c0,float(stp),mn0))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *W, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { W = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (W)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(W),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (W)" << endl; return 1; }
    if (ov::window_univar_c(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),W,dim,c0,float(stp),mn0))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] Y;
}

//Finish

