//Includes
#include "frame_univar.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
int dim;
double fs, fr, stp;
int L, c0, T;

//Description
string descr;
descr += "Takes univariate X and produces a series of (overlapping) frames.\n";
descr += "The output Y has size LxT or TxL, where L is the length of each frame, \n";
descr += "and T is the total number of frames (output time points).\n";
descr += "\n";
descr += "Use -d (--dim) to give the output dimension of frames [default=0].\n";
descr += "If d=0, then Y is LxT [default], but if d=1, then Y is TxL.\n";
descr += "\n";
descr += "Use -l (--winlength) to give L, the length of each frame.\n";
descr += "The default is the closest odd integer to 25-ms frame length.\n";
descr += "However, only rely on the default if fs is in Hz.\n";
descr += "\n";
descr += "Use -r (--fr) to give the frame rate in Hz [default=100.0].\n";
descr += "This actually has to be in the same units as fs (sample rate),\n";
descr += "but this is usually in Hz (if fs is in kHz, then use kHz for fr).\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is the sample number (within X) at the center of the first frame.\n";
descr += "\n";
descr += "Output (Y) has the same data type and file format as X,\n";
descr += "with data efficiently copied from X into Y.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ frame_univar -l255 X -o Y \n";
descr += "$ frame_univar -l255 -d1 X > Y \n";
descr += "$ cat X | frame_univar -l127 -r200 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"length in samps of each frame");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"time dimension of Y [default=0]");
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

//Get L
if (a_wl->count==0)
{
    L = int(ceil(0.025*fs));
    if (L%2==0) { L = int(floor(0.025*fs)); }
    if (L%2==0) { L++; }
}
else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be positive" << endl; return 1; }
else { L = a_wl->ival[0]; }
if (L<1 || L>int(i1.N())-1) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be positive" << endl; return 1; }

//Checks
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (c0>=int(i1.N())) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < length of X" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = 1 + int(floor((int(i1.N())-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(L) : uint32_t(T);
o1.C = (dim==0) ? uint32_t(T) : uint32_t(L);
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::frame_univar_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,float(stp)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::frame_univar_c(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,float(stp)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

