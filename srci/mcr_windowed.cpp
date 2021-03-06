//Includes
#include "mcr_windowed.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 2, O = 1;
int dim, g, c0, T;
double fs, fr, stp;

//Description
string descr;
descr += "Gets mean crossing rate (MCR) of X along dim.\n";
descr += "Output (Y) has a moving average of MCs using the window W.\n";
descr += "W is the output of a window function, e.g. hamming.\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate of X in Hz [default=10000].\n";
descr += "This can also be entered in kHz, as long as fr (frame rate) is also in kHz.\n";
descr += "\n";
descr += "Use -r (--fr) to give the output frame rate in Hz [default=100.0].\n";
descr += "But if fs is entered in kHz, then also use kHz here.\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is the sample number (within X) at the center of the first frame.\n";
descr += "\n";
descr += "Use -g (--going) to specify positive- or negative-going MCs.\n";
descr += "Use -g0 to detect positive- and negative-going MCs [default].\n";
descr += "Use -g1 to detect only positive-going MCs.\n";
descr += "Use -g-1 to detect only negative-going MCs.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Y has the same size, data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ mcr_windowed X -o Y \n";
descr += "$ mcr_windowed -d1 X > Y \n";
descr += "$ cat X | mcr_windowed -g1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,W)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int    *a_g = arg_intn("g","going","<uint>",0,1,"if using positive- or negative-going MCs [default=0]");
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
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get c0
if (a_c0->count==0) { c0 = 0; }
else if (a_c0->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be nonnegative" << endl; return 1; }
else { c0 = a_c0->ival[0]; }

//Get g
g = (a_g->count>0) ? a_g->ival[0] : 0;
if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be 1D or 2D" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) must be a vector" << endl; return 1; }
if (i2.N()>=i1.N()) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < length of X" << endl; return 1; }
if (dim==0 && i1.R<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
if (dim==1 && i1.C<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = (dim==0) ? 1 + int(floor((int(i1.R)-1-c0)/stp)) : 1 + int(floor((int(i1.C)-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(T) : i1.R;
o1.C = (dim==1) ? uint32_t(T) : i1.C;
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
    if (ov::mcr_windowed_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),W,int(i2.N()),dim,c0,float(stp),g))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] Y;
}

//Finish

