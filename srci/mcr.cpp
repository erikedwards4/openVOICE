//Includes
#include "mcr.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
int dim, g, c0, T, L;
double fs, fr, stp;

//Description
string descr;
descr += "Gets the mean crossing rate (MCR) of X along dim.\n";
descr += "This uses the usual definition (i.e., rectangular window of length L).\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate of X in Hz [default=10000.0].\n";
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
descr += "Y has the data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ mcr X -o Y \n";
descr += "$ mcr -d1 X > Y \n";
descr += "$ cat X | mcr -g1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int    *a_l = arg_intn("l","L","<uint>",0,1,"winlength over which to average MCs [default=1]");
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

//Get L
if (a_l->count==0) { L = 1; }
else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L must be positive" << endl; return 1; }
else { L = a_l->ival[0]; }

//Get g
g = (a_g->count>0) ? a_g->ival[0] : 0;
if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be 1D or 2D" << endl; return 1; }
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
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::mcr_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),L,dim,c0,float(stp),g))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

