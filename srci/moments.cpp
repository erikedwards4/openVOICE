//Includes
#include "moments.c"

//Declarations
const valarray<uint8_t> oktypes = {1};
const size_t I = 1, O = 1;
int dim;
double fs, fr, stp;
int c0, T;

//Description
string descr;
descr += "Gets the first M standardized statistical moments of each row or col of X.\n";
descr += "These are the mean, variance, skewness, kurtosis, etc.\n";
descr += "This operates on each frame of length L individually.\n";
descr += "\n";
descr += "Use -d (--dim) to give the main time dimension along which to operate.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate of X in Hz [default=100.0].\n";
descr += "Use -r (--fr) to give the frame rate of Y in Hz [default=5.0].\n";
descr += "\n";
descr += "Use -l (--winlength) to give L, the length of each frame.\n";
descr += "The default is the closest odd integer to 0.5 sec frame length.\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is an integer in units of samples (within X).\n";
descr += "\n";
descr += "Use -m (--M) to give the number of statistical moments to compute [default=4].\n";
descr += "\n";
descr += "If d=0, then X has size NxP; and if d=1, then X has size PxN,\n";
descr += "where N is the number of samples, and P is the number of features.\n";
descr += "If d=0, then Y has size Tx(P*M); and if d=1, then Y has size (P*M)xT,\n";
descr += "where T is the number of frames (as determined by fs, fr and c0).\n";
descr += "\n";
descr += "Within Y, the first block of size PxT or TxP has the first functional\n";
descr += "for each of the P original features. The second block of size PxT or TxP\n";
descr += "has the second functional, and so on.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ moments -r10 -l51 X -o Y \n";
descr += "$ moments -d1 X > Y \n";
descr += "$ cat W | moments -r10 X > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_m = arg_intn("m","M","<uint>",0,1,"number of moments to compute [default=4]");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X in Hz [default=100]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y in Hz [default=5]");
struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"length of each frame [default is ~0.5 s]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"main time dimension of X and Y [default=0]");
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

//Get m
if (a_m->count==0) { M = 4; }
else if (a_m->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "m (num moments) must be positive" << endl; return 1; }
else { M = a_m->ival[0]; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a matrix" << endl; return 1; }
if (dim==0 && L>=int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < num samps in X" << endl; return 1; }
if (dim==1 && L>=int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < num samps in X" << endl; return 1; }
if (dim==0 && c0>=int(i1.R) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < num samps of X" << endl; return 1; }
if (dim==1 && c0>=int(i1.C) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < num samps of X" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = 1 + int(floor((int(i1.N())-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(T) : i1.R*M;
o1.C = (dim==1) ? uint32_t(T) : i1.C*M;
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
    if (ov::moments_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,float(stp)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

