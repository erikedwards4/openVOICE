//Includes
#include <float.h>
#include "ac2mvdr.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
int dim, F;
double preg;

//Description
string descr;
descr += "Gets multivariate distortionless response (MVDR),\n";
descr += "starting from the autocorrelation (AC) function at each frame.\n";
descr += "\n";
descr += "Use -f (--F) to specify the number of output frequencies at which to compute.\n";
descr += "This determines the nfft used for the FFT of the AC at each frame.\n";
descr += "Since nfft must be power of 2, and F=nfft/2+1, F must be 1 plus a power of 2.\n";
descr += "Also F must be > the number of lags at each AC frame.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then Y has size B x C.\n";
descr += "If dim==1, then Y has size R x B.\n";
descr += "\n";
descr += "A small regularizing power is added to the raw power before log compression.\n";
descr += "Use -p (--preg) to specify this constant [default=10*EPS].\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ ac2mvdr -b16 X -o Y \n";
descr += "$ ac2mvdr -d1 -b16 X -o Y \n";
descr += "$ cat X | ac2mvdr -b20 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_f = arg_intn("f","F","<uint>",0,1,"number of output frequencies to compute [default from X]");
struct arg_dbl *a_preg = arg_dbln("p","preg","<dbl>",0,1,"power regularization constant [default=10*EPS]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get F
if (a_f->count==0)
{
    F = 1;
    if (dim==0) { while (F<int(i1.R)) { F *= 2; } }
    else { while (F<int(i1.C)) { F *= 2; } }
    F++;
}
else if (a_f->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "B (num output bands) must be positive" << endl; return 1; }
else { F = a_f->ival[0]; }

//Get preg
preg = (a_preg->count>0) ? a_preg->dval[0] : 10.0*DBL_EPSILON;
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0 && F<int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "F must be >= nrows X for dim==0" << endl; return 1; }
if (dim==1 && F<int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "F must be >= ncols X for dim==0" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(F) : i1.R;
o1.C = (dim==1) ? uint32_t(F) : i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::ac2mvdr_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,F,float(preg))) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

