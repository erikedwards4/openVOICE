//Includes
#include "lev_durb.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 2;
int dim, L;

//Description
string descr;
descr += "Does Levinson-Durbin recursion of each row or col of X.\n";
descr += "where X has the autocovariance (AC) functions,\n";
descr += "and output Y holds the linear prediction (LP) coefficients.\n";
descr += "The 2nd output, V, holds the noise variances for each row or col.\n";
descr += "\n";
descr += "Use -l (--L) to specify the number of AC lags used to compute.\n";
descr += "This is the same as P, the number of resulting LP coefficients.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then Y has size L x C.\n";
descr += "If dim==1, then Y has size R x L.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ lev_durb -l127 X -o Y -o V \n";
descr += "$ lev_durb -d1 -l127 X -o Y -o V \n";
descr += "$ cat X | lev_durb -l127 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_l = arg_intn("l","L","<uint>",0,1,"number of lags to use [default uses full size of X]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (Y,V)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get L
if (a_l->count==0) { L = (dim==0) ? int(i1.R) : int(i1.C); }
else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L (nlags) must be positive" << endl; return 1; }
else { L = a_l->ival[0]; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = o2.F = i1.F;
o1.T = o2.T = i1.T;
o1.R = (dim==0) ? uint32_t(L) : i1.R;
o1.C = (dim==1) ? uint32_t(L) : i1.C;
o2.R = (dim==0) ? 1u : i1.R;
o2.C = (dim==1) ? 1u : i1.C;
o1.S = o2.S = i1.S;
o1.H = o2.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y, *V;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y)" << endl; return 1; }
    try { V = new float[o2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::lev_durb_s(Y,V,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,L)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (V)" << endl; return 1; }
    }
    delete[] X; delete[] Y; delete[] V;
}

//Finish

