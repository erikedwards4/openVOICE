//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
int dim, ndct;

//Description
string descr;
descr += "Does 1D DCT-II along rows or cols of RxC input matrix X.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to transform.\n";
descr += "Use -d0 to operate along cols, and -d1 to operate along rows.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -n (--ndct) to specify transform length [default is R or C].\n";
descr += "X is zero-padded as necessary to match ndct.\n";
descr += "\n";
descr += "The output (Y) is real-valued with size: \n";
descr += "d=0 :   ndct x C \n";
descr += "d=1 :   R x ndct \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ dct -n256 X -o Y \n";
descr += "$ dct -n256 -d1 X > Y \n";
descr += "$ cat X | dct -n256 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to take DCT [default=0]");
struct arg_int    *a_n = arg_intn("n","ndct","<uint>",0,1,"transform length [default is R or C]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get ndct
if (a_n->count==0) { ndct = (dim==0) ? int(i1.R) : int(i1.C); }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be positive" << endl; return 1; }
else { ndct = a_n->ival[0]; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0 && ndct<int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be > nrows of X" << endl; return 1; }
if (dim==1 && ndct<int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be > ncols of X" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(ndct) : i1.R;
o1.C = (dim==1) ? uint32_t(ndct) : i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (ov::dct_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,ndct)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

