//Includes
#include "iir.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 2, O = 1;
uint32_t n;
int dim;

//Description
string descr;
descr += "Mth-order causal IIR filter:\n";
descr += "Y(t) = X(t) - a1*Y(t-1) -a2*Y(t-2) ... \n";
descr += "\n";
descr += "The first element of A should be a0==1.\n";
descr += "\n";
descr += "Output (Y) has the same size and data type as X.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is along cols, unless X is a row vector.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ iir X A -o Y \n";
descr += "$ iir -d1 A > Y \n";
descr += "$ cat X | iir A > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,A)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Checks
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (A) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be 1D or 2D" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (A) must be a vector" << endl; return 1; }
if (dim==0 && i1.R<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
if (dim==1 && i1.C<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *A;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { A = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (A)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(A),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (A)" << endl; return 1; }
    if (A[0]!=1.0f) { cerr << progstr+": " << __LINE__ << warstr << "A[0] (1st element of A coeffs) is not 1" << endl; }
    for (n=0u; n<i2.N()/2u; n++) { swap(A[n],A[i2.N()-1u-n]); }
    if (ov::iir_s(X,i1.iscolmajor(),int(i1.R),int(i1.C),A,int(i2.N()),dim)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] A;
}
else if (i1.T==101)
{
    float *X, *A;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { A = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (A)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(A),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (A)" << endl; return 1; }
    if (A[0]!=1.0f) { cerr << progstr+": " << __LINE__ << warstr << "A[0] (1st element of A coeffs) is not 1" << endl; }
    for (n=0u; n<i2.N()/2u; n++) { swap(A[2u*n],A[2u*i2.N()-2u-2u*n]); swap(A[2u*n+1u],A[2u*i2.N()-1u-2u*n]); }
    if (ov::iir_c(X,i1.iscolmajor(),int(i1.R),int(i1.C),A,int(i2.N()),dim)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] A;
}

//Finish

