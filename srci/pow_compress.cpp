//Includes
#include "pow_compress.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
double p, preg;

//Description
string descr;
descr += "Applies power compression to each element of X:\n";
descr += "For p>0 : Y = pow(X+preg,p) = X.^p \n";
descr += "For p==0: Y = log(X+preg) \n";
descr += "where p is in [0.0 1.0].\n";
descr += "where preg is a small regularization constant [default=0.0].\n";
descr += "\n";
descr += "X is therefore assumed to be nonnegative and real-valued.\n";
descr += "The output (Y) has the same size, data type, and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ pow_compress -p0.5 X -o Y \n";
descr += "$ pow_compress -p0.1 X > Y \n";
descr += "$ cat X | pow_compress -r1e-6 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"power exponent [default=0.0 -> log]");
struct arg_dbl *a_preg = arg_dbln("r","preg","<dbl>",0,1,"power regularization constant [default=0.0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get p
p = (a_p->count>0) ? a_p->dval[0] : 0.0;
if (p<0.0 || p>1.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be in [0.0 1.0]" << endl; return 1; }

//Get preg
preg = (a_preg->count>0) ? a_preg->dval[0] : 0.0;
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (ov::pow_compress_s(X,int(i1.N()),float(p),float(preg)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish

