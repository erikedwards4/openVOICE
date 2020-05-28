//Includes
#include "get_ccs.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
int dim, ndct, K;
double Q;

//Description
string descr;
descr += "Gets cepstral coefficients (CCs) of RxC input matrix X,\n";
descr += "where X is usually a spectrogram.\n";
descr += "This does 1D DCT-II along rows or cols, and then applies lifter.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to transform.\n";
descr += "Use -d0 to operate along cols, and -d1 to operate along rows.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -n (--ndct) to specify transform length [default is R or C].\n";
descr += "X is zero-padded as necessary to match ndct.\n";
descr += "\n";
descr += "Use -q (--Q) to specify the lifter \"bandwidth\".\n";
descr += "The default [Q=22.0] is that of HTK, Kaldi, and D. Ellis.\n";
descr += "Set Q to 0 to skip the lifter.\n";
descr += "\n";
descr += "Use -k (--K) to specify how many CCs to keep at the end.\n";
descr += "The default is to keep all [K=ndct], but K=13 is a typical choice.\n";
descr += "\n";
descr += "The output (Y) is real-valued with size: \n";
descr += "d=0 :   K x C \n";
descr += "d=1 :   R x K \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ get_ccs -n256 X -o Y \n";
descr += "$ get_ccs -n256 -d1 -k13 X > Y \n";
descr += "$ cat X | get_ccs -n256 -k13 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to take DCT [default=0]");
struct arg_int    *a_n = arg_intn("n","ndct","<uint>",0,1,"transform length [default is R or C]");
struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num CCs to keep [default=ndct]");
struct arg_dbl    *a_q = arg_dbln("q","Q","<dbl>",0,1,"lifter bandwidth parameter [default=22.0]");
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

//Get K
if (a_k->count==0) { K = ndct; }
else if (a_k->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positive" << endl; return 1; }
else if (a_k->ival[0]>ndct) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= ndct" << endl; return 1; }
else { K = a_k->ival[0]; }

//Get Q
Q = (a_q->count>0) ? a_q->dval[0] : 22.0;
if (Q<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "Q must be positive" << endl; return 1; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0 && ndct<int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be > nrows of X" << endl; return 1; }
if (dim==1 && ndct<int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be > ncols of X" << endl; return 1; }
if (K>ndct) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= ndct" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(K) : i1.R;
o1.C = (dim==1) ? uint32_t(K) : i1.C;
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
    if (ov::get_ccs_s(Y,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,ndct,float(Q),K)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

