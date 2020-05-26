//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 0, O = 1;
int L, dim;
char normalize;

//Description
string descr;
descr += "Gets 1D Blackman-Harris window.\n";
descr += "\n";
descr += "Use -l (--winlength) to give the window length in number of taps (sample points).\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "\n";
descr += "Include -n (--normalize) to normalize the output such that it sums to 1.\n";
descr += "In this case, it represents a moving-average window.\n";
descr += "\n";
descr += "Since this is a generating function (no inputs),\n";
descr += "the output data type and file format can be specified by\n";
descr += "-t and -f, respectively (these are the usual CMLI opts).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ blackmanharris -l255 -o Y \n";
descr += "$ blackmanharris -l255 > Y \n";
descr += "$ blackmanharris -l127 -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"window length [default=7]");
struct arg_lit  *a_nrm = arg_litn("n","normalize",0,1,"normalize output to sum to 1");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 102; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = uint8_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 2; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = uint8_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get L
if (a_wl->count==0) { L = 7; }
else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "winlength must be a positive int" << endl; return 1; }
else { L = a_wl->ival[0]; }

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get normalize
normalize = (a_nrm->count>0);

//Checks

//Set output header info
o1.R = (dim==0) ? uint32_t(L) : 1u;
o1.C = (dim==1) ? uint32_t(L) : 1u;
o1.S = (dim==2) ? uint32_t(L) : 1u;
o1.H = (dim==3) ? uint32_t(L) : 1u;

//Other prep

//Process
if (o1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (ov::blackmanharris_s(Y,L,normalize)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}
else if (o1.T==101)
{
    float *Y;
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (ov::blackmanharris_c(Y,L,normalize)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish

