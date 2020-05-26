//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
uint32_t r, c;
int N, dim;

//Description
string descr;
descr += "Gets 1st order differences (deltas) of X along dim,\n";
descr += "and appends them to X, so that Y has twice the size of X.\n";
descr += "\n";
descr += "The method is as used typically in ASR:\n";
descr += "a linear ramp FIR from -N to N, with edges treated\n";
descr += "as if using numpy.pad(A,N,'edge') (repeat edge sample).\n";
descr += "\n";
descr += "Use -n (--N) to give the number of samps before/after current samp.\n";
descr += "That is, the FIR has length 2*N+1, with values -N to N.\n";
descr += "The default is 2, as used in Kaldi, but 4 also appears typical.\n";
descr += "\n";
descr += "Output (Y) has the same data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ add_deltas X -o Y \n";
descr += "$ add_deltas X > Y \n";
descr += "$ cat X | add_deltas -n4 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"num samps such that FIR length is 2*N+1 [default=2]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get N
if (a_n->count==0) { N = 2; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N must be positive" << endl; return 1; }
else { N = a_n->ival[0]; }

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a matrix" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (dim==0 && N>int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "N must be <= nrows X for dim==0" << endl; return 1; }
if (dim==1 && N>int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "N must be <= ncols X for dim==1" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? i1.R : 2*i1.R;
o1.C = (dim==1) ? i1.C : 2*i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if ((o1.iscolmajor() && dim==0) || (o1.isrowmajor() && dim==1))
    {
        try { ifs1.read(reinterpret_cast<char*>(Y),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    }
    else if (dim==0)
    {
        for (r=0u; r<i1.R; r++)
        {
            try { ifs1.read(reinterpret_cast<char*>(&Y[r*o1.C]),i1.sz()*i1.C); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at row " << r << endl; return 1; }
        }
    }
    else
    {
        for (c=0u; c<i1.C; c++)
        {
            try { ifs1.read(reinterpret_cast<char*>(&Y[c*o1.R]),i1.sz()*i1.R); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at col " << c << endl; return 1; }
        }
    }
    if (ov::add_deltas_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),dim,N)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish

