//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 1, O = 1;
int dim;
double fs, tau, spl, maxspl;

//Description
string descr;
descr += "Scales each row or column of X to a target RMS (root-mean-squared) amplitude.\n";
descr += "This works under a basic assumption of conversational speech with reasonable SNR.\n";
descr += "\n";
descr += "Output (Y) has the same size, data type and file format as X,\n";
descr += "but has units of Pascals and is at ~70 dB SPL.\n";
descr += "\n";
descr += "Use -s (--spl) to give the target dB SPL (sound pressure level) [default=70.0].\n";
descr += "Use -m (--maxspl) to give the max dB SPL (sound pressure level) [default=85.0].\n";
descr += "The later is only used if some point in X exceeds this value.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is along cols, unless X is a row vector.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ rms_scale X -o Y \n";
descr += "$ rms_scale -d1 X > Y \n";
descr += "$ cat X | rms_scale -p0.98 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_dbl  *a_tau = arg_dbln("t","tau","<dbl>",0,1,"time constant tau for running RMS estimate [default=0.2]");
struct arg_dbl  *a_spl = arg_dbln("p","spl","<dbl>",0,1,"target dB SPL for output [default=70]");
struct arg_dbl  *a_max = arg_dbln("m","maxspl","<dbl>",0,1,"max dB SPL for output [default=85]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get fs
fs = (a_fs->count>0) ? a_fs->dval[0] : 10000.0;
if (fs<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fs must be positive" << endl; return 1; }

//Get tau
tau = (a_tau->count>0) ? a_tau->dval[0] : 0.2;
if (tau<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "tau must be positive" << endl; return 1; }

//Get spl
spl = (a_spl->count>0) ? a_spl->dval[0] : 70.0;
if (spl<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "spl must be positive" << endl; return 1; }

//Get maxspl
maxspl = (a_max->count>0) ? a_max->dval[0] : 85.0;
if (maxspl<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "maxspl must be positive" << endl; return 1; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0 && i1.R<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
if (dim==1 && i1.C<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (ov::rms_scale_s(X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,float(fs),float(tau),float(spl),float(maxspl)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish

