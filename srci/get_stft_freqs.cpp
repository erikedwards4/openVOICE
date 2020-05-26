//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 0, O = 1;
int dim, nfft, F;
double fs;

//Description
string descr;
descr += "Gets the linearly-spaced frequencies for an STFT analysis,\n";
descr += "with sample rate fs and FFT length nfft.\n";
descr += "These go from DC (0 Hz) to Nyquist (fs/2 Hz),\n";
descr += "and then through the negative frequencies by default.\n";
descr += "\n";
descr += "Include -p (--positive) to output only the floor(nfft/2) positive freqs.\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate in Hz.\n";
descr += "The default is 2*pi, in which case the output is in radians.\n";
descr += "\n";
descr += "Use -n (--nfft) to give the FFT length.\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "\n";
descr += "Since this is a generating function (no inputs),\n";
descr += "the output data type and file format can be specified by\n";
descr += "-t and -f, respectively (these are the usual CMLI opts).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ get_stft_freqs -n512 -o Y \n";
descr += "$ get_stft_freqs -n512 -s10000 -p > Y \n";
descr += "$ get_stft_freqs -n128 -s10000 -p -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_dbl   *a_fs = arg_dbln("s","fs","<uint>",0,1,"sample rate in Hz [default=2*pi]");
struct arg_int   *a_nf = arg_intn("n","nfft","<uint>",0,1,"FFT length [default=0]");
struct arg_lit   *a_po = arg_litn("p","positive",0,1,"include to output only the positive freqs");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get fs
if (a_fs->count==0) { fs = 2.0*M_PI; }
else if (a_fs->dval[0]<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fs (sample rate) must be positive" << endl; return 1; }
else { fs = a_fs->dval[0]; }

//Get nfft
if (a_nf->count==0) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be specified with -n (--nfft) opt" << endl; return 1; }
else if (a_nf->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = a_nf->ival[0]; }

//Get positive and F
F = (a_po->count>0) ? nfft/2+1 : nfft;

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

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

//Checks

//Set output header info
o1.R = (dim==0) ? uint32_t(F) : 1u;
o1.C = (dim==1) ? uint32_t(F) : 1u;
o1.S = (dim==2) ? uint32_t(F) : 1u;
o1.H = (dim==3) ? uint32_t(F) : 1u;

//Other prep

//Process
if (o1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (ov::get_stft_freqs_s(Y,int(o1.N()),float(fs),nfft)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish

