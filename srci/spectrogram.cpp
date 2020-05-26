//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 3, O = 1;
int dim, nfft, c0, T;
double fs, fr, stp, p, preg;
char mn0;

//Description
string descr;
descr += "Gets spectrogram of univariate X, which is the STFT power\n";
descr += "transformed to another frequency scale and then compressed.\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate of X in Hz [default=10000].\n";
descr += "This can also be entered in kHz, as long as fr (frame rate) is also in kHz.\n";
descr += "\n";
descr += "Use -r (--fr) to give the frame rate in Hz [default=100.0].\n";
descr += "But if fs is entered in kHz, then also use kHz here.\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is the sample number (within X) at the center of the first frame.\n";
descr += "\n";
descr += "A window W of length L must also be input (e.g. from hamming).\n";
descr += "\n";
descr += "A freq-scale transform matrix, H, must also be input (e.g. from get_spectrogram_T_mat).\n";
descr += "\n";
descr += "Use -n (--nfft) to specify transform length [default=nextpow2(L)].\n";
descr += "Frames of X will be zero-padded as necessary to match nfft.\n";
descr += "\n";
descr += "Use -d (--dim) to give the output FFT dimension [default=0].\n";
descr += "If d=0, then Y is BxT [default]; if d=1, then Y is TxB,\n";
descr += "where: \n";
descr += "B is the number of center frequencies (cfs) = num rows in H;\n";
descr += "and T is the total number of frames (output time points).\n";
descr += "\n";
descr += "Output (Y) has the same data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ spectrogram X W H -o Y \n";
descr += "$ spectrogram -n256 -g1e-6 X W H > Y \n";
descr += "$ cat H | spectrogram -n256 -g1e-6 -r200 X W > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,W,H)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_lit  *a_mn0 = arg_litn("m","mean0",0,1,"include to subtract mean from each frame");
struct arg_int *a_nfft = arg_intn("n","nfft","<uint>",0,1,"transform length [default=nextpow2(L)]");
struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"power compression exponent [default=0.0 -> log]");
struct arg_dbl *a_preg = arg_dbln("g","preg","<dbl>",0,1,"power regularization constant [default=0.0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"frequency dimension of Y [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get fs
fs = (a_fs->count>0) ? a_fs->dval[0] : 10000.0;
if (fs<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fs (sample rate) must be positive" << endl; return 1; }

//Get fr
fr = (a_fr->count>0) ? a_fr->dval[0] : 100.0;
if (fr<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fr (frame rate) must be positive" << endl; return 1; }
if (fr>fs) { cerr << progstr+": " << __LINE__ << errstr << "fr (frame rate) must be <= fs (sample rate)" << endl; return 1; }

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = a_d->ival[0]; }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get c0
if (a_c0->count==0) { c0 = 0; }
else if (a_c0->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be nonnegative" << endl; return 1; }
else { c0 = a_c0->ival[0]; }

//Get mn0
mn0 = (a_mn0->count>0);

//Get nfft
if (a_nfft->count==0) { nfft = 1; while (nfft<int(i2.N())) { nfft *= 2; } } //fastest nextpow2 (I tested)
else if (a_nfft->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = a_nfft->ival[0]; }

//Get p
p = (a_p->count>0) ? a_p->dval[0] : 0.0;
if (p<0.0 || p>1.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be in [0.0 1.0]" << endl; return 1; }

//Get preg
preg = (a_preg->count>0) ? a_preg->dval[0] : 0.0;
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) must be a vector" << endl; return 1; }
if (i2.N()>=i1.N()) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < length of X" << endl; return 1; }
if (c0>=int(i1.N())) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < length of X" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = 1 + int(floor((int(i1.N())-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? i3.R : uint32_t(T);
o1.C = (dim==1) ? i3.R : uint32_t(T);
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *W, *H, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { W = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (W)" << endl; return 1; }
    try { H = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (H)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(W),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (W)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(H),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (H)" << endl; return 1; }
    if (ov::spectrogram_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),W,int(i2.N()),H,dim,c0,float(stp),mn0,nfft,float(p),float(preg)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] H; delete[] Y;
}

//Finish

