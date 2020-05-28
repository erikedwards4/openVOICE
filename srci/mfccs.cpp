//Includes
#include "mfccs.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 3, O = 1;
int c0, T, dim, nfft, ndct, K;
double fs, fr, stp, Q, preg;
char mn0;

//Description
string descr;
descr += "Gets MFCCs (mel frequency cepstral coeffs) of univariate X.\n";
descr += "\n";
descr += "Use -s (--fs) to give the sample rate of X in Hz [default=10000].\n";
descr += "\n";
descr += "Use -r (--fr) to give the final frame rate in Hz [default=100.0].\n";
descr += "\n";
descr += "Use -c (--c0) to give the first center sample [default=0].\n";
descr += "This is the sample number (within X) at the center of the first frame.\n";
descr += "\n";
descr += "The framing of X is controlled by c0, frame rate, and\n";
descr += "input W, which is a window of length L (e.g. from hamming).\n";
descr += "Each frame is transformed by an nfft-point FFT (nfft>=L).\n";
descr += "\n";
descr += "Use -f (--nfft) to specify transform length [default=nextpow2(L)].\n";
descr += "Frames of X will be zero-padded as necessary to match nfft.\n";
descr += "\n";
descr += "A small regularizing power is added to the raw FFT power before log compression.\n";
descr += "Use -p (--preg) to specify this constant [default=10*EPS].\n";
descr += "\n";
descr += "The FFT power is transformed by the third input, H, to a mel\n";
descr += "frequency scale with B bands, as would appear in a mel spectrogram..\n";
descr += "\n";
descr += "The mel spectrogram is subjected to an ndct-point DCT-II transform.\n";
descr += "Use -n (--ndct) to specifty the DCT transform length [default=B].\n";
descr += "\n";
descr += "Use -q (--Q) to specify the lifter \"bandwidth\".\n";
descr += "The default [Q=22.0] is that of HTK, Kaldi, and D. Ellis.\n";
descr += "Set Q to 0 to skip the lifter.\n";
descr += "\n";
descr += "Use -k (--K) to specify how many CCs to keep at the end [default=13].\n";
descr += "\n";
descr += "Use -d (--dim) to give the output cepstral dimension [default=0].\n";
descr += "If d=0, then Y is KxT [default]; if d=1, then Y is TxK,\n";
descr += "\n";
descr += "Output (Y) has the same data type and file format as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ mfccs X W H -o Y \n";
descr += "$ mfccs -n256 X W H > Y \n";
descr += "$ cat H | mfccs -n256 -p1e-6 -r200 X W > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,W,H)");
struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
struct arg_lit  *a_mn0 = arg_litn("m","mean0",0,1,"include to subtract mean from each frame");
struct arg_int *a_nfft = arg_intn("f","nfft","<uint>",0,1,"FFT transform length [default=nextpow2(L)]");
struct arg_dbl *a_preg = arg_dbln("p","preg","<dbl>",0,1,"power regularization constant [default=10*EPS]");
struct arg_int *a_ndct = arg_intn("n","ndct","<uint>",0,1,"DCT transform length [default is B]");
struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num CCs to keep [default=13]");
struct arg_dbl    *a_q = arg_dbln("q","Q","<dbl>",0,1,"lifter bandwidth parameter [default=22.0]");
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
if (a_nfft->count==0) { nfft = 1; while (nfft<int(i2.N())) { nfft *= 2; } }
else if (a_nfft->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = a_nfft->ival[0]; }

//Get preg
preg = (a_preg->count>0) ? a_preg->dval[0] : 10.0*DBL_EPSILON;
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Get ndct
if (a_ndct->count==0) { ndct = int(i3.R); }
else if (a_ndct->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be positive" << endl; return 1; }
else { ndct = a_ndct->ival[0]; }

//Get K
if (a_k->count==0) { K = 13; }
else if (a_k->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positive" << endl; return 1; }
else if (a_k->ival[0]>ndct) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= ndct" << endl; return 1; }
else { K = a_k->ival[0]; }

//Get Q
Q = (a_q->count>0) ? a_q->dval[0] : 22.0;
if (Q<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "Q must be positive" << endl; return 1; }

//Checks
if (i1.T!=i2.T || i1.T!=i3.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) must be a vector" << endl; return 1; }
if (!i3.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (H) must be a matrix" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) found to be empty" << endl; return 1; }
if (i3.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (H) found to be empty" << endl; return 1; }
if (i2.N()>=i1.N()) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be < length of X" << endl; return 1; }
if (c0>=int(i1.N())) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < length of X" << endl; return 1; }
if (ndct>nfft) { cerr << progstr+": " << __LINE__ << errstr << "ndct must be <= nfft" << endl; return 1; }
if (K>ndct) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= ndct" << endl; return 1; }

//Set output header info
stp = fs/fr;
T = 1 + int(floor((int(i1.N())-1-c0)/stp));
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? uint32_t(K) : uint32_t(T);
o1.C = (dim==1) ? uint32_t(K) : uint32_t(T);
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
    if (ov::mfccs_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),W,int(i2.N()),H,int(i3.R),dim,c0,float(stp),mn0,nfft,float(preg),ndct,float(Q),K))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] H; delete[] Y;
}

//Finish

