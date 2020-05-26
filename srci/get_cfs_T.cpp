//Includes
#include "/home/erik/codee/openvoice/openvoice.h"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 0, O = 1;
int dim, B;
string fscale;
string::size_type c;
double fs;

//Description
string descr;
descr += "Gets cfs for use in making T, the spectrogram transform matrix.\n";
descr += "These place the lofreq lower edge at DC (0 Hz),\n";
descr += "and the hifreq upper edge at Nyquist (fs/2 Hz).\n";
descr += "\n";
descr += "Thus, this is identical to get_cfs with the -e option,\n";
descr += "and with lofreq=0 and hifreq=fs/2.\n";
descr += "Also here the default fscale is 'mel'.\n";
descr += "\n";
descr += "Use -b (--ncfs) to give the number of cfs (B).\n";
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
descr += "$ get_cfs_T -m'mel' -b40 -o Y \n";
descr += "$ get_cfs_T -m'mel' -r20000 -b40 > Y \n";
descr += "$ get_cfs_T -m'bark' -b20 -r20000 -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_dbl   *a_fs = arg_dbln("s","fs","<uint>",0,1,"sample rate in Hz [default=10000.0]");
struct arg_int    *a_b = arg_intn("b","ncfs","<uint>",0,1,"number of frequency bands (cfs) [default=23]");
struct arg_str   *a_sc = arg_strn("m","fscale","<str>",0,1,"frequency scale [default='mel']");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get fscale
if (a_sc->count==0) { fscale = "mel"; }
else
{
	try { fscale = string(a_sc->sval[0]); }
	catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem getting string for frequency scale" << endl; return 1; }
}
for (c=0; c<fscale.size(); c++) { fscale[c] = char(tolower(fscale[c])); }

//Get fs
if (a_fs->count==0) { fs = 10000.0; }
else if (a_fs->dval[0]<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "fs (sample rate) must be positive" << endl; return 1; }
else { fs = a_fs->dval[0]; }

//Get B
if (a_b->count==0) { B = 23; }
else if (a_b->ival[0]<2) { cerr << progstr+": " << __LINE__ << errstr << "B must be a positive int > 1" << endl; return 1; }
else { B = a_b->ival[0]; }

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
o1.R = (dim==0) ? uint32_t(B) : 1u;
o1.C = (dim==1) ? uint32_t(B) : 1u;
o1.S = (dim==2) ? uint32_t(B) : 1u;
o1.H = (dim==3) ? uint32_t(B) : 1u;

//Other prep

//Process
if (o1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (ov::get_cfs_T_s(Y,B,float(fs),fscale.c_str())) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish

