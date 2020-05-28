//@author Erik Edwards
//@date 2019-2020


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <complex>
#include <unordered_map>
#include <argtable2.h>
#include "/home/erik/codee/cmli/cmli.hpp"
#include "frame_univar.c"

#ifdef I
#undef I
#endif


int main(int argc, char *argv[])
{
    using namespace std;


    //Declarations
    int ret = 0;
    const string errstr = ": \033[1;31merror:\033[0m ";
    const string warstr = ": \033[1;35mwarning:\033[0m ";
    const string progstr(__FILE__,string(__FILE__).find_last_of("/")+1,strlen(__FILE__)-string(__FILE__).find_last_of("/")-5);
    const valarray<uint8_t> oktypes = {1,2,101,102};
    const size_t I = 1, O = 1;
    ifstream ifs1; ofstream ofs1;
    int8_t stdi1, stdo1, wo1;
    ioinfo i1, o1;
    int dim;
    double fs, fr, stp;
    int L, c0, T;


    //Description
    string descr;
    descr += "Takes univariate X and produces a series of (overlapping) frames.\n";
    descr += "The output Y has size LxT or TxL, where L is the length of each frame, \n";
    descr += "and T is the total number of frames (output time points).\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the output dimension of frames [default=0].\n";
    descr += "If d=0, then Y is LxT [default], but if d=1, then Y is TxL.\n";
    descr += "\n";
    descr += "Use -l (--winlength) to give L, the length of each frame.\n";
    descr += "The default is the closest odd integer to 25-ms frame length.\n";
    descr += "However, only rely on the default if fs is in Hz.\n";
    descr += "\n";
    descr += "Use -r (--fr) to give the frame rate in Hz [default=100.0].\n";
    descr += "This actually has to be in the same units as fs (sample rate),\n";
    descr += "but this is usually in Hz (if fs is in kHz, then use kHz for fr).\n";
    descr += "\n";
    descr += "Use -c (--c0) to give the first center sample [default=0].\n";
    descr += "This is the sample number (within X) at the center of the first frame.\n";
    descr += "\n";
    descr += "Output (Y) has the same data type and file format as X,\n";
    descr += "with data efficiently copied from X into Y.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ frame_univar -l255 X -o Y \n";
    descr += "$ frame_univar -l255 -d1 X > Y \n";
    descr += "$ cat X | frame_univar -l127 -r200 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_dbl   *a_fs = arg_dbln("s","fs","<dbl>",0,1,"sample rate of X [default=10000]");
    struct arg_dbl   *a_fr = arg_dbln("r","fr","<dbl>",0,1,"frame rate of Y [default=100]");
    struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"length in samps of each frame");
    struct arg_int   *a_c0 = arg_intn("c","c0","<uint>",0,1,"first center sample [default=0]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"time dimension of Y [default=0]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_fs, a_fr, a_wl, a_c0, a_d, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdin
    stdi1 = (a_fi->count==0 || strlen(a_fi->filename[0])==0 || strcmp(a_fi->filename[0],"-")==0);
    if (stdi1>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Open input
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file" << endl; return 1; }


    //Read input header
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }


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

    //Get L
    if (a_wl->count==0)
    {
        L = int(ceil(0.025*fs));
        if (L%2==0) { L = int(floor(0.025*fs)); }
        if (L%2==0) { L++; }
    }
    else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be positive" << endl; return 1; }
    else { L = a_wl->ival[0]; }
    if (L<1 || L>int(i1.N())-1) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be positive" << endl; return 1; }


    //Checks
    if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
    if (c0>=int(i1.N())) { cerr << progstr+": " << __LINE__ << errstr << "c0 must be < length of X" << endl; return 1; }


    //Set output header info
    stp = fs/fr;
    T = 1 + int(floor((int(i1.N())-1-c0)/stp));
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (dim==0) ? uint32_t(L) : uint32_t(T);
    o1.C = (dim==0) ? uint32_t(T) : uint32_t(L);
    o1.S = i1.S; o1.H = i1.H;


    //Open output
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }


    //Write output header
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }


    //Other prep


    //Process
    if (i1.T==1)
    {
        float *X, *Y;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::frame_univar_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,float(stp)))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==2)
    {
        double *X, *Y;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::frame_univar_d(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,double(stp)))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==101)
    {
        float *X, *Y;
        try { X = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::frame_univar_c(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,float(stp)))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==102)
    {
        double *X, *Y;
        try { X = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::frame_univar_z(Y,o1.iscolmajor(),int(o1.R),int(o1.C),X,int(i1.N()),dim,c0,double(stp)))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

